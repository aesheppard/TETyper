#!/usr/bin/env python
"""
TETyper: a tool for typing transposable elements from whole-genome sequencing data
Copyright (C) 2018 Anna Sheppard

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import sys
import os
import logging
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
import pysam
import vcf
from collections import Counter

VERSION = '1.1'


class ProfileError(Exception):
    def __init__(self, value):
        self.value = value


class ProfileMatcher:
    def __init__(self, profile_file, nonestring, delim, sep='\t', start_pos=None, end_pos=None):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.sep = sep
        self.delim = delim
        self.nonestring = nonestring
        self.profile_dict = {}
        with open(profile_file) as profile_handle:
            for line in profile_handle:
                if len(line.strip()) > 0:  # If a line is only white space, simply ignore it
                    strippedline = line.rstrip('\n')
                    if self.sep not in strippedline:
                        raise ProfileError('"{0}" does not contain separate columns for name and profile'.format(strippedline))
                    profile_name, profile = strippedline.split(self.sep, 1)
                    if profile in self.profile_dict:
                        raise ProfileError('Duplicate entries for profile "{0}"'.format(profile))
                    self.validate_profile(profile)
                    self.profile_dict[profile] = profile_name

    def get_profile(self, profile, default=None):
        # returns corresponding profile name if there's a match, otherwise None/default
        if profile in self.profile_dict:
            return self.profile_dict[profile]
        return default

    def validate_profile(self, profile):
        pass
        

class StructProfileMatcher(ProfileMatcher):
    def validate_profile(self, profile):
        if profile == self.nonestring:
            return
        # each element should be: start-end, where start and end are (positive) integers and within range (if specified)
        elements = profile.split(self.delim)
        prev_start, prev_end = None, None
        for element in elements:
            try:
                start, end = map(int, element.split('-', 1))
            except ValueError:  # covers splitting with inappropriate number of elements and integer conversion
                raise ProfileError('"{0}" is not a valid range'.format(element))
            if start > end:
                raise ProfileError('"{0}" is not a valid range'.format(element))
            if (self.start_pos is not None and start < self.start_pos) or (self.end_pos is not None and end > self.end_pos):
                raise ProfileError('"{0}" is not contained within the allowable range of {1}'.format(element, str(self.start_pos) + '-' + str(self.end_pos)))
            # elements should be ordered by start and ranges should be non-adjacent and non-overlapping (i.e. only most simplified version is valid)
            if prev_end is not None and start <= prev_end + 1:
                raise ProfileError('Incorrectly ordered or overlapping ranges: "{0}" and "{1}"'.format(str(prev_start) + '-' + str(prev_end), str(start) + '-' + str(end)))
            prev_start, prev_end = start, end
        
                                                
class SNPProfileMatcher(ProfileMatcher):
    def validate_profile(self, profile):
        try:
            var_element_string, N_element_string = profile.split(self.sep)
        except ValueError:
            raise ProfileError('SNP profile "{0}" does not contain the required two columns'.format(profile))
        if var_element_string != self.nonestring:
            var_elements = var_element_string.split(self.delim)
            self.validate_snps(var_elements)
        if N_element_string != self.nonestring:
            N_elements = N_element_string.split(self.delim)
            self.validate_snps(N_elements, Nsite=True)

    def validate_snps(self, snps, Nsite=False):
        prev_pos = None
        for snp in snps:
            try:
                ref, pos, alt = snp[0], int(snp[1:-1]), snp[-1]
            except (IndexError, ValueError):
                raise ProfileError('"{0}" is not a valid snp'.format(snp))
            if pos < self.start_pos or pos > self.end_pos:
                raise ProfileError('Error parsing SNP "{0}": Position "{1}" is not contained within the allowable range of {2}'.format(snp, pos, str(self.start_pos) + '-' + str(self.end_pos)))
            bases = ['A','C','G','T']
            if ref not in bases:
                raise ProfileError('Error parsing SNP "{0}": "{1}" is not a valid reference base. Allowed values are: A,C,G,T'.format(snp, ref))
            if Nsite == False and alt not in bases:
                raise ProfileError('Error parsing SNP "{0}": "{1}" is not a valid homozygous SNP call. Allowed values are: A,C,G,T'.format(snp, alt))
            if Nsite == True and alt not in ['M','R','W','S','Y','K']:
                raise ProfileError('Error parsing SNP "{0}": "{1}" is not a valid heterozygous SNP call. Allowed values are: M,R,W,S,Y,K'.format(snp, alt))
            if prev_pos is not None and pos <= prev_pos:
                raise ProfileError('Incorrect ordering for "{0}" and "{1}". SNPs should be ordered by position.'.format(prev_pos, pos))



class TETyper:

    DELIM = '|'
    NONESTRING = 'none'
    UNKSTRING = 'unknown'

    def __init__(self, args):
        loglevels = {1:logging.ERROR, 2:logging.WARNING, 3:logging.INFO, 4:logging.DEBUG}
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        stderrhandler = logging.StreamHandler()
        stderrhandler.setLevel(loglevels[args.verbosity])
        stderrhandler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(stderrhandler)
        self.outprefix = args.outprefix
        logfile = self.outprefix + '.log'
        self.loghandle = None
        if args.no_overwrite and os.path.isfile(logfile):
            self.exitonerror('File {0} already exists. Exiting.'.format(logfile))       
        try:
            self.loghandle = open(logfile, 'w')
        except IOError as e:
            self.exitonerror('Error opening file {0}: {1}'.format(logfile, e))

        logfilehandler = logging.StreamHandler(self.loghandle)
        logfilehandler.setLevel(logging.DEBUG)
        logfilehandler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
        logger.addHandler(logfilehandler)

        logging.info('TETyper version {0}'.format(VERSION))
        logging.info('TETyper command: {0}'.format(' '.join(sys.argv)))

        self.ref = args.ref
        try:
            contigs = {contig.name: len(contig) for contig in SeqIO.parse(self.ref, 'fasta')}
        except IOError as e:
            self.exitonerror('Error reading reference file {0}: {1}'.format(self.ref, e))
        if len(contigs) != 1:
            self.exitonerror('Error reading reference file {0}: Expected 1 contig, found {1}'.format(self.ref, len(contigs)))
        self.ref_start = 1
        self.ref_contig, self.ref_end = list(contigs.items())[0]  
        logging.info('Successfully read in reference of length {0}'.format(self.ref_end))

        # if blast database is explicitly specified, make sure it actually exists
        if args.refdb:
            self.refdb = args.refdb
            self.check_file(self.refdb + '.nin')
        # otherwise check for blast database with default naming and create one if it doesn't already exist
        else:
            self.refdb = self.ref
            if not os.path.isfile(self.refdb + '.nin'):
                logging.info('Blast database not found. Creating database automatically.')
                self.run_external_call([['makeblastdb', '-dbtype', 'nucl', '-in', self.ref]], 'blast database creation')

        if args.bam:
            self.bam = args.bam
            self.bam_provided = True
        else:
            self.fq_files = [fq for fq in [args.fq1, args.fq2, args.fq12] if fq is not None]
            self.interleaved = True if args.fq12 else False
            if not ((len(self.fq_files) == 2 and args.fq1 and args.fq2) or (len(self.fq_files) == 1 and self.interleaved == True)):
                self.exitonerror('Invalid input. Exactly one of the following must be provided: (--fq1 AND --fq2) OR --fq12 OR --bam')
            self.bam_provided = False

        if not self.bam_provided and not os.path.isfile(self.ref + '.bwt'):
            logging.info('Bwa index files not found. Creating index automatically.')
            self.run_external_call([['bwa', 'index', self.ref]], 'bwa indexing')

        if args.assembly:
            self.assembly = args.assembly
            self.assembly_provided = True
        else:
            self.assembly_provided = False

        self.struct_profiles = args.struct_profiles
        self.snp_profiles = args.snp_profiles
        self.flank_len = args.flank_len
        self.min_reads = args.min_reads
        self.min_each_strand = args.min_each_strand
        self.min_mapped_len = args.min_mapped_len
        self.min_qual = args.min_qual

        if args.show_region:
            self.show_region = args.show_region
            try:
                self.show_region_start, self.show_region_end = map(int, args.show_region.split('-'))
            except ValueError:
                self.exitonerror('--show_region must be in the format "startpos-endpos"')
            if self.show_region_start < self.ref_start or self.show_region_end > self.ref_end or self.show_region_start > self.show_region_end:
                self.exitonerror('--show_region must be a valid range within the given reference')
        else:
            self.show_region = None
        
        if args.threads <= 0:
            self.exitonerror('--threads must be a positive integer.')
        self.threads = args.threads


    def exitonerror(self, errorstring):
        logging.error(errorstring)
        self.cleanup()
        sys.exit(1)


    def cleanup(self):
        if self.loghandle is not None:
            self.loghandle.close()
        logging.shutdown()


    def check_file(self, filename, checkexist='error', checkempty='warning'):
        if checkexist is not None and not os.path.isfile(filename):
            existmsg = 'File "{0}" does not exist'.format(filename)
            if checkexist == 'warning':
                logging.warning(existmsg)
            elif checkexist == 'error':
                self.exitonerror(existmsg)
            else:
                self.exitonerror('Unexpected value for checkexist: {0}'.format(checkexist))
            return False
        elif checkempty is not None and os.path.getsize(filename) == 0:
            emptymsg = 'File "{0}" is empty'.format(filename)
            if checkempty == 'warning':
                logging.warning(emptymsg)
            elif checkempty == 'error':
                self.exitonerror(emptymsg)
            else:
                self.exitonerror('Unexpected value for checkempty: {0}'.format(checkempty))
            return False
        return True


    def check_return_code(self, program_name, retcode):
        if retcode != 0:
            self.exitonerror('Command {0} returned non-zero exit status {1}'.format(program_name, retcode))


    def run_external_call(self, proc_args, stage_name, outfile=None):
        logging.debug('Running {0} using the following command: "{1}"'.format(stage_name, ' | '.join([' '.join(arglist) for arglist in proc_args])))
        try:
            procs = []
            if len(proc_args) == 1:
                procs.append(Popen(proc_args[0], stdout=self.loghandle, stderr=STDOUT))
            else:
                procs.append(Popen(proc_args[0], stdout=PIPE, stderr=self.loghandle))  # first in pipe
                for args in proc_args[1:-1]:  # everything in the middle
                    procs.append(Popen(args, stdin=procs[-1].stdout, stdout=PIPE, stderr=self.loghandle))
                procs.append(Popen(proc_args[-1], stdin=procs[-1].stdout, stdout=self.loghandle, stderr=STDOUT))  # last in pipe
                for proc in procs[:-1]:
                    proc.stdout.close()
            for proc,args in zip(procs[::-1], proc_args[::-1]):
                proc.wait()
                progname = args[0] if args[1][0] == '-' else ' '.join(args[0:2])
                self.check_return_code(progname, proc.returncode)
        except OSError:
            proglist = list(set([args[0] for args in proc_args]))
            if len(proglist) == 1:
                progs = proglist[0] + ' is'
            else:
                progs = ', '.join(proglist[:-1]) + ' and ' + proglist[-1] + ' are'
            self.exitonerror('Error executing {0}. Check that {1} available via $PATH.'.format(stage_name, progs))

        outputmessage = '' if outfile is None else 'Output stored in {0}'.format(outfile)
        logging.debug('{0} completed successfully. {1}'.format(stage_name[:1].upper() + stage_name[1:], outputmessage))


    def do_map(self):
        if self.bam_provided:
            logging.info('Skipping mapping step and using bam file provided: ' + self.bam)
            self.check_file(self.bam)
            return

        # Build up bwa mem command
        self.bam = self.outprefix + '.bam'
        bwa_mem_args = ['bwa', 'mem', '-t', str(self.threads)]
        if self.interleaved:
            bwa_mem_args.append('-p')
        bwa_mem_args.append(self.ref)
        for fq_file in self.fq_files:
            self.check_file(fq_file)
            bwa_mem_args.append(fq_file)

        samtools_view_args = ['samtools', 'view', '-bu', '-F2048', '-G12', '--threads', str(self.threads), '-']
        samtools_sort_args = ['samtools', 'sort', '-o', self.bam , '--threads', str(self.threads), '-']

        self.run_external_call([bwa_mem_args, samtools_view_args, samtools_sort_args], 'mapping', self.bam)


    def do_assembly(self):
        if self.assembly_provided:
            logging.info('Skipping assembly step and using assembly file provided: {0}'.format(self.assembly))
            return

        # Convert mapped reads to fastq format
        samtools_namesort_args = ['samtools', 'sort', '-n', '--threads', str(self.threads), self.bam]
        samtools_fastq_args = ['samtools', 'fastq', '--threads', str(self.threads)]
        fq_prefix = self.outprefix + '_mappedreads'
        fq_args = ['-1', fq_prefix + '_1.fq', '-2', fq_prefix + '_2.fq']  # Same for both samtools fastq and spades
        samtools_fastq_args.extend(fq_args)
        samtools_fastq_args.append('-')
        fq_files = fq_args[1::2]
        self.run_external_call([samtools_namesort_args, samtools_fastq_args], 'fastq conversion', ' and '.join(fq_files))
        for fq_file in fq_files:
            self.check_file(fq_file)

        # Run assembly
        spadesdir = self.outprefix + '_spades'
        spades_args = ['spades.py', '-t', str(self.threads), '-o', spadesdir] + fq_args + args.spades_params.split()
        self.run_external_call([spades_args], 'spades assembly', spadesdir + '/')
        self.assembly = spadesdir + '/contigs.fasta'


    def do_blast(self):
        if not self.check_file(self.assembly, checkexist='warning'):
            logging.warning('Assuming zero length assembly and proceeding anyway')
            logging.info('Skipping blast step')
            self.blastfile = None
            return

        self.blastfile = self.outprefix + '_blast.txt'
        blastn_args = ['blastn', '-db', self.refdb, '-query', self.assembly, '-outfmt', '6', '-out', self.blastfile, '-num_threads', str(self.threads)]
        self.run_external_call([blastn_args], 'blastn', self.blastfile) 


    def parse_blast(self):
        hitranges = []  # list of tuples containing ranges covered by blast hits
        if self.blastfile is not None:
            logging.debug('Parsing blast output...')
            with open(self.blastfile) as blastfilehandle:
                for blasthit in blastfilehandle:
                    blasthitfields = blasthit.strip().split()
                    sstart, send = int(blasthitfields[8]), int(blasthitfields[9])
                    hitranges.append((min(sstart, send), max(sstart, send))) 
        hitranges.sort()  # sort by hit start position

        # merge hit overlaps
        mergedhitranges = []
        startpos, endpos = None, None
        for (currstart, currend) in hitranges:
            if startpos is None:
                startpos, endpos = currstart, currend
            else:
                if (currstart > endpos + 1):
                    mergedhitranges.append((startpos, endpos))
                    startpos, endpos = currstart, currend
                else:
                    assert currstart >= startpos
                    endpos = max(endpos, currend)
        if startpos is not None:
            mergedhitranges.append((startpos, endpos))

        self.deleted_regions = []
        del_start, del_end = self.ref_start, self.ref_end
        for (startpos, endpos) in mergedhitranges:
            if startpos > del_start:
                self.deleted_regions.append((del_start, startpos - 1))
            del_start = endpos + 1
        if del_start <= del_end:
            self.deleted_regions.append((del_start, del_end))
        self.deletion_string = self.DELIM.join(['{0}-{1}'.format(start, end) for (start, end) in self.deleted_regions]) if len(self.deleted_regions) > 0 else self.NONESTRING

        logging.info('Deletions identified: {0}'.format(self.deletion_string))

        if self.show_region is not None:
            if len([(start,end) for (start,end) in mergedhitranges if start <= self.show_region_start and end >= self.show_region_end]) > 0:
                self.region_present = 1
                logging.info('Status of region {0}: present'.format(self.show_region))
            else:
                self.region_present = 0
                logging.info('Status of region {0}: absent'.format(self.show_region))


    def get_struct_profile(self):
        self.struct_profile = self.UNKSTRING
        if self.struct_profiles:
            logging.debug('Structural variant profile file provided. Searching for a matching profile...')
            self.check_file(self.struct_profiles)
            try:
                struct_matcher = StructProfileMatcher(self.struct_profiles, delim=self.DELIM, nonestring=self.NONESTRING, start_pos=self.ref_start, end_pos=self.ref_end)
            except IOError:
                self.exitonerror('IOError reading structural profile file {0}'.format(self.struct_profiles))
            except ProfileError as e:
                self.exitonerror('Error reading structural profile file: {0}'.format(e.value))
            self.struct_profile = struct_matcher.get_profile(self.deletion_string, default=self.UNKSTRING)
            logging.info('Structural variant profile identified: {0}'.format(self.struct_profile))


    def call_struct(self):
        self.do_assembly()
        self.do_blast()
        self.parse_blast()
        self.get_struct_profile()


    def generate_vcf(self):
        self.snpfile = self.outprefix + '.vcf'
        samtools_mpileup_args = ['samtools', 'mpileup', '-uAIf', self.ref, self.bam]
        bcftools_call_args = ['bcftools', 'call', '-mv', '-o', self.snpfile]
        self.run_external_call([samtools_mpileup_args, bcftools_call_args], 'SNP calling', self.snpfile)
        self.check_file(self.snpfile)


    def get_ambiguous_code(self, bases):
        sorted_bases = sorted(map(str,bases))
        if sorted_bases == ['A','C']:
            return 'M'
        elif sorted_bases == ['A','G']:
            return 'R'
        elif sorted_bases == ['A','T']:
            return 'W'
        elif sorted_bases == ['C','G']:
            return 'S'
        elif sorted_bases == ['C','T']:
            return 'Y'
        elif sorted_bases == ['G','T']:
            return 'K'
        elif sorted_bases == ['A','C','G']:
            return 'V'
        elif sorted_bases == ['A','C','T']:
            return 'H'
        elif sorted_bases == ['A','G','T']:
            return 'D'
        elif sorted_bases == ['C','G','T']:
            return 'B'
        else:
            return 'N'


    def list2delimsep(self, inlist):
        if len(inlist) == 0:
            return self.NONESTRING
        return self.DELIM.join(map(str, inlist))


    def parse_vcf(self):
        logging.debug('Parsing SNP calling output...')
        varsites = [] 
        Nsites = []
        Nsite_counts = []
        del_varsites = []
        del_Nsites = []
        with open(self.snpfile) as snpfilehandle:
            vcf_reader = vcf.Reader(snpfilehandle)
            for record in vcf_reader:
                # use this record if it's not in a deleted region and there's at least one forward and one reverse read supporting the alternative call. If it's in a deleted region, just report in log file
                if record.INFO['DP4'][2] > 0 and record.INFO['DP4'][3] > 0:
                    if len([(start, end) for (start, end) in self.deleted_regions if record.POS >= start and record.POS <= end]) == 0:
                        if record.samples[0]['GT'] == '1/1':
                            varsites.append('{0}{1}{2}'.format(record.REF, record.POS, record.ALT[0]))
                        else:
                            Nsites.append('{0}{1}{2}'.format(record.REF, record.POS, self.get_ambiguous_code([record.REF, record.ALT[0]])))
                            refcounts = record.INFO['DP4'][0] + record.INFO['DP4'][1]
                            altcounts = record.INFO['DP4'][2] + record.INFO['DP4'][3]
                            Nsite_counts.append('{0}{1},{2}{3}'.format(record.REF, refcounts, record.ALT[0], altcounts))
                    else:
                        if record.samples[0]['GT'] == '1/1':
                            del_varsites.append('{0}{1}{2}'.format(record.REF, record.POS, record.ALT[0]))
                        else:
                            del_Nsites.append('{0}{1}{2}'.format(record.REF, record.POS, self.get_ambiguous_code([record.REF, record.ALT[0]])))

        self.varsite_string = self.DELIM.join(varsites) if len(varsites) > 0 else self.NONESTRING
        self.Nsite_string = self.DELIM.join(Nsites) if len(Nsites) > 0 else self.NONESTRING
        self.Nsite_count_string = self.list2delimsep(Nsite_counts)
        logging.info('Homozygous SNPs identified: {0}'.format(self.varsite_string))
        logging.info('Heterozygous SNPs identified: {0}'.format(self.Nsite_string))
        if len(del_varsites) > 0:
            logging.info('Homozygous SNPs identified in deleted regions (not reported in summary file): {0}'.format(self.DELIM.join(del_varsites)))
        if len(del_Nsites) > 0:
            logging.info('Heterozygous SNPs identified in deleted regions (not reported in summary file: {0}'.format(self.DELIM.join(del_Nsites)))


    def get_snp_profile(self):
        self.snp_profile = self.UNKSTRING
        if self.snp_profiles:
            logging.debug('SNP profile file provided. Searching for a matching profile...')
            self.check_file(self.snp_profiles)
            try:
                snp_matcher = SNPProfileMatcher(self.snp_profiles, delim=self.DELIM, nonestring=self.NONESTRING, start_pos=self.ref_start, end_pos=self.ref_end)
            except IOError:
                self.exitonerror('IOError reading SNP profile file {0}.'.format(self.snp_profiles))
            except ProfileError as e:
                self.exitonerror('Error reading SNP profile file: {0}'.format(e.value))
            self.snp_profile = snp_matcher.get_profile('\t'.join([self.varsite_string, self.Nsite_string]), default=self.UNKSTRING)
            logging.info('SNP profile identified: {0}'.format(self.snp_profile)) 


    def call_snps(self):
        self.generate_vcf()
        self.parse_vcf()
        self.get_snp_profile()


    def extract_flanks(self):
        samtools_index_args = ['samtools', 'index', self.bam]
        self.run_external_call([samtools_index_args], 'bam file indexing')

        logging.debug('Identifying flanking sequences...')

        samfile = pysam.AlignmentFile(self.bam, 'rb')

        lflanks_highqual = []  # list of tuples with left flank seqs (i.e. flanking the start position) and strand, for reads that pass quality threshold.
        rflanks_highqual = []  # As above, for right flanking seqs (i.e. flanking the end position)

        # Fetch all reads that overlap the start position
        for a in samfile.fetch(self.ref_contig, self.ref_start-1, self.ref_start):  # 0-based coordinates
            ref_posns = a.get_reference_positions(full_length=True)  # List of same length as read length, containing reference positions that each base is aligned to (None if a base is unaligned).
            # First test ensures that the start position is included in the alignment (i.e. if there's a deletion in the read at exactly this position, ignore as it's not possible to determine where the flanking sequence should start). Second test checks that a minimum length within the region of interest has been mapped to (avoiding spurious mapping). Note 0-based coordinates.
            if self.ref_start-1 in ref_posns and len([pos for pos in ref_posns if pos is not None and pos >= self.ref_start-1 + (self.min_mapped_len - 1)]) > 0:
                q_start = ref_posns.index(self.ref_start-1)  # Position in the read corresponding to the start position
                lflank = a.query_sequence[q_start-self.flank_len:q_start]  # Flanking sequence for this read
                lqual = a.query_qualities[q_start-self.flank_len:q_start]  # Corresponding quality values
                # Use read if there is sufficient flanking sequence to extract the full region, and minimum quality value within flanking sequence is above threshold
                if len(lflank) == self.flank_len and min(lqual) >= self.min_qual:
                    strand = -1 if a.is_reverse else 1
                    lflanks_highqual.append((lflank, strand))

        # Fetch all reads with mapping that overlaps the end position
        for a in samfile.fetch(self.ref_contig, self.ref_end-1, self.ref_end):
            ref_posns = a.get_reference_positions(full_length=True)
            if self.ref_end-1 in ref_posns and len([pos for pos in ref_posns if pos is not None and pos <= self.ref_end-1 - (self.min_mapped_len - 1)]) > 0:
                q_end = ref_posns.index(self.ref_end-1)
                rflank = a.query_sequence[q_end+1:q_end+self.flank_len+1]
                rqual = a.query_qualities[q_end+1:q_end+self.flank_len+1]
                if len(rflank) == self.flank_len and min(rqual) >= self.min_qual:
                    strand = -1 if a.is_reverse else 1
                    rflanks_highqual.append((rflank, strand))

        # Count the number of reads containing each unique flanking sequence
        lcounts_highqual = Counter([flank for (flank,strand) in lflanks_highqual])
        rcounts_highqual = Counter([flank for (flank,strand) in rflanks_highqual])

        # Deduplicated, sorted list containing each flanking sequence that passes quality threshold + minimum read threshold + strand threshold
        lseqlist_highqual = sorted([seq for seq in lcounts_highqual if lcounts_highqual[seq] >= self.min_reads and lflanks_highqual.count((seq, 1)) >= self.min_each_strand and lflanks_highqual.count((seq, -1)) >= self.min_each_strand])
        rseqlist_highqual = sorted([seq for seq in rcounts_highqual if rcounts_highqual[seq] >= self.min_reads and rflanks_highqual.count((seq, 1)) >= self.min_each_strand and rflanks_highqual.count((seq, -1)) >= self.min_each_strand])

        self.lflanks, self.rflanks, self.lcounts, self.rcounts = map(self.list2delimsep, [lseqlist_highqual, rseqlist_highqual, [lcounts_highqual[lseq] for lseq in lseqlist_highqual], [rcounts_highqual[rseq] for rseq in rseqlist_highqual]])
    
        logging.info('Flanking sequences identified (left flank, right flank, left counts, right counts): {0}, {1}, {2}, {3}'.format(self.lflanks, self.rflanks, self.lcounts, self.rcounts))


    def write_output(self):
        header = ['Deletions', 'Structural_variant', 'SNPs_homozygous', 'SNPs_heterozygous', 'Heterozygous_SNP_counts', 'SNP_variant', 'Combined_variant', 'Left_flanks', 'Right_flanks', 'Left_flank_counts', 'Right_flank_counts']
        contents = [self.deletion_string, self.struct_profile, self.varsite_string, self.Nsite_string, self.Nsite_count_string, self.snp_profile, self.struct_profile + '-' + self.snp_profile, self.lflanks, self.rflanks, self.lcounts, self.rcounts]
        if self.show_region:
            header.append('{0}_presence'.format(self.show_region))
            contents.append(str(self.region_present))
        outfile = self.outprefix + '_summary.txt'
        try:
            with open(outfile, 'w') as outhandle:
                outhandle.write('\t'.join(header) + '\n')
                outhandle.write('\t'.join(map(str, contents)) + '\n')
        except IOError as e:
            self.exitonerror('Error writing file {0}'.format(outfile))
        logging.info('Final output written to: {0}'.format(outfile))


    def run_typing(self):
        self.do_map()
        self.call_struct()
        self.call_snps()
        self.extract_flanks()
        self.write_output()
        self.cleanup()



def get_argsparser():
    parser = argparse.ArgumentParser(description = 'TETyper version {0}. Given a set of input reads and a reference, TETyper performs typing to identify: 1. deletions and SNP variation relative to the reference, and 2. the immediate (up to ~20bp) sequence(s) flanking the reference.'.format(VERSION))
    parser.add_argument('--outprefix', help = 'Prefix to use for output files. Required.', required = True)
    parser.add_argument('--ref', help = 'Reference sequence in fasta format. If not already indexed with bwa, this will be created automatically. A blast database is also required, again this will be created automatically if it does not already exist. Required.', required = True)
    parser.add_argument('--refdb', help = 'Blast database corresponding to reference file (this argument is only needed if the blast database was created with a different name).')
    parser.add_argument('--fq1', help = 'Forward reads. Can be gzipped.')
    parser.add_argument('--fq2', help = 'Reverse reads. Can be gzipped.')
    parser.add_argument('--fq12', help = 'Interleaved forward and reverse reads.')
    parser.add_argument('--bam', help = 'Bam file containing reads mapped to the given reference. If the reads have already been mapped, this option saves time compared to specifying the reads in fastq format. If this option is specified then --fq* are ignored.')
    parser.add_argument('--assembly', help = 'Use this assembly (fasta format) for detecting structural variants instead of generating a new one. This option saves time if an assembly is already available.')
    parser.add_argument('--spades_params', help = 'Additional parameters for running spades assembly. Enclose in quotes and precede with a space. Default: " --cov-cutoff auto --disable-rr". Ignored if --assembly is specified.', default=' --cov-cutoff auto --disable-rr')
    parser.add_argument('--struct_profiles', help = 'File containing known structural variants. Tab separated format with two columns. First column is variant name. Second column contains a list of sequence ranges representing deletions relative to the reference, or "none" for no deletions. Each range should be written as "startpos-endpos", with multiple ranges ordered by start position and separated by a "|" with no extra whitespace.')
    parser.add_argument('--snp_profiles', help = 'File containing known SNP profiles. Tab separated format with three columns. First column: variant name, second column: homozygous SNPs, third column: heterozygous SNPs. SNPs should be written as "refPOSalt", where "POS" is the position of that SNP within the reference, "ref" is the reference base at that position (A, C, G or T), and "alt" is the variant base at that position (A, C, G or T for a homozygous SNP; M, R, W, S, Y or K for a heterozygous SNP). Multiple SNPs of the same type (homozygous or heterozygous) should be ordered by position and separated by a "|". "none" indicates no SNPs of the given type.')
    parser.add_argument('--flank_len', help = 'Length of flanking region to extract. Required.', type = int, required = True)
    parser.add_argument('--min_reads', help = 'Minimum read number for including a specific flanking sequence. Default 10.', type = int, default = 10)
    parser.add_argument('--min_each_strand', help = 'Minimum read number for each strand for including a specific flanking sequence. Default 1.', type = int, default = 1)
    parser.add_argument('--min_mapped_len', help = 'Minimum length of mapping for a read to be used in determining flanking sequences. Higher values are more robust to spurious mapping. Lower values will recover more reads. Default 30.', type = int, default = 30)
    parser.add_argument('--min_qual', help = 'Minimum quality value across extracted flanking sequence. Default 10.', type = int, default = 10)
    parser.add_argument('--show_region', help = 'Display presence/absence for a specific region of interest within the reference (e.g. to display blaKPC presence/absence with the Tn4401b-1 reference, use "7202-8083")')
    parser.add_argument('--threads', help = 'Number of threads to use for mapping and assembly steps. Default: 1', type = int, default = 1)
    parser.add_argument('-v', '--verbosity', help = 'Verbosity level for logging to stderr. 1 = ERROR, 2 = WARNING, 3 = INFO, 4 = DUBUG. Default: 3.', type = int, choices = [1,2,3,4], default = 3)
    parser.add_argument('--no_overwrite', help = 'Flag to prevent accidental overwriting of previous output files. In this mode, the pipeline checks for a log file named according to the given output prefix. If it exists then the pipeline exits without modifying any files.', action = 'store_true')
    return(parser)



if __name__ == '__main__':
    parser = get_argsparser()
    args = parser.parse_args()
    typer = TETyper(args)
    typer.run_typing()
