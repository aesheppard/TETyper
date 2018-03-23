# TETyper

TETyper is a command line tool designed for typing a specific transposable element (TE) of interest from paired-end sequencing data. It determines single nucleotide variants (SNVs) and deletions within the TE, as well as flanking sequences surrounding the TE.


## Requirements

- python 2.7 (with Biopython, pysam, vcf)
- samtools (tested on version 1.4.1)
- bwa (tested on version 0.7.12)
- spades (tested on version 3.10.1)
- BLAST+ (tested on version 2.2.25)


## Running TETyper

### Basic usage

```
TETyper.py --ref REFERENCE.fasta --fq1 FORWARD_READS.fq.gz --fq2 REVERSE_READS.fq.gz --outprefix OUTPREFIX --flank_len FLANK_LENGTH
```

REFERENCE.fasta should be indexed with bwa, and a blast database generated prior to running TEtyper. This can be achieved with the following commands:
```
bwa index REFERENCE.fasta
makeblastdb -in REFERENCE.fasta -dbtype nucl
```

### Output

TETyper produces the following output files:
- **OUTPREFIX_summary.txt**: Summary file containing all typing outputs, described in detail below
- **OUTPREFIX.log**: Log file containing details of individual steps performed, as well as any errors
- **OUTPREFIX.bam**: Reads mapped to REFERENCE.bam
- **OUTPREFIX_mappedreads_1.fq, OUTPREFIX_mappedreads_2.fq**: Mapped reads converted to fastq format
- **OUTPREFIX_spades**: Directory containing spades assembly from mapped reads
- **OUTPREFIX_blast.txt**: Tabular blastn results comparing the spades assembly to REFERENCE.fasta
- **OUTPREFIX.vcf**: SNVs identified


The summary file contains 11-12 columns. These are:
- **Deletions**: A list of sequence ranges corresponding to regions of the reference classified as deletions for this sample, or "none" for no deletions.
- **Structural_variant**: If --struct_profiles is specified and the pattern of deletions above corresponds to one of these profiles, then the profile name is given, otherwise "unknown".
- **SNPs_homozygous**: A list of homozygous SNPs identified, or "none".
- **SNPs_heterozygous**: A list of heterozygous SNPs identified, or "none".
- **Heterozygous_SNP_counts**: For each heterozygous SNP, the number of reads supporting the reference and alternative calls, or "none" if there are no heterozygous SNPs.
- **SNP_variant**: If --snp_profiles is specified and the pattern of homozygous and heterozygous SNPs corresponds to one of these profiles, then the profile name is given. Otherwise "unknown".
- **Combined_variant**: Single name combining Structural_variant and SNP_variant, separate by "-".
- **Left_flanks**: A list of distinct sequences passing quality filters that flank the start position of the reference. 
- **Right_flanks**: A list of distinct sequences passing quality filters that flank the end position of the reference.
- **Left_flank_counts**: The number of high quality reads supporting each of the left flanking sequences.
- **Right_flank_counts**: The number of high quality reads supporting each of the right flanking sequences.
- **X_Y_presence**: If --show_region is specified as --show_region X-Y, this column shows 1 if the entirety of that region is classified as present (i.e. no overlap with deleted regions), or 0 otherwise. If --show_region is unspecified, this column is omitted.


### Advanced usage (all options):
-  -h, --help:            show this help message and exit
-  --outprefix OUTPREFIX:
                        Prefix to use for output files. Required.
-  --ref REF:             Reference sequence in fasta format. Must be indexed
                        with bwa index. Required.
-  --refdb REFDB:         Blast database corresponding to reference file. If
                        unspecified, blast database is assumed to have the
                        same name as the reference file.
-  --fq1 FQ1:             Forward reads. Can be gzipped.
-  --fq2 FQ2:             Reverse reads. Can be gzipped.
-  --fq12 FQ12:           Interleaved forward and reverse reads.
-  --bam BAM:             Bam file containing reads mapped to the given
                        reference. If the reads have already been mapped, this
                        option saves time compared to specifying the reads in
                        fastq format. If this option is specified then --fq*
                        are ignored.
-  --assembly ASSEMBLY:   Use this assembly (fasta format) for detecting
                        structural variants instead of generating a new one.
                        This option saves time if an assembly is already
                        available.
-  --struct_profiles STRUCT_PROFILES:
                        File containing known structural variants. Tab
                        separated format with two columns. First column is
                        variant name. Second column contains a list of
                        sequence ranges representing deletions relative to the
                        reference, or "none" for no deletions. Each range
                        should be written as "startpos-endpos", with multiple
                        ranges ordered by start position and separated by a
                        "|" with no extra whitespace.
-  --snp_profiles SNP_PROFILES:
                        File containing known SNP profiles. Tab separated
                        format with three columns. First column: variant name,
                        second column: homozygous SNPs, third column:
                        heterozygous SNPs. SNPs should be written as
                        "refPOSalt", where "POS" is the position of that SNP
                        within the reference, "ref" is the reference base at
                        that position (A, C, G or T), and "alt" is the variant
                        base at that position (A, C, G or T for a homozygous
                        SNP; M, R, W, S, Y or K for a heterozygous SNP).
                        Multiple SNPs of the same type (homozygous or
                        heterozygous) should be ordered by position and
                        separated by a "|". "none" indicates no SNPs of the
                        given type.
-  --flank_len FLANK_LEN:
                        Length of flanking region to extract. Required.
-  --min_reads MIN_READS:
                        Minimum read number for including a specific flanking
                        sequence. Default 10.
-  --min_each_strand MIN_EACH_STRAND:
                        Minimum read number for each strand for including a
                        specific flanking sequence. Default 1.
-  --min_mapped_len MIN_MAPPED_LEN:
                        Minimum length of mapping for a read to be used in
                        determining flanking sequences. Higher values are more
                        robust to spurious mapping. Lower values will recover
                        more reads. Default 30.
-  --min_qual MIN_QUAL:   Minimum quality value across extracted flanking
                        sequence. Default 10.
-  --show_region SHOW_REGION:
                        Display presence/absence for a specific region of
                        interest within the reference (e.g. to display blaKPC
                        presence/absence with the Tn4401b-1 reference, use
                        "7202-8083")
-  --threads THREADS:     Number of threads to use for mapping and assembly
                        steps. Default: 1
-  -v {1,2,3,4}, --verbosity {1,2,3,4}:
                        Verbosity level for logging to stderr. 1 = ERROR, 2 =
                        WARNING, 3 = INFO, 4 = DUBUG. Default: 3.



### Suggested usage for KPC isolates

TETyper was designed with the blaKPC transposon Tn4401 in mind. A Tn4401b reference sequence is provided with TETyper (Tn4401b-1.fasta), as well as example profile definitions for SNVs / deletions with respect to this reference (struct_profiles.txt and snp_profiles.txt). To use these additional resources, TEtyper can be run as follows:
```
TETyper.py --ref Tn4401b-1.fasta --fq1 FORWARD_READS.fq.gz --fq2 REVERSE_READS.fq.gz --outprefix OUTPREFIX --flank_len FLANK_LENGTH --struct_profiles struct_profiles.txt --snp_profiles snp_profiles.txt --show_region 7202-8083
```


### Re-running samples

TETyper provides options for specifying the mapped bam file and/or assembly file in order to save processing time if the same samples are rerun with different parameters.

For example, newly discovered profiles can be manually appended to the profile files. TETyper can then be rerun with the modified profile files, without redoing all the processing steps, by specifiying the mapped bam file and spades assembly as parameters instead of the original reads. E.g.:
```
TETyper.py --ref Tn4401b-1.fasta --outprefix RERUN --bam OUTPREFIX.bam --assembly OUTPREFIX_spades/contigs.fasta --flank_len FLANK_LENGTH --struct_profiles STRUCT_PROFILES_MODIFIED.txt --snp_profiles SNP_PROFILES_MODIFIED.txt --show_region 7202-8083
```

The --bam and --assembly options can also be useful for re-running samples with different parameters for flanking sequence extraction (e.g. a different flank length).

