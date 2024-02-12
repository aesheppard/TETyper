# TETyper

TETyper is a command line tool designed for typing a specific transposable element (TE) of interest from paired-end sequencing data. It determines single nucleotide variants (SNVs) and deletions within the TE, as well as flanking sequences surrounding the TE.

TETyper can be cited as follows:

[Sheppard et al bioRxiv 288001; doi: https://doi.org/10.1101/288001](https://www.biorxiv.org/content/early/2018/03/23/288001)


## Installation

Requirements:

- python 3 (<3.11) (with Biopython, pysam, pyvcf3)
- [samtools, bcftools](http://www.htslib.org/) (tested on version 1.18)
- [bwa](http://bio-bwa.sourceforge.net/) (tested on version 0.7.17)
- [spades](http://cab.spbu.ru/software/spades/) (tested on version 3.11.1)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (tested on version 2.15.0)
- [snakemake](https://snakemake.github.io/) (tested on version 7.32.4)

*Still to be implemented*:
TETyper is conda-installable, enabling automatic installation of python dependencies:

```
conda install TETyper
```

To check whether TETyper is running correctly, a small test dataset has been provided. If the following command is executed:

```
TETyper.py --ref Tn4401b-1.fasta --fq1 SRR1582895_sm_1.fq.gz --fq2 SRR1582895_sm_2.fq.gz --outprefix test --flank_len 5 --struct_profiles struct_profiles.txt --snp_profiles snp_profiles.txt --show_region 7202-8083
```

then the file test_summary.txt should look something like this:

Deletions | Structural_variant | SNPs_homozygous | SNPs_heterozygous | Heterozygous_SNP_counts | SNP_variant | Combined_variant | Left_flanks | Right_flanks | Left_flank_counts | Right_flank_counts | 7202-8083_presence
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
none | Tn4401b | C8015T | none | none | 2 | Tn4401b-2 | AGATA\|GTTCT | AGATA\|GTTCT | 56\|151 | 52\|101 | 1

Note that the exact values in the Left_flank_counts and Right_flank_counts columns may vary slightly depending on the bwa version used.


## Running TETyper

### Basic usage

When running one sample, you can run TETyper directly using this command:

'''
TETyper.py --ref REFERENCE.fasta --fq1 FORWARD_READS.fq.gz --fq2 REVERSE_READS.fq.gz --outprefix OUTPREFIX --flank_len FLANK_LENGTH
'''
If you want to run multiple samples, TETyper can be run via snakemake. For example, the test dataset can be run with:
''' 
snakemake --cores 1
'''
Parameters to run the test dataset are already provided in the config.json file.
Parameters can be modified by directly editing the **config.json** file. If a parameter is not needed, write **null**.
You can also provide a tab separated file called **sample_data.txt** (e.g. exported from Excel) in the following format for some of the parameters:

| outprefix | fq1                      | fq2                      | bam            |
| --------- | ------------------------ | ------------------------ | -------------- |
| sample_1  | sample_1_fwd_reads.fq.gz | sample_1_rev_reads.fq.gz |                |
| sample_2  | sample_2_fwd_reads.fq.gz | sample_2_rev_reads.fq.gz | sample_2.bam   |
| sample_3  |                          |                          | sample_3.bam   |

Please note that if the sample_data.txt file is present, it will be used. Rename or delete it if you would like to use arguments from config.json.
If you are new to TETyper, ignore the bam column (but still include it in the header line). See Advanced Usage below for an explanation of how to use it. 
If both bam and fq files are provided, the bam file will be used.

```
snakemake --cores *n*
```
n is the number of cores you want to use. For example, *snakemake --cores 4* uses 4 cores.

### Output

TETyper produces the following output files:
- **OUTPREFIX_summary.txt**: Summary file containing all typing outputs, described in detail below
- **OUTPREFIX.log**: Log file containing details of individual steps performed, as well as any errors
- **OUTPREFIX.bam**: Reads mapped to REFERENCE.bam
- **OUTPREFIX_mappedreads_1.fq, OUTPREFIX_mappedreads_2.fq**: Mapped reads converted to fastq format
- **OUTPREFIX_spades**: Directory containing spades assembly from mapped reads. By default, only contigs.fasta and spades.log are kept. See *--keep-spades* option below to keep all files produced by spades.
- **OUTPREFIX_blast.txt**: Tabular blastn results comparing the spades assembly to REFERENCE.fasta
- **OUTPREFIX.vcf**: SNVs identified
- **all_summary.txt**: Summary of all results generated in various **OUTPREFIX_summary.txt** files.


The summary file contains up columns when running mode "all". These are:
- **Deletions**: A list of sequence ranges corresponding to regions of the reference classified as deletions for this sample, or "none" for no deletions.
- **Structural_variant**: If --struct_profiles is specified and the pattern of deletions above corresponds to one of these profiles, then the profile name is given, otherwise "unknown".
- **SNPs_homozygous**: A list of homozygous SNPs identified, or "none".
- **SNPs_heterozygous**: A list of heterozygous SNPs identified, or "none".
- **Heterozygous_SNP_counts**: For each heterozygous SNP, the number of reads supporting the reference and alternative calls, or "none" if there are no heterozygous SNPs.
- **SNP_variant**: If --snp_profiles is specified and the pattern of homozygous and heterozygous SNPs corresponds to one of these profiles, then the profile name is given. Otherwise "unknown".
- **Combined_variant**: Single name combining Structural_variant and SNP_variant, separated by "-".
- **Left_flanks**: A list of distinct sequences passing quality filters that flank the start position of the reference. 
- **Right_flanks**: A list of distinct sequences passing quality filters that flank the end position of the reference.
- **Left_flank_counts**: The number of high quality reads supporting each of the left flanking sequences.
- **Right_flank_counts**: The number of high quality reads supporting each of the right flanking sequences.
- **X_Y_presence**: If --show_region is specified as --show_region X-Y, this column shows 1 if the entirety of that region is classified as present (i.e. no overlap with deleted regions), or 0 otherwise. If --show_region is unspecified, this column is omitted.

When running mode "variants", only columns 1-7 are produced. When running mode "flanks", only columns 8-11 are produced. Column 12 is present in mode "all" and "variants" if the option *--show-region* (see Advanced usage for more information) is enabled.

### Advanced usage (all options):
-  -h, --help:            show this help message and exit
-  --outprefix OUTPREFIX:
                        Prefix to use for output files. Required.
-  --ref REF:             Reference sequence in fasta format. If not already
                        indexed with bwa, this will be created automatically.
                        A blast database is also required, again this will be
                        created automatically if it does not already exist.
                        Required.
-  --refdb REFDB:         Blast database corresponding to reference file (this
                        argument is only needed if the blast database was
                        created with a different name).
-  --fq1 FQ1:             Forward reads. Can be gzipped.
-  --fq2 FQ2:             Reverse reads. Can be gzipped.
-  --bam BAM:             Bam file containing reads mapped to the given
                        reference. If the reads have already been mapped, this
                        option saves time compared to specifying the reads in
                        fastq format. If this option is specified then --fq*
                        are ignored.
-  --assembly ASSEMBLY:   Use this assembly (fasta format) for detecting
                        structural variants instead of generating a new one.
                        This option saves time if an assembly is already
                        available.
-  --spades_params SPADES_PARAMS:
                        Additional parameters for running spades assembly.
                        Enclose in quotes and precede with a space. Default: "
                        --cov-cutoff auto --disable-rr". Ignored if
                        --assembly is specified.
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
-  -v, --verbosity:
                        Verbosity level for logging to stderr. 1 = ERROR, 2 =
                        WARNING, 3 = INFO, 4 = DEBUG. Default: 3.
-  --no_overwrite:        Flag to prevent accidental overwriting of previous
                        output files. In this mode, the pipeline checks for a
                        log file named according to the given output prefix.
                        If it exists then the pipeline exits without modifying
                        any files.
-   --keep-spades:      Keeps spades folder intact. Default: off (spades.log 
                        and params.txt are always preserved.)
-   --mode:             Run only variants (deletions, homo/heterozygous SNVs), only flank extraction, or whole
                        pipeline. Options: "all", "flanks", "variants". Default: "all".                


### A note on the profile files

The arguments --struct_profiles and --snp_profiles enable the user to specify a naming system for identified variants. This simply provides a convenient way of generating output that is more human-readable than a list of deletions or SNVs. This is probably not helpful when running TETyper for the first time with a new reference, but may be useful once relevant variants have been identified.


### Suggested usage for KPC isolates
TETyper was designed with the blaKPC transposon Tn4401 in mind. A Tn4401b reference sequence is provided with TETyper (Tn4401b-1.fasta), as well as example profile definitions for SNVs / deletions with respect to this reference (struct_profiles.txt and snp_profiles.txt). To use these, TETyper can be run as follows:
```
TETyper.py --ref Tn4401b-1.fasta --fq1 FORWARD_READS.fq.gz --fq2 REVERSE_READS.fq.gz --outprefix OUTPREFIX --flank_len FLANK_LENGTH --struct_profiles struct_profiles.txt --snp_profiles snp_profiles.txt --show_region 7202-8083
```


### Re-running samples

TETyper provides options for specifying the mapped bam file and/or assembly file in order to save processing time if the same samples are rerun with different parameters. For example, newly discovered profiles can be manually appended to the profile files. TETyper can then be rerun with the modified profile files, without redoing all the processing steps, by specifiying the mapped bam file and spades assembly as parameters instead of the original reads. E.g.:


The --bam and --assembly options can also be useful for re-running samples with different parameters for flanking sequence extraction (e.g. a different flank length).


### Changes in version 1.2
- Removed support for reverse/forward interleaved fq files.
- Updated dependencies to the latest versions (where possible).
- Added options for running variant calling (deletions, homo/heterozygous SNVs) and flank extraction separately.
- Implemented cleanup of files in the spades folder (with the option to retain them).
- Introduced ability to run multiple samples concurrently using Snakemake.
- Included a summary file encompassing all samples for easier processing.

## A note on spades
Before you can run TETyper, you must execute the following command one time:
```
source spades_corrector.sh
```
 However, you can also do the following things manually if you encounter any errors:
- Find the folder where your base environment is located by typing
```conda env list ``` (the path next to the environment)
- Open *<PATH>/share/spades-3.11.1*/share/spades/pyyaml3*
- In line 126, change *if not isinstance(key, collections.Hashable)* to *if not isinstance(key, collections.abc.Hashable):* and save your changes.
See https://github.com/ablab/spades/issues/873 for more information.
