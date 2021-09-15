# Osteo-WGS-SNV-Workflow
Workflow for identifying single-nucleotide-variants (SNVs) in osteosarcoma Whole Genome Sequencing (WGS) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the Burrows-Wheeler Aligner.  The version used was bwa-0.7.17
- GATK version 3.6 (Picard included)
- All Python scripts were run on Python version 2.7.13 unless otherwise noted.  

# Workflow
**Notes:**
- Each sample has it's own directory for output files.  Each individual directory is labeled by the "Sample ID"

**Step 1) Downloading the hg38 genome files:** The following files were downloaded from the GATK resource bundle on 9/18/2018 and saved to a directory named "hg38_osteo".  The .fasta.gz file was then unzipped.  
- Homo_sapiens_assembly38.fasta.gz
- Homo_sapiens_assembly38.dict
- Homo_sapiens_assembly38.fasta.fai

```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
```
**Step 2) Create a BWA genome index files:** The BWA "index" command is used to create the genome index files from the previously downloaded hg38 .fasta genome sequence.  These files will then be used in the upcoming alignment step.  BWA places the output files into the same foloder as the genome .fasta file was called from.  The BWA index files were then moved to a directory names "BWA_indexes" for later use.

The -a bwtsw flag will specify that we want to use the indexing algorithm that is capable of handling the whole human genome.

```
fasta_dir=<path to hg38 .fasta directory>

bwa-0.7.17/bwa index -a bwtsw $fasta_dir/Homo_sapiens_assembly38.fasta
```
**Step 3) Revert .bam files:** The downloaded osteo sequence files have already been processed and aligned.  In order to use your own analysis workflow, the files first must be reverted back to their pre-processed and unmapped form.  To accomplish this, we first use the program, Picard, with the "RevertSam" function.  This takes an aligned .bam file and removes alignment information in order to generate an unmapped BAM (uBAM).  Details on this part of the workflow can be found in the GATK Tutorial #6484: (How to) Generate an unmapped BAM from FASTQ or aligned BAM. https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top.  It removes alignment information and any recalibrated base quality information.  It makes it possible to re-analyze the file using your own pipeline.  Separate revert steps were carried out for the tumor and matched normal samples.  

Standard tags cleared by default are NM, UQ, PG, MD, MQ, SA, MC and AS.  Additionally, the OQ tag is removed by the default "RESTORE_ORIGINAL_QUALITIES" parameter.  Any nonstandard tags should be removed.  To list all tags within a BAM, use the command: 
`samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}'` You should leave the RG tag.

The "SANITIZE" option removes reads that cause problems for certain tools such as MarkIlluminaAdapters. It removeds paired reads with missing mates, duplicate records and records with mismatches in length of bases and qualities.

For paired read files: because each read in a pair has the same query name, they are listed consecutively.  We make sure to alter the previous sort order.  Coordinate sorted reads result in the aligner incorrectly estimating insert size from blocks of paired reads as they are not randomly distributed.

I use the "OUTPUT_BY_READGROUP=true" option in order to create a separate file for each readgroup.  This is neccessary because there are two different MarkDuplicate steps (one before merging each readgroup and one after merging).  This ensures that both optical and PCR duplicates are identified (see GATK tutorial 6483)

The "MAX_DISCARD_FRACTION" option is informational only.  It does nto affect processing.

The "SORT_ORDER=queryname", "RESTORE_ORIGINAL_QUALITIES=true", REMOVE_DUPLICATE_INFORMATION=true" and "REMOVE_ALIGNMENT_INFORMATION=true" options are all default but I kept them in the code anyway.

```
ALIGNMENT_RUN=<Sample ID>
INPUT=<path to input directory>"/"$ALIGNMENT_RUN".grch38_NoAltAnalysisSet_bwa.bam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar RevertSam \
    I=$INPUT \
    O=$OUTPUT_DIR \
    OUTPUT_BY_READGROUP=true \
    SANITIZE=true \
    MAX_DISCARD_FRACTION=0.005 \
    ATTRIBUTE_TO_CLEAR=XS \
    ATTRIBUTE_TO_CLEAR=XA \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=true \
    REMOVE_DUPLICATE_INFORMATION=true \
    REMOVE_ALIGNMENT_INFORMATION=true \
    TMP_DIR=<path to temp directory>/working_temp_rsn
```
**Step 4) Mark adapters:** Using Picard and the "MarkIlluminaAdapters" function, this step takes an uBAM file from the previous step and rewrites it with new adapter-trimming tags.  Per tutorial 6483 on the GATK website: https://software.broadinstitute.org/gatk/documentation/article?id=6483 This is needed so that the sequenced adapters contribute minimally to the alignments.  The tool adds an "XT" tag to a read record to mark the 5' start position of the specified adapter sequence.  It also produces a metrics file.  Note: each .bam from Step 3 may produce multiple uBAM files depending on how many read groups are present.  Each individual uBAM must undergo "MarkIlluminaAdapters".  

```
ALIGNMENT_RUN=<Sample ID>
INPUT=<path to uBAM file>
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar MarkIlluminaAdapters \
    I=$INPUT \
    O=$OUTPUT_DIR/uBAM_markedAdapters.bam \
    M=$OUTPUT_DIR/adapter_metrics.txt \
    TMP_DIR=<path to temp directory>/working_temp
```
**Step 5) Convert uBAM to FASTQ:** Now we use Picard and the "SamToFastq" function to take the uBAM files, which have been modified so that Illumina adapter sequences are marked with the "XT" tag, and converts them to .fastq files for further processing.  It produces a .fastq file in which all extant meta data (read group info, alignment info, flags and tags) are purged.  What remains are the read query names prefaced with the @ symbol, read sequences and read base quality scores.  The meta data will be added back later in the MergeBam step.  See GATK Tutorial 6483 for details: https://software.broadinstitute.org/gatk/documentation/article?id=6483

By setting the CLIPPING_ATTRIBUTE to "XT" and by setting the CLIPPING_ACTION to 2, the program effectively removes the previously marked adapter sequences by changing their quality scores to two.  This makes it so they don't contribute to downstream read alignment or alignment scoring metrics.

This program will produce an interleaved .fastq file (paired reads are retained and marked appropriately in the same .fastq file).  This can then be fed into BWA_mem with the "-p" option.

The NON_PF=true option is set per the GATK tutorial.  This tells the program to retain reads marked with the 0x200 flag bit that denotes reads that do not pass quality controls (reads failing platform or vendor quality checks).
```
ALIGNMENT_RUN=<Sample ID>
INPUT=<path to input directory>"/"$ALIGNMENT_RUN"/uBAM_markedAdapters.bam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar SamToFastq \
    I=$INPUT \
    FASTQ=$OUTPUT_DIR/interleaved.fastq \
    CLIPPING_ATTRIBUTE=XT \
    CLIPPING_ACTION=2 \
    INTERLEAVE=true \
    NON_PF=true \
    TMP_DIR=<path to temp directory>/working_temp
```
