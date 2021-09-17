# Osteo-WGS-SNV-Workflow
Workflow for identifying single-nucleotide-variants (SNVs) in osteosarcoma Whole Genome Sequencing (WGS) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the Burrows-Wheeler Aligner.  The version used was bwa-0.7.17
- GATK version 4.0.12.0 (Picard included)
- All Python scripts were run on Python version 2.7.13 unless otherwise noted.  

# Workflow
**Notes:**
- Each sample has it's own directory for output files.  Each individual directory is labeled by the "Sample ID"

**Step 1) Downloading the hg38 genome files:** The following files were downloaded from the GATK resource bundle on 9/18/2018 and saved to a directory named "hg38_osteo".  The .fasta.gz file was then unzipped.  
- Homo_sapiens_assembly38.fasta.gz
- Homo_sapiens_assembly38.dict
- Homo_sapiens_assembly38.fasta.fai
- Homo_sapiens_assembly38.fasta.64.alt

```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.64.alt
```
- Note: The .fasta.64.alt file is used with BWA-MEM for ALT-aware alignment (Step ).  Once downloaded it was moved to the same directory the holds the BWA genome index files created in Step 2.  

**Step 2) Create a BWA genome index files:** The BWA "index" command is used to create the genome index files from the previously downloaded hg38 .fasta genome sequence.  These files will then be used in the upcoming alignment step.  BWA places the output files into the same foloder as the genome .fasta file was called from.  The BWA index files were then moved to a directory named "BWA_indexes/GATK_hg38" for later use.

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
**Step 4) Mark adapters:** Using Picard and the "MarkIlluminaAdapters" function, this step takes an uBAM file from the previous step and rewrites it with new adapter-trimming tags.  Per tutorial 6483 on the GATK website: https://software.broadinstitute.org/gatk/documentation/article?id=6483 This is needed so that the sequenced adapters contribute minimally to the alignments.  The tool adds an "XT" tag to a read record to mark the 5' start position of the specified adapter sequence.  It also produces a metrics file.  Note: each .bam from Step 3 may produce multiple uBAM files depending on how many read groups are present.  Each individual uBAM must undergo "MarkIlluminaAdapters".  Separate MarkIlluminaAdapter steps were carried out for the tumor and matched normal samples.

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
**Step 5) Convert uBAM to FASTQ:** Now we use Picard and the "SamToFastq" function to take the uBAM files, which have been modified so that Illumina adapter sequences are marked with the "XT" tag, and converts them to .fastq files for further processing.  It produces a .fastq file in which all extant meta data (read group info, alignment info, flags and tags) are purged.  What remains are the read query names prefaced with the @ symbol, read sequences and read base quality scores.  The meta data will be added back later in the MergeBam step.  See GATK Tutorial 6483 for details: https://software.broadinstitute.org/gatk/documentation/article?id=6483.  Separate SamToFastq steps were carried out for the tumor and matched normal samples.

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
**Step 6) Alignment of sequencing reads to the hg38 genome:** The program BWA is used with the "mem" function to align sequencing reads to the hg38 reference genome.  Sequencing reads from each lane are aligned individually.  They will later be combined into a single .sam file for each sample.  Separate alignments are carried out for the tumor and matched normal samples.  See GATK Tutorial 6483 for details: https://software.broadinstitute.org/gatk/documentation/article?id=6483

The "-M" flag tells the program to mark shorter split hits as secondary (for Picard compatibility)

The "-t" flag tells the program how many threads to use.

the "-p" flag tells the program that the input contains interleaved paired reads.
```
fasta_dir=<path to BWA index files created in Step 2>
ALIGNMENT_RUN=<Sample ID>
input=<path to input directory>"/"$ALIGNMENT_RUN"/lane_"${SLURM_ARRAY_TASK_ID}"_interleaved.fastq"

bwa-0.7.17/bwa mem -M -t 32 -p $fasta_dir/Homo_sapiens_assembly38.fasta $input
```
**Step 7) Post-Alt Processing:** Separate Post_Alt Processing steps were carried out for the tumor and matched normal samples.
```
ALIGNMENT_RUN=<Sample ID>
k8_DIR=<path to BWA directory>"/bwakit/bwa.kit"
POST_ALT_DIR=<path to BWA directory>"/bwa-0.7.17/bwakit"
ALT_INDEX=<path to BWA index files created in Step 2>"/Homo_sapiens_assembly38.fasta.64.alt"
ALIGNED=<path to input directory>"/"$ALIGNMENT_RUN"/aligned_lane_.sam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

$k8_DIR/k8 $POST_ALT_DIR/bwa-postalt.js $ALT_INDEX $ALIGNED > $OUTPUT_DIR/postalt_aligned_lane_.sam
```
- Note that the ".64" in the ALT_INDEX file indicates that this index file was generated with the version 0.6 or later of BWA and is the 64-bit index (as opposed to files generated by earlier versions, which were 32-bit).  
- Notes on ALT-aware alignment are described in this tutorial: https://software.broadinstitute.org/gatk/documentation/article.php?id=8017
- The ".fasta.64.alt" file is the ALT index file.  BWA-MEM uses this file to prioritize primary assembly alignments for reads that can map to both the primary assembly and an alternate contig.  The ALT index also contains decoy contig records as unmappped SAM records.  This is relevant to the postalt-processing step.  See https://software.broadinstitute.org/gatk/documentation/article.php?id=8017 for more details.
- Note that for this to work, the .fasta.64.alt basename needs to be the same as the other index and .dict files ("Homo_sapiens_assembly38" in this case).

**Step 8) MergeBamAlignments:** Here, we use Picard with the "MergeBamAlignment" function to merge defined information from the unmapped BAM with that of the aligned BAM to conserve read data (original read information and base quality scores).  The tool also generates additional meta information based on the information generated by the aligner, which may alter aligner-generated designations (mate info and secondary alignment flags). The tool then makes adjustments so that all meta information is congruent (e.g. read and mate strand info based on proper mate designations).  For details on this, see the GATK Tutorial 6483 at:
https://software.broadinstitute.org/gatk/documentation/article?id=6483.  Separate MergeBamAlignment steps were carried out for the tumor and matched normal samples.

The aligned BAM generated from the BWA_mem_align step lacks read group info and certain tags (UQ, MC, MQ as examples). It has hard-clipped sequences from split reads and altered base qualities.  MergeBamAlignment adjusts the read and read mate strand orientations for reads in a proper pair.  Finally, the alignment records are sorted by query name. MergeBamAlignment applies read group info from the uBAM and retains the program group info from the aligned BAM. In restoring original sequences, the tool adjusts CIGAR strings from hard-clipped to soft-clipped.  If the alignment file is missing reads present in the unaligned file, then these are retained as unmapped records. Additionally, the tool evaluates primary alignment designations according to a user-specified stragegy (e.g. for optimal mate paie mapping, and changes secondary alignment and mate unmapped flags based on its calculations).

The "PRIMARY_ALIGNMENT_STRATEGY=MostDistant" option marks primary alignments based on the alignment pair with the largest insert size.  This strategy is based on the premise that of chimeric sections of a read aligning to consecutive regions, the alignment giving the largest insert size with the mate gives the most information.

The "MAX_INSERTIONS_OR_DELETIONS=-1" option retains reads irregardless of the number of insertons and deletions.

The "ALIGNER_PROPER_PAIR_FLAGS" option is set to its default.  Thus, MergeBamAlignment will reassess and reassign proper pair designations made by the aligner.

The "ATTRIBUTES_TO_RETAIN=XS" option tells the program to keep reads flagged by BWA-mem's suboptimal alignment scores.

The "CLIP_ADAPTERS=false" option leaves reads unclipped.
```
ALIGNMENT_RUN=<Sample ID>
INPUT_ALIGNED=<path to input directory>"/"$ALIGNMENT_RUN"/postalt_aligned_lane_.sam"
INPUT_uBAM=<path to input directory>"/"$ALIGNMENT_RUN"/"<name of uBAM file matching the lane assignment of the ALIGNED file>
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN
FASTA_DIR=<path to directory containing the hg38 genome files downloaded in Step 1>

java -Xmx16G -jar picard.jar MergeBamAlignment \
    ALIGNED_BAM=$INPUT_ALIGNED \
    UNMAPPED_BAM=$INPUT_uBAM \
    O=$OUTPUT_DIR/MergedBamAlignment_postalt_lane_.bam \
    R=$FASTA_DIR/Homo_sapiens_assembly38.fasta \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=<path to temp directory>/working_temp_mban
```
**Step 9) Merge BAM files:** Here, we use Picard and the "MergeSamFiles" function to take the individual alignemnts from each of the WGS lanes and merge them into one complete .bam file

I have also included the "SORT_ORDER=coordinate" option in order to ensure that the .bam file is sorted correctly for downstread applications
```
#for n lanes

ALIGNMENT_RUN=<Sample ID>
INPUT_DIR=<path to input directory>"/"$ALIGNMENT_RUN
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN
#
file_1=$INPUT_DIR"/MergedBamAlignment_postalt_lane_1.bam"
file_2=$INPUT_DIR"/MergedBamAlignment_postalt_lane_2.bam"
#           .
#           .
#           .
file_<n>=$INPUT_DIR"/MergedBamAlignment_postalt_lane_<n>.bam"

java -Xmx8G -jar picard.jar MergeSamFiles \
    I=$file_1 \
    I=$file_2 \
#           .
#           .
#           .
    I=$file_<n> \
    SORT_ORDER=coordinate \
    O=$OUTPUT_DIR/aligned.bam \
    TMP_DIR=<path to temp directory>/working_temp_msn
```
**Step 10) Mark Duplicates:** Here, we use Picard with the "MarkDuplicates" function to mark any duplicate reads from a sequence alignment file.  Separate MarkDuplicate steps were carried out for the tumor and matched normal samples.

Per the GATK pipeline, I've included the CREATE_INDEX=true command which will create a index for the outputed .bam file.  This is needed
for downstream GATK tools

Also per the GATK pipeline, I've included the VALIDATION_STRINGENCY=SILENT command.
```
ALIGNMENT_RUN=<Sample ID>
INPUT_FILE=<path to input directory>"/"$ALIGNMENT_RUN"/aligned.bam"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

java -Xmx8G -jar picard.jar MarkDuplicates \
    I=$INPUT_FILE \
    O=$OUTPUT_DIR/rg_added_aligned_MarkedDup.bam \
    CREATE_INDEX=true \
    TAGGING_POLICY=All \
    VALIDATION_STRINGENCY=SILENT \
    M=$OUTPUT_DIR/Markdup_metrics.txt \
    TMP_DIR=<path to temp directory>/working_temp_mdn
```
**Step 11) Download hg38 known sites of variation:** The following file was downloaded from the GATK resource bundle on 8/17/2018 and saved to a directory named "hg38_osteo".  The file was then unzipped.
```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
gunzip dbsnp_146.hg38.vcf.gz
```
**Step 12) Index Feature File:** The .vcf file downloaded in Step 11 is indexed using the GATK tool, "IndexFeatureFile.  This will allow for querying features by a genomic interval.
```
INPUT_FILE=<path to directory containing the .vcf file downloaded in Step 11>"/dbsnp_146.hg38.vcf"

gatk --java-options "-Xmx4g" IndexFeatureFile -F $INPUT_FILE
```
**Step 13) Base Quality Score Recalibration:** The GATK tool, "BaseRecalibrator" is used to take a .bam file, the GATK reference for hg38 as well as a file containing the known sites of variation in hg38 according to dbsnp (downloaded from GATK site).  It produces a recal_data.table as the output which consists of several tables:
- The list of arguments
- The quantized qualities talbe
- The recalibration table by read group
- The recalibration talbe by quality score
- The recalibration table for all the optional covariates
```
ALIGNMENT_RUN=<Sample ID>
REF=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens_assembly38.fasta"
INPUT_FILE=<path to input directory>"/"$ALIGNMENT_RUN"/rg_added_aligned_MarkedDup.bam"
KNOWN_SITES=<path to directory containing the .vcf file downloaded in Step 11>"/dbsnp_146.hg38.vcf"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

gatk --java-options "-Xmx4g" BaseRecalibrator -R $REF -I $INPUT_FILE --known-sites $KNOWN_SITES -O $OUTPUT_DIR/recal_data.table
```
**Step 14) 2nd Pass Base Quality Score Recalibration:** The GATK tool, "BaseRecalibrator" is used to take the .bam file that was produced by the first past BQSR step and run a second pass using the GATK reference for hg38 as well as a file containing the known sites of variation in hg38 according to dbsnp (downloaded from GATK site).  It produces a second recal_data.table as the output that can be used in the "AnalyzeCovariates" step.
```
ALIGNMENT_RUN=<Sample ID>
REF=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens_assembly38.fasta"
INPUT_FILE=<path to input directory>"/"$ALIGNMENT_RUN"/recal_reads.bam"
KNOWN_SITES=<path to directory containing the .vcf file downloaded in Step 11>"/dbsnp_146.hg38.vcf"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

gatk --java-options "-Xmx4g" BaseRecalibrator -R $REF -I $INPUT_FILE --known-sites $KNOWN_SITES -O $OUTPUT_DIR/2ndP_recal_data.table
```
**Step 15) Analyze Covariates:** Here, we use the GATK tool, AnalyzeCovariates, to generate a document called "recalibration_plots.pdf" that contains plots that show how the reported base qualities match up to the empirical qualities calculated by the BaseRecalibrator. This is a method of checking the effect of the base recalibration process that will be applied to the sequence data.  For details see "Base Quality Score Recalibration (BQSR)" at
https://software.broadinstitute.org/gatk/documentation/article?id=11081

Takes as input a genome reference file, the before recal_data.table and the after post-recal_data.table.
```
ALIGNMENT_RUN=<Sample ID>
REF=<path to directory containing the hg38 genome files downloaded in Step 1>"/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR=<path to output directory>"/"$ALIGNMENT_RUN

gatk --java-options "-Xmx4g" AnalyzeCovariates \
    -before $OUTPUT_DIR/recal_data.table \
    -after $OUTPUT_DIR/2ndP_recal_data.table \
    -plots $OUTPUT_DIR/recalibration_plots.pdf
```
