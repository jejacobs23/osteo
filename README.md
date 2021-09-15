# Osteo-WGS-SNV-Workflow
Workflow for identifying single-nucleotide-variants (SNVs) in osteosarcoma Whole Genome Sequencing (WGS) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
- Alignment of sequencing reads was accomplished using the Burrows-Wheeler Aligner.  The version used was bwa-0.7.17
- GATK version 3.6 (Picard included)
- All Python scripts were run on Python version 2.7.13 unless otherwise noted.  

# Workflow
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
