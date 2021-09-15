# Osteo-WGS-SNV-Workflow
Workflow for identifying single-nucleotide-variants (SNVs) in osteosarcoma Whole Genome Sequencing (WGS) samples

# Version Notes
These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.7.1908 unless otherwise noted
Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts.
Alignment of sequencing reads was accomplished using the Burrows-Wheeler Aligner.  The version used was bwa-0.7.17
GATK version 3.6 (Picard included)
All Python scripts were run on Python version 2.7.13 unless otherwise noted.  
#Workflow
Step 1) 
