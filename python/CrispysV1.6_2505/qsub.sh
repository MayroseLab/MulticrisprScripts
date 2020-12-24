#!/bin/tcsh
#$ -N MultiCRISPR_1452010925
#$ -S /bin/tcsh
#$ -cwd
#$ -l bioseq
#$ -e /bioseq/data/results/multicrispr/1452010925/$JOB_NAME.$JOB_ID.ER
#$ -o /bioseq/data/results/multicrispr/1452010925/$JOB_NAME.$JOB_ID.OU
cd /bioseq/data/results/multicrispr/1452010925
echo "/bioseq/data/results/multicrispr/1452010925/phylip_file.ph\nF\n/bioseq/data/results/multicrispr/1452010925/protdist_file.ph\nY" | protdist
