#!/bin/bash
# Author: Hannes Reinwald

#  This script, executed in your project folder will navigate to the raw_reads folder and will generate 
# a FastQC report for each file in it. At the End MultiQC will summarize the results in a single html file.

########################################
# got to raw_reads
cd ./raw_reads

# create fastqc outdir
mkdir ../fastQC
touch ../fastQC/fastqc.progress

### FastQC --------------------------------------------------
# Activate Conda Env
source activate multiqc #need the older fastQC version because of java runtime env problems
CPUs=$[$(nproc)-2] #threads to use for multithreading

fastqc --outdir ../fastQC/ \
--threads $CPUs \
--noextract \
--format fastq \
--kmers 7 $(find *.fastq.gz) 2>&1 | tee -a ../fastQC/fastqc.progress

conda deactivate
cd ..

### MultiQc ------------------------------------------------
source activate rnaseq #newer version of
mkdir multiQC   # output dir
NAME="`basename "$(pwd)"`"

# Run multiqc (for all parameters check multiqc --help)
multiqc --zip-data-dir \
    --outdir ./multiQC \
    --filename ${NAME}"_multiQC"\
    -v -f .  #the "." specificies the current wd! ~ to "./"

# Deactivate the conda env
conda deactivate

exit 0
### END of Script ###