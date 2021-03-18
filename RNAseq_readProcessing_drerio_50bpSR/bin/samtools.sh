#!/bin/bash

# Author: Hannes Reinwald
# Execute this script within the STAR folder containing the sorted BAM files after mapping.

################
#   Samtools   #
################

# Conda environment
#conda activate rnaseq
source activate rnaseq

cd STAR
# Running samtools with flagstat, index and idxstats on SORTED! BAM file 
for BAM in *.out.bam
do 
    echo " 
 `date` - Starting samtool analysis with ${BAM/.Aligned.sortedByCoord.out.bam/}"
    samtools flagstat $BAM > ${BAM/.Aligned.sortedByCoord.out.bam/}".txt"
    echo " Indexing ${BAM/.Aligned.sortedByCoord.out.bam/} ..."
    samtools index $BAM
    samtools idxstats $BAM > ${BAM/.Aligned.sortedByCoord.out.bam/}"_idxstats.txt"
    echo " `date` - Finished samtool analysis with ${BAM/.Aligned.sortedByCoord.out.bam/}"
done 2>&1 | tee -a samtools.report

# the upper code section is not very efficient yet. Uses only one core at a time. Try to spread it for multiple cores
# https://www.biostars.org/p/259509/
#ls *.bam | parallel --progress --eta -j 8 'samtools sort -@ 7  {} > {.}_sorted.bam'

# This could work ...
# ls *.out.bam | parallel --progress --eta -j 8 'samtools flagstat {} > ${/.Aligned.sortedByCoord.out.bam/}".txt"'


# Organize samtools output
mkdir idxstats
mv *idxstats.txt ./idxstats

mkdir flagstat
mv *.txt ./flagstat

#Go back to project folder
cd ..
conda deactivate
exit 0
### END ###
