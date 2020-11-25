#!bin/bash
# Author: Hannes Reinwald

##################
## FastQ Screen ##
##################
home=$(pwd)
CPUs=$[$(nproc)-1]

# Specify path to config file
# Config file specifies genomes to map against. Look at the file for details. 
# All bowtie2 indexed ref genomes to map against can be found under:
# /srv/attract/ref_genome/FastQ_Screen_Genomes
configFile=/srv/attract/ref_genome/FastQ_Screen_Genomes/fastQscreen_mod.conf

# go to raw reads
cd raw_reads

# output dir
mkdir ../fastQ_screen
touch ../fastQ_screen/fastQscreen.log

# Check if provided file exists
if [ -f "$configFile" ]; then
	echo "$configFile exists." 2>&1 | tee -a ../fastQ_screen/fastQscreen.log
else
	echo "$configFile does not exists. 
	Exiting script" 2>&1 | tee -a ../fastQ_screen/fastQscreen.log
	cd $home
	exit 0
fi

# Run fastq_screen
source activate rnaseq
fastq_screen --threads $CPUs \
    --conf $configFile \
    --aligner bowtie2 \
    --subset 100000 \
    --force \
    --outdir ../fastQ_screen $(ls *.fastq.gz) 2>&1 | tee -a ../fastQ_screen/fastQscreen.log

cd $home
exit 0
##### END OF SCRIPT ###### 