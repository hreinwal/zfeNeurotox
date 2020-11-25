#!/bin/bash

# Author: Hannes Reinwald
# Navigate to your project folder (/srv/attract/seq_files/project_folder) then execute this script.

# This script will:
# - map raw reads via STAR and will generate a mapped reads output file in BAM format
# - generate a Count List (ReadsPerGene.out.tab) similar to htseq count with default parameters
#       col1: gene ID   col2: counts of unstranded RNAseq col3: counts for 1st read col4: counts for 2nd read

#############################################################################################

#conda activate rnaseq
source activate rnaseq
home=$(pwd)

############
#   STAR   #
############

## STAR parameters 
# Specify path to your refference genome index files directory
#refgenome=/srv/attract/ref_genome/saccharomyces_cerevisiae/star_genome_index
refgenome=/srv/attract/ref_genome/danio_rerio.GRCz11.100/genome_index_files_50bp
CPUs=$[$(nproc)-2] # 14
OUT=BAM
TYPE=SortedByCoordinate #Unsorted (Sorting requiered for Samtools)
MAXmultimap=20  #Max number of multiple alignments allowed for a single read (Set to 1 for uniquely aligned reads); DEFAULT = 20
BAMsortRAM=20000000000  # > 10Gb!

echo "
 `date` - Start STAR mapping with the following settings:
    - MaxMultimap  : $MAXmultimap
    - Output format: $OUT $TYPE
    - nThreads     : $CPUs
    - BAM sort RAM : $[BAMsortRAM/1000000000] GB
    - Ref. Genome  : ${refgenome/"/star_genome_index"/}
"

# Create a STAR output folder
mkdir STAR

# Run STAR mapping with RAW READS (Untrimmed!!!)
cd ./raw_reads

#STAR --runThreadN $CPUs \
#--genomeDir $refgenome \
#--genomeLoad LoadAndExit
STAR

for SID in *.gz
do
    echo " 
    Processing sample file $SID in:
    $pwd"
	STAR --runThreadN $CPUs \
    --seedPerWindowNmax 49 \
    --runMode alignReads \
    --genomeDir $refgenome \
    --genomeLoad LoadAndKeep \
    --readFilesCommand gunzip -c \
    --readFilesIn $SID \
    --outFileNamePrefix ../STAR/${SID/fastq.gz/} \
    --outFilterMultimapNmax $MAXmultimap \
    --outSAMtype $OUT $TYPE \
    --limitBAMsortRAM $BAMsortRAM \
    --quantMode GeneCounts 
done

# Further potential options for STAR mapping run. 
# --readMapNumber 1000000 \
#--twopassMode Basic #this command can be only run when --genomeLoad NoSharedMemory

# remove loaded genome from shared memory storage
echo "
 Removing loaded genome from shared memory."
$StarPATH --runThreadN $CPUs \
--genomeDir $refgenome \
--genomeLoad Remove
# remove STAR generated folders as sideproduct of shared memory clearing
rm Log.progress.out Log.out Aligned.out.sam #SJ.out.tab Log.final.out
rmdir _STARtmp/
echo "
 `date` - Finished STAR mapping with Max Multimap = $MAXmultimap
"

# --------------------------------------------------------

# Now organize the output in ./STAR
cd ../STAR

#Compile Log.final.out files
mkdir log_files
mv -v *final.out ./log_files

#Compile Log.progress.out & Log.out files
mkdir log_files/progress_log
mv -v *progress.out ./log_files/progress_log

mkdir log_files/full_report
mv -v *Log.out ./log_files/full_report

#Compile SJ.out.tab files
mkdir splicing_counts
mv -v *SJ.out.tab ./splicing_counts

#Compile Reads per gene counts to a single folder
mkdir reads_per_gene
mv -v *ReadsPerGene.out.tab ./reads_per_gene

# --------------------------------------------------------

### Extract gene counts ###
# The following one liner code is designed to filter STAR's --quantMode GeneCounts output (FILENAME.ReadsPerGene.out.tab).
# The filtered output txt tables can then be directly used as input for CountMatrixGenerator.R script.
# The counts from STAR coincide with those produced by htseq-count with default parameters.
# STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
#column 1: gene ID
#column 2: counts for unstranded RNA-seq => this is the column we are interested for single read sequencing
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

# CODE: -------------------------------------
# cut: filter 1st and 2nd row of input file (specified by --fields (-f) 1-2)
# tail: remove first four rows (+NUM to output starting with line NUM)
# cut -f 1-2 *out.tab | tail -n +5 > GeneCounts.txt

# Start in your project dir, creating a GeneCounts folder
cd $home
mkdir GeneCounts

# navigating to STARs output dir
cd ./STAR/reads_per_gene/ 

for VAR in *.ReadsPerGene.out.tab; do
    echo 'Filtering STARs GeneCount output for:' $VAR
    cut -f 1-2 $VAR | tail -n +5 > ../../GeneCounts/${VAR/.ReadsPerGene.out.tab/_GeneCounts.txt}
done

# clip file names down to sample ID
cd ../../GeneCounts
for VAR in *_GeneCounts.txt; do
    [ -f "$VAR" ] || continue
    mv -vnT "$VAR" "${VAR/*R/R}" #matches everything from p to R in string and replaces by R
done

cd $home
exit 0

#### END OF SCRIPT ####