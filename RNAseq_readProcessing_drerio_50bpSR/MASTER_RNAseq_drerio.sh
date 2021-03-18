#!/bin/bash
# Author: Hannes Reinwald

# Full RNAseq pipeline for: DANIO RERIO

# MASTER SCRIPT for Danio rerio (dre) RNAseq processing
# This script will run:
# - Quality check on raw reads (fastQC)
# - Read mapping and feature gene counts (STAR)
# - 2nd QC via Samtools's flagstat function (samtools)
# - Summarize all QC reports in a MultiQC report (multiqc)

# To run this script go to your project folder containing a 'raw_reads' folder with
# fastq files and execute this script via:
# > bash /some/path/Drerio_RNAseq_readProcessing/dre_MASTER_RNAseq.sh

##############################################################################
# 
echo "Executing: `realpath $0`" >> runtime.log

# Activate Conda env
source activate rnaseq

# Specify path to script repo
ScriptPATH="`dirname $(realpath $0)`/bin" #more advanced and more flexible

# Runtime log 
SECONDS=0
printf "`date` - START RNAseq pipeline \n" >> runtime.log


### FastQ Screen on raw reads ----------------------------------------------------
printf "\n#### FastQ Screen ####" 2>&1 | tee -a runtime.log
echo "
 `date` ### START FastQ Screen on raw reads - Checking for contaminants ###
" 2>&1 | tee -a MasterRun.report
{ time bash $ScriptPATH/fastQscreen.sh 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log
echo "
 `date` ### END FastQ Screen on raw reads ###
" 2>&1 | tee -a MasterRun.report


### FastQC on raw reads & MultiQC report -----------------------------------------
printf "\n#### FastQC & MultiQC ####" 2>&1 | tee -a runtime.log
echo "
 `date` - START FastQC on raw reads & MultiQC summary
" 2>&1 | tee -a MasterRun.report
{ time bash $ScriptPATH/fastQC_multiQC.sh 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log
echo "
 `date` - END FastQC on raw reads & MultiQC summary
" 2>&1 | tee -a MasterRun.report


### STAR - Raw read (primer clipped no other trimming) mapping & gene count ------
printf "\n#### STAR mapping ####" 2>&1 | tee -a runtime.log
echo "
 `date` ### START STAR mapping, gene counts ###
" 2>&1 | tee -a MasterRun.report
{ time bash $ScriptPATH/STAR_alignGTF_drerio.sh 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log
echo "
 `date` ### END STAR mapping, gene counts ###
" 2>&1 | tee -a MasterRun.report


### Samtools ---------------------------------------------------------------------
printf "\n#### Samtools ####" 2>&1 | tee -a runtime.log
echo "
 `date` ### START Samtools - QC on aligned BAM files ###
" 2>&1 | tee -a MasterRun.report
{ time bash $ScriptPATH/samtools.sh 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log
echo "
 `date` ### END Samtools - QC on aligned BAM files ###
" 2>&1 | tee -a MasterRun.report


### In silico RNA integrity check!
#
#
#


### Final MultiQC report ---------------------------------------------------------
printf "\n#### MultiQC summary ####" 2>&1 | tee -a runtime.log
echo "
 `date` ### Start MultiQC summary ###
" 2>&1 | tee -a MasterRun.report
{ time bash $ScriptPATH/multiQC.sh 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log
echo "
 `date` ### End MultiQC summary ###" 2>&1 | tee -a MasterRun.report


### Merging GeneCount files to single CountMatrix via R --------------------------
printf "\n#### CountMatrix generation (R) ####" 2>&1 | tee -a runtime.log | tee -a MasterRun.report
printf "\n" >> MasterRun.report
{ time Rscript $ScriptPATH/CountMatrix_generator.R --verbose 2>&1 | tee -a MasterRun.report ; } 2>> runtime.log


### Analysis finished
echo "
 `date` - FINISHED SUCCESFULLY :) Yay!" 2>&1 | tee -a MasterRun.report

conda deactivate

# Runtime log
duration=$SECONDS
echo "
`date` - END RNAseq pipeline
$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." 2>&1 | tee -a runtime.log

exit 0
### END OF SCRIPT