#!/bin/bash
# Author: Hannes Reinwald

# This script, executed in your project folder will run MultiQC to summarize QC results in a single html file.

########################################
# Conda environment
#conda activate rnaseq
source activate rnaseq

### MultiQc ###
mkdir multiQC   # output dir
NAME="`basename "$(pwd)"`"

# Run multiqc (for all parameters check multiqc --help)
multiqc --zip-data-dir \
    --outdir ./multiQC \
    --filename ${NAME}"_multiQC"\
    -v -f .  #the "." specificies the current wd! ~ to "./"

# Deactivate the conda multiqc env again
conda deactivate
exit 0
### END of Script ###