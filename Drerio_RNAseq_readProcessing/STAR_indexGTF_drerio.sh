#!/bin/bash

# This script will: 
# 1) Generate a user defined folder for GTF and toplevel ref genome files. 
# 2) Will run STAR and generate ref. index files and store in ./genome_index_files

# It is kind of tricky to make this script scalable for any available genome. So please consider
# revising the binbits parameter for index generation as this scales after Nbrs contigs in your 
# genome as well as genome size. A rough estimate is computed as follows:
# binbits = log2(genome size in bytes (unzipped!) / Number of contigs)

# You can check the number of contigs in your fasta genome file using:
# awk '/>/ {a++} END {print "number of sequences in this file: " a}' genome.fa
# or
# grep -c "^>" genome.fa

# For the zebrafish genome this would be: 
# 1.5 GB = 1.5e9 bytes => log2(1.5e9/993) ~ 20

###################################################################################
# Conda environment
#conda activate rnaseq
source activate rnaseq
cd /srv/attract/ref_genome

#  Define variables 
echo "
Hello there friend! :)

My name is J.A.R.V.I.S. and I am programmed to help you with generating genome index files via STAR aligner.

I will generate a directory (folder) for you in which you are supposed to store your fasta genome file and respective gtf file. 
These files will be then used to generate genome index files which are requiered for read mapping via STAR. 

Lets start by creating your reference genome directory in which we will store everything.

Please type in your species scientifc name. 
Only use small letters and separate genus and species with _. 
i.e. 'danio_rerio'"
read species

echo "
Please type in the reference genome version and release version. 
Do not use empty space to separate objects. Use '.'' instead.
i.e. GRCz11.100"
read version

#  Create folders named after defined variables
mkdir ./$species"."$version
mkdir ./$species"."$version/genome_index_files #output dir
touch ./$species"."$version/Genome_INFO.txt

# Inform user about the created env and tell them to copy toplevel and gtf file into /srv/attract/ref_genome/$species.$version""
echo "
Your reference genome directory $species.$version was generated. 
Please copy the zipped (*.fa.gz) reference genome and gtf (*.gtf.gz) file into this directory:
 /srv/attract/ref_genome/$species.$version

You can use the following comand for that:
 cp /copy/file/from/filename /srv/attract/ref_genome/$species.$version


If you don't know where to look for these files you can check the ENSEMBL repo online:
 ftp://ftp.ensembl.org/pub/release-100
If you wish to directly download the files run the following commands (exemplarly for D.rerio):
 cd /srv/attract/ref_genome/$species.$version
 wget ftp://ftp.ensembl.org/pub/release-100/gtf/danio_rerio/Danio_rerio.GRCz11.100.gtf.gz
 wget ftp://ftp.ensembl.org/pub/release-100/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz


An empty text file (Genome_INFO.txt) was generated for you in $species.$version 
Please manually annotate the follwoing information in the text file.
    - Species & Genome version 
    - Source (i.e. URL code)
    - Additional info like:
        - Long or short reads? 
        - Paired end or single reads? 
        - Which sequencing method (Illumina,...)

The more Info the better ;) Feel free to add as much information about the reference genome as possible.
"

# Check with user whether files were transfered into the folder
read -p "
Did you copy the ref.genome file (*.fa.gz) and the GTF file (*.gtf.gz) into $species.$version? 
Do NOT continue with this script until you are done copying.

Do you want to continue (y/n)?" -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
   echo "
Cool! B-)"
else
    echo "
Sorry to hear that. Please copy the files into /srv/attract/ref_genome/$species.$version and rerun the script"
    exit 0
fi

# Define STAR index generating parameters
echo "
OK Great! :) Before we tell STAR to generate our genome index files we need to set a few parameters first.
What is your sequence read length? For 50 bp just type: 50"
read nbp

echo "
What genomeChrBinNbits size should be used? (For D.rerio.GRCz11 20 is fine)
If you are not familiar with this parameter please check the STAR's manual."
read binbits

CPUs=$[$(nproc)-2]  # Number of available CPUs - 2
overhang=$[$nbp-1]  # read length-1 (i.e. 50 bp; overhang= 50-1 = 49) 
#binbits=20           # Scalable factor to limit ram usage of STAR. For further Info please check STAR manual. 

# Now execute STAR
read -p "
OK now everything is set. Do you really want to continue with STAR genome index generation(y/n)?" -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "
Great we will continue then. Star is now generating your genome index files.
Depending on your genome size this might take some time."
    cd ./$species"."$version
    #unpigz -tv *.gtf.gz *.fa.gz test if files corrupted
    unpigz -kf *.gtf.gz *.fa.gz

    # output dir for genome index files
    outdir='genome_index_files_${nbp}bp'
    #Running STAR
    STAR --runThreadN $CPUs \
    --runMode genomeGenerate \
    --genomeDir ./$outdir \
    --genomeFastaFiles ./*dna.primary_assembly.fa \
    --sjdbGTFfile ./*.gtf \
    --sjdbOverhang $overhang \
    --genomeChrBinNbits $binbits
else
    echo "Boooooo!"
    exit 0
fi

echo "If everything worked fine, your genome index files are located in:

/srv/attract/ref_genome/$species.$version/genome_index_files

Use this path when running the read mapping job with STAR.
J.A.R.V.I.S. over and out."
exit 0
##################### END #######################################