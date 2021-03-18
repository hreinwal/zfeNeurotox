##############################
## COUNT MATRIX GENERATION ###
##############################

### README ####
# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# This script will take the multiple output files from STAR ~GeneCounts or HTseq count
# (Each file = Gene count list for the respective sample)
# and will arrange the count numbers in a single data array called "CountMatrix.csv"
# and save the file in the same directory as the count files. 

setwd("./GeneCounts") #set working directory (where *_geneCounts.txt are located); for RNAseq pipeline

# Count Matrix generation
message(" Importing GeneCount files ...")
myFiles <- list.files(pattern="*[Gg]ene[Cc]ounts.t[xs][tv]") #import all file names matching the string
count.files <- lapply(myFiles, FUN=read.delim, header=FALSE) #Creates a list of datasets col1=ENSEMBLgenID col2=counts

### PROOFREADING ####
# Sort gene lists (so they are all equally sorted)
#head(count.files[[1]],10)
#head(count.files[[1]][order(count.files[[1]]$V1),],10) #great!
genesort <- function(x){
  x <- x[order(x$V1),]
}
sorted <- lapply(count.files, FUN = genesort)
# Check if gene count order is identical throughout all files
for (i in c(2:length(sorted))){
  stopifnot(all(sorted[[1]]$V1==sorted[[i]]$V1))
}
count.files <- sorted
rm(sorted)
####################

colNbr <- 2 # specifies the 2nd col! 
# For STAR's GeneCounts output colNbr = 2, 
# depending on the strandedness of your library number may be 2,3 or 4. 
count.matrix <- as.data.frame(sapply(count.files, function(x) x[,colNbr])) 
ID <- gsub("^.*sR","R", myFiles)
ID <- gsub("_.*","", ID)     # creates a vector list of unique sample identifiers (ID)
colnames(count.matrix) <- ID # assigns ID to respective columns 
row.names(count.matrix) <- count.files[[1]]$V1 # assigns gene names as row names in count.data
count.files <- NULL  # Empty this large data set
message(" GeneCount files imported. Saving CountMatrix.csv ...")

# Export CountMatrix
# write.csv2(count.matrix, file = paste0("CountMatrix.csv"))
# better to write to txt file as txt format is requiered for ArrayExpress upload
setwd('../') # this line important for RNAseq pipeline to get project dir name into file name
fna = paste0(gsub("_.*","",basename(getwd())),"_CountMatrix.txt") #file name
write.table(count.matrix, file = fna, sep = "\t", row.names = T, col.names = T, dec = ".")

message(paste0("\n That's it!\n\n ",fna," saved in: ",getwd(),"\n Well done! :)"))
##### END OF SCRIPT #####
