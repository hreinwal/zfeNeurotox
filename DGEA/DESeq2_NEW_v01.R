########################################################################
### DESeq2 RNA Seq Analysis Script for multiple factors / treatments ###
########################################################################

# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

### README ### -----------------------------------------------------------
# This script is desgined to run DESeq2 normalization and statistical testing on RNAseq experiments
# in Ecotox testing for multiple concentrations. This script will automatically adopt to the k numbers
# of samples and N numbers of conditions (High,Mid,Low or Treat1, Treat2, Treat3, ... , Control) and 
# will generate Log2FC values with respect to control groups.

# Requiered Input: 
# - CountMatrix (txt or csv is fine, if csv read.csv2 is used)
# - coldata.csv

# Main Analysis Steps are: 
# 1) RLE normalization (DESeq2) across all samples element of CountMatrix
# 2) Calc LFC (with respect to control) and LFcutOff (upper 90% quantile of abs(LFC values))
# 3) apeglm shrinking on LFC

# 4) Multiple t-testing with Benjamin-Hochberg correction (padj < 0.05) 
#    and independet hypothesis weighing (IHW) to identify DEGs
#    for LFC and apeglm(LFC) values for H0: LFC = 0; apeglm(LFC) = 0
# 4.1) Annotation of Genelists via org.Dr.eg.db [https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html]
# 4.2) DEGs with LFC / apeglm(LFC) > < abs(LFcutOff) are kept as potential molecular marker features

# 5) Plotting
# 5.1 MA- & Vulcano plots for N-1 conditions for LFC [LFcut & padj cut]
# 5.2 MA- & Vulcano plots for N-1 conditions for apeglm(LFC) [LFcut & padj cut]
# 5.3 DataQC
#     - Correlation bio. replicates (rlog(meanCounts))
#     - MeanSD plots (rlog(meanCounts))
#     - RLE_normalisation
#     - normCounts_rlogTransformation

# 6) Output: 
# 6.1 DESeq2 Results (xcel; still working on html output though ... )
#     - N-1 DESeq2 result files for LFC
#     - N-1 DESeq2 result files for apeglm(LFC)
# 6.2 DEGs (padj < 0.05 and LFcut for: 
#     - N-1 files for LFC
#     - N-1 files for apeglm(LFC)
##############

### LOAD PACKAGES ### -------------------------------------------------
require(DESeq2) #installed
require(IHW)    #installed
require(apeglm) #installed
require(qvalue)
#require(locfdr)
#require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(hexbin)
#####################


### DATA IMPORT & Color settings ### ---------------------------------------------------
message("\nImporting CountMatrix and coldata files ...")

## coldata import
tmp = list.files(full.names = T, pattern = "[Cc]oldata")
if(grepl(".csv",tmp)==T){
  coldata = read.csv2(file = tmp, row.names = 1, header = T)
  # In case csv file was not differently exported run read.csv
  if (ncol(coldata) <= 1) {
    coldata = read.csv(file = tmp, row.names = 1, header = T)
  }
} else {
  # if not ending with csv import with read.delim2
  coldata = read.delim2(file = tmp, row.names = 1, header = T)
  if (ncol(coldata) <= 1) {
    # if file was differently formated use read.delim
    coldata = read.delim(file = tmp, row.names = 1, header = T)
  }
}
coldata <- droplevels(coldata[which(coldata$Condition != ""), # filtering for non empty rows
                              which(!(grepl("X.",colnames(coldata))) == T)]) # filtering cols
coldata$Condition <- gsub(" ","",coldata$Condition) #removing any spaces in Condition levels
coldata$Tank <- gsub("-","_",coldata$Tank)

# set factor levels
coldata$Tank      <- factor(coldata$Tank)
coldata$Condition <- factor(coldata$Condition)
coldata$Substance <- factor(coldata$Substance)
substance         <- levels(coldata$Substance) #Extract the name(s) of the tested substance(s)

# relevel factors in Condition in right order:
condition = levels(coldata$Condition)
tmp = as.character(coldata$Condition[grep("[Cc]ontrol",coldata$Condition)][1]) #extract control factor level
if(all(grepl("[Ee]xposure",condition[2:4]))==T & length(condition)==4) {# Control, LE, ME, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,4,2)])
} else if (all(grepl("[Ee]xposure",condition[2:3]))==T & length(condition)==3) { # Control, LE, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,2)])
} else { # just specifying the reference level --> Control
  coldata$Condition <- relevel(coldata$Condition, ref = tmp)
}
message(paste0("Reference level:\t",levels(coldata$Condition)[1]),"\nSorting ",length(condition)-1," Treatments:\t",paste(levels(coldata$Condition)[-1], collapse = " "))
condition <- levels(coldata$Condition)

# CountMatrix import
tmp = list.files(full.names = T, pattern = "[Cc]ount[Mm]atrix")
if(any(grepl(".txt",tmp))==T) { # import txt
  tmp1 = tmp[grepl(".txt",tmp)]
  count.matrix = read.table(file = tmp1, row.names = 1, header = T)
  if(ncol(count.matrix) <= 1) {
    count.matrix = read.delim2(file = tmp1, row.names = 1, header = T)
    if(ncol(count.matrix) <= 1) {
      count.matrix = read.delim(file = tmp1, row.names = 1, header = T)
    }
  }
} else { # import csv
  tmp1 = tmp[grepl(".csv",tmp)]
  count.matrix = read.csv2(tmp1, header = T, row.names = 1, fill = F)
  if(ncol(count.matrix) <= 1) {
    count.matrix = read.csv(tmp1, header = T, row.names = 1, fill = F)
  }
}
rm(tmp,tmp1)

# subset cols in count.matrix after rows in coldata. So there is only the need to trim the coldata file manually
count.matrix <- subset(count.matrix, select = row.names(coldata))

# Sort & check correct data import:
count.matrix <- count.matrix[,rownames(coldata)] # reduce the to analyse samples in mtx by the samples listed in coldata
stopifnot(all(rownames(coldata) == colnames(count.matrix))) # Must be TRUE!
message("Done!\n")

### Set annotation colors for downstream plotting ###
col.Cond = colorRampPalette(RColorBrewer::brewer.pal(n=8, name="YlOrRd"))(length(levels(coldata$Condition)))
col.Cond[1] <- "#56B4E9" #Defines color for control
col.Tank = colorRampPalette(c("gray95","gray50"))(length(levels(coldata$Tank)))
ann_colors = list(Tank = setNames(col.Tank, levels(coldata$Tank)),
                  Condition = setNames(col.Cond, levels(coldata$Condition)))
#scales::show_col(col.Cond)

####################################

### Setup output environment ###
dir.create("DESeq2_Pairwise", showWarnings = F)
setwd("DESeq2_Pairwise")
home = getwd()

##################
###   DESeq2   ###
##################
## Create DESeq2 object -------------------------------------------------------
p = 0.05 # padj cutoff
message(paste0("\n Starting DESeq2 Analysis - Pairwise Wald's t-test and IHW \n padj < ",p,"\n Tested Substance: \t",substance,"\n"))

# Import entire count matrix into a DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = coldata,
                              design = ~ Tank + Condition) # Multifactor model controlling for variability in Tanks

## Count Matrix filtering -----------------------------------------------------
message(paste0("\nRemoving low abundant counts from count matrix.",
               "\nMinimum number of gene counts per row: \t\t",ncol(count.matrix)))

# Filter count matrix - removing low count genes
keep = row.names(counts(dds)[rowSums(count.matrix) >= length(colnames(count.matrix)), ])
dds  <- dds[keep,]

# Get a list of different disp estimates for later QC plotting
fitType = c("parametric", "local") #additonally "mean" possible
estDisp.ls <- list()
for(i in fitType){
  message(paste("\nEstimate Dispersions for fitType: ",i))
  x <- estimateSizeFactors(dds)
  x <- estimateDispersions(x, fitType = i)
  estDisp.ls[[i]] <- x
}
rm(x,i,fitType,keep)

## RUN DESeq() ----------------------------------------------------------------
dds <- DESeq(dds, test = "Wald")                   # Pairwise comparison
#dds <- DESeq(dds, test = "LRT", reduced = ~ Tank) # ANOVA-like approach

## Extract norm. & transformed count.matrix ------------------------------------
normMtx =list()
normMtx[["norm"]]   = round(counts(dds, normalized = T),3) # [[1]]  normalized read counts
normMtx[["ntd"]]    = round(assay(normTransform(dds)),3)   # [[2]] (n+1)log2 transformed norm counts
normMtx[["nrl"]]    = round(assay(rlog(dds, blind = F)),3) # [[3]] rlog transformed mean read counts
normMtx[["nrl.bl"]] = round(assay(rlog(dds, blind = T)),3) # [[4]] rlog transformed mean read counts; BLIND!

message("\nExporting DESeq2 normalized count matrix. This might take a while...\nSaving to txt table ...")
write.table(normMtx[[1]], paste0(substance,"_DESeqNormCounts.txt")) # Extract DESEq2 normalized count matrix
#write.table(df.rld, paste0(substance,"_rlogDESeqNormCounts.txt")) # Extract rlog(norm counts) matrix
message("Done!\n")

## Extract DESeq2 result tables to list objects --------------------------------
resNames = resultsNames(dds)[grepl("^[Cc]ondition_",resultsNames(dds))]
## Non-shrunk results
res.ls <- list() #create an empty list to store objects in
for(i in resNames) {
  x = gsub("[Cc]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, non-shrunk Log2FC] ..."))
  res.ls[[gsub("_.+","",x)]] <- results(dds, name = i, filterFun = ihw, alpha = p)
  message("Done!\n")
}
## apeglm shunk results
resLfs.ls <- list() # save log2FC shrunk results in different list
for(i in resNames) {
  x = gsub("[cC]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, apeglm shrunk Log2FC] ..."))
  resLfs.ls[[gsub("_.+","",x)]] <- lfcShrink(dds, coef = i, type = "apeglm", 
                                             res = res.ls[[gsub("_.+","",x)]])
  message("Done!\n")
}
rm(x,i)

## Append qvalue to res.ls / resLfs.ls -----------------------------------------
message("\nCalculating and appending qvalues to DESeq2 results ...")
# res.ls
res.ls = lapply(res.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
# resLfs.ls
resLfs.ls = lapply(resLfs.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
message("Done!\n")
###################


###################
###   biomaRt   ###
###################
## Create annotation object for DESeq2 result tables ---------------------------
message("\nAnnotating ENSEMBL Gene IDs in result tables using biomaRt ...")
# if rerio mart object found locally load it, create new mart object from host 
if(file.exists("~/biomaRt/drerio_mart.Robj")) {
  message("\nDanio rerio mart object found locally. \nLoading 'drerio_mart.Robj' into R session ...")
  load("~/biomaRt/drerio_mart.Robj")
} else if (file.exists("S:/data/biomaRt/drerio_mart.Robj")){
  message("\nDanio rerio mart object found in S:/data/biomaRt/ \nLoading 'drerio_mart.Robj' into R session ...")
  load("S:/data/biomaRt/drerio_mart.Robj")
  } else {
    message("\nCould not find 'drerio_mart.Robj'. Creting new mart object from 'www.ensembl.org'",
            "\nMake sure you have a working interent connection. Otherwise this will fail!\nConnecting to server ...")
    rerio = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="https://www.ensembl.org")
}

## Build gene annotation object & merge ----------------------------------------
id     = rownames(res.ls[[1]])
idType = "ensembl_gene_id"
attr   = c("ensembl_gene_id","external_gene_name","description","gene_biotype",
           "entrezgene_id")
message("\nAnnotating Ensembl gene IDs ...")
GeneAnno = biomaRt::getBM(attributes = attr, mart = rerio, uniqueRows = T,
                          filters = idType, values = id)
# !!! In a few cases there can be more than one entrez id for a single ensembl gene id !!!
# Here we remove duplicated entries (~ 360 out of ~ 25 000, so it's really a minor fraction)
# As it turns out higher integer values in the entrez ID correspond to a deeper level of organization
# i.e. ID 999 -> Protein A, ID 1000348 -> Protein A subunit a1. For our purpose it is good enough 
# to retrieve the overall gene / protein information.
# => We sort GeneAnno by entrezgene_id and then remove the duplicates with !duplicate().
# That way, in case multiple entrez ids are assigned to a single ensembl id the lower entrez id
# will be kept for downstream analysis. 
GeneAnno = dplyr::arrange(GeneAnno, entrezgene_id)
GeneAnno = GeneAnno[!duplicated(GeneAnno$ensembl_gene_id),] #keep only non duplicated
stopifnot(length(id) == nrow(GeneAnno))
row.names(GeneAnno) = GeneAnno$ensembl_gene_id
# Finally, to ensure compatability with downstream GSEA scripts rename external_gene_name & entrezgene_id
colnames(GeneAnno) = c("ensembl_gene_id","SYMBOL","description","biotype","ENTREZID")

# Now merge GeneAnno with DESeq2 result tables
AnnoFun = function(x){
  tmp = merge(as.data.frame(x), GeneAnno, by=0)[,-1]
  row.names(tmp) = tmp$ensembl_gene_id
  tmp[,-which(colnames(tmp) %in% "ensembl_gene_id")]
}
res.ls    = lapply(res.ls, FUN = AnnoFun)
resLfs.ls = lapply(resLfs.ls, FUN = AnnoFun)
# Finally sort res.ls / resLfs.ls after the padj value
res.ls    = lapply(res.ls, function(x){dplyr::arrange(x, padj)})
resLfs.ls = lapply(resLfs.ls, function(x){dplyr::arrange(x, padj)})
rm(id, rerio)
message("Done!\n")
###################


###  EXPORT DESeq2 RESULTS  ### -------------------------------------------------------
message("Exporting DESeq2 result tables. This might take a while :)")
dir.create("Results", showWarnings = F)
for(k in names(res.ls)){
  n = gsub("[Ee]xposure","",k)
  message(paste("Saving DESeq2 result table",k,"[IHW, non-shrunk Log2FC] to csv ..."))
  write.csv2(res.ls[[k]], file = paste0("Results/",substance,"_res_",n,".csv"))
  message(paste("Saving DESeq2 result table",k,"[IHW, apeglm shrunk Log2FC] to csv ..."))
  write.csv2(resLfs.ls[[k]], file = paste0("Results/",substance,"_reslfs_",n,".csv"))
}
message("Done!\n")

## Select only DEGs for export and downstream plotting ##
deg.ls    <- lapply(res.ls, function(x){subset(x, padj <= p)}) # pcut
degCut.ls <- lapply(res.ls, function(X){                       # pcut & log2FC cutoff 
  lfcut = quantile(abs(X$log2FoldChange), 0.9) #(90% quantile of abs(log2FC))
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  x = subset(X, padj <= p)
  subset(x, abs(log2FoldChange) >= lfcut)
})
# Selecting only DEGs with padj < p & abs(apeglm shrunk Log2FC) > lfcut (90% quantile from abs(non shrunk log2FC))
degLfs.ls = list()
for(k in names(resLfs.ls)){
  # determine lfc cut off based on non-shrunk values
  lfcut = quantile(abs(res.ls[[k]]$log2FoldChange), 0.9)
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  message(paste("Log2FC cut off for:\t",k,"\t----->\t",round(lfcut,2)))
  # Select padj < p from resLfs.ls and subset by log2FC
  degLfs.ls[[k]] <- subset(subset(resLfs.ls[[k]], padj <= p), abs(log2FoldChange) >= lfcut)
}

## Exporting DEGs in separate result csv tables ##
dir.create("DEGs", showWarnings = F)
out = c("unshrink","unshrinkFCcut","lfsFCcut") # output dir. for different levels of filtering
for(k in out){dir.create(paste0("DEGs/",k))}
# Export to respective dir
for(k in names(deg.ls)) {
  message(paste("Saving selection of DEGs for",k,"to csv ..."))
  n = gsub("[Ee]xposure","",k)
  write.csv2(deg.ls[[k]],paste0("DEGs/",out[1],"/",substance,"_pcut_",n,".csv"))   # unshrink
  write.csv2(degCut.ls[[k]],paste0("DEGs/",out[2],"/",substance,"_pcut_LFcut_",n,".csv"))# unshrink + FCcut
  write.csv2(degLfs.ls[[k]],paste0("DEGs/",out[3],"/",substance,"_lfs_pcut_LFcut_",n,".csv"))# shrink + FCcut
}
message("Done!\n")

rm(n,k,lfcut,out)
###############################


######################################   PLOTTING   ##########################################
## Output directories for plots
dir.create(paste0(home,"/Plots"), showWarnings = F)
dir.create(paste0(home,"/Plots/DataQC"), showWarnings = F)
message("\nStarting Data Composition and QC plotting ...")

### QC - DESeq2 norm. & count distr. & countMtx filtering  ### -----------
message("DESeq2 normalization & gene count distribution ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_CountNorm_mtxFiltering.pdf"),
    width = 9, height = 11.7, onefile = T, bg = "transparent", fg ="black")
par(mfrow=c(4,2))

## Normalization bar plots ##
barplot(colSums(counts(dds)),
        ylab= "Total gene counts",
        main= "Raw counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
barplot(colSums(counts(dds, normalized=T)),  #plot mean normalized counts
        ylab= "Total gene counts",
        main= "DESeq2 norm. counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
boxplot(log10(counts(dds)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "Raw counts")
boxplot(log10(counts(dds, normalized=T)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "DESeq2 norm. counts")

## Gene count distr. & countMtx filtering ##
cut = ncol(count.matrix) # min counts per row
# quantile distribution of gene count matrix row sums
tmp = apply(count.matrix, 1, sum)
tmp = quantile(log10(tmp+1), probs = seq(0, 1, .05))
tmp = data.frame(x=seq(0,1,.05)*100, y=tmp)
plot(tmp, xlab = "Quantile-%", ylab = "log10(gene row sum + 1)", main = "Gene row sum quantile distribution")
abline(h = log10(cut+1), lwd=1.5, lty=2, col = "firebrick")
text(x = 0, y = log10(cut+1)*1.3, pos = 4, offset = 0.5,labels = paste("Cutoff:",cut))

# rowsum cutoff vs remaining number of genes in count matrix
cutOffPlot <- function(countMtx, cut) {
  n = ncol(countMtx)
  if(missing(cut)){cut = n}
  if(n < 6) {stop("Number of count matrix columns < 6")}
  X = seq(n-6,n+50,1) 
  X[X == 0] <- 1
  X = sort(union(X, c(1:5)))
  
  rNbr = c() # empty vector to store Nbr of genes per cutoff in
  for(x in X) {
    tmp  = nrow(countMtx[which(rowSums(countMtx) >= x), ])
    rNbr = append(rNbr, tmp)
  }
  df <- data.frame(x = X, y = rNbr)
  plot(df, xlab = "Row sum cut off", ylab = "Nbr of genes",
       main = "Genes in count matrix")
  abline(v=cut, lwd=1.5, lty=2, col = "firebrick")
  #abline(h=df[cut,2], lwd=1.5, lty=2, col = "blue")
  text(x = cut, y = df[cut,2], pos = 4, offset = 1.5,
       labels = paste("Cutoff:",cut,"=",df[cut,2],"genes"))
}
cutOffPlot(count.matrix, cut)

# gene count distribution before AND after row sum cutoff (based on on row means)
tmp  = apply(count.matrix,1,mean)
tmp1 = apply(count.matrix[which(rowSums(count.matrix) >= cut), ], 1, mean)
hist(log2(tmp+1), breaks = 40, main = "countMtx", 
     ylim = c(0,2500), xlim = c(0,20),
     xlab = "log2(Row mean count +1)")
hist(log2(tmp1 +1), breaks = 40, main = paste("filtered countMtx - min. rowSum:",cut),
     ylim = c(0,2500), xlim = c(0,20),
     xlab ="log2(Row mean count +1)") #histogram of filtered counts

dev.off()
rm(cut,tmp, tmp1)
##############################################################

### QC - DESeq2's dispersion estimate models ### ---------------------------------------------
message("DESeq2's dispersion estimate models (estimateDispersions) ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_DispEstimates.pdf"),
    width = 7, height = 5.9, onefile = T, bg = "transparent")
for(i in names(estDisp.ls)){
  plotDispEsts(estDisp.ls[[i]], main = paste(substance,'- DESeq2 fit type:',i))
}
dev.off()
rm(estDisp.ls, i)
gc() #free some memory 
################################################

### QC - pvalue and LFC distribution ### -----------------------------------------------------
message("pvalue and LFC distribution ...")
while (!is.null(dev.list()))  dev.off()
n = length(condition)-1
pdf(file = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.pdf"),
    width = 3*n, height = 10, onefile = T, bg = "transparent")
#png(filename = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.png"),
#    width = 6*n, height = 18, units = "cm", bg = "white", pointsize = 7, res = 450)
par(mfrow=c(3,n))
## pval distr
for(i in names(res.ls)){
  hist(res.ls[[i]]$pvalue, main= paste(substance,i),col = "gray50", 
       border = "gray50", xlab = "p-value", breaks = 40)
}
## pval vs padj & qval
pVSpadjPlot <- function(res, pval=.05, q=.05, p=.05, title = ""){
  plot(res$pvalue, res$padj, xlab="pvalue", ylab="conversion of pvalue",
       xlim= c(0,0.25), ylim = c(0,1), main= paste(substance,title),
       col = "dodgerblue1",pch = 20)
  points(res$pvalue, res$qval, col="springgreen3", pch = 20)
  abline(h=p, lwd=1, lty=2)
  abline(v=pval, lwd=1, lty=2)
  xpos <- 0.1
  text(x = .25, y = .2, pos = 2, offset = 0, col = "black",
       labels = paste("Genes with pvalue <",pval,":\t",sum(res$pvalue <= pval, na.rm = T)))
  text(x = .25, y = .15, pos = 2, offset = 0, col = "springgreen3",
       labels = paste("Genes with Qvalue <",q,":\t",sum(res$qval <= q, na.rm = T)))
  text(x = .25, y = .1, pos = 2, offset = 0, col = "dodgerblue",
       labels = paste("Genes with padj (BH) <",p,":\t",sum(res$padj <= p, na.rm = T)))
}
for(i in names(res.ls)) { pVSpadjPlot(res.ls[[i]], title = i) }
## Log2FC distr.
lfcDistPlot <- function(res, LFcut=quantile(abs(res$log2FoldChange),.9), title="") {
  hist(res$log2FoldChange, breaks= 500, main= paste(substance,title), col= "gray50", 
       border= "gray50", xlab= "Log2(fc)", xlim= c(-3,3), ylim= c(0,300))
  abline(v= c(LFcut,-LFcut), col= "dodgerblue", lty= 2, lwd= 1)
  legend(x="topright", text.col = "dodgerblue", bty = "n",
         legend = paste0("LFcut: +/-",round(LFcut, digits = 2)))
}
for(i in names(res.ls)) { lfcDistPlot(res.ls[[i]], title = i) }
dev.off()
########################################

### QC - Norm. gene count transformation ### -------------------------------------------------
message("Norm. count matrix transformation ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_NormCount_Transformation.pdf"),
    width = nrow(coldata)-2, height = 4, onefile = T, bg = "transparent")
par(mfrow=c(1,3))
boxplot(normMtx$norm , notch = TRUE,las = 3, #rotating sample labels 90
        main = "Normalized read counts", ylab = "norm. read counts")
boxplot(normMtx$ntd, notch = TRUE, las = 3, #rotating sample labels 90
        main = "log2 Transformation", ylab = "log2 (norm. read counts + 1)")
boxplot(normMtx$nrl.bl, notch = TRUE, las = 3, #rotating sample labels 90
        main = "rlog Transformation", ylab = "rlog (norm. read counts)")
dev.off()
############################################

### QC - PearsonCor of biol. replicates ### ----------------------------------------------------
message("Biological replicate gene count correlation ...")
id.ls = list()
for(i in condition) {id.ls[[i]] = rownames(coldata[coldata$Condition %in% i,])}
# get all possible combinations of samples for each condition
comb = lapply(id.ls, function(x){gtools::combinations(length(x),2,x)})
# plot dim.
w = nrow(comb[[1]])   #width
h = length(condition) #height

# plotting function
corplot.1 = function(df, tmp, Lab="", title="") {
  d  = densCols(df[,tmp[1]], df[,tmp[2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  ggplot(df) + geom_point(aes(df[,tmp[1]], df[,tmp[2]], col = d), size = .8, alpha = .4) +
    labs(x=paste(Lab,tmp[1]),y=paste(Lab,tmp[2]),
         subtitle = (paste0(title," [",coldata[tmp[1],"Tank"]," vs ",coldata[tmp[2],"Tank"],"]"))) +
    annotate("text", label = paste0("Pearson: ",round(cor(df[,tmp[1]], df[,tmp[2]]),2)),
             x = (max(df[,tmp[1]])), y = (min(df[,tmp[2]])), hjust=1, vjust=0) + 
    annotate("text", label = paste0("R2: ",round((cor(df[,tmp[1]], df[,tmp[2]]))^2,2)),
             x = (min(df[,tmp[1]])), y = (max(df[,tmp[2]])), hjust=0, vjust=1) +
    scale_color_identity() + coord_equal(ratio=1) + theme_bw()
}
multiCorPlot = function(mtx, label="") {
  gg.ls = list()
  for(k in names(comb)){
    x = comb[[k]]
    for(i in 1:nrow(x)){
      tmp = x[i,] #IDs to plot with
      df = as.data.frame(mtx[ ,tmp])
      gg.ls[[paste0(k,i)]] <- corplot.1(df, tmp, Lab = label,
                                        title = gsub("[Ee]xposure","",k))
    }
  }
  gg.ls
}

## log2 ##
gg = multiCorPlot(normMtx$ntd, label = "log2(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_log2.png"),
    width = 7.33*w, height = 7.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 500)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()

## rlog ##
gg = multiCorPlot(normMtx$nrl.bl, label = "rlog(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_rlog.png"),
    width = 7.33*w, height = 7.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 500)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()

# clear env
rm(w,h,i,gg,comb)
gc()
###########################################



## Session Information & save Rdata  ## ------------
sink(paste0(home,"/DESeq2_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()
message("Saving RData object. This might take a while ...")
save.image(paste0(home,"/DESeq2_pairwise.RData"))

message(paste0("\n\nFinished Wald's pairwise testing DESeq2's DEG analysis & data plotting.",
               "\nAll outputs were saved in:\n",home,
               "\n\nWhat a ride! Finally END of SCRIPT! :)\nJ.A.R.V.I.S. over and out"))
###    END OF SCRIPT   ####
