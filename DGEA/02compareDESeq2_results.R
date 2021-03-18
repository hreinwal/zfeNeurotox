# Author:   Hannes Reinwald
# Contact:  hannes.reinwald@ime.fraunhofer.de
# Type:     R script
# source("S:/manuscripts/neurotox/Rscripts/compareDESeq2.R")

### PACKAGES ### --------------
require(DESeq2)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(data.table)
require(plotrix)
require(pcaMethods)
require(Rtsne)
require(biomaRt)
################

# Prepare output environment
home = "S:/manuscripts/neurotox/data/"
out = "compareDESeq2"
dir.create(paste0(home,out), showWarnings = F)
setwd(paste0(home,out))

# File location for import
pathToData = "S:/data/RNA_seq_gene_counts/"
dataSet=c("Abamectin","Chlorpyrifos_2","Fipronil_2","Imidacloprid","Carbaryl","Methoxychlor")

### Import Data ### ---------------
message(paste0("\nImporting DESeq2 results for:\n",paste0(dataSet, collapse = " ")))
count.ls <- list()
coldat.ls <- list()
for(k in dataSet) {
  message(paste0("Importing data for:\t",k))
  setwd(paste0(pathToData,k))
  tmp <- list.files(pattern = "CountMatrix.txt")
  mtx <- read.delim(file = tmp, header = T, row.names = NULL)
  colnames(mtx)[1] <- "ID"
  mtx$ID <- as.factor(mtx$ID)
  count.ls[[k]] <- mtx
  tmp <- list.files(pattern = "coldata.csv")
  coldat.ls[[k]] <- read.csv2(tmp, header = T, row.names = 1)
  setwd(paste0(home,out))
}
rm(k,tmp,mtx)

# filter and merge coldata files
coldat.ls <- lapply(coldat.ls, function(x){
  #x <- coldat.ls[[1]]
  df <- subset(x, select = c("Condition","Tank","Substance","Group","MoA_Group","MoA",
                             "IRAC_Classification","SamplingDate","Conc_ug.L","SendForSeq","Solvent",
                             "Coagulated_24hpf","Dead_96hpf","Hatched_48hpf","Hatched_96hpf","RIN"))
  df$ID <- row.names(df)
  # Cut out empty rows!
  df[which(row.names(df) != ""),]
})
coldat <- as.data.frame(data.table::rbindlist(coldat.ls))
rownames(coldat) <- as.character(coldat$ID)

coldat$Condition <- as.factor(as.character(coldat$Condition))
coldat$Tank <- as.factor(as.character(coldat$Tank))
coldat$Substance <- as.factor(as.character(coldat$Substance))

#shapiro.test(coldat$RIN) # W = 0.87971, p-value = 2.59e-05 => Data not normal distributed! 

# Re-level the factor Condition!
lev <- levels(coldat$Condition)
if(length(lev) == 4 & all(grepl("Exposure", lev[2:length(lev)]))){
  message("\nSorting 3 Treatment levels as:\tLow, Mid & High Exposure \nControl set as reference level.\n")
  lev = lev[c(1,3,4,2)]
  coldat$Condition <- factor(coldat$Condition, levels = lev)
  
} else if(length(lev) == 3 & all(grepl("Exposure", lev[2:length(lev)]))) {
  message("\nSorting 2 Treatment levels as:\tLow & High Exposure \nControl set as reference level.")
  lev = lev[c(1,3,2)]
  coldat$Condition <- factor(coldat$Condition, levels = lev)
  
} else {
  message("\nControl set as reference level for multiple treatments.")
  coldat$Condition <- relevel(coldat$Condition, ref = "Control")
}

# cbind count.ls to countMtx
myMerge <- function(df1,df2){merge(df1,df2, by="ID", all = T)}
countMtx <- Reduce(myMerge, count.ls)
row.names(countMtx) <- countMtx$ID
countMtx <- as.matrix(subset(countMtx, select = -ID)) # Great! :) Now we have our final countMtx

# Check if imported data adds up
stopifnot(all(colnames(countMtx) %in% rownames(coldat)))
stopifnot(sort(row.names(coldat)) == sort(colnames(countMtx)))

# Reorder countMtx & coldat file
# Create ID order by sorting ID to its substance and condition to reorder countMtx
id.order <- c()
for(i in levels(coldat$Substance)){
  for(k in levels(coldat$Condition)){
    id <- rownames(coldat[coldat$Condition %in% k & coldat$Substance %in% i,])
    id.order <- append(id.order,id)
  }
}
coldat <- coldat[id.order,]
countMtx <- countMtx[, id.order]
stopifnot(colnames(countMtx) == row.names(coldat))

# remove ls objects
rm(coldat.ls, count.ls)
###################


### Some minor data exploration ### -----
pdf(file = "QC_CountFrequency.pdf", width = 15, height = 15, onefile = T, bg = "transparent", fg ="black")
par(mfrow=c(2,2))
# plot row sums
tmp <- apply(log10(countMtx + 1), 1, sum) #get the sum of each row (log10(N+1))
hist(tmp, breaks = 100, xlab = "log10(N+1) row counts", main = "Histogram of gene row counts",
     ylim = c(0,1300))

q <- 0.02 #Quantile steps
tmp <- quantile(apply(log10(countMtx + 1), 1, sum), probs = seq(0,1,q), type = 8)
df <- data.frame(x=seq(0,1,q)*100,y=tmp)
plot(df, xlab = "Quantile-%", ylab = "log10(N+1) row counts", main = paste0(q*100,"% Quantiles"))

# plot row means
tmp <- apply(log10(countMtx + 1), 1, mean) #get the mean of each row (log10(N+1))
hist(tmp, breaks = 100, xlab = "log10(N+1) row mean counts", main = "Histogram of row mean counts",
     ylim = c(0,800))

q <- 0.02 #Quantile steps
tmp <- quantile(apply(log10(countMtx + 1), 1, mean), probs = seq(0,1,q), type = 8)
df <- data.frame(x=seq(0,1,q)*100,y=tmp)
plot(df, xlab = "Quantile-%", ylab = "log10(N+1) row mean counts", main = paste0(q*100,"% Quantiles"))
dev.off()

rm(df,tmp,q,k,i,id)
#15%-quantile in this case is ~64, so setting the minimum count cutoff to rowsum >= ncol(mtx) seems reasonable
###################################


### DESeq2 data norm ### --------------------------
# Idea: As the idea is to mainly compare control to treated samples we could add a feature to 
# coldat separating samples in two groups: Control - Treated. Hence everything which is not control is
# then considered as treated. The design could look like this then: design = Feature
# This could potentially improve the vst model to clearly separate between induced effects
# and variance among control groups.
f <- c() #empty feature vector
for(i in coldat$Condition) {if(i == "Control"){f <- append(f,"Control")}else{f <- append(f,"Treated")}}
stopifnot(nrow(coldat)==length(f))
coldat$feature <- as.factor(f)
coldat$batch <- as.factor(paste(coldat$Tank, coldat$Substance, sep="."))
rm(f,i)

# Import entire count matrix into a DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData=countMtx, colData=coldat, design = ~ batch + Condition) #Multifactor Level
#dds <- DESeqDataSetFromMatrix(countData=countMtx, colData=coldat, design = ~ Substance)
#dds <- DESeqDataSetFromMatrix(countData=countMtx, colData=coldat, design = ~ feature)

dds <- dds[rowSums(countMtx) >= ncol(countMtx),] # Filter DDS object
dds <- estimateSizeFactors(dds) # Deseq2's median of ratios normalization
normMtx <- round(counts(dds, normalized=TRUE),2)

## Use following section to check manually which fitType is best suited
skip = T
if(skip != T) {
  ## Get a list of different disp estimates for later QC plotting
  estDisp.ls <- list()
  fitType = c("parametric", "local")
  for(i in fitType){
    message(paste0("\nDispersion estimates for fitType:\t",i))
    x <- estimateSizeFactors(dds) # have to create a new DESeq Object (x) to store different fitType disp. estimate in it.
    x <- estimateDispersions(x, fitType = i)
    estDisp.ls[[i]] <- x
  }
  ## Plot dispEstimates
  pdf(file = "QC_DispEstimates.pdf", width = 7, height = 5.9, onefile = T, bg = "transparent")
  for(i in fitType){
    plotDispEsts(estDisp.ls[[i]], main = paste("DESeq2 dispEstimate fit type:",i))
  }
  dev.off()
  rm(estDisp.ls, x, i, fitType, skip)
} else {rm(skip)}
##

## Run vst with best fitting disp estimate model
bestFit <- "parametric"
dds     <- estimateDispersions(dds, fitType = bestFit)
vst.bl  <- assay(vst(dds, blind = T, fitType = bestFit)) # blind=T; use for QC. Assesses data unbiased by Condition or Tank
vst     <- assay(vst(dds, blind = F, fitType = bestFit)) # blind=F; use for downstream analysis if QC ok. Takes full use of data information.


### Plot RLE Norm. ### 
pdf(file = "QC_DESeq2_RLE_Norm.pdf", width = 20, height = 11, onefile = T,  #multiple figures in one file
    bg = "transparent", fg ="black")
par(mfrow=c(2,2))
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
dev.off()

# Clear Env
#rm(dds)
#gc()
########################


### Color Settings ### ---------------------------------------------------------------
stopifnot(colnames(normMtx) == row.names(coldat)) #check
Substance = levels(as.factor(coldat$Substance))
Condition = levels(as.factor(coldat$Condition))
Tank = levels(as.factor(coldat$Tank))
# color blind friendly palette [https://www.r-bloggers.com/2013/10/creating-colorblind-friendly-figures/]
cbl = c("#00FF4D", "#E69F00", "#56B4E9","#F0E442","#CC79A7","#0072B2","#009E73","#D55E00","#005DFF")
if(length(cbl) > length(Substance)){cbl = cbl[1:length(Substance)]}else{cbl=cbl}
col.Subs = colorRampPalette(cbl)(length(Substance))
col.Cond = colorRampPalette(c("white","gray20"))(length(Condition))
col.Tank = colorRampPalette(c("gray95","gray50"))(length(Tank))
ann_colors = list(Tank = setNames(col.Tank, Tank),
                  Substance = setNames(col.Subs, Substance),
                  Condition = setNames(col.Cond, Condition))
rm(lev)
######################


### Dissim. Mtx ### ----------------------------------------------------------------
myDismtx <- function(mtx, method, top, title) {
  if(missing(method)){method = "euclidean"}
  if(missing(title)){title = ""}
  if(missing(top)){top = 2000}
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting Sample Dist for Var.Top:\t\t",top))
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # transpose input, calculate sample distance and create a distance matrix
  sampleDist <- dist(X, method = method) #must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distMtx <- as.matrix(sampleDist)
  #rownames(distMtx) <- paste(coldat$Condition, coldat$Tank, sep="-")
  
  # create annotation object for heatmap
  ann_col <- subset(coldat, select = c("Condition","Substance"))
  ann_row <- subset(coldat, select = "Tank")
  #rownames(ann_row) <- paste(coldat$Condition, coldat$Tank, sep="-")
  
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(15)
  pheatmap::pheatmap(
    distMtx,
    angle_col = "90",
    #treeheight_col = 40, #default 50
    #fontsize = 12, #default 10
    #fontsize_number = 0.7*12,
    #display_numbers = T,
    #cellwidth = 26,
    #cellheight = 26,
    drop_levels = T,
    number_format = if(method == "manhattan"){"%1.0f"}else{"%.1f"},
    main = paste0(title," - SampleDist [",method,"] - Var.Top: ",top),
    clustering_distance_rows = sampleDist,
    clustering_distance_cols = sampleDist,
    annotation_col = ann_col, # Condition!
    annotation_row = ann_row, # Tank!
    annotation_colors = ann_colors,
    col = colors)
}

while (!is.null(dev.list()))  dev.off()
pdf(file = "SampleDistMtx_bl.pdf", width = 10, height = 8.3, onefile = T,  #multiple figures in one file
    bg = "transparent", fg ="black")
for(method in c("euclidean","canberra","manhattan","maximum")){
  message(paste0("\nComputing Sample Dist with dist measure:\t",method))
  if(method == "maximum") {
    myDismtx(vst.bl, method, top = nrow(vst.bl), title = "blind")
  } else {
    for(top in c(1000,5000,10000,nrow(vst.bl))){
      myDismtx(vst.bl, method, top = top, title = "blind")
    }
  }
}
dev.off()

pdf(file = "SampleDistMtx.pdf", width = 10, height = 8.3, onefile = T,  #multiple figures in one file
    bg = "transparent", fg ="black")
for(method in c("euclidean","canberra","manhattan","maximum")){
  message(paste0("\nComputing Sample Dist with dist measure:\t",method))
  if(method == "maximum") {
    myDismtx(vst, method, top = nrow(vst))
  } else {
    for(top in c(1000,5000,10000,nrow(vst))){
      myDismtx(vst, method, top = top)
    }
  }
}
dev.off()

rm(method, top)
gc()
###################


### PCA - t-SNE ### ----------------------------------------------------------------
gg.bl = list() #for blind vst transformed data - list to store ggplot objects in
gg = list()
topVar = c(1000,5000,10000) #selecting top N genes with maximum variance!

## Colors 
col.Subs = colorRampPalette(cbl)(length(Substance))
col.Cond = colorRampPalette(brewer.pal(n=8, name="YlOrRd"))(length(Condition))
col.Cond[1] <- "#56B4E9" #Defines color for control
names(col.Cond) <- levels(coldat$Condition)

## PCA ##
myPCA <- function(mtx, pcaM ="svd", top = 2000, title = "") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting PCA for Var.Top:\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  pc <- pca(X, method = pcaM, center = T, nPcs=2)
  pcaDf <- merge(coldat, scores(pc), by=0)
  ggplot(pcaDf, aes(PC1, PC2, colour = Condition, shape = Substance)) +
    geom_point(size = 3, alpha = .65) + scale_colour_manual(values = col.Cond) +
    ggtitle(paste0(title," - expl.Var[%]: ",round(pc@R2cum[2]*100,1)," - ",pcaM," - Top:",pc@nVar)) +
    xlab(paste0("PC1: ",round((pc@R2)[1]*100,1),"% variance")) +
    ylab(paste0("PC2: ",round((pc@R2)[2]*100,1),"% variance")) + #stat_ellipse() +
    theme_bw() + theme(aspect.ratio = 1)
}

for(pcaM in listPcaMethods()[1]){
  message(paste0("\nRunning PCA method:\t\t",pcaM))
  for(top in topVar){
    name <- paste0("pca.",top)
    gg.bl[[name]] <- myPCA(vst.bl, pcaM, top, title = "blind")
    gg[[name]] <- myPCA(vst, pcaM, top)
  }
}

## t-SNE ##
message(paste0("\nRunning t-SNE"))
mytSNE <- function(mtx,top,title) {
  if(missing(title)){title = ""}
  if(missing(top)){top = 2000}
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting t-SNE for Var.Top:\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # perform Rtsne calc.
  perplex <- (nrow(X)-1)/3
  set.seed(42) #to make results reproducable! 
  tsne <- Rtsne::Rtsne(X, dims=2, initial_dims=nrow(X), check_duplicates=F, num_threads=2,
                       pca = TRUE, perplexity = perplex, #(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
                       theta = 0.0, #Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
                       partial_pca = F, #(requires the irlba package). This is faster for large input matrices (default: FALSE)
                       max_iter = 10000, #number of iterations
                       is_distance = FALSE, pca_center = TRUE, pca_scale = FALSE, verbose = F,
                       normalize = F) #Default True; Set to F as DESeq's RLE norm was performed prior!
  row.names(tsne$Y) <- row.names(X)
  colnames(tsne$Y) <- c("tSNE_1","tSNE_2")
  
  # ggplot
  ggtsne <- merge(coldat, tsne$Y, by=0)
  ggplot(ggtsne) + geom_point(aes(x=tSNE_1, y=tSNE_2, color = Condition, shape = Substance), size=3, alpha = .7) +
    ggtitle(paste0("t-SNE on vst[norm counts] ",title," - Top: ",top)) + xlab('t-SNE 1') + ylab('t-SNE 2') +
    scale_colour_manual(values = col.Cond) + theme_bw() + theme(aspect.ratio = 1)
}
for(top in topVar){
  name <- paste0("tSNE.",top)
  gg.bl[[name]] <- mytSNE(vst.bl, top, title = "blind")
  gg[[name]] <- mytSNE(vst, top)
}

## Plot ##
while (!is.null(dev.list()))  dev.off()
pdf(file = "PCA_tSNE.pdf", width = 5.4*length(topVar), height = 5.2*2, 
    onefile = T, bg = "transparent", fg ="black")
print(ggpubr::ggarrange(plotlist = gg.bl, ncol=length(topVar), nrow=2))
print(ggpubr::ggarrange(plotlist = gg, ncol=length(topVar), nrow=2))
dev.off()


## Individual PCA & t-SNE plots for each Substance
myPCA2 <- function(mtx, pcaM = "svd", top = 2000, title = "") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting PCA for ",title," Var.Top:\t\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  pc <- pca(X, method = pcaM, center = T, nPcs=2)
  pcaDf <- merge(coldat, scores(pc), by=0)
  
  ggplot(pcaDf, aes(PC1, PC2, colour = Condition, shape=Tank)) +
    geom_point(size = 3, alpha = .65) + scale_colour_manual(values = col.Cond) +
    ggtitle(paste0(title," - expl.Var[%]: ",round(pc@R2cum[2]*100,1)," - ",pcaM," - Top:",pc@nVar)) +
    xlab(paste0("PC1: ",round((pc@R2)[1]*100,1),"% variance")) +
    ylab(paste0("PC2: ",round((pc@R2)[2]*100,1),"% variance")) + #stat_ellipse() +
    theme_bw() + theme(aspect.ratio = 1) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
}
mytSNE2 <- function(mtx, top = 2000, title = "") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting t-SNE for ",title," Var.Top:\t\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # perform Rtsne calc.
  perplex <- (nrow(X)-1)/3
  set.seed(42) #to make results reproducable! 
  tsne <- Rtsne::Rtsne(X, dims=2, initial_dims=nrow(X), check_duplicates=F, num_threads=2,
                       pca = TRUE, perplexity = perplex, #(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
                       theta = 0.0, #Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
                       partial_pca = F, #(requires the irlba package). This is faster for large input matrices (default: FALSE)
                       max_iter = 10000, #number of iterations
                       is_distance = FALSE, pca_center = TRUE, pca_scale = FALSE, verbose = F,
                       normalize = F) #Default True; Set to F as DESeq's RLE norm was performed prior!
  row.names(tsne$Y) <- row.names(X)
  colnames(tsne$Y) <- c("tSNE_1","tSNE_2")
  
  # ggplot
  ggtsne <- merge(coldat, tsne$Y, by=0)
  ggplot(ggtsne) + geom_point(aes(x=tSNE_1, y=tSNE_2, color = Condition, shape = Tank), size=3, alpha = .7) +
    ggtitle(paste0("t-SNE on vst[norm counts] ",title," - Top: ",top)) + xlab('t-SNE 1') + ylab('t-SNE 2') +
    scale_colour_manual(values = col.Cond) + theme_bw() + theme(aspect.ratio = 1)
}

# PCA
gg <- list()
gg.bl <- list()
for(k in levels(coldat$Substance)) {
  # subset vst & vst.bl for PCA
  id = rownames(coldat[coldat$Substance %in% k,])
  gg[[k]] <- myPCA2(vst[,id], title = k, top = 5000)
  gg.bl[[k]] <- myPCA2(vst.bl[,id], title = paste0(k,".blind"), top = 5000)
  
}

N <- ceiling(length(levels(coldat$Substance))/3)
pdf(file = "Indiv_PCA.pdf", width = 16, height = N*4, 
    onefile = T, bg = "transparent", fg ="black")
print(ggpubr::ggarrange(plotlist = gg,    ncol=3, nrow=N))
print(ggpubr::ggarrange(plotlist = gg.bl, ncol=3, nrow=N))
dev.off()

# t-SNE
gg <- list()
gg.bl <- list()
for(k in levels(coldat$Substance)) {
  # subset vst & vst.bl for PCA
  id = rownames(coldat[coldat$Substance %in% k,])
  gg[[k]] <- mytSNE2(vst[,id], title = k, top = 5000)
  gg.bl[[k]] <- mytSNE2(vst.bl[,id], title = paste0(k,".blind"), top = 5000)
  
}

N <- ceiling(length(levels(coldat$Substance))/3)
pdf(file = "Indiv_tSNE_individual.pdf", width = 16, height = N*4, 
    onefile = T, bg = "transparent", fg ="black")
print(ggpubr::ggarrange(plotlist = gg,    ncol=3, nrow=N))
print(ggpubr::ggarrange(plotlist = gg.bl, ncol=3, nrow=N))
dev.off()


# Clear env
rm(top,topVar,name,pcaM,gg,gg.bl) #k, id, N
gc()
###################


### Coldata data visualization ###
## RIN values ----------------------------------------------------------------------------------
pdf(file = "QC_RINscores.pdf", width = 7.3, height = 7, onefile = T, bg = "transparent", fg ="black")
ggboxplot(coldat, x = "Substance", y = "RIN", color = "Substance", notch = F, error.plot = "sd",
          add = "jitter", legend = "none") + rotate_x_text(angle = 35)+ ylim(0,10.2)+xlab(NULL)+
  scale_color_manual(values = ann_colors$Substance)+
  geom_hline(yintercept = mean(coldat$RIN), linetype = 3)+ # Add horizontal line at base mean
  annotate("text", x=.5, y=0.5, hjust = "inward", label=paste(
    "RIN mean: ",round(mean(coldat$RIN),2),"\nRIN range: ",min(coldat$RIN)," - ",max(coldat$RIN)))+
  geom_hline(yintercept = 8, linetype = 2, color = "firebrick1")+ # Add min RIN quality cut off value
  stat_compare_means(method = "kruskal.test", label.y = 1.2, hjust = "inward")+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", label.y = 10.1,
                     ref.group = ".all.", hide.ns = F)
dev.off()


## Survival rate -------------------------------------------------------------------------------
coldat$Hatched_48hpf <- as.numeric(coldat$Hatched_48hpf)
coldat$Hatched_96hpf <- as.numeric(coldat$Hatched_96hpf)
tmp = coldat
# Compute survival 24 and 96 hpf
tmp$survival24h <- 100 - (tmp$Coagulated_24hpf/15)*100 #15 is the number of individuals used in the experiment
tmp$survival96h <- 100 - ((tmp$Coagulated_24hpf + tmp$Dead_96hpf)/15)*100
tmp$Hatched_48hpf <- tmp$Hatched_48hpf*100
tmp$Hatched_96hpf <- tmp$Hatched_96hpf*100

pdf(file = "QC_SurvivalRate.pdf", width = 7, height = 7, onefile = T, bg = "transparent", fg ="black")
#24 hpf
ggboxplot(tmp, x = "Condition", y = "survival24h", color = "Condition", notch = F, error.plot = "sd",
          add = "jitter", facet.by = "Substance", short.panel.labs = T) +
  rotate_x_text(angle = 35) + ylim(-.1,105) + xlab(NULL) + ylab("Survival rate 24 hpf [%]") +
  scale_color_manual(values = col.Cond) + geom_hline(yintercept = 80, linetype = 3) +
  stat_compare_means(method = "kruskal.test", label.y = 5, hjust = "inward") + # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Control", hide.ns = F)
#96 hpf
ggboxplot(tmp, x = "Condition", y = "survival96h", color = "Condition", notch = F, error.plot = "sd",
          add = "jitter", facet.by = "Substance", short.panel.labs = T) +
  rotate_x_text(angle = 35) + ylim(-.1,105) + xlab(NULL) + ylab("Survival rate 96 hpf [%]") +
  scale_color_manual(values = col.Cond) + geom_hline(yintercept = 80, linetype = 3) +
  stat_compare_means(method = "kruskal.test", label.y = 5, hjust = "inward") + # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Control", hide.ns = F)
dev.off()


## Hatching rate -------------------------------------------------------------------------------
pdf(file = "QC_HatchRate.pdf", width = 7, height = 7, onefile = T, bg = "transparent", fg ="black")
#24 hpf
ggboxplot(tmp, x = "Condition", y = "Hatched_48hpf", color = "Condition", notch = F, error.plot = "sd",
          add = "jitter", facet.by = "Substance", short.panel.labs = T) +
  rotate_x_text(angle = 35) + ylim(-.1,105) + xlab(NULL) + ylab("Hatching rate 48 hpf [%]") +
  scale_color_manual(values = col.Cond) +
  stat_compare_means(method = "kruskal.test", label.y = 99, hjust = "inward") + # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Control", hide.ns = F, label.y = 90)
#96 hpf
ggboxplot(tmp, x = "Condition", y = "Hatched_96hpf", color = "Condition", notch = F, error.plot = "sd",
          add = "jitter", facet.by = "Substance", short.panel.labs = T) +
  rotate_x_text(angle = 35) + ylim(-.1,105) + xlab(NULL) + ylab("Hatching rate 96 hpf [%]") +
  scale_color_manual(values = col.Cond) +
  stat_compare_means(method = "kruskal.test", label.y = 50, hjust = "inward") + # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Control", hide.ns = F, label.y = 103)
dev.off()
###################################


# -------------------------------------------------------------------------------------
###  COMPARING DEGs  ###
# -------------------------------------------------------------------------------------

### Data Import ### --------------------------
# File location for import
pathToData = "S:/data/RNA_seq_gene_counts/"
dataSet=c("Abamectin","Chlorpyrifos_2","Fipronil_2","Imidacloprid","Carbaryl","Methoxychlor")

# Objects to store DESeq2 result tables in 
res.ls <- list()
reslfs.ls <- list()
# Import DESeq2 result tables
message("\nStart comparing DEG lists")
for(k in dataSet) {
  message(paste0("Importing DESeq2 results for:\t",k))
  setwd(paste0(pathToData,k,"/DESeq2_Pairwise/Results"))
  # res - non-shrunk lfc results
  tmp <- list.files(pattern = "_res_")
  for(i in tmp) {
    n = gsub(".csv","",i)
    n = gsub(".+?_","",n)
    n = paste0(k,".",n)
    res.ls[[n]] <- read.csv2(i, header = T, row.names = 1)
  }
  # reslfs - lfc shrunk results
  tmp <- list.files(pattern = "_reslfs_")
  for(i in tmp) {
    n = gsub(".csv","",i)
    n = gsub(".+?_","",n)
    n = paste0(k,".",n)
    reslfs.ls[[n]] <- read.csv2(i, header = T, row.names = 1)
  }
  setwd(paste0(home,out)) # back home again
}
message("Done importing DESeq2 result files.\n")
rm(tmp,k,i,n)

# Select differentially expressed genes across all conditions.
# In this case we are picking the genes which have a padj <= 0.05 &
# with an apeglm(lfc) >= LFcut.
# LFcut is the 90% qunatile of the abs(non-shrunk lfc) values
lfcut.ls <- lapply(res.ls, function(x) {quantile(abs(x$log2FoldChange),.9)})
deg.ls <- list()
degPcut.ls <- list() #for non LFcut filtered DEGs
for(k in names(reslfs.ls)) {
  df <- subset(reslfs.ls[[k]], padj <= .05, select = c(baseMean:padj, SYMBOL, ENTREZID)) #select by padj value
  degPcut.ls[[k]] <- df
  df <- df[which(abs(df$log2FoldChange) >= lfcut.ls[[k]]),] #select by lfcut
  deg.ls[[k]] <- df
}
rm(k,df, lfcut.ls)

commonGenes <- Reduce(intersect, lapply(res.ls, row.names))
 #24676 common observed genes across all experiments

tmp = lapply(res.ls, nrow)
tmp = lapply(tmp, function(x) length(commonGenes)/x)
tmp = t(as.data.frame(tmp))
max(tmp) # 0.9857388 ~ 98.6 %
min(tmp) # 0.9520796 ~ 95.2 %
###################

dir.create("DEGcomparison", showWarnings = F)
setwd("DEGcomparison/")

### BAR PLOT of DEGs ### -------------------------------
ggformat <- function(x){
  df <- as.data.frame(unlist(lapply(x, nrow)))
  colnames(df)[1] <- "x"
  df$Condition <- as.factor(gsub(".+?[.]","",row.names(df)))
  if(length(levels(df$Condition)) == 2){lev = c(2,1)}
  if(length(levels(df$Condition)) == 3){lev = c(2,3,1)}
  if(length(levels(df$Condition)) > 3){lev = c(1:length(levels(df$Condition)))}
  df$Condition <- factor(df$Condition, levels = df$Condition[lev])
  tmp <- gsub("[.].+$","",row.names(df))
  df$Substance <- as.factor(gsub("_.+$","",tmp))
  #colnames(df) <- c("Counts",colnames(df)[2:3])
  df
}
mygg <- function(df, subtitle = ""){
  ggplot(df, aes(x=factor(Substance), y=x, fill = Condition)) + theme_light() + theme(panel.grid.minor = element_blank()) +
    labs(subtitle = subtitle, x = "", y = 'DEG counts') + rotate_x_text(angle = 40) +
    theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) +
    geom_bar(stat="identity", alpha=.8, width = .8, position = position_dodge(.89)) +
    geom_text(aes(x=factor(Substance), y=x+(max(x)*.03), label=x), angle = "0", position = position_dodge(.89)) +
    scale_fill_manual(values = as.character(col.Cond[-1]))
}

pdf(file = "DEG_count.pdf", width = 8, height = 5.6, onefile = T, bg = "transparent", fg ="black")
mygg(ggformat(deg.ls), subtitle = "padj < 0.05 & LFcut") ## pcut & LFcut 
mygg(ggformat(degPcut.ls), subtitle = "padj < 0.05") ## pcut only
dev.off()
########################


### Venn Diagram & ID selection (for heatmap) ### ---------------------------
# plotting function
multiVenn <- function(id.ls, Subs, shape = "ellipse", ...) {
  name <- names(id.ls)[grep(Subs, names(id.ls))]
  # set order & colors
  #pal <- colorRampPalette(c("cadetblue1","deepskyblue1","royalblue4"))(length(name))
  pal <- colorRampPalette(ann_colors$Condition[-1])(length(name))
  
  if(length(name) == 2) {
    name <- name[c(2,1)]
    anno.col <- setNames(c(pal,"firebrick1"), c(name,"Common"))
  }
  if(length(name) == 3) {
    name <- name[c(2,3,1)]
    anno.col <- setNames(c(pal), c(name))
  }
  
  # Venn fun
  set.seed(42)
  venn <- eulerr::euler(id.ls[name], shape = shape, ...)
  s <- round(venn$stress,3)
  e <- round(venn$diagError,3)
  
  # plot
  plot(venn,
       fills = list(fill = anno.col, alpha = .6), #labels = list(col = "black", font = 4),
       legend = list(col = "black", font = 4),
       main = paste0("DEGs [Stress:",s," ; Diag.Er:",e,"]"), #main = paste0(fn,": sign. terms [padj < ",pcut,"]"),
       quantities = TRUE, shape = shape, lty = 1) #lty=0 for transparent
}
intSelect <- function(id.ls, Subs) {
  name <- names(id.ls)[grep(Subs, names(id.ls))]
  Reduce(intersect, id.ls[name])
}
intVenn <- function(int.ls, shape = "ellipse"){
  # Venn fun
  venn <- eulerr::euler(int.ls, shape = "ellipse")
  s <- round(venn$stress,3)
  e <- round(venn$diagError,3)
  # plot
  plot(venn,
       fills = list(fill = col.Subs, alpha = .8), #labels = list(col = "black", font = 4),
       legend = list(col = "black", font = 4),
       main = paste0("Inters.DEGs (LFcut) [Stress:",s," ; Diag.Er:",e,"]"), #main = paste0(fn,": sign. terms [padj < ",pcut,"]"),
       quantities = TRUE, shape = "ellipse", lty = 0)
}

# ---> Needed for heat map row annotation later on!
DEGlfc <- list()
DEG    <- list()
for(k in Substance) {DEGlfc[[k]] <- intSelect(lapply(deg.ls, row.names), k)}
for(k in Substance) {DEG[[k]] <- intSelect(lapply(degPcut.ls, row.names), k)}
univ.DEGlfc <- Reduce(union, DEGlfc) # 232 - with lfcut
univ.DEG <- Reduce(union, DEG) # 1006 - without lfcut

## intersecting DEGs Venn -----------------------------
gg  <- intVenn(DEGlfc) #Lfcut & pcut
gg1 <- intVenn(DEG)    #only pcut
pdf(file = "DEGinters_Venn.pdf", width = 8, height = 8, onefile = T, bg = "transparent", fg ="black")
ggpubr::ggarrange(gg, gg1, ncol=1, nrow=1)
dev.off()

## Venn with pcut & LFcut -----------------------------
id.ls <- lapply(deg.ls, row.names)
gg.ls <- list()
for (k in Substance) {gg.ls[[k]] <-  multiVenn(id.ls, k)}
pdf(file = "DEG_multiVenn_LFcut.pdf", width = 25, height = 14, onefile = T, bg = "transparent", fg ="black")
print(ggpubr::ggarrange(plotlist = gg.ls, ncol=3, nrow=2))
dev.off()

## Venn with pcut & LFcut -----------------------------
id.ls <- lapply(degPcut.ls, row.names)
gg.ls <- list()
for (k in Substance) {gg.ls[[k]] <-  multiVenn(id.ls, k)}
pdf(file = "DEG_multiVenn.pdf", width = 25, height = 14, onefile = T, bg = "transparent", fg ="black")
print(ggpubr::ggarrange(plotlist = gg.ls, ncol=3, nrow=2))
dev.off()

## Create row annotation df ----------------------------
#require(org.Dr.eg.db)
message("\nConnecting to BioMart database ... ")
if(file.exists("~/biomaRt/drerio_mart.Robj")) {
  message("Danio rerio mart object found locally. \nLoading 'drerio_mart.Robj' into R session ...")
  load("~/biomaRt/drerio_mart.Robj")
} else if (file.exists("S:/data/biomaRt/drerio_mart.Robj")){
  message("Danio rerio mart object found in S:/data/biomaRt/ \nLoading 'drerio_mart.Robj' into R session ...")
  load("S:/data/biomaRt/drerio_mart.Robj")
} else {
  message("Could not find 'drerio_mart.Robj'. Creting new mart object from 'www.ensembl.org'",
          "Make sure you have a working interent connection. Otherwise this will fail!\nConnecting to server ...")
  rerio = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="https://www.ensembl.org")
}
message("BioMart Db connected!\n")

mkRowAnnoDf <- function(univ, degList, key = "ENSEMBL", annoPkg="biomart"){
  ## annotate ids with gene SYMBOL ---------------------------------
  # ! ENSDARG00000115405 is falsely annotated with hbbe1.2 !!! => better to use biomaRt for annotation!!!
  if(annoPkg == "biomart") { #stick with biomaRt as it covers more external gene names! Alligns better with ENSEMBL entries
    message("\nAnnotating ENSEMBL IDs with biomaRt ...")
    attr <- c("ensembl_gene_id","external_gene_name","gene_biotype") # list of attributes to extract from biomaRt
    tmp  <- biomaRt::getBM(attributes = attr, mart = rerio, uniqueRows = T,
                 #filters = c("ensembl_gene_id","biotype"), values = list(univ,"protein_coding"), # to extract only protein coding
                 filters = c("ensembl_gene_id"), values = univ)
    colnames(tmp) <- c(key,"SYMBOL","gene_biotype")
  } else {
    message("\nAnnotating ENSEMBL IDs with AnnotationDbi ...")
    AnnoFun <- function(X, key, colInfo){
      x <- as.data.frame(X)
      colnames(x) <- key
      anno <- AnnotationDbi::mapIds(org.Dr.eg.db::org.Dr.eg.db , keys = x[,key], 
                                    keytype = key, column = colInfo, multiVals = "first")
      x[,paste0(colInfo)] <- anno
      x
    }
    tmp <- AnnoFun(univ, key = key, colInfo = "SYMBOL")
  }
  
  ## return T / F vector for matching ids --------------------------
  for(k in names(degList)){
    tmp[,paste0(k)] <- tmp[,1] %in% degList[[k]]
  }
  
  ## replace NAs in SYMBOL column with ENSEMBL gene id if there are missing values -------------
  if(length(which(is.na(tmp$SYMBOL))) > 0){
    #1) Add missing factor levels
    Id <- tmp[which(is.na(tmp$SYMBOL)),1] # ENSEMBL IDs of missing SYMBOLs
    levels <- levels(as.factor(tmp$SYMBOL))
    for (i in c(1:length(Id))){
      levels[length(levels) + 1] <- Id[i]
    }
    tmp$SYMBOL <- factor(tmp$SYMBOL, levels = levels)
    #2) replace NAs in ref$SYMBOL with ENSEMBL ID
    tmp[is.na(tmp$SYMBOL),"SYMBOL"] <- tmp[which(is.na(tmp$SYMBOL)),1]
  } else {message("All ENSEMBL IDs annotated with external gene name")}
  
  ## add Gene IDs from degList which are not element of univ 
  if(length(univ) < length(Reduce(union, degList))){
    #warning("\n'univ' argument contains fewer IDs then need to be annotated in 'degList'!")
    # Annotate degList IDs which are not element of univ with biomart
    attr <- c("ensembl_gene_id","external_gene_name","gene_biotype") # list of attributes to extract from biomaRt
    tmp1  <- biomaRt::getBM(attributes = attr, mart = rerio, uniqueRows = T,
                           filters = c("ensembl_gene_id"), values = setdiff(Reduce(union, degList),univ))
    colnames(tmp1) <- c(key,"SYMBOL","gene_biotype")
    
    for(k in names(degList)){
      tmp1[,paste0(k)] <- "FALSE"
    }
    tmp <- rbind(tmp, tmp1)
    }
  
  # Print output
  rownames(tmp) <- tmp$ENSEMBL
  tmp
}
rowAnLfc <- mkRowAnnoDf(univ.DEGlfc, DEGlfc) #row annotation for Lfc - cutoff and pcut filtered
rowAnPc <- mkRowAnnoDf(univ.DEG, DEG) #row annotation for ONLY pcut filtered 
rm(gg.ls, id.ls, k, gg, gg1)


## ggplot for gene_biotype ##
ggbar <- function(x, subtitle="", xlab = ""){
  # remove rows that only contain FALSE; this is needed as rowALfc was extendet with gene IDs from input list DEG
  rmv <- apply(x[,-c(1:3)],1,function(x){all(x == F)})
  if(any(rmv==T) == T) {X <- droplevels(x[-c(which(rmv==T)),])} else {X <- x}
  # create df for plotting
  df <- as.data.frame(table(X$gene_biotype))
  df$Var1 <- factor(df$Var1, levels = names(sort(table(X$gene_biotype), decreasing = T)))
  S <- sum(df$Freq)
  df$perc <- round((df$Freq/S)*100,2)
  colnames(df)[1] <- "Biotype"
  # plot
  ggplot(df, aes(x = Biotype, y = Freq, fill = Biotype)) + theme_light() +
    theme(panel.grid.minor = element_blank(), legend.position = "none") +
    labs(subtitle = subtitle, x = xlab, y = 'Counts') + rotate_x_text(angle = 25) +
    theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) +
    geom_bar(stat="identity", alpha=.8, width = .8) +
    geom_text(aes(x=Biotype, y=Freq+(max(Freq)*.03), label= paste0(perc," %")), angle = "0")
}

pdf(file = "DEG_biotypes.pdf", width = 8.2, height = 6, onefile = T, bg = "transparent", fg ="black")
ggbar(rowAnLfc, "padj < 0.05 & LFcut")
ggbar(rowAnPc, "padj < 0.05")
dev.off()
#################################################


### TOP SELECTION ### ------------------------------------------------------------------
# Now we are having two major ENSEMBL Gene ID lists, univ.DEG & univ.DEGlfcut
# univ.DEG contains ALL genes that were identified as deferentially regulated
# with a padj <= 0.05 (log2-fc value are apeglm shrunk)
# univ.DEGlfcut contains DEG based on padj cut off AND log2-FC cut off.
# However as there are more common sets for Abamectin and Imidacloprid in the univ.DEG list
# I would suggest to pick the top 25 from each list and run a heatmap with that. 
# The question remains, HOW do we define the top DEGs? By log2fc? padj? baseMean?
# Or: I could run a LDA to separate different MoA groups and see which are most
# powerful to separate the groups from each other (kind of loading score or RandomForest)
# As log2fc values are apeglm shrunk, large log2-fc values with small baseMean might be ignored
# Hence, let's go for log2fc selection now. 

#selecting top N by: log2FoldChange
#df.ls  <- degPcut.ls       # df list with padj, log2-fc & baseMean values of DEGs
#sortBy <- "log2FoldChange" # one of: "baseMean" "log2FoldChange" "padj"
#top    <- 25               # top N genes
#pcut   <- .05 # Not needed in this case, but might be helpful if input ls is not prior filtered

## top selection function ##
topSelection <- function(df.ls, sortBy = "padj", top = 25, pcut = .05) {
  
  # Extract Substance names from list object 
  Subs <- gsub("[.].+$","",names(df.ls))
  Subs <- as.factor(gsub("_.+$","",Subs))
  Subs <- levels(Subs)
  # To store exports in
  top.ls <- list()
  # Merging unction
  myMerge <-  function(df1,df2){merge(df1,df2, by="ENSEMBL" , all = F)}
  
  for(k in Subs) {
    name <- names(df.ls)[grep(k, names(df.ls))]
    tmp <- df.ls[name]
    # Append ENSEMBL ID column to each object in list. Needed for selection later on
    tmp <- lapply(tmp, function(x){
      df <- subset(x, padj <= pcut)
      df$ENSEMBL <- row.names(df)
      df
    })
    # bind tmp list to single df of sortBy value (lfc, padj or baseMean) to compute the mean value and sort by it
    tmp <- Reduce(myMerge, tmp)
    row.names(tmp) <- tmp$ENSEMBL
    # subset tmp to a sortBy df
    tmp <- subset(tmp, select = colnames(tmp)[grep(sortBy, colnames(tmp))])
    # compute row's mean
    tmp <- apply(tmp,1,mean)
    # order by absolute values
    if(sortBy == "padj") {
      # smallest to largest
      X <- sort(tmp, decreasing = F)
    } else {
      # largest to smallest
      X <- sort(abs(tmp), decreasing = T)
    }
    if(length(X) == 0) {stop(paste0("\nNo common DEGs were found for ",k," with padj < ",pcut))}
    if(top > length(X)) {
      warning(paste0("For ",k," 'top' parameter is larger than the set of common DEGs\n'top' ",top," was set to ",length(X)))
      top <- length(X)
    }
    # Final out
    top.ls[[k]] <- names(X)[1:top]
  }
  top.ls
}

## Run topSelection function - works like a charm :)
top <- 15
topLfc = Reduce(union, topSelection(degPcut.ls, sortBy = "log2FoldChange", top = top))
topP   = Reduce(union, topSelection(degPcut.ls, sortBy = "padj", top = top))
topM   = Reduce(union, topSelection(degPcut.ls, sortBy = "baseMean", top = top))
#####################


### HEATMAP with TOP DEG selection ### --------------------------------------
message("\n Start heatmap plotting ...")
## magic center function ##
magiCentFun <- function(mtx, coldat) {
  # this function is designed to center the vst values of each substance / experiment to its
  # respective control. This makes sense as there could be differences in general gene expression
  # values due to differences in age among different experiments.
  
  # Step 1: Subset the input mtx into each experiment
  id.order <- c()
  for(i in levels(coldat$Substance)){
    for(k in levels(coldat$Condition)){
      id <- rownames(coldat[coldat$Condition %in% k & coldat$Substance %in% i,])
      id.order <- append(id.order,id)
    }
  }
  sub.ls <- list()
  for(k in levels(coldat$Substance)){sub.ls[[k]] <- mtx[,which(coldat[id.order,"Substance"] %in% k)]}
  
  # Step 2: Calc the control's mean for each experiment and center the values around it
  centerForNC <- function(mtx, coldat){
    cold <- coldat[colnames(mtx),]
    id   <- row.names(cold[cold$Condition %in% "Control",])
    ctrM <- apply(mtx[,id], 1, FUN = mean) #calcs the mean from the control samples
    mtx-ctrM
  }
  sub.ls <- lapply(sub.ls, function(x){centerForNC(mtx = x, coldat)})
  
  # Step 3: merge mtx back together again & compute SD of each row to scale values
  X <- rlist::list.cbind(sub.ls)
  Sd <- apply(X, 1, FUN = sd) # calcs Sd of each row
  stopifnot(length(Sd) == nrow(X))
  X/Sd # scales for overall Sd of the row
}
vst.cent <- magiCentFun(vst, coldat)

## function to compute colGaps for heatmap ##
colGapsFun <- function(coldat, id.order){
  if(missing(id.order)){
    id.order <- c()
    for(i in levels(coldat$Substance)){
      for(k in levels(coldat$Condition)){
        id <- rownames(coldat[coldat$Condition %in% k & coldat$Substance %in% i,])
        id.order <- append(id.order,id)
      }
    }
  }
  Subs <- levels(coldat$Substance)
  colGaps <- c()
  for(k in Subs[1:(length(Subs)-1)]){
    x <- which(coldat[id.order,'Substance'] %in% k)
    x <- tail(x,1)
    colGaps <- append(colGaps,x)
  }
  colGaps
}

## heatmap function ##
magicHeat <- function(mtx, coldat, rowN = F, colclust = F, rowclust = T, rowdat, symbol = T,
                      clustM = "ward.D2", distM = "euclidean", title = "", ...) {
  # Define Column gaps per Substance 
  colGaps <- colGapsFun(coldat)
  
  # set anno colors & breaks 
  x <- 10 #length of colors above and below zero values
  x <- 2*x + 1
  col <- colorRampPalette(c("mediumblue","white","red2"))(x) #defines color palette in heatmap
  s <- sd(mtx[,rownames(coldat[which(coldat$Condition %in% 'Control'),])])
  m <- mean(mtx[,rownames(coldat[which(coldat$Condition %in% 'Control'),])])
  myBreaks <- c(seq(min(mtx), m-s, length.out=ceiling(x/2)), 
                seq(m+s, max(mtx), length.out=floor(x/2)))
  #scales::show_col(col)
  
  # Prepare row data annotation
  rowData <- rowdat[rownames(mtx),] #cuts rowdat down to the common set and sets order identical to mtx order
  stopifnot(all(rownames(mtx)==rownames(rowData)))
  
  # subset rowData for row annotation & extend ann_colors with respective colors for row annotation
  rowanno <- rowData[,-c(1:3)]
  rowanno <- rowanno[,sort(colnames(rowanno), decreasing = T)] #resort for nicer order in the heatmap
  rowanno[colnames(rowanno)] <- lapply(rowanno[colnames(rowanno)], factor)  ## as.factor() could also be used
  
  anncol <- ann_colors
  for(k in names(ann_colors$Substance)){
    anncol[[k]] <- setNames(c(ann_colors$Substance[[k]],"white"), c("TRUE","FALSE"))
  }
  
  # plot heatmap
  pheatmap::pheatmap(mtx, angle_col = "45", treeheight_col = 30, drop_levels = T, border_color = "gray80",
                     show_rownames = rowN, show_colnames = F, breaks = myBreaks, color = col,
                     #legend = if(symbol==T){F}else{T}, annotation_legend = if(symbol==T){F}else{T},
                     cluster_rows= rowclust,
                     cluster_cols= colclust, gaps_col = colGaps,
                     clustering_distance_rows = distM,
                     clustering_distance_cols = distM,
                     clustering_method = clustM,
                     annotation_col = coldat[,c('Condition','Substance')],
                     #labels_row = rownames(rowData),
                     labels_row = if(symbol == T){rowData$SYMBOL}else{rownames(rowData)},
                     annotation_row = rowanno,
                     annotation_colors = anncol,
                     main = paste0(nrow(mtx)," DEGs ",title," [",distM,", ",clustM,"]")#, ...)
  )
}

# INPUT FOR HEATMAP - mtx object
while (!is.null(dev.list()))  dev.off()
pdf(file = "DEG_heatmaps.pdf", width = 12, height = 24, onefile = T,
    bg = "transparent", fg ="black")

for(k in c("univ.DEG","topLfc","topP","topM")) {
  message(paste0("Heatmap plotting for:\t",k))
  # Subset to only protein coding transcripts!
  x <- intersect(get(k), rowAnPc[rowAnPc$gene_biotype %in% "protein_coding",1])
  mtx <- vst.cent[x,]
  magicHeat(mtx, coldat, title = paste(k,"Prot.coding"), rowN = T, rowdat = rowAnPc)
  magicHeat(mtx, coldat, title = paste(k,"Prot.coding"), rowdat = rowAnPc, rowN = T, colclust = T)
}

for(k in c("univ.DEGlfc")) {
  message(paste0("Heatmap plotting for:\t",k))
  # Subset to only protein coding transcripts!
  x <- intersect(get(k), rowAnPc[rowAnPc$gene_biotype %in% "protein_coding",1])
  mtx <- vst.cent[x,]
  magicHeat(mtx, coldat, title = paste(k,"Prot.coding"), rowN = T, rowdat = rowAnLfc)
  magicHeat(mtx, coldat, title = paste(k,"Prot.coding"), rowdat = rowAnLfc, rowN = T, colclust = T)
  message("Done!")
}
dev.off()
rm(k,x,mtx)


### Additional heatmap for SPECIAL CUSTOM selection! ### -------
run = T
if(run == T){
  # make heatmap with potential biomarker genes overlapping among any other substance + 8 Abamectin DEGs
  # get the custom selection from DEGlfc list object; !Venn only supports up to 5 sets!!!
  # remove Abamectin from list object
  tmp <- DEGlfc[names(DEGlfc)[-1]]
  #tmp <- eulerr::venn(tmp)
  x <- gplots::venn(tmp, show.plot = F)
  x <- attr(x, "intersections")
  x <- x[grep(":",names(x))]
  cust1 <- Reduce(union, x) #custom selection of DEGs without Abamectin
  x[names(DEGlfc)[1]] <- DEGlfc[names(DEGlfc)[1]] # skip this section if you wish to exclude Abamectin
  cust2 <- Reduce(union, x) #custom selection of DEGs with Abamectin
  
  while (!is.null(dev.list()))  dev.off()
  pdf(file = "DEG_heatmaps_CUSTOM.pdf", width = 8.6, height = 6.2, onefile = T,
      bg = "transparent", fg ="black")
  for (k in c("cust2","cust1")){
    x <- intersect(get(k), rowAnPc[rowAnPc$gene_biotype %in% "protein_coding",1])
    mtx <- vst.cent[x,]
    magicHeat(mtx, coldat, title = "Custom selection - Prot.coding", rowdat = rowAnLfc, rowN = T)
    magicHeat(mtx, coldat, title = "Custom selection - Prot.coding", rowdat = rowAnLfc, rowN = T, colclust = T)
  }
  dev.off()
  rm(tmp,x, mtx)
}
rm(run)

### DEG base mean plot - sorted by base mean ### ---------------------------------
ggform <- function(Mtx,select){
  # Mtx: input matrix = normalized gene count matrix
  # select: vector list of ENSEMBL gene IDs which you like to plot
  # use id param to compute mean and sd ONLY from control groups 
  #id <- row.names(coldat[coldat$Condition %in% "Control",]) 
  #mtx <- Mtx[select,id] #params: inputMtx and DEG list
  
  # As i observed that in many cases the compound treatment induces the expression of genes, using only 
  # the control groups mean might underrepresent the actual number of relevant DEGs
  # => run mean and sd calc on ALL gene counts although this will mean a much greater sd 
  # especially for genes like cyp1a which are often an early and heavy responder to env.stress
  mtx <- Mtx[select,]
  # format for ggplot
  m <- cbind(apply(mtx,1,mean),apply(mtx,1,sd))
  colnames(m) <- c("m","sd")
  m <- merge(as.data.frame(m),rowAnPc[,1:3], by=0)[,-1]
  m <- dplyr::arrange(m, desc(m))
  #m$SYMBOL <- factor(m$SYMBOL, levels = unique(m$SYMBOL))
  m
}
myggbar <- function(m){
  ggplot(m, aes(x=reorder(SYMBOL, -m), y=m)) + ylab("mean normalized counts") + theme_light() +  xlab("") + rotate_x_text(angle = 40) +
    theme(panel.grid.minor = element_blank()) + ggtitle(paste0("Custom prot.coding DEG set: ",nrow(m))) +
    geom_errorbar(aes(x=reorder(SYMBOL, -m), ymin=m, ymax=m+sd), colour="black", alpha=0.8, size=.3, width=.4) +
    geom_bar(stat="identity", fill="skyblue", alpha=1) +
    geom_hline(yintercept = lim, linetype = 2, color="firebrick")
}
myggbar2 <- function(m){ # with different title type for base mean cut off 
  ggplot(m, aes(x=reorder(SYMBOL, -m), y=m)) + ylab("mean normalized counts") + theme_light() +  xlab("") + rotate_x_text(angle = 40) +
    theme(panel.grid.minor = element_blank()) + ggtitle(paste0("Custom prot.coding DEG set: ",nrow(m)," - min. mean count: ",lim)) +
    geom_errorbar(aes(x=reorder(SYMBOL, -m), ymin=m, ymax=m+sd), colour="black", alpha=0.8, size=.3, width=.4) +
    geom_bar(stat="identity", fill="skyblue", alpha=1, width=.5) +
    geom_hline(yintercept = lim, linetype = 2, color="firebrick")
}


## Plot mean counts of DEG genes
lim = 600 #minimum mean count number
# Small custom set of overlaps
pdf(file = "DEG_controlsMeanCount_CUSTOM.pdf", width = 7, height = 5, onefile = T, bg = "transparent", fg ="black")
myggbar(ggform(normMtx, intersect(cust2, rowAnPc[rowAnPc$gene_biotype %in% "protein_coding",1])))
dev.off()

# ALL DEGs from overlapping conditions
pdf(file = "DEG_controlsMeanCount.pdf", width = 12, height = 7, onefile = T, bg = "transparent", fg ="black")
m <- ggform(normMtx, intersect(univ.DEGlfc, rowAnPc[rowAnPc$gene_biotype %in% "protein_coding",1]))
myggbar(m)
# Now with mean count cut off (lim parameter)
for(lim in c(200,400,500,600)){
  print(myggbar2(subset(m, m > lim)))
}
dev.off()
customGeneSet <- m

# Compute the quantile% for a particular baseMean cutoff
ecdf_fun <- function(x,val) ecdf(x)(val) #x = vector of values; val = value for which to compute the %quantile for
ecdf_fun(apply(normMtx[commonGenes, ],1,mean) ,500) # top 67% 
ecdf_fun(apply(normMtx[customGeneSet$ENSEMBL,],1,mean) ,500) #60.6%

### Create heatmap with custom gene set after baseMean cutoff (protein coding only!!!) ### -------------------
baseMeanQFun = function(x, Q = 65, title = "", ylab = "baseMean"){
  if(Q > 85){Q = 85} # top-% quantile
  q = quantile(x, seq(0, 1, 0.01))[Q+1]
  plot(quantile(x, seq(0, 1, 0.01)), main = paste("Average baseMean gene count quantile distr.",title),
       ylab = ylab, xlab = "%-Quantile")
  abline(v = Q, lwd=1.5, lty=2, col = "firebrick")
  abline(h = q, lwd=1.5, lty=2, col = "firebrick")
  legend(x="topleft", text.col = "black", bty = "n",
         legend = paste0(Q-15,"% quantile = ",round(quantile(x, (Q-15)*.01),1),"\n",
                         Q,"% quantile = ",round(q,1),"\n",
                         Q+15,"% quantile = ",round(quantile(x, (Q+15)*.01),1),"\n"))
}
while (!is.null(dev.list()))  dev.off()
pdf("DEG_baseMean_quantile_distr.pdf", width = 8 , height = 7, onefile = T, bg = "transparent")
baseMeanQFun(m$m, 65)
baseMeanQFun(log2(m$m), 65, title = "log2(baseMean)", ylab = "log2(baseMean)")
dev.off()

m = customGeneSet
while (!is.null(dev.list()))  dev.off()
pdf(file = "DEG_heatmaps_CUSTOM_baseMeanCut.pdf", width = 9.8, height = 10.8, onefile = T, bg = "transparent")
for(lim in c(500,600,1000)){# baseMean cutoff for genes in heatmap
  mtx <- vst.cent[subset(m, m > lim)[,"ENSEMBL"],]
  #magicHeat(mtx, coldat, title = paste("Cust.select - Prot.coding baseMeanCut:",lim), rowdat = rowAnLfc, rowN = T, symbol = F)
  magicHeat(mtx, coldat, title = paste("- Prot.coding - baseMeanCut:",lim), rowdat = rowAnLfc, rowN = T, symbol = T)
  magicHeat(mtx, coldat, title = paste("- Prot.coding - baseMeanCut:",lim), rowdat = rowAnLfc, rowN = F, symbol = T) #to get borders
  magicHeat(mtx, coldat, title = paste("- Prot.coding - baseMeanCut:",lim), rowdat = rowAnLfc, rowN = T, colclust = T)
  magicHeat(mtx, coldat, title = paste("- Prot.coding - baseMeanCut:",lim), rowdat = rowAnLfc, rowN = F, colclust = T)
}
dev.off()

rm(m,lim,mtx)
######################################


## Session Information ------------
setwd(paste0(home,out))
sink(paste0("SessionInfo_compareDESeq2.txt"))
print(date())
print(devtools::session_info())
sink()

## Clear Env and save Rdata -------------
rm(dds, DEG, DEGlfc, res.ls, reslfs.ls,rerio, id, k, N)
message(paste0("\nSaving RData to:\t",home,out,"\nThis might take a while ..."))
save.image(paste0(home,out,"/compareDESeq2.RData"))
message("Done!")
setwd(home)
message("\nPuh! What a ride! All Done!\nJ.A.R.V.I.S over and out! :)")
####  END OF SCRIPT  ####