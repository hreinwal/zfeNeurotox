# Author:   Hannes Reinwald
# Contact:  hannes.reinwald@ime.fraunhofer.de
# Type:     R script

### README: -------------------------------------------------------
# This script is designed to summarize and compare the GSEA results obtained from clusterProfiler.
# Input format are the generated csv output files located in:
# ./clusterProfiler/GSEA/gse[GO][KEGG][Reactome]
# The idea is to first import these tables for a specified set of substances/projects.

# ! Copy ALL the csv files from different outputs you wish to compare into a new dir !
# This should be the structure:
#
#     compareGSEA/
#               gseKEGG/
#               gseReactome/
#               gseGO/
#                     BP/
#                     MF/
#                     CC/
# 
# The section below will automatically import the respective files for RNAseq experiments.
# (Located within the DESeq2 output folder) If you have a different dir or file structure please
# create the compareGSEA dir structure yourself. 
###############

### PACKAGES ----------------------------------------------------------------------------------------
require(dplyr)
require(data.table)
require(ggpubr)
require(ggplot2)
require(ggpmisc)
require(corrplot)
require(RColorBrewer)
require(R.utils)
require(plyr)
require(viridis)
################

### Output environment for compareGSEA ###
home = "S:/manuscripts/neurotox/data/"
dir.create(paste0(home,"/compareGSEA"), showWarnings = F)
setwd(paste0(home,"/compareGSEA"))


### Data Import from Transcriptomics -----------------------------------------------------------------
# Specify experiments to compare
pathToData = "S:/data/RNA_seq_gene_counts/SUBSTANCE/DESeq2_Pairwise/clusterProfiler/GSEA"
#dir("S:/data/RNA_seq_gene_counts/") #to see all the Substances available
#dataSet=c("Abamectin","Carbaryl","Chlorpyrifos_2","Fipronil","Imidacloprid","Methoxychlor") #FULL SET
dataSet=c("Abamectin","Chlorpyrifos_2","Fipronil_2","Imidacloprid","Carbaryl","Methoxychlor") #for now, only this

# empty list objects to store imported data in
kegg.ls <- list() # gseKEGG results
reac.ls <- list() # gseReactome results
goBP.ls <- list() # gseGO - biol. processes results
goMF.ls <- list() # gseGO - molec. function results
goCC.ls <- list() # gseGO - cel. component results
gseaLists = objects(pattern = "[gcPFC].ls")

# loop to navigate into project folder and import respective csv table
for(subs in dataSet){
  for(gse in list.dirs(path = paste0(gsub("SUBSTANCE",subs,pathToData)),full.names = F, recursive = F)){
    if(gse == "gseGO"){
      # import from subdir
      for(godir in c("BP","CC","MF")){
        # Importing gseGO results
        message(paste(subs,"- Importing",gse,godir,"results"))
        # get file paths to subdir
        files = list.files(path = paste0(gsub("SUBSTANCE",subs,pathToData),"/",gse,"/",godir), 
                           recursive = F, full.names = T, pattern = paste0("*_",gse,".",godir,".csv"))
        for(f in files){
          name = gsub(paste0("^.*",gse,"/",godir,"/"),"",f)
          name = gsub("reslfs_","",name)
          name = gsub(".csv","",name)
          # import files to respective GO list object
          if(godir == "BP"){
            # import BP
            goBP.ls[[name]] = read.csv2(f, header = T)[,-1]
          }else if(godir == "CC"){
            # import CC
            goCC.ls[[name]] = read.csv2(f, header = T)[,-1]
          }else if(godir == "MF"){
            # import MF
            goMF.ls[[name]] = read.csv2(f, header = T)[,-1]
          }
        }
      }
    }else{
      # Importing gseKEGG / gseReactome results
      message(paste(subs,"- Importing",gse,"results"))
      files = list.files(path = paste0(gsub("SUBSTANCE",subs,pathToData),"/",gse), 
                         recursive = F, full.names = T, pattern = paste0("*_",gse,".csv"))
      for(f in files){
        name = gsub(paste0("^.*",gse,"/"),"",f)
        name = gsub("reslfs_","",name)
        name = gsub(".csv","",name)
        if(gse == "gseKEGG"){
          # gseKEGG import
          kegg.ls[[name]] = read.csv2(f, header = T)[,-1]
        }else{
          # gseReactome import
          reac.ls[[name]] = read.csv2(f, header = T)[,-1]
        }
      }
    }
  }
}
rm(files,gse,name,subs,f,godir)


# To organize some of the downstream plots a little bit nicer, let's rename the sample feature Low Mid High to 
# treat1, treat2, treat3 (treat = treatment level) to make downstream plotting nicer with alphabetically sorting working: 
for(n in gseaLists){
  message(paste("\nRenaming list objects in:\t",n))
  N = gsub("_[Ll]ow", "_treat1",names(get(n)))
  N = gsub("_[Mm]id", "_treat2",N)
  N = gsub("_[Hh]igh","_treat3",N)
  x = get(n)
  names(x) = N
  message(paste("Sorting list objects in:\t",n))
  assign(n,x[sort(names(x))])
}
rm(N,n,x)


### Define UNIVERSE ---------
# = common set of possible comparable terms / pathways
# Compute universe parameters for fisher's exact t-test (Background of possible common terms/pathways)
for(fn in gseaLists) {
  ls <- get(fn)
  tmp <- lapply(ls, function(x){as.factor(x$ID)})
  assign(paste0(gsub(".ls",".Univ",fn)), Reduce(intersect, tmp))
}
rm(tmp,ls,fn)
########################################


### Merge imported GSEA results -----------------------------------------------------------------------
# The enrichment score depends on gene set size. Consequently, unless we are in an unlikely scenario
# where all the gene sets we test are of equal size, the enrichment scores will not be identically 
# distributed and cannot be directly compared. This precludes the calculation of an empirical false
# discovery rate. GSEA solves this problem by applying a transformation to calculated 
# enrichment scores (ES) such that they lie on a comparable scale.
# This normalized enrichment score (NES) is the ES divided by the expected value (i.e. average) of 
# the corresponding null distribution. Hence, the NES is the relevant factor to compare with.
# See: https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/

# Merged datasets are reduced to the common set of possible observable terms / pathways (Universe)

## Multiple file integration function
fileIntegration = function(ls, filter.by, pcut = .05, universe){
  # shorten / filter data
  short.ls <- list() #object list to store shortened and column name extended dataframes in
  for(k in names(ls)){
    #ls = reac.ls
    #k = names(ls)[1]
    x = ls[[k]][,c("ID","Description","NES","pvalue","p.adjust")]
    if(nrow(x) > 0){
          colnames(x) = c("ID","Description","NES","p","padj")
          if(missing(filter.by)){filter.by = 'None'}
          
          # filter signif. results by p or padj value
          if(filter.by == 'p' | filter.by == 'padj'){
            message(paste(filter.by,"<=",pcut,"applied on:\t",k))
            if(filter.by == 'p'){x = subset(x, p <= pcut)}else{x = subset(x, padj <= pcut)}
          }else{message(paste("No pval/padj filter applied on:\t",k))}
          if(nrow(x) == 0){warning(paste("\nNo sign. enriched terms found in",k,"for",filter.by,"<=",pcut))}
          # colname extension
          ext = paste(colnames(x), gsub("_gse.*","",k), sep = ".")
          ext[1] = "ID" # set first character to common ID name
          ext[2] = "Description"
          colnames(x) = ext
          
          # Set 'ID' and 'Description' as factor!
          x$ID = as.factor(x$ID)
          x$Description = as.factor(x$Description)
          short.ls[[k]] = x
    }else{
      warning(paste0("\n",k," is an empty dataframe!"))
    }
  }
  # merge data frames in short.ls into wide df -> for heatmap
  myMerge = function(df1,df2){
    merge(df1,df2, by=c('ID','Description') , all = T)
  }
  res <- Reduce(myMerge, short.ls)
  res[which(res$ID %in% universe),] #filter by universe
}

## Run fileIntegration in loop for the *.ls objects
for(k in gseaLists){
  x = get(k)
  univ = get(gsub(".ls",".Univ",k))
  n = gsub(".ls","",k)
  message(paste0("\nIntegrating data from:\t",k," \t---> \t",paste0(n,".df")))
  assign(paste0(n,".df"), fileIntegration(x, filter.by = 'padj', universe = univ))
  message(paste0("\nIntegrating data from:\t",k," \t---> \t",paste0(n,".ALL")))
  assign(paste0(n,".ALL"), fileIntegration(x, universe = univ))
}

## Export wide df
for(k in objects(pattern = "[gcPFC].ALL")){
  if(k == 'reac.ALL'){n = 'reactome'}else{n = gsub('.ALL','',k)}
  message(paste0("\nSaving merged data tables to csv for:\t",n))
  dir.create(n, showWarnings = F)
  write.csv2(get(k), file = paste0(n,"/",n,".gse_mergedRes.csv"), row.names = F)
  message("Done!")
}


# Ignore this section for now ...
# Thought I might need them later but I didn't. But keep code section for any case
skip = T
if(skip == F){
  ## Merge multiple files from ls objects to long df 
mkLongDf = function(ls){
  # merge ls into long df -> ggplot2 / correlation plots
  tmp.ls <- list()
  for(k in names(ls)){
    x = ls[[k]][,c("ID","Description","NES","pvalue","p.adjust")]
    n = gsub(".ls",".df",k)
    colnames(x) = c("ID","Description","NES","p","padj")
    if(nrow(x) > 0){
      n = gsub("_gse.*","",n)
      x$Type = as.factor(c(n))
      x$ID = as.factor(x$ID)
      x$Description = as.factor(x$Description)
      tmp.ls[[n]] = x
    }else{warning(paste0(n," is an empty dataframe!"))}
  }
  rbindlist(tmp.ls)
}
## Run mkLongDf in loop for the *.ls objects to create long df for ggplot
for(k in gseaLists){
  x = get(k)
  n = gsub(".ls",".ldf",k) #ldf = long data frame
  message(paste0("Integrating data from:\t",k," \t---> \t",n))
  assign(n, mkLongDf(x))
}
}
###

rm(x,n,k,skip, univ)
###################################


### PLOTTING ### #####################################################################

# Explain what happens here ...
# QC
# Correlation plots & Correlation analysis
# Dist Mtx
# Heatmap


### QC Plots ### --------------------------------------------------------------------------------------------
message("\nQC Plots ... ")

## overall sign. enriched terms (barplot)
myBarplot <- function(n){ # n = character vector!!!
  if(n == "reac.df"){fn = "reactome"}else{fn = gsub(".df","",n)} #set title for out file name (fn)
  x = get(n)
  df = x[,grep("NES.",colnames(x))]
  row.names(df) <- x$ID
  res = apply(df,2,function(x){length(na.omit(x))}) # counts N of sign. enriched terms
  res = setNames(res,gsub("NES.","",names(res))) ## shorten row and col names of M
  
  # widen the margin of plot (b,l,t,r) DEFAULT: 5 4 4 1
  par(mar = c(8,4,4,1)+.1)
  res <- res[sort(names(res))] #sort alphabetically
  xx = barplot(res, las = 3, #rotating sample labels 90
               #xlab = "Treatment",
               ylim = c(0,1.2*max(res)),
               ylab = paste0("Sign. enriched terms: ",fn),
               main = fn
  )
  text(x=xx, y=res, labels = res, pos = 3, cex = 1, col = "blue", srt=0)
}

while (!is.null(dev.list()))  dev.off()
pdf(file = "SignEnrichedTerms.pdf", width = 12.5, height = 7, onefile = T, bg = "transparent")
par(mfrow=c(2,3))
obj = objects(pattern='*[.]df')#[c(1,4,2,5,3)]
for(n in obj){myBarplot(n)}

par(mfrow=c(2,3))
obj = objects(pattern='*[.]ALL')#[c(1,4,2,5,3)]
for(n in obj){myBarplot(n)}
dev.off()

## pval & NES distr
for(n in gseaLists){
  # QC plotting for a ls object containing clusterProfiler results.
  if(n == "reac.ls"){fn = "reactome"}else{fn = gsub(".ls","",n)} #set title for out file name (fn)
  message(paste0("\nStart QC plotting for:\t", fn," results"))
  # Output for QC plots
  dir.create(fn, showWarnings = F)
  pdf(file = paste0(fn,"/",fn,"_pval_NES_distr.pdf"), #print directly to pdf
      width = 7.8, height = 7, onefile = T,
      bg = "transparent", #Background color
      fg ="black"         #Foreground color
  )
  ls = get(n)
  for(k in names(ls)){
    df = ls[[k]]
    if(nrow(df) < 1){
      warning(paste(k, "is an empty data frame. No QC plots created!"))
    }else{
      par(mfrow=c(2,2))
      # pval hist
      hist(df$pvalue,12, xlab = "p value", main = k)
      # pval vs padj
      plot(df$pvalue,df$p.adjust,
           xlab = "p value", ylab = "BH corrected p value (padj)", main = k)
      pval = .05 
      padj = .05
      abline(h=padj, lwd=1, lty=2)
      abline(v=pval, lwd=1, lty=2)
      legend("bottomright", inset = c(.02,0.15),
             text.col = "springgreen3",
             bty = "n",
             legend = paste("p <",pval,":",sum(df$pvalue <= pval, na.rm = T))
      )
      legend("bottomright", inset = .02,
             text.col = "dodgerblue",
             bty = "n",
             legend = paste("padj <",padj,":",sum(df$p.adjust <= padj, na.rm = T))
      )
      # ES hist
      hist(df$enrichmentScore,50, xlab = "Enrichment score (ES)", main = k)
      # NES hist
      hist(df$NES,50, xlab = "Normalized enrichment score (NES)", main = k)
    }
  }
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  message("Done!")
  rm(ls,fn,n,k,df,pval,padj)
}
################


### Treatmen QCR (Quadrant count ratio) plot ###
### Functions ### -----------------------------
## QC-ratio matrix function ##
# requires wide *.df as input, which is prefiltered for padj <= 0.05; can deal with NAs
qcrMtx  <- function(mat, rmNA) {
  if(missing(rmNA)){rmNA = T}
  
  # read in data
  mat <- as.matrix(mat)
  n <- ncol(mat)
  
  # Create output mtx
  matQ <- matN <- matrix(NA, n, n)
  colnames(matQ) <- row.names(matQ) <- colnames(matN) <- row.names(matN) <- colnames(mat)
  
  # get diag values for n.mat (Number of common sign. enriched terms)
  nVal <- c()
  for (i in 1:n) {nVal <- append(x = nVal, values = length(na.omit(mat[,i])))}
  diag(matN) <- nVal
  diag(matQ) <- 1
  
  ### quadrant count function ###
  quadrantCount <- function(x,y){
    df <- cbind(x,y)
    quadrant <- apply(df, 1, function(x){
      if(x[1]>0 & x[2]>0){1}else if(x[1]>0 & x[2]<0){2}else if(x[1]<0 & x[2]<0){3}else if(x[1]<0 & x[2]>0){4} 
    })
    as.data.frame(table(quadrant))
  }
  ### qcr computation function ###
  QCR <- function(x,y){
    res <- quadrantCount(x,y)
    # Compute QCR (Quadrant count ratio)
    UP = sum(res[which(res$quadrant %in% c(1,3)),"Freq"])
    DOWN = sum(res[which(res$quadrant %in% c(2,4)),"Freq"])
    ALL = UP+DOWN
    qcr = (UP/ALL)-(DOWN/ALL)
    qcr
  }
  
  # run loop
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp = cbind(mat[,i],mat[,j])
      tmp = as.data.frame(na.omit(tmp))
      if(nrow(tmp > 0)) {
        qcr <- QCR(x = tmp[, 1], y = tmp[, 2])
        matQ[i, j] <- matQ[j, i] <- qcr
        matN[i, j] <- matN[j, i] <- nrow(tmp)
      } else {
        warning(paste("\nNo common observations found for:",colnames(mat)[i],"vs",colnames(mat)[j]))
        matQ[i, j] <- matQ[j, i] <- NA
        matN[i, j] <- matN[j, i] <- nrow(tmp)
      }
    }
  }
  
  if(rmNA == T){
    matQ[is.na(matQ)] <- 0
    list(qcr = matQ, counts = matN)
    # Not quite what I wanted ... 
    # check how many columns have no common sets
    # 1) Create tmp matrix to store one side of the mat object in
    #tmp <- matrix(1,ncol(matQ),ncol(matQ))
    #for (i in 1:(n - 1)) {
    #  for (j in (i + 1):n) {
    #    tmp[j,i] <- matQ[i,j]
    #  }
    #}
    # Check which rows are not complete, those are to be removed
    #rm <- which(complete.cases(tmp) != T)
    #stopifnot(nrow(matQ[-rm,-rm]) == ncol(matQ[-rm,-rm]) & nrow(matN[-rm,-rm]) == ncol(matN[-rm,-rm]))
    #stopifnot(nrow(matQ[-rm,-rm]) == nrow(matN[-rm,-rm])) # Check!
    #list(qcr = matQ[-rm,-rm], counts = matN[-rm,-rm])
    
  } else {
    list(qcr = matQ, counts = matN)
  }
}

### Fisher"s exact t-test ###
conTable <- function(x,m,k,N,rowCat,colCat){
  # x = overlap between m & k
  # m = set size of m (fixed marginal value)
  # k = set size of k (fixed marginal value)
  # N = total set size background (universe); fixed marginal value
  
  # Make sure input is correctly formated
  if(missing(x) | missing(m) | missing(k) | missing(N)){stop("x,m,k & N must be provided!") }
  if(N < k | N < m | N < x){stop("N (universe) must be largest integer!")}
  if(x > k | x > m){stop("x can not be greater than k or m !")}
  if(((m+k-x) > N) == T){stop("N can not be smaller than ( m + k - x ) !")}
  
  # Check if all inputs are integers
  a <- c(x,m,k,N)
  stopifnot(length(a) == 4)
  if(all.equal(a, as.integer(a)) != T){stop("All input values must be integers!")}
  
  if(missing(rowCat)){rowCat=c("In m","Not m")}
  if(missing(colCat)){colCat=c("In k","Not k")}
  
  res <- as.table(matrix( c(x, m-x, k-x, N-m-k+x), nrow = 2, 
                          dimnames = list(rowCat,
                                          colCat)))
  stopifnot(sum(as.data.frame(res)$Freq) == N)
  res
}
fisher.mtest <- function(mat, N, alternative, padjM) {
  # mat:    is an equal dim. count matrix (cont.table) displaying the count of common sign. enriched terms
  #         observed between the variables. Diagonals = full set of terms in variable.
  # N:      Universe; entire background of possible detectable GO / Pathway terms
  #         Equal to the number of all described terms.
  # for reactome ~ 814
  # for kegg     ~ 153
  # for GO BP    ~ 866
  #        CC    ~ 3148
  #        MF    ~ 514
  #
  # padjM:  one of p.adjust(method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  #         default is "bonferroni"
  
  n <- nrow(mat) # dimensions for output mtx
  d <- diag(mat) # diagonal values = total sets
  stopifnot(n == length(d))
  matP <- matConf1 <- matConf2 <- matEst <- matrix(0, n, n) # output mtx
  
  ### contigency table function (N = universe) ###
  # Creates contigency table from the following input:
  conTable <- function(x,m,k,N,rowCat,colCat){
    # x = overlap between m & k
    # m = set size of m (fixed marginal value)
    # k = set size of k (fixed marginal value)
    # N = total set size background (universe); fixed marginal value
    
    # Make sure input is correctly formated
    if(missing(x) | missing(m) | missing(k) | missing(N)){stop("x,m,k & N must be provided!") }
    if(N < k | N < m | N < x){stop("N (universe) must be largest integer!")}
    if(x > k | x > m){stop("x can not be greater than k or m !")}
    if(((m+k-x) > N) == T){stop("N can not be smaller than ( m + k - x ) !")}
    
    # Check if all inputs are integers
    a <- c(x,m,k,N)
    stopifnot(length(a) == 4)
    if(all.equal(a, as.integer(a)) != T){stop("All input values must be integers!")}
    
    if(missing(rowCat)){rowCat=c("In m","Not m")}
    if(missing(colCat)){colCat=c("In k","Not k")}
    
    res <- as.table(matrix( c(x, m-x, k-x, N-m-k+x), nrow = 2, 
                            dimnames = list(rowCat,
                                            colCat)))
    stopifnot(sum(as.data.frame(res)$Freq) == N)
    res
  }
  
  # Running Fisher"s Exact Test for Count Data
  if(missing(alternative)){alternative = "greater"}
  for(i in 1:(n-1)){
    m = as.integer(d[i])
    for(j in (i+1):n){
      #for m set by i, respective k set by j
      k = as.integer(d[j])
      x = mat[j,i]
      ft = fisher.test(conTable(x, m, k, N), alternative = alternative)
      matP[i, j] <- matP[j, i] <- ft$p.value
      matEst[i, j] <- matEst[j, i] <- ft$estimate
      matConf1[i, j] <- matConf1[j, i] <- ft$conf.int[1]
      matConf2[i, j] <- matConf2[j, i] <- ft$conf.int[2]
    }
  }
  
  ## pvalue correction ##
  if(missing(padjM)){padjM = "bonferroni"}
  if(padjM != "none"){
    matPadj <- matrix(0,ncol(mat),ncol(mat))
    
    # Extract pval from matrix
    pval <- c()
    for(i in 1:(ncol(mat)-1)){
      for(j in (i+1):ncol(mat)){
        pval <- append(pval,matP[j,i])
      }
    }
    
    # Run pval correction
    padj <- p.adjust(pval, method = padjM)
    
    # Sort corrected pval back in matrix
    tmp <- padj
    for(i in 1:(ncol(mat)-1)){
      for(j in (i+1):ncol(mat)){
        matPadj[i,j] <- matPadj[j,i] <- tmp[1]
        tmp <- tmp[-1]
      }
    }
    # Bind results in list WITH corrected pvalues
    fisherMtx <- list(p.adjust = matPadj,
                      p.value = matP, 
                      estimate = matEst,
                      lowCI = matConf1, 
                      uppCI = matConf2)
  } else {
    # Bind results in list without corrected pvalues
    fisherMtx <- list(p.value = matP, 
                      estimate = matEst,
                      lowCI = matConf1, 
                      uppCI = matConf2)
  }
  
  # Append row and col names
  lapply(fisherMtx, function(x){
    row.names(x) <- colnames(mat)
    colnames(x) <- colnames(mat)
    x
  })
  # Done! :) #
}
#ft <- fisher.mtest(qcrMtx(mtx, rmNA = T)$counts, 153, padjM = "bonferroni")
#round(ft.res$p.adjust,3)

### Plotting function based on corrplot ###
QCRplot <- function(fn, plotType, method, order, clust, sig.level, rmNA, diag = T, universe, padjM) {
  
  # fn = "kegg.df" <- df object name
  
  # if you wish to compute fisher"s exact t-test for overlap you must provide the 'universe' parameter!!!
  # for reactome ~ 814
  # for kegg     ~ 153
  # for GO BP    ~ 866
  #        CC    ~ 3148
  #        MF    ~ 514
  if(missing(plotType)){plotType = "qcr"} #or: "counts" , "rmvNonSig"
  if(plotType == "rmvNonSig" & rmNA == F) {stop()}
  
  # set title for out file name (n)
  if(fn == "reac.df"){n = "reactome"}else{n = gsub(".df","",fn)}
  
  # read data
  x <- get(fn)
  df <- x[,grep("NES.",colnames(x))] #reduce input df to NES columns
  row.names(df) <- x$ID
  colnames(df) <- gsub("NES.","",colnames(df))
  
  # corrplot params
  if(missing(method)){method = "circle"}
  if(missing(order)){
    order = "alphabet"
    clust = NULL
  }
  
  if(order != "hclust") {
    # QCR matrix with NAs 
    M <- qcrMtx(df, rmNA = rmNA)
    title = paste0(n," NES qc-ratios for common sign. enriched terms")
  }
  
  if(order == "hclust") {
    # QCR matrix without NAs 
    M <- qcrMtx(df, rmNA = T)
    if(missing(clust)){clust = "average"}
    title = paste0(n," NES qc-ratios ",clust," clust")
  }
  
  ### Corrplot ###
  
  # widen the margin of plot (b,l,t,r) DEFAULT: 5 4 4 1
  par(mar = c(5,4,6,1)+.1)
  
  # set colors 
  col=rev(brewer.pal(n=8, name="RdYlBu"))
  
  # plot
  if(plotType == "qcr") {
    if(missing(sig.level)){sig.level = -1}
    corrplot(M$qcr, type = "upper", tl.col="black", tl.srt=35, diag = diag, addgrid.col = "grey90", na.label = " ",
             col = col, cl.align.text = "l",
             method = method,
             order = order,
             hclust.method = clust,
             title = title,
             p.mat = M$counts, sig.level = sig.level, insig = "p-value", #add counts
             mar=c(0,0,2,0)
    )
  }
  
  if(plotType == "rmvNonSig" & rmNA == T){ #this works only for matrix input without NAs!!!
    if(missing(sig.level)){sig.level = .05}
    if(missing(padjM)){padjM = "bonferroni"}
    
    ft <- fisher.mtest(M$counts, N = universe, padjM = padjM)
    corrplot(M$qcr, type = "upper", tl.col="black", tl.srt=35, diag = diag, addgrid.col = "grey90", na.label = " ",
             col = col, cl.align.text = "l",
             method = method,
             order = order,
             hclust.method = clust,
             title = paste(title,padjM,"pcut:",sig.level),
             p.mat = ft[[1]], sig.level = sig.level, insig = "blank", #remove circles below sig.level
             mar=c(0,0,2,0)
    )
  }
  
  if(plotType == "counts") {
    vir <- viridis(14)[8:14]
    col <- c()
    col <- append(col,c(rev(vir),rev(vir)))
    
    mtx <- log2(M$counts)+0.1 #to make values visible that are 0 (log2(1)=0)
    mtx[mtx == -Inf] <- 0
    corrplot(mtx, type = "upper", tl.col="black", tl.srt=35, diag = diag, addgrid.col = "grey90", na.label = "nc",
             col = col,cl.align.text = "l",
             method = method,
             order = order,
             hclust.method = clust,
             title = gsub("NES qc-ratios","Common sign counts",title),
             p.mat = M$qcr, sig.level = -1.1, insig = "p-value", #add qcr
             #p.mat = M$counts, sig.level = sig.level, insig = "p-value", #add counts
             is.corr = FALSE, cl.lim = c(0,ceiling(max(mtx, na.rm = T))),
             mar=c(0,0,2,0)
    )
  }
  
  if(plotType != "counts" & plotType != "qcr" & plotType != "rmvNonSig"){
    stop("\nDefined plotType unknown! plotType must be 'qcr', 'counts' or 'rmvNonSig'")
  }
}
################# 

message("\nTreatment QCR - Comparative plot ... ")
obj = objects(pattern='*[.]df')[c(4,5,1,2,3)]# resort for nicer arranged output

while (!is.null(dev.list()))  dev.off()
pdf(file = "NES_QCR_sign.pdf", width = 14, height = 20, onefile = T, bg = "transparent")
par(mfrow=c(5,2))
for(fn in obj) {
  univ <- get(gsub("[.]df",".Univ",fn))
  QCRplot(fn, rmNA = T, order="alphabet", clust = NULL, method = 'circle')
  QCRplot(fn, rmNA = T, order="alphabet", clust = NULL, method = 'circle', 
          plotType = "rmvNonSig", universe = length(univ), padjM = 'BH')
}
for(k in c("average","ward.D2")) {
  for(fn in obj) {
    univ <- get(gsub("[.]df",".Univ",fn))
    QCRplot(fn, rmNA = T, order="hclust", clust = k, method = 'circle')
    QCRplot(fn, rmNA = T, order="hclust", clust = k, method = 'circle', 
            plotType = "rmvNonSig", universe = length(univ), padjM = 'BH')
  }
}
dev.off()

while (!is.null(dev.list()))  dev.off()
pdf(file = "NES_QCR_sign_NOzero.pdf", width = 14, height = 20, onefile = T, bg = "transparent")
par(mfrow=c(5,2))
for(fn in obj) {
  univ <- get(gsub("[.]df",".Univ",fn))
  # create tmp list file in which empty df are removed
  tmp = get(fn)
  tmp = tmp[ ,grep("NES.",colnames(tmp))]
  keep = apply(tmp, 2, function(x){all(is.na(x))})
  tmp = tmp[,keep == F]
  # plot
  QCRplot("tmp", rmNA = T, order="alphabet", clust = NULL, method = 'circle')
  QCRplot("tmp", rmNA = T, order="alphabet", clust = NULL, method = 'circle', 
          plotType = "rmvNonSig", universe = length(univ), padjM = 'BH')
}
for(k in c("average","ward.D2")) {
  for(fn in obj) {
    univ <- get(gsub("[.]df",".Univ",fn))
    # create tmp list file in which empty df are removed
    tmp = get(fn)
    tmp = tmp[ ,grep("NES.",colnames(tmp))]
    keep = apply(tmp, 2, function(x){all(is.na(x))})
    tmp = tmp[,keep == F]
    # plot
    QCRplot("tmp", rmNA = T, order="hclust", clust = k, method = 'circle')
    QCRplot("tmp", rmNA = T, order="hclust", clust = k, method = 'circle', 
            plotType = "rmvNonSig", universe = length(univ), padjM = 'BH')
  }
}
dev.off()

message("Done!")
rm(keep,tmp,k,univ,fn)
################################################


### Treatment Correlation Plots -----------------------------------------------------------------------------------
# https://towardsdatascience.com/beautiful-correlation-plots-in-r-a-new-approach-d3b93d9c77be
# http://www.sthda.com/english/wiki/correlation-analyses-in-r

# Plot GO, KEGG & Reactome corrplots all togehter in one plot
# 1) Without NA replacement (circle, alphabetically sorted)
# 2) NA replacement by 0 with clustering after: 
#     - average linkage
#     - single "
#     - complete "
#     - ward.D2

### Functions ### --------------------
# Wrote my own cor.mtest function to cope with missing values and with groups that have < 2 common terms
# from which no pvalue can be computed via the cor.mtest function.
# Further I integrated the p.adjust() in the function. So obtained pvalues can be directly corrected for 
# multiple comparisons. Default padjM = "bonferroni". See ?p.adjust for details.
### MY.cor.mtest ### 
MY.cor.mtest <- function(mat, padjM, method, ...){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  padj.mat <- p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(padj.mat) <- diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  
  ## cor.test loop ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      df <- na.omit(data.frame(x = mat[, i], y = mat[, j]))
      if(nrow(df) > 2) { 
        
        # common set greater 2 -> compute cor value!
        tmp <- cor.test(x = df$x, y = df$y, method = method, ...)
        #tmp <- cor.test(x = df$x, y = df$y)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        if (!is.null(tmp$conf.int)) {
          lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
          uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
        }
        
      } else {
        
        # assign custom stats param 
        warning(paste("Not enough common observations for:",colnames(mat)[i],"vs",colnames(mat)[j],"to compute a p.value.\np.value set to 1\n"))
        p.mat[i, j] <- p.mat[j, i] <- 1
        lowCI.mat[i, j] <- lowCI.mat[j, i] <- NA
        uppCI.mat[i, j] <- uppCI.mat[j, i] <- NA
        
      }
    }
  }
  
  ## pvalue correction ##
  if(missing(padjM)){padjM = "bonferroni"}
  if(padjM != "none"){
    # Extract pval from matrix
    pval <- c()
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        pval <- append(pval, p.mat[j,i])
      }
    }
    # Run pval correction
    padj <- p.adjust(pval, method = padjM, ...)
    #padj <- p.adjust(pval, method = padjM)
    
    # Sort corrected pval back in matrix
    tmp <- padj
    for(i in 1:(ncol(mat)-1)){
      for(j in (i+1):ncol(mat)){
        padj.mat[i,j] <- padj.mat[j,i] <- tmp[1]
        tmp <- tmp[-1]
      }
    }
    # Bind results in list WITH corrected pvalues
    corMtx <- list(p.adj = padj.mat,
                   p = p.mat, 
                   lowCI = lowCI.mat, 
                   uppCI = uppCI.mat)
  } else {
    corMtx <- list(p = p.mat, 
                   lowCI = lowCI.mat, 
                   uppCI = uppCI.mat)
  }
  
  ## append row and col names
  lapply(corMtx, function(x){
    row.names(x) <- colnames(mat)
    colnames(x) <- colnames(mat)
    x})
  # Done! :)
}

## function of n, where n is character string of the *.df object (Output of fileIntegration function)
# ! Keep in mind that *.df are prefiltered !!! Default padj < 0.05
myCorrPlot1 <- function(n, method ,cor , order, clust, padjM, Counts, rmvNonSig, sig.level, ...){
  # n = "kegg.df"
  # set title for out file name (fn)
  if(n == "reac.df"){fn = "reactome"}else{fn = gsub(".df","",n)}
  if(missing(rmvNonSig) & missing(Counts)){
    Counts = T # display number of observed common counts in corrplot
    rmvNonSig = F # If T: remove non-sign correlation values from corrplot, threshold set by sig.level
  }
  if(rmvNonSig == T){Counts = F}else{Counts = T} # can either show one or the other unfortunately ... problem of corrplot
  
  # read data
  x <- get(n)
  df <-  x[,grep("NES.",colnames(x))] #reduce input df to NES columns
  row.names(df) <- x$ID
  colnames(df) <- gsub("NES.","",colnames(df))
  
  # Cor. Matrix & countMtx
  if(missing(cor)){cor = 'spearman'} #pearson is default
  if(missing(padjM)){padjM = "BH"}
  M <- cor(df, use = "pairwise.complete.obs", #"pairwise.complete.obs": correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables.
           method = cor) # try 'spearman'(on ranks) or pearson (default)
  P <- MY.cor.mtest(df, padjM = padjM, method = cor)[[1]] #store only padj / pvalues in P
  C <- qcrMtx(df)$counts
  
  # set colors 
  col=rev(brewer.pal(n=8, name="RdYlBu"))
  
  # corrplot params
  if(missing(method)){method = "circle"}
  if(missing(order)){
    order = "alphabet"
    clust = NULL
  }
  if(order != "hclust"){title = paste0(fn," NES ",cor,".cor of common enriched terms")}
  if(order == "hclust" & missing(clust)){clust = "average"}
  if(order == "hclust"){
    title = paste0(fn," NES ",cor,".cor: ",clust," clustered")
    M[is.na(M)] <- 0 #replace NAs in cor matrix with 0
  }
  
  ## Corrplot
  # widen the margin of plot (b,l,t,r) DEFAULT: 5 4 4 1
  par(mar = c(5,4,6,1)+.1)
  if(Counts == T){
    # plot with counts!
    corrplot(M, type = "upper", tl.col="black", tl.srt=35, diag = T, addgrid.col = NA, na.label = " ",
           col = col,
           method = method,
           order = order,
           hclust.method = clust,
           title = title,
           p.mat = C, sig.level = -1, insig = "p-value", #add counts
           mar=c(0,0,2,0),
           ...
           )
  } else {
    # plot non-sign. removed!
    if(missing(sig.level)){sig.level = .05}
    corrplot(M, type = "upper", tl.col="black", tl.srt=35, diag = T, addgrid.col = NA, na.label = "nc",
             col = col,
             method = method,
             order = order,
             hclust.method = clust,
             title = paste(title,padjM,sig.level),
             p.mat = P, sig.level = sig.level, insig = "blank", #remove circles below sig.level
             mar=c(0,0,2,0),
             ...
    )
  }
}
################# 

## run plotting loop for sign enriched NES values
message("\nTreatment Correlations - Comparative plot ... ")
obj = objects(pattern='*[.]df')[c(4,5,1,2,3)]# resort for nicer arranged output

pdf(file = "NES_corr_sign.pdf", width = 14, height = 20, onefile = T, bg = "transparent")
par(mfrow=c(5,2))
# non clustered 
for(n in obj){
  myCorrPlot1(n, cor = 'spearman')
  myCorrPlot1(n, cor = 'spearman', rmvNonSig = T, padjM = "BH", sig.level = .05)
}
# clustered 
for(k in c("ward.D2","average","single","complete")){
  par(mfrow=c(5,2))
  for(n in obj){
    myCorrPlot1(n, cor = 'spearman', order = 'hclust', clust = k)
    myCorrPlot1(n, cor = 'spearman', order = 'hclust', clust = k,
                rmvNonSig = T, padjM = "BH", sig.level = .05)
  }
}
dev.off()
message("Done!")

## run plotting loop for ALL NES values (padj <= filter)
obj = objects(pattern='*[.]ALL')[c(1,4,2,5,3)]# resort for nicer arranged output

pdf(file = "NES_corr_unfiltered.pdf", width = 14, height = 20, onefile = T, bg = "transparent")
par(mfrow=c(5,2))
# non clustered 
for(n in obj){
  myCorrPlot1(n, cor = 'spearman')
  myCorrPlot1(n, cor = 'spearman', rmvNonSig = T, padjM = "BH", sig.level = .05)
}
# clustered 
for(k in c("ward.D2","average","single","complete")){
  par(mfrow=c(5,2))
  for(n in obj){
    myCorrPlot1(n, cor = 'spearman', order = 'hclust', clust = k)
    myCorrPlot1(n, cor = 'spearman', order = 'hclust', clust = k,
                rmvNonSig = T, padjM = "BH", sig.level = .05)
  }
}
dev.off()
message("Done!")

### This section is not needed anymore ... for now. Take it out of the pipeline
skip=T
if(skip != T){
  ## Plot function for corrplot analysis on *.df objects
myCorrPlot <- function(x,db,corM){
  if(missing(db)){db = NULL}
  df = x[,grep("NES.",colnames(x))] #reduce input df
  row.names(df) = x$ID
  # Cor. Matrix
  if(missing(corM)){corM = "pearson"}
  M <- cor(df, use = "pairwise.complete.obs", #"pairwise.complete.obs": correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables.
           method = corM)
  # shorten row and col names of M 
  row.names(M) <- gsub("NES.","",row.names(M))
  colnames(M) <- gsub("NES.","",colnames(M))
  
  # set colors 
  col=rev(brewer.pal(n=8, name="RdYlBu"))
  
  # Unclustered - entire M without NA removal
  par(mfrow=c(2,1))
  par(mar = c(5,4,6,1)+.1)
  corrplot(M, method="circle", type = "upper", tl.col="black", tl.srt=35, diag = F, addgrid.col = NA, 
           order = 'alphabet', col=col, na.label = "nc",
           title = paste0(gsub(".df","",db)," - ",corM,".cor of common sign. NES"),
           mar=c(0,0,2,0))
  
  corrplot(M, method="ellipse", type = "upper", tl.col="black", tl.srt=35, diag = F, addgrid.col = NA, 
           order = 'alphabet', col=col, na.label = "nc",
           #title = paste0(gsub(".df","",db)," - Pearson Correlation of common sign. NES"),
           mar=c(0,0,2,0))
  
  # Replace NAs with 0 for clustering
  m = M
  m[is.na(m)] <- 0
  # for loop with clustering methods
  clustMethod = c("average","single","complete","ward.D2")
  par(mfrow=c(2,2))
  par(mar = c(5,4,6,1)+.1)
  for(k in clustMethod){
    corrplot(m, type="upper", method = "circle", order="hclust", hclust.method= k,
             tl.col="black", tl.srt=35, diag = F, addgrid.col = NA, col=col,
             title = paste0(gsub(".df","",db)," - ",corM,".cor. - ",k," clustered"),
             mar=c(1,0,2,0)
    )
  }
}
## Run plotting loop
message("\nTreatment Correlations ... ")

obj = objects(pattern='*[.]df')
for(n in obj){
  if(n == "reac.df"){fn = "reactome"}else{fn = gsub(".df","",n)} #set title for out file name (fn)
  message(paste0("\nStart treatment correlation plotting for:\t", fn," results"))
  x = get(n)
  dir.create(fn, showWarnings = F)
  pdf(file = paste0(fn,"/",fn,"_Correlation_sign.pdf"), #print directly to pdf
      width = 9.6, height = 7.3, onefile = T,
      bg = "transparent", #Background color
      fg ="black"         #Foreground color
  )
  myCorrPlot(x,fn, corM = "spearman")
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  message("Done!")
}

obj = objects(pattern='*[.]ALL')
for(n in obj){
  if(n == "reac.ALL"){fn = "reactome"}else{fn = gsub(".ALL","",n)} #set title for out file name (fn)
  message(paste0("\nStart treatment correlation plotting for:\t", fn," results"))
  x = get(n)
  dir.create(fn, showWarnings = F)
  pdf(file = paste0(fn,"/",fn,"_Correlation_ALL.pdf"), #print directly to pdf
      width = 9.6, height = 7.3, onefile = T,
      bg = "transparent", #Background color
      fg ="black"         #Foreground color
  )
  myCorrPlot(x,fn, corM = "spearman")
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  message("Done!")
}
}

rm(n, obj, skip)
###############################


### Individual Condition Correlation plots ----------------------------------------------------
# Plot common set of observation with padj < 0.05
# Color red dots which are common AND padj < 0.05

# Indiv. Corrplot/Corrtest function 
magicFun = function(fn, pcut = .05, filter = 1, method = "spearman", alternative = "greater", ...){
  ls <- get(fn)
  # filter input list to common set of terms (universe = univ)
  univ <- Reduce(intersect, lapply(ls, function(x){as.factor(x$ID)}))
  ls <- lapply(ls, function(x){x[which(x$ID %in% univ),]})
  comb = combn(names(ls),2,simplify = F) # Get all possible combinations for corrplots
  
  for(i in c(1:length(comb))){
    x=comb[[i]]
    # names for axis
    x1 = gsub("_gse.*","",x[1])
    x2 = gsub("_gse.*","",x[2])
    # data import / filtering
    df1 = subset(ls[[x[1]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= filter)
    #row.names(df1) = df1$ID
    df2 = subset(ls[[x[2]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= filter)
    #row.names(df2) = df2$ID
    
    # continue only if there are common sets
    if((length(intersect(df1$ID,df2$ID)) > 0) == T){
      # filter only significant enriched terms (padj <= pcut)
      tmp = list(df1 = subset(df1, p.adjust <= pcut, select = c(ID,NES,p.adjust)),
                 df2 = subset(df2, p.adjust <= pcut, select = c(ID,NES,p.adjust)))
      names(tmp) = gsub("_gse.*","",x)
      int = intersect(tmp[[1]]['ID'],tmp[[2]]['ID'])
      
      # Run venn & NES corr. plot for combinations that show common sign. enriched terms
      if(length(int$ID) > 0){
        message(paste0("Common enriched terms for:\n",x[1]," vs ",x[2]," with padj < ",pcut,":\t",length(int$ID)))
        
        # set colors
        anno.col = setNames(c("deepskyblue1","mediumblue","red2"),
                            c(gsub("_gse.*","",x[1]),gsub("_gse.*","",x[2]),"Common"))
        
        ## Venn plot
        tmp1 = lapply(tmp, function(x){x$ID})
        gg = plot(eulerr::euler(tmp1),
                  fills = list(fill = anno.col, alpha = 0.4),
                  #labels = list(col = "black", font = 4),
                  legend = list(col = "black", font = 4),
                  #main = paste0(fn,": sign. terms [padj < ",pcut,"]"),
                  main = paste0(fn," sign. terms"),
                  quantities = TRUE,
                  shape = "circle",
                  lty = 0)
        
        ## Corr plot (using df1 & df2)
        # format df for ggplot
        tmp2 = merge(df1,df2, by=c('ID','Description') , all = F)
        anno = apply(tmp2[,c('p.adjust.x','p.adjust.y','ID')], 1, function(X){
          if(X[1] > pcut & X[2] > pcut){'n.s.'}
          else if(X[1] <= pcut & X[2] <= pcut){'Common'}
          else if(X[1] <= pcut){gsub("_gse.*","",x[1])}
          else if(X[2] <= pcut){gsub("_gse.*","",x[2])}
        })
        stopifnot(length(anno) == nrow(tmp2)) #check
        tmp2$Type <- c(anno)
        df = tmp2[which(!(tmp2$Type %in% 'n.s.')),]  # remove n.s. types
        
        # cor test - common (padj < 0.05)
        if(nrow(df[which(df$Type %in% 'Common'),]) > 2){
          # run cor test & cor value computation
          pt = cor.test(df[which(df$Type %in% 'Common'),'NES.x'], # only common sign. terms
                        df[which(df$Type %in% 'Common'),'NES.y'],
                        alternative = alternative,
                        method = method,
                        ...)
          sl = round(pt$p.value, 8) # sign. level
          if(sl < .0001){m = '***'} # m = mark for plot annotation
          else if(sl < .001){m = '**'}
          else if(sl <= .05){m = '*'}
          else if(sl > .05){m = 'n.s.'}
        }else{
          # asign empty pt object so plotting does not return an error
          sl = 'NA'
          m = 'NA'
          pt = list(statistic = 0,
                    parameter = 0)
          # cannot run cor test due to too few common sign data points
          warning(paste("Less than 2 common significant enriched terms for",x1,"vs",x2,"\nNo cor.test possible."))
        }
        
        # run corrplot function on df object
        gg1 = ggplot2::ggplot(df, aes(x=NES.x, y=NES.y, color = Type)) +
          geom_point(size = 3, alpha = .4, shape = 16) +
          geom_point(data = df[which(df$Type %in% 'Common'),], #2nd layer
                     aes(x=NES.x, y=NES.y),
                     shape = 16, size = 3, alpha = .5) +
          scale_colour_manual(values = anno.col) +
          labs(title = paste0(fn," - NES cor. for padj < ",pcut," [",m,"]"), 
               #subtitle = paste(x1,"vs",x2),
               subtitle = paste0(method," ",alternative,": t=",round(pt$statistic,1),"  df=",pt$parameter,"  pval=",sl),
               x = paste0("NES [",x1,"]"), 
               y = paste0("NES [",x2,"]")) +
          geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
          theme_bw() + theme(aspect.ratio = 1, legend.position = 'bottom') +
          stat_quadrant_counts(colour = "black", xintercept = 0, yintercept = 0) +
          annotate("text",
                   label = paste0(method," = ",round(cor(df[df$Type %in% "Common","NES.x"],df[df$Type %in% "Common","NES.y"]),2),
                                  "\nR2 = ",round(cor(df$NES.x,df$NES.y)^2,2)),
                   color = "red2",
                   x = (min(df$NES.x)), 
                   y = 0,
                   vjust = "inward", hjust = "inward") +
          annotate("text",
                   label = paste0(method," = ",round(cor(df$NES.x,df$NES.y),2),
                                  "\nR2 = ",round(cor(df[df$Type %in% "Common","NES.x"],df[df$Type %in% "Common","NES.y"])^2,2)),
                   x = (max(df$NES.x)), 
                   y = 0,
                   vjust = "inward", hjust = "inward")
        
        # PRINT Venn & corrplot
        print(ggpubr::ggarrange(gg1,gg, ncol = 2, nrow = 1))
        
      }else{
        # message that no venn plot was created
        warning(paste0("No common sign. enriched terms found for:\n",x[1]," vs ",x[2]," with padj < ",pcut,"\nNo Venn plot created."))
      }
    }else{
      warning(paste0("No common enriched terms found for:\n",x[1]," vs ",x[2]," with padj < ",filter,"\nNo plot created."))
    }
    # END of lapply
  }
  # END OF magicFUN  
}

# Indiv. Corrplot/Corrtest ploting 
for(k in gseaLists){
  if(k == 'reac.ls'){n = 'reactome'}else{n = gsub('.ls','',k)}
  message(paste0("\nStart venn & pairwise correlation plotting for:\t",n))
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste0(n,"/",n,"_VennAndPairwCorrelation.pdf"), width = 10.6, height = 5.9, onefile = T,
      bg = "transparent", fg ="black")
  magicFun(k)
  dev.off()
  message("\nDONE!\n")
  gc()
}
rm(k,n)
##############################################


### Heatmap: Selection and plotting -----------------------------------
# Heatmap function
myHeatFun = function(fn, pcut = .05,top,preFilter,rowclust,colclust,distM){
  ls = get(fn)
  # Data set selection --------
  if(fn == 'reac.ls'){n = 'reactome'}else{n = gsub('.ls','',fn)}
  # Get ids of terms that show overlap between pairwise comparison in each condition
  comb = combn(names(ls),2,simplify = F) 
  tmp = lapply(comb, function(x){
    df1 = subset(ls[[x[1]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= pcut)
    df2 = subset(ls[[x[2]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= pcut)
    intersect(df1$ID,df2$ID)
  })
  ids = unique(unlist(tmp)) #546
  
  # Now selecting those ids from the wide df objects
  # *.df is preFilter with padj < 0.05 / *.ALL with padj < 0.2
  if(missing(preFilter)){preFilter = T}
  if(preFilter == T){
    wdf = get(gsub('.ls','.df',fn))
  }else{
    wdf = get(gsub('.ls','.ALL',fn))
  }
  
  row.names(wdf) = wdf$ID
  # reduce to mtx
  mtx = wdf[ids, grep("NES.",colnames(wdf))]
  colnames(mtx) = gsub('NES.','',colnames(mtx)) #shorten col names
  # replace NAs with 0
  mtx[is.na(mtx)] <- 0 
  # Idea: compute the row sum of the abs(NES) values. 
  # Rows with larger values have fewer zeroes and/or greater NES scores. 
  # This can be used to sort the rows in mtx accordingly. 
  rank = apply(mtx,1,function(x){sum(abs(x))})
  rank = sort(rank, decreasing = T)
  stopifnot(nrow(mtx) == length(rank)) #check
  mtx = mtx[names(rank),]
  
  # Heatmap plotting ----------------------------------------------------------------
  if(missing(rowclust)){rowclust = T}
  if(missing(colclust)){colclust = T}
  if(missing(distM)){distM = 'euclidean'} # euclidean, maximum, manhattan, ...
  if(missing(top)){top = nrow(mtx)}
  if(top == 'max'){top = nrow(mtx)}
  if(top > nrow(mtx)){top = nrow(mtx)}
  
  # set anno colors & breaks
  x = 80 #length of color palette
  col = colorRampPalette(c("#4575B4","white","#D73027"))(x)
  df = as.matrix(mtx)
  s = sd(df)*.25
  m = 0
  myBreaks = c(seq(min(df), m-s, length.out=ceiling(x/2) + 1), 
               seq(m+s, max(df), length.out=floor(x/2)))
  
  # create annotation df
  coldat = data.frame(ID = colnames(mtx),
                      Substance = as.factor(gsub('_.*','',colnames(mtx))),
                      Condition = as.factor(gsub('.+?_','',colnames(mtx))),
                      row.names = 1)
  stopifnot(colnames(mtx) == row.names(coldat)) #check
  Substance = levels(as.factor(coldat$Substance))
  Condition = levels(as.factor(coldat$Condition))
  if(length(Condition)==3){Condition = Condition[c(2,3,1)]} #set to low, mid,high
  
  # color blind friendly palette [https://www.r-bloggers.com/2013/10/creating-colorblind-friendly-figures/]
  cbl = c("#00FF4D", "#E69F00", "#56B4E9","#F0E442","#CC79A7","#0072B2","#009E73","#D55E00","#005DFF")
  if(length(cbl) > length(Substance)){cbl = cbl[1:length(Substance)]}else{cbl=cbl}
  col.Subs = colorRampPalette(cbl)(length(Substance))
  col.Cond = colorRampPalette(c("gray95","gray50"))(length(Condition))
  ann_colors = list(Substance = setNames(col.Subs, Substance),
                    Condition = setNames(col.Cond, Condition))
  
  # heatmap
  M = na.omit(mtx[c(1:top),])
  print(
    pheatmap::pheatmap(M,
                       clustering_method = "average",
                       cluster_rows = rowclust, cluster_cols = colclust, 
                       clustering_distance_rows = distM, clustering_distance_cols = distM,
                       show_rownames= T, show_colnames = T, border_color = NA, color = col, 
                       breaks = myBreaks, #treeheight_col = 30, #default 50
                       annotation_col = coldat[,c(2,1)], #resorted
                       annotation_colors = ann_colors,
                       main = if(preFilter == T){
                         paste0("NES of sign. (padj<",pcut,") ",n," Top: ",nrow(M))
                       }else{
                         paste0("NES of ",n," Top: ",nrow(M))
                       }
    )
  )
  # -----------------
}

# Heatmap plotting
for(fn in gseaLists){
  if(fn == 'reac.ls'){n = 'reactome'}else{n = gsub('.ls','',fn)}
  message(paste0("\nStart heatmap plotting for:\t",n))
  
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste0(n,"/",n,"_NES_Heatmaps.pdf"),
      width = 10, height = 8, onefile = T,
      bg = "transparent", #Background color
      fg ="black"         #Foreground color
  )
  
  for(i in c('max',250,100,50,25)){
    message(paste0("NES Top:\t",i))
    myHeatFun(fn, top = i, preFilter = T)
    myHeatFun(fn, top = i, preFilter = F)
    myHeatFun(fn, top = i, preFilter = T, colclust = F)
    myHeatFun(fn, top = i, preFilter = F, colclust = F)
  }
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  message("DONE!\n")
}
rm(fn,i,n)

# Export NES matrix
mtxExp = function(fn, pcut = .05, preFilter = T){
  ls = get(fn)
  if(fn == 'reac.ls'){n = 'reactome'}else{n = gsub('.ls','',fn)}
  # Get ids of terms that show overlap between pairwise comparison in each condition
  comb = combn(names(ls),2,simplify = F) 
  tmp = lapply(comb, function(x){
    df1 = subset(ls[[x[1]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= pcut)
    df2 = subset(ls[[x[2]]][,c("ID","Description","NES","pvalue","p.adjust")], p.adjust <= pcut)
    intersect(df1$ID,df2$ID)
  })
  ids = unique(unlist(tmp))
  
  # Now selecting those ids from the wide df objects
  # *.df is preFilter with padj < 0.05 / *.ALL with padj < 0.2
  if(preFilter == T){
    wdf = get(gsub('.ls','.df',fn)) #*.df is prefiltered with padj < 0.05
  }else{
    wdf = get(gsub('.ls','.ALL',fn)) #*.df is prefiltered with padj < 0.2
  }
  
  row.names(wdf) = wdf$ID
  # reduce to mtx
  mtx = wdf[ids, c(1,2,grep("NES.",colnames(wdf)))]
  colnames(mtx) = gsub('NES.','',colnames(mtx)) #shorten col names
  # replace NAs with 0
  mtx[is.na(mtx)] <- 0
  mtx
}
mtx.ls = list() #Can be used for PCA, t-SNE, LDA, ... 
for(fn in gseaLists){
  if(fn == 'reac.ls'){n = 'reactome'}else{n = gsub('.ls','',fn)}
  pcut = 0.05
  # ".df" => only contains the statistically sign. enriched pathways
  tmp = mtxExp(fn, preFilter = T, pcut = pcut)
  mtx.ls[[paste0(n,".df")]] = tmp
  write.csv2(tmp, row.names = F, file = paste0(n,'/',n,'_NESmtx_prefilt',pcut,'.csv')) #*.df
  
  # ".ALL" => all values
  #tmp = mtxExp(fn, preFilter = F, pcut = pcut)
  #mtx.ls[[paste0(n,".ALL")]] = tmp
  #write.csv2(tmp, row.names = F, file = paste0(n,'/',n,'_NESmtx_ALL.csv')) #*.ALL
}
rm(pcut,tmp,fn,n)
#######################################


### MCA, Sim.Mtx on CATEGORIES --------------------------------------------------------------------
# filter mtx.ls; Ready to be used for LDA,t-SNE,PCA, ...
#mtx.ls = lapply(mtx.ls, function(x){x[,grep('_',colnames(x))]})
#require(FactoMineR)
# CAT-PCA #https://www.ibm.com/support/knowledgecenter/en/SSLVMB_23.0.0/spss/categories/idh_cpca.html
# Idea: Run a category PCA on the three categories:
# +1 = GSEA term sign. upregulated (NES > 0)
# -1 = GSEA term sign. downregulated (NES < 0)
# 0  = GSEA term not-sign regulated
# => Run a Multiple Correspondence Analysis (MCA) on this categorical matrix
# https://stats.stackexchange.com/questions/159705/would-pca-work-for-boolean-binary-data-types
##################################


### END of SCRIPT -----------
message(paste0("\nFinished compareGSEA script comparing:\n",paste(dataSet, collapse = " "),"\nSaving compareGSEA.RData ..."))
save.image("compareGSEA.RData")
message(paste0("Done!\t compareGSEA.RData saved under:\n",getwd(),"/compareGSEA.RData"))

## Session Information ##
sink(paste0("SessionInfo_compareGSEA.txt"))
print(date())
print(devtools::session_info())
sink()

setwd(home)
message("\nPuh! What a ride! All Done!\nJ.A.R.V.I.S over and out! :)")
#####################