######################################################
###  ENRICHED PATHWAY ANALYSIS for DESeq2 Results  ###
######################################################

# Author:   Hannes Reinwald
# Contact:  hannes.reinwald@ime.fraunhofer.de

# README: -------------------------------------------------------------------------------
# This script performs Over Represenation Analysis (ORA) with a set of predefined diff. expressed
# genes / proteins against the background of all commonly detected genes/proteins &
# Gene Set Enrichment Analysis (GSEA) based on the sorted log2-FCs. 
# Here we use ClusterProfiler & ReactomePA for this purpose. For details please see: 
# https://yulab-smu.github.io/clusterProfiler-book/index.html
# To increase the power of GSEA we recommend to use apeglm shrunk Effect sizes (log2-FC) for the analysis

# More details about the principals of GSEA:
# https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/

# Input ORA:
# Data table list of DEGs/DEPs with ENSEMBL/PROT IDs, apeglm(LFC), ENTREZ, padj 
# Universe - set of background genes (all detected genes in the experiment) => This background can be used for GSEA then

# Input GSEA: 
# Data table with ENSEMBL / PROT IDs, apeglm(LFC), ENTREZ, padj 
# Transcriptomics - Ideally the result table from the DESeq2 with apeglm shrunk LFC (named with _reslfs) and annotated ENTREZ IDs
# Proteomics - I recommend to shrinking here as well on the LFC tables.

# For now I am developing this code mostly for the DESeq2 output files (Check path dependency).
# The idea is to put all the results files, which shall be compared into a single folder and then run the script. 
# Navigate to substance folder and execute this script (i.e "S:/data/RNA_seq_gene_counts/6-PTU_2/")
#############

# Load packages --------------------------------------------------------------------------
require(ReactomePA)
require(clusterProfiler)
require(AnnotationDbi)
require(org.Dr.eg.db)
require(ggplot2)
require(DOSE)
#BiocManager::install(c("ReactomePA","clusterProfiler","AnnotationDbi","org.Dr.eg.db","DOSE"), update = T, ask=F)
###################

# Navigate to your DESeq2 output dir (i.e. "./DESeq2_Pairwise")
# The dir must contain two sub dir called:
# DEG: Filtered DESeq2 output for the DEGs (for now this script uses the lfs as input, change that if needed)
# Results: Entire DESeq2 result tables
setwd("DESeq2_Pairwise")
HOME = getwd()


### DATA IMPORT ###
# Import DEG tables (apgelm shrunk lfcut & pcut) for ORA ----
files <- list.files(path = "DEGs/lfsFCcut/", pattern = ".csv")
if(length(files) == 2){files <- files[c(2,1)]}
if(length(files) == 3){files <- files[c(2,3,1)]}
inputORA <- list()
for(i in files){
  x <- read.csv2(paste0("DEGs/lfsFCcut/",i),header = T, row.names = 1)
  name <- gsub(".csv","",i)
  name <- gsub("_lfs_pcut_LFcut","",name)
  inputORA[[name]] <- x
}

# Import complete res tables (apgelm shrunk lfcut & pcut) for GSEA & universe definition ----
files <- list.files(path = "Results/", pattern = "_reslfs_")
if(length(files) == 2){files <- files[c(2,1)]}
if(length(files) == 3){files <- files[c(2,3,1)]}
inputGSEA <- list()
for(i in files){
  x <- read.csv2(paste0("Results/",i),header = T, row.names = 1)
  x <- x[which(!(is.na(x$pvalue))),]# rmv genes with no computed pvalue
  name <- gsub(".csv","",i)
  name <- gsub("_lfs_pcut_LFcut","",name)
  inputGSEA[[name]] <- x
}
rm(name,x,i)

# Unify tables from protein and transcript output ---------------------------------------
# Replace "log2FoldChange" with "log2FC"
if(any(colnames(inputGSEA[[1]]) %in% "log2FoldChange")){
  inputGSEA <- lapply(inputGSEA, function(x){
    colnames(x)[which(colnames(inputGSEA[[1]]) %in% "log2FoldChange")] <- "log2FC"
    x
  })
  inputORA <- lapply(inputORA, function(x){
    colnames(x)[which(colnames(inputORA[[1]]) %in% "log2FoldChange")] <- "log2FC"
    x
  })
}
# Replace "adj.pvalue" with "padj"
if(any(colnames(inputGSEA[[1]]) %in% "adj.pvalue")){
  inputGSEA <- lapply(inputGSEA, function(x){
    colnames(x)[which(colnames(inputGSEA[[1]]) %in% "adj.pvalue")] <- "padj"
    x
  })
  inputORA <- lapply(inputORA, function(x){
    colnames(x)[which(colnames(inputORA[[1]]) %in% "log2FoldChange")] <- "padj"
    x
  })
}

# Check type of Input (Ensembl prot or gene ID): Set KEY parameter -----------------------
x = unlist(lapply(inputORA, row.names))
if(all(grepl("ENSDARG", x))){
  message("\nInput: ENSEMBL Gene IDs - All good :)\n")
  KEY = "ENSEMBL"
}else if(all(grepl("ENSDARP", x))){
  message("\nInput: ENSEMBL Protein IDs - All good :)\n")
  KEY = "ENSEMBLPROT"
}else if(any(grepl("ENSDARG", x))){
  stop("\nNot all Genes in your datafile are annotated with an ENSEMBL Gene ID.\nPlease check if there might be some other IDs than ENSEMBL Gene IDs in your list.")
}else if(any(grepl("ENSDARP", x))){
  stop("\nNot all Proteins in your datafile are annotated with an ENSEMBL Protein ID.\nPlease check if there might be some other IDs than ENSEMBL Protein IDs in your list.")
}else{
  stop("\nCould not identify the ID TYPE provided. Please check your input files again (inputORA).\nThis script only supports ENSEMBL Gene or Protein Identifiers")
}

# Check input - If ENTREZID missing, annotate! ------------------------------------------
if(any(unlist(lapply(inputGSEA, colnames)) == "ENTREZID") == T){
  message("ENTREZID found in input - All good :) Keep going ...\n")
}else{
  message("No ENTREZID found in input - Appending ENTREZIDs ...\n")
  # Map EntrezIDs to ENSEMBL Prot IDs for enrichPathway & enrichKEGG
  #columns(org.Dr.eg.db) #ENSEMBLPROT; ENTREZID
  AnnoFun <- function(X, key, colInfo){
    x <- X
    anno <- mapIds(org.Dr.eg.db,
                   keys=row.names(x),
                   keytype= key,
                   column= colInfo,
                   multiVals= "first")
    x[,paste0(colInfo)] <- anno
    x
  }
  inputORA <- lapply(inputORA, FUN = AnnoFun, key=KEY, colInfo = "ENTREZID")
  inputGSEA <- lapply(inputGSEA, FUN = AnnoFun, key=KEY, colInfo = "ENTREZID")
}

# Define common gene sets = universe (background for ORA) -------------------------------

# Common gene set definition only realy needed when comparing results from differen experiments
# as identified proteins / genes could differ between experiments
# But we will still do it anyways to make this script appliciable for a wider usage.
ls <- lapply(inputGSEA, function(x){row.names(x)})
univ <- Reduce(intersect, ls) #univ (universe) is the common set of genes/proteins from all inputs

# Still need a univ.entrez (universe with ENTREZID)
ls <- lapply(inputGSEA, function(x){x[which(!is.na(x$ENTREZID)),"ENTREZID"]})
univ.entrez <- Reduce(intersect, ls)
rm(ls)

# Create FINAL INPUT tables GSEA / ORA --------------------------------------------------
# (with common gene sets (univ)) for GSEA & ORA (name & log2FC)

# Sort input after Log2FC 
inputGSEA <- lapply(inputGSEA, function(x){x[order(x$log2FC, decreasing = T),]})
inputORA <- lapply(inputORA, function(x){x[order(x$log2FC, decreasing = T),]})

# ENSEMBL / ENSEMBLPROT for GO
GSEA <- lapply(inputGSEA, function(x){
  x <- x[row.names(x) %in% univ,] # select only common set of genes
  y <- na.omit(x$log2FC)
  names(y) <- as.character(row.names(x))
  y
})

# creat union list with all identified DEGs & compare to univ. Then only pick those,
# which are element of univ for the final ORA input list.
ls <- lapply(inputORA, function(x){row.names(x)})
tmp <- Reduce(union, ls)
tmp <- intersect(univ,tmp)
ORA <- lapply(inputORA, function(x){
  x <- x[row.names(x) %in% tmp,]
  y <- na.omit(x$log2FC)
  names(y) <- as.character(row.names(x))
  y
})

# ENTREZID (rmv ENTREZID NAs) for Reactome & KEGG (& DOSE, not implemented yet)
GSEA.entrez <- lapply(inputGSEA, function(x){
  x <- x[x$ENTREZID %in% univ.entrez,] # select only common set of genes
  x <- x[which(!is.na(x$ENTREZID)),] # rmv entries without ENTREZID, multiple copies of the same ENTREZID present!!!
  y <- na.omit(x$log2FC)
  names(y) <- as.character(x$ENTREZID)
  y
})

ls <- lapply(inputORA, function(x){x$ENTREZID[which(!is.na(x$ENTREZID))]})
tmp <- Reduce(union, ls)
tmp <- intersect(univ.entrez,tmp)
ORA.entrez <- lapply(inputORA, function(x){
  x <- x[x$ENTREZID %in% tmp,]
  y <- na.omit(x$log2FC)
  names(y) <- as.character(x$ENTREZID)
  y
})
rm(ls,tmp)

# Create Output folders in "clusterProfiler"
dir.create("clusterProfiler")
dir.create("clusterProfiler/ORA")
dir.create("clusterProfiler/GSEA")

# extract substance name/s! Could be used for between substance comparison as well
subs <- unique(sub("\\_.*", "", files))
subs <- paste(subs, collapse = "_") # merge to single string with _ seperator
###################


#### ORA - only DE genes /proteins + universe! #############################################
## ORA Parameters ------------
qcut = 1.1       # qvalueCutoff -> to export ALL results!
pcut = 1.1      # pvalueCutoff -> to export ALL results!
padjm = "BH"    # pAdjustMethod: "none", "BH" - benjamini&hochber, ... see ?enrichPathway for details
MIN = 10        # minGSSize
MAX = 500       # maxGSSize

## CompareCluster Function -----------------------------------------------------------------
compClustFUN <- function(res, substance, top, db){
  # For debugging:
  #res=kegg
  #top=15
  #db="KEGG"
  #substance=subs
  
  # Parameters
  outdir = "./clusterProfiler/ORA/"
  if(missing(substance)){substance = ""} # substance name which was tested
  if(missing(top)){top = 15} # top enriched sets to be displayed in dotplot
  if(missing(db)){stop("\nPlease provide info which enrichment function was used.\nUse one of the following for the db parameter:\nReactome, KEGG, KEGGM, GO")}
  # run function
  message(paste0("Saving compareCluster results to csv in:\n",getwd(),outdir))
  write.csv2(as.data.frame(res), paste0(outdir,db,"_",substance,"_compareCluster.csv"))
  
  if(nrow(as.data.frame(res)) > 0){
    message(paste0("Printing dotplots to pdf ... "))
    pdf(paste0(outdir,db,"_",substance,"_compareCluster.pdf"), 
        width = 12, height = top/2, onefile = T, bg = "transparent", fg ="black"
    )
    # if all padjuste != NA, print plot
    if(any(!is.na(as.data.frame(res)[,"p.adjust"]))){
      print(
        dotplot(res, x=~Count, showCategory = top, color = "p.adjust",
                title = paste(db,"-",substance)) + facet_grid(~Cluster) + 
          aes(Count,reorder(stringr::str_wrap(Description, 55),Count)) + ylab(NULL)
      )
    }else{message("p.adjusted values not a continous scale.\nNo plot created!")}
    # if all qvalues != NA, print plot
    if(any(!is.na(as.data.frame(res)[,"qvalue"]))){
      print(
        dotplot(res, x=~Count, showCategory = top, color = "qvalue",
                title = paste(db,"-", substance)) + facet_grid(~Cluster) +
          aes(Count,reorder(stringr::str_wrap(Description, 55),Count)) + ylab(NULL)
      )
    }else{message("qvalues not a continous scale.\nNo plot created!")}
    dev.off()
  }else{
    message(paste("No significant enriched terms found for",substance,"in -",db))
  }
} # END of function

## REACTOME & KEGG (using ENTREZID) --------------------------------------------------------
ls.ora <- lapply(ORA.entrez, function(x){names(x)})
# Reactome
message("\n\nORA for: Reactome Pathways (compareCluster())\nComputing ... ")
reactome <- compareCluster(ls.ora, organism = "zebrafish",
                           fun = "enrichPathway", #One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "?enrichPathway"
                           universe = univ.entrez,
                           pvalueCutoff = pcut, qvalueCutoff = qcut, pAdjustMethod = padjm, 
                           minGSSize = MIN, maxGSSize = MAX, readable = T)
# Export results
compClustFUN(reactome, substance = subs, top = 15, db = "Reactome")
message("Done!\n")

# KEGG
message("\nORA for: KEGG Pathways (compareCluster())\nComputing ... ")
kegg <- compareCluster(ls.ora, organism = "dre", keyType = "ncbi-geneid",
                       fun = "enrichKEGG", 
                       universe = univ.entrez, 
                       pvalueCutoff = pcut, qvalueCutoff = qcut, pAdjustMethod = padjm,
                       minGSSize = MIN, maxGSSize = MAX, use_internal_data = F)
# Export results
compClustFUN(kegg, substance = subs, top = 15, db = "KEGG")
message("Done!\n")

#mkegg <- compareCluster(ls.ora, organism = "dre", keyType = "ncbi-geneid",
#fun = "enrichMKEGG", 
#universe = univ.entrez, 
#pvalueCutoff = pcut, qvalueCutoff = qcut, pAdjustMethod = padjm,
#minGSSize = MIN, maxGSSize = MAX)
#compClustFUN(mkegg, substance = subs, top = 15, db = "MKEGG")

## GO (unique ENSEMBL /ENSEMBLPROT IDs) -----------------------------------------------------
ls.ora <- lapply(ORA, function(x){names(x)})
for(ont in c("BP","MF","CC")){
  message(paste0("\nORA for: GO(",ont,")\nComputing ... "))
  GO <- compareCluster(ls.ora, OrgDb = org.Dr.eg.db, keyType = KEY, ont = ont,
                       fun = "enrichGO", universe = univ, 
                       pvalueCutoff = pcut, qvalueCutoff = qcut, pAdjustMethod = padjm,
                       minGSSize = MIN, maxGSSize = MAX, readable = T, pool = T)
  compClustFUN(GO, substance = subs, top = 15, db = paste0("GO.",ont))
  message("Done!\n")
}

rm(GO,kegg,reactome)
gc() # Clear memory 
############################################################################################


### Gene Set Enrichment Analysis (GSEA) via gseGO, gsePathway, gseKEGG #####################
# Reactome -----------------------------------------------------------------
#?gsePathway
MyGsePathway <- function(gse.ls, title, top, pcut, padjm, MIN, MAX){
  gseFun = "gseReactome"
  if(missing(title)){title = ""}
  if(missing(top)){top = 15}
  if(missing(pcut)){pcut = .1}
  if(missing(padjm)){padjm = "none"} #padj-method: "none", "BH" - benjamini&hochber see ?enrichPathway for details
  if(missing(MIN)){MIN=10}
  if(missing(MAX)){MAX=500}
  
  # gsePathway()
  gse <- gsePathway(gse.ls, organism = "zebrafish",
                    minGSSize = MIN, maxGSSize = MAX, pvalueCutoff = pcut, nPerm = 100,
                    pAdjustMethod = padjm, verbose = TRUE, seed = FALSE, by = "fgsea")
  gseR <- setReadable(gse, OrgDb = org.Dr.eg.db, keyType="ENTREZID") # convert ENTREZ IDs to readable files
  
  # export results -----
  write.csv2(as.data.frame(gseR), paste0(title,"_",gseFun,".csv"))
  
  if(nrow(as.data.frame(gse)) > 1){
    pdf(paste0(title,"_",gseFun,".pdf"),
        width = 13, height = top/2, onefile = T, bg = "transparent", fg ="black")
    # dotplot ----
    print(
      dotplot(object = gse, showCategory = top,
              x = "GeneRatio", #"GeneRatio", #"Count"
              color = "p.adjust", #"p.adjust" "qvalue" "pvalue", 
              title = paste0(gseFun," - ",title)) + 
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    print(
      dotplot(object = gse, showCategory = top, split = ".sign", #separate result by "category" variable
              x = "GeneRatio",#"GeneRatio", #"Count"
              color = "p.adjust", 
              title = paste0(gseFun," - ",title)) + 
        facet_grid(.~.sign) +
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    # Network ---- 
    print(emapplot(gseR, color = "p.adjust", showCategory = top))
    print(ridgeplot(gse, showCategory = top, fill = "p.adjust", core_enrichment = T))
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F,
               node_label="category")
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, colorEdge = T)
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T,node_label="category")
    )
    
    dev.off()
    
    # gsea model -----
    pdf(paste0(title,"_",gseFun,"_gseaModel.pdf"),
        width = 13, height = 8, onefile = T, bg = "transparent", fg ="black")
    
    if(nrow(as.data.frame(gse)) > 0){
      if(nrow(as.data.frame(gse)) > top){x = top}else{x = nrow(as.data.frame(gse))}
      for(i in c(1:x)){
        print(gseaplot(gse, geneSetID = i, title = gse$Description[i]))
      }
    }
    dev.off()
  }else{
    message("No plotting performed")
  }
}

# KEGG ---------------------------------------------------------------------
#?gseKEGG
MyGseKEGG <- function(gse.ls, title, top, pcut, padjm, MIN, MAX){
  gseFun = "gseKEGG"
  if(missing(title)){title = ""}
  if(missing(top)){top = 15}
  if(missing(pcut)){pcut = .1}
  if(missing(padjm)){padjm = "none"} #padj-method: "none", "BH" - benjamini&hochber see ?enrichPathway for details
  if(missing(MIN)){MIN=10}
  if(missing(MAX)){MAX=500}
  
  # gseKEGG()
  gse <- gseKEGG(gse.ls, organism = "dre", keyType = "ncbi-geneid", nPerm = 10000,
                 minGSSize = MIN, maxGSSize = MAX, pvalueCutoff = pcut, use_internal_data = F,
                 pAdjustMethod = padjm, verbose = TRUE, seed = FALSE, by = "fgse")
  gseR <- setReadable(gse, OrgDb = org.Dr.eg.db, keyType="ENTREZID") # convert ENTREZ IDs to readable files
  # export results -----
  write.csv2(as.data.frame(gseR), paste0(title,"_",gseFun,".csv"))
  
  if(nrow(as.data.frame(gse)) > 1){
    pdf(paste0(title,"_",gseFun,".pdf"),
        width = 13, height = top/2, onefile = T, bg = "transparent", fg ="black")
    # dotplot ----
    print(
      dotplot(object = gse, showCategory = top,
              x = "GeneRatio", #"GeneRatio", #"Count"
              color = "p.adjust", #"p.adjust" "qvalue" "pvalue", 
              title = paste0(gseFun," - ",title)) + 
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    print(
      dotplot(object = gse, showCategory = top, split = ".sign", #separate result by "category" variable
              x = "GeneRatio",#"GeneRatio", #"Count"
              color = "p.adjust", 
              title = paste0(gseFun," - ",title)) + 
        facet_grid(.~.sign) +
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    # Network ---- 
    print(emapplot(gseR, color = "p.adjust", showCategory = top))
    print(ridgeplot(gse, showCategory = top, fill = "p.adjust", core_enrichment = T))
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F,
               node_label="category")
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, colorEdge = T)
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T,node_label="category")
    )
    
    dev.off()
    
    # gsea model -----
    pdf(paste0(title,"_",gseFun,"_gseaModel.pdf"),
        width = 13, height = 8, onefile = T, bg = "transparent", fg ="black")
    
    if(nrow(as.data.frame(gse)) > 0){
      if(nrow(as.data.frame(gse)) > top){x = top}else{x = nrow(as.data.frame(gse))}
      for(i in c(1:x)){
        print(gseaplot(gse, geneSetID = i, title = gse$Description[i]))
      }
    }
    dev.off()
  }else{
    message("No plotting performed")
  }
}

# GO -----------------------------------------------------------------------
#?gseGO
MyGseGO <- function(gse.ls, ont, key, title, top, pcut, padjm, MIN, MAX){
  gseFun = paste0("gseGO.",ont)
  if(missing(title)){title = ""}
  if(missing(top)){top = 15}
  if(missing(pcut)){pcut = .1}
  if(missing(padjm)){padjm = "none"} #padj-method: "none", "BH" - benjamini&hochber see ?enrichPathway for details
  if(missing(MIN)){MIN=10}
  if(missing(MAX)){MAX=500}
  
  # gseGO()
  gse <- gseGO(gse.ls, OrgDb = org.Dr.eg.db, keyType = key, ont = ont, nPerm = 10000,
               minGSSize = MIN, maxGSSize = MAX, pvalueCutoff = pcut,
               pAdjustMethod = padjm, verbose = TRUE, seed = FALSE, by = "fgse")
  gseR <- setReadable(gse, OrgDb = org.Dr.eg.db, keyType=key) # convert ENTREZ IDs to readable files
  
  # export results -----
  write.csv2(as.data.frame(gseR), paste0(title,"_",gseFun,".csv"))
  
  if(nrow(as.data.frame(gse)) > 1){
    pdf(paste0(title,"_",gseFun,".pdf"),
        width = 13, height = top/2, onefile = T, bg = "transparent", fg ="black")
    # dotplot ----
    print(
      dotplot(object = gse, showCategory = top,
              x = "GeneRatio", #"GeneRatio", #"Count"
              color = "p.adjust", #"p.adjust" "qvalue" "pvalue", 
              title = paste0(gseFun," - ",title)) + 
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    print(
      dotplot(object = gse, showCategory = top, split = ".sign", #separate result by "category" variable
              x = "GeneRatio",#"GeneRatio", #"Count"
              color = "p.adjust", 
              title = paste0(gseFun," - ",title)) + 
        facet_grid(.~.sign) +
        aes(GeneRatio,reorder(stringr::str_wrap(Description, 55),GeneRatio)) + ylab(NULL)
    )
    
    # Network ---- 
    print(emapplot(gseR, color = "p.adjust", showCategory = top))
    print(ridgeplot(gse, showCategory = top, fill = "p.adjust", core_enrichment = T))
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = top, colorEdge = F,
               node_label="category")
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, colorEdge = T)
    )
    
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T)
    )
    print(
      cnetplot(gseR, foldChange = gse.ls, categorySize="p.adjust", showCategory = 5, 
               colorEdge = T, circular = T,node_label="category")
    )
    
    dev.off()
    
    # gsea model -----
    pdf(paste0(title,"_",gseFun,"_gseaModel.pdf"),
        width = 13, height = 8, onefile = T, bg = "transparent", fg ="black")
    
    if(nrow(as.data.frame(gse)) > 0){
      if(nrow(as.data.frame(gse)) > top){x = top}else{x = nrow(as.data.frame(gse))}
      for(i in c(1:x)){
        print(gseaplot(gse, geneSetID = i, title = gse$Description[i]))
      }
    }
    dev.off()
  }else{
    message("No plotting performed")
  }
}

# FULL GSEAnalysis LOOP -------------------------------------------------------------------
pcut = 1.1
# Reactome
dir.create(paste0(HOME,"/clusterProfiler/GSEA/gseReactome"), showWarnings = F)
setwd(paste0(HOME,"/clusterProfiler/GSEA/gseReactome/"))
for(i in names(GSEA.entrez)){
  message(paste0("\nStarting Reactome GSEA for:\t",i))
  gse.ls <- GSEA.entrez[[i]]
  MyGsePathway(gse.ls, title = i, top = 15, pcut = pcut, padjm = "BH", MIN = 10, MAX = 500)
  message("Done!\n")
}

# KEGG
dir.create(paste0(HOME,"/clusterProfiler/GSEA/gseKEGG"), showWarnings = F)
setwd(paste0(HOME,"/clusterProfiler/GSEA/gseKEGG/"))
for(i in names(GSEA.entrez)){
  message(paste0("\nStarting KEGG GSEA for:\t",i))
  gse.ls <- GSEA.entrez[[i]]
  MyGseKEGG(gse.ls, title = i, top = 15, pcut = pcut, padjm = "BH", MIN = 10, MAX = 500)
  message("Done!\n")
}

# GO
dir.create(paste0(HOME,"/clusterProfiler/GSEA/gseGO"), showWarnings = F)
setwd(paste0(HOME,"/clusterProfiler/GSEA/gseGO/"))
for(i in names(GSEA)){
  gse.ls <- GSEA[[i]]
  for(ont in c("MF","BP","CC")){
    message(paste0("\nStarting GO(",ont,") GSEA for:\t",i))
    dir.create(ont, showWarnings = F)
    setwd(ont)
    MyGseGO(gse.ls, title = i, ont = ont, key = KEY,
            top = 15, pcut = pcut, padjm = "BH", MIN = 10, MAX = 500)
    message("Done!\n")
    setwd(paste0(HOME,"/clusterProfiler/GSEA/gseGO/"))
  }
}
rm(ont,i,qcut,pcut,gse.ls)
gc()
############################################################################################


### Save R Image ### --------------------------
message(paste0("Saving R image under:\n",HOME,"/clusterProfiler/Pathway.RData\nThis might take a while ..."))
save.image(paste0(HOME,"/clusterProfiler/Pathway.RData"))

### R Session Info ### -------------------
sink(paste0(HOME,"/clusterProfiler/SessionInfo_clusterProfiler.txt"))
print(date())
print(devtools::session_info())
sink()

setwd(HOME)
setwd("../") #Back to where we started
message("\nAll Done! END OF ORA & GSEA SCRIPT\nJ.A.R.V.I.S over and out! :)\n")
####### END OF SCRIPT #############