### merge GSEA results to a single df ###
load("S:/manuscripts/neurotox/data/compareGSEA/compareGSEA.RData")

# the object gseaLists contains the names of the imported GSEA results lists
gseaLists
ls = get(gseaLists[4])

##########################################
###   AWESOME FUNCTIONS - GSEAcombinR  ###
##########################################
filterGSEA    <- function(ls, filter= "p.adjust", pcut= .05, top, topBy= "pvalue", minCount= 10) {
  # ls    = list with GSEA result tables from clusterProfiler package
  # topBy = value to select top by. One of: "NES","pvalue","GeneRatio","Count"
  # top   = integer to select the n - top results from each object in ls
  
  # (1) Subset the df in lists for only sign. enriched terms. 'filterBy' variable will 
  # determine for which statistical value to filter for (one of: "pvalue", "p.adjust", "qvalues")
  pCutOff <- function(x, filterBy = "p.adjust", cut = .05){
    message(paste("Applied statistical cut off:\t",filterBy)," < ",cut)
    if(filterBy == "p.adjust"){ subset(x, p.adjust <= cut)
    } else if (filterBy == "pvalue") {
      subset(x, pvalue   <= cut)
    } else if(filterBy == "qvalues") {
      subset(x, qvalues  <= cut)
    } else {stop("'filter' parameter missing!\nPlease provide one of: 'pvalue', 'p.adjust', 'qvalues'")}
  }
  ext.ls = lapply(ls, pCutOff, filterBy = filter, cut = pcut)
  
  # (2) Compute gene ratio & counts (Nbr of enriched genes in cluster)
  geneRcount <- function(x){
    # compute counts - count strings separated with / in core_enrichment column
    x$Count <- stringr::str_count(x$core_enrichment, "/") + 1 # Nbr of genes is the number of "/" + 1
    # compute gene ratio - divide Count by setSize
    x$GeneRatio <- x$Count / x$setSize
    # Resort df
    col = c(1:9,12,13,10,11)
    x <- x[,col]
    dplyr::arrange(x, desc(Count)) # arrange by Counts
  }
  ext.ls = lapply(ext.ls, geneRcount)
  
  # (3) Remove terms with less counts than minCount (DEFAULT = 10)
  ext.ls = lapply(ext.ls, function(x) { subset(x, Count >= minCount) })
  
  # (4) Remove empty list objects
  ext.ls <- ext.ls[which(lapply(ext.ls, nrow) != 0)]
  
  # (5) add a "Cluster" col to each df in list. 'Cluster' explains which sample the values correspond to.
  for(i in names(ext.ls)){
    message(paste("Adding Cluster column to:\t", gsub("_gse.*","",i)))
    ext.ls[[i]]$Cluster <- gsub("_gse.*","",i)
  }
  
  # (6) Select topmost significant results - specified by 'top' parameter
  # Only run this when 'top' is specified!
  if(!missing(top)){ 
    # select function
    topSelect <- function(x, topNrow = nrow(x), BY = "pvalue"){
      if(topNrow > nrow(x)) { N = nrow(x) } else { N = topNrow }
      if(BY == "NES") {
        x <- dplyr::arrange(x, desc(abs(NES)))
      } else if (BY == "pvalue") {
        x <- dplyr::arrange(x, pvalue)
      } else if (BY == "GeneRatio") {
        x <- dplyr::arrange(x, desc(GeneRatio))
      } else if (BY == "Count") {
        x <- dplyr::arrange(x, desc(Count))
      } else {
        stop(paste("topBy <-",BY,"not recognized!\n'topBy' must be one of:\t'NES','pvalue','GeneRatio','Count'"))
      }
      x[1:N,"ID"] # get only IDs!
    }
    # select topSet - union of the top N IDs across all df in list 
    topSet <- Reduce(union, lapply(ext.ls, topSelect, topNrow = top, BY = topBy))
    Reduce(union, lapply(ext.ls, topSelect, topNrow = top, BY = "Count"))
    
    # subset ext.ls
    ext.ls = lapply(ext.ls, function(x){ x[x$ID %in% topSet, ] })
  }
  
  ext.ls
}
mergeGSEAres  <- function(ls, filter= "p.adjust", pcut= .05, top, topBy= "pvalue", minCount= 10) {
  # ls    = list with GSEA result tables from clusterProfiler package
  # topBy = value to select top by. One of: "NES","pvalue","GeneRatio","Count"
  # top   = integer to select the n - top results from each object in ls
  
  # (1) Subset the df in lists for only sign. enriched terms. 'filterBy' variable will 
  # determine for which statistical value to filter for (one of: "pvalue", "p.adjust", "qvalues")
  pCutOff <- function(x, filterBy = "p.adjust", cut = .05){
    message(paste("Applied statistical cut off:\t",filterBy)," < ",cut)
    if(filterBy == "p.adjust"){ subset(x, p.adjust <= cut)
    } else if (filterBy == "pvalue") {
      subset(x, pvalue   <= cut)
    } else if(filterBy == "qvalues") {
      subset(x, qvalues  <= cut)
    } else {stop("'filter' parameter missing!\nPlease provide one of: 'pvalue', 'p.adjust', 'qvalues'")}
  }
  ext.ls = lapply(ls, pCutOff, filterBy = filter, cut = pcut)
  
  # (2) Compute gene ratio & counts (Nbr of enriched genes in cluster)
  geneRcount <- function(x){
    # compute counts - count strings separated with / in core_enrichment column
    x$Count <- stringr::str_count(x$core_enrichment, "/") + 1 # Nbr of genes is the number of "/" + 1
    # compute gene ratio - divide Count by setSize
    x$GeneRatio <- x$Count / x$setSize
    # Resort df
    col = c(1:9,12,13,10,11)
    x <- x[,col]
    dplyr::arrange(x, desc(Count)) # arrange by Counts
  }
  ext.ls = lapply(ext.ls, geneRcount)
  
  # (3) Remove terms with less counts than minCount (DEFAULT = 10)
  ext.ls = lapply(ext.ls, function(x) { subset(x, Count >= minCount) })
  
  # (4) Remove empty list objects
  ext.ls <- ext.ls[which(lapply(ext.ls, nrow) != 0)]
  
  # (5) add a "Cluster" col to each df in list. 'Cluster' explains which sample the values correspond to.
  for(i in names(ext.ls)){
    message(paste("Adding Cluster column to:\t", gsub("_gse.*","",i)))
    ext.ls[[i]]$Cluster <- gsub("_gse.*","",i)
  }
  
  # (6) Select topmost significant results - specified by 'top' parameter
  # Only run this when 'top' is specified!
  if(!missing(top)){ 
    # select function
    topSelect <- function(x, topNrow = nrow(x), BY = "pvalue"){
      if(topNrow > nrow(x)) { N = nrow(x) } else { N = topNrow }
      if(BY == "NES") {
        x <- dplyr::arrange(x, desc(abs(NES)))
      } else if (BY == "pvalue") {
        x <- dplyr::arrange(x, pvalue)
      } else if (BY == "GeneRatio") {
        x <- dplyr::arrange(x, desc(GeneRatio))
      } else if (BY == "Count") {
        x <- dplyr::arrange(x, desc(Count))
      } else {
        stop(paste("topBy <-",BY,"not recognized!\n'topBy' must be one of:\t'NES','pvalue','GeneRatio','Count'"))
      }
      x[1:N,"ID"] # get only IDs!
    }
    # select topSet - union of the top N IDs across all df in list 
    topSet <- Reduce(union, lapply(ext.ls, topSelect, topNrow = top, BY = topBy))
    Reduce(union, lapply(ext.ls, topSelect, topNrow = top, BY = "Count"))
    
    # subset ext.ls
    ext.ls = lapply(ext.ls, function(x){ x[x$ID %in% topSet, ] })
  }
  
  # (7) finally cbind the individual Lists together
  df = rlist::list.rbind(ext.ls)
  row.names(df) <- NULL
  df
}
multiGSEAplot <- function(ls, title = "", top, topBy= "pvalue", pcut= .05, ... ) {
  # run merge FUN
  res <- mergeGSEAres(ls, pcut=pcut, top=top, topBy=topBy, ...)
  # plot
  if(missing(top)) {top <- NULL}
  ggplot(res, aes(x= GeneRatio, y= Description, color= NES, size= log2(Count) )) + geom_point() +
    scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
    facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste("GSEA -",title,":",length(unique(res$ID))), 
         subtitle = paste("p.adj <",pcut,"/ Top:",top,"/ Top select by:",topBy)) + 
    aes(GeneRatio, reorder(stringr::str_wrap(Description, 80), GeneRatio)) + ylab(NULL)
}

## Function to select common enriched terms among treatments of a tested substance 
# for a list of clusterProfiler result tables
commonGS <- function(ls, unite = T) {
  int <- list()
  for(sub in unique(gsub("_.*","",names(ls)))) {
    tmp = ls[grepl(sub, names(ls))]
    # get common set of enriched terms
    tmp <- lapply(tmp, function(x) x$ID)
    enr <- Reduce(intersect, tmp)
    int[[sub]] <- enr
  }
  if(unite == T) { Reduce(union, int) } else { int }
}

# commonGS(filterGSEA(goBP.ls, top = 3))
#############################

setwd(paste0(home,"/compareGSEA"))
# plot top selection
gg = list()
## Top 5 ##
for(i in gseaLists) { gg[[i]] <- multiGSEAplot(get(i), title = gsub(".ls","",i), top = 5) }
pdf("multiGSEA_Top5.pdf", width = 20, height = 8, bg = "transparent")
gg
dev.off()
## Top 7 ##
for(i in gseaLists) { gg[[i]] <- multiGSEAplot(get(i), title = gsub(".ls","",i), top = 7) }
pdf("multiGSEA_Top7.pdf", width = 20.5, height = 10, bg = "transparent")
gg
dev.off()
## Top 10 ##
for(i in gseaLists) { gg[[i]] <- multiGSEAplot(get(i), title = gsub(".ls","",i), top = 10) }
pdf("multiGSEA_Top10.pdf", width = 20.5, height = 12, bg = "transparent")
gg
dev.off()
## ALL ##
for(i in gseaLists) { gg[[i]] <- multiGSEAplot(get(i), title = gsub(".ls","",i)) }
pdf("multiGSEA_ALL.pdf", width = 20, height = 20, bg = "transparent")
gg
dev.off()

## Overlapping set with top selection ##
# get list of overlaps with commonGS function
TOP = 5
gg = list()
for(i in gseaLists) {
  sel <- commonGS(filterGSEA(get(i), top = TOP))
  ls  <- lapply(get(i), function(x){ x[x$ID %in% sel,] })
  gg[[i]] <- multiGSEAplot(ls, top = TOP, title = gsub(".ls","",i))
}
pdf("multiGSEA_commonTop5.pdf", width = 16, height = 8, bg = "transparent")
gg
dev.off()

TOP = 3
gg = list()
for(i in gseaLists) {
  sel <- commonGS(filterGSEA(get(i), top = TOP))
  ls  <- lapply(get(i), function(x){ x[x$ID %in% sel,] })
  gg[[i]] <- multiGSEAplot(ls, top = TOP, title = gsub(".ls","",i))
}
pdf("multiGSEA_commonTop3.pdf", width = 16, height = 7, bg = "transparent")
gg
dev.off()

rm(i,gg,TOP)


## EXPORT RESULTS ## --------------------------
dir.create(home,"/compareGSEA/multiGSEAres", showWarnings = F)
multiGSEA.ls <- list()
for(i in gseaLists){
  df <- mergeGSEAres(get(i), top = 5)
  multiGSEA.ls[[i]] <- df
  write.csv2(df, file = paste0("multiGSEAres/multiGSEAres_top5_",gsub(".ls","",i),".csv"),row.names = F)
}

multiGSEA.ls <- list()
for(i in gseaLists){
  df <- mergeGSEAres(get(i))
  multiGSEA.ls[[i]] <- df
  write.csv2(df, file = paste0("multiGSEAres/multiGSEAres_",gsub(".ls","",i),".csv"),row.names = F)
}

rm(df,i)

## TOPMOST ## ----------------------
TOP  = 5 # selecting the n topmost sign. results from each condition
BY   = "NES" # NES pvalue Count GeneRatio
PCUT = .05

topGSEA.ls <- list()
for(i in gseaLists) { topGSEA.ls[[i]] <- mergeGSEAres(get(i), top = TOP, topBy = BY, pcut = PCUT) }
#rm(TOP,PCUT,BY,i)
# Nbr of unique terms in df
lapply(topGSEA.ls, function(x) {length(unique(x$ID))})

### Now let's try to plot this mess ... 
res = topGSEA.ls$reac.ls
#View(res)
table(res$Cluster)
length(unique(res$ID))
#write.csv2(res, file = "~/DATA/multiGSEA_reactome.csv", row.names = F) #For Daniel 
##################


title = "Reactome"
BY = "pvalue"
TOP  = 7
PCUT = .05

res = mergeGSEAres(reac.ls)
res = mergeGSEAres(reac.ls, top = TOP, topBy = BY, pcut = PCUT)
# plot
ggplot(res, aes(x= GeneRatio, y= Description, color= NES, size= log2(Count) )) + geom_point() +
  scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = paste("GSEA -",DB,":",length(unique(res$ID))), 
       subtitle = paste("p.adj <",PCUT,"/ Top:",TOP,"/ Top select by:",BY)) + 
  aes(GeneRatio, reorder(stringr::str_wrap(Description, 80), GeneRatio)) + ylab(NULL)

#ggplot(res, aes(x= Count, y= Description, color= NES, size= GeneRatio)) + geom_point() +
#  scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
#  facet_grid(~Cluster) + labs(title = paste("GSEA -",DB,":",length(unique(res$ID))),
#                              subtitle = paste("p.adj <",PCUT,"/ Top:",TOP,"/ Top select by:",BY)) + 
#  theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#  aes(Count, reorder(stringr::str_wrap(Description, 80),Count)) + ylab(NULL)

#ggplot(res, aes(x= GeneRatio, y= Description, color= NES, size= Count)) + geom_point() +
#  scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
#  facet_grid(~Cluster)+ labs(title = paste("GSEA -",DB,":",length(unique(res$ID))),
#                             subtitle = paste("p.adj <",PCUT,"/ Top:",TOP,"/ Top select by:",topBy)) + 
#  theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#  aes(GeneRatio, reorder(stringr::str_wrap(Description, 80), Count)) + ylab(NULL)

