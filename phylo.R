
# A function to perform phylogenetic analysis using SNP data in a specified genomic region.
# Change to the directory of MaizeSNPDB using the setwd function of R.
# Usage: type the next three lines in R Console without the leading #
# source("Global.R")
# phy.plot <- phylo(chr="chr9", start=37800, end=41400, accession=NULL, mutType=NULL, snpSites = NULL)
# print(phy.plot)
# Then the NJ tree would be displayed in a plotting device.
# For more info, please check the Phylogenetic menu of the MaizeSNPDB database.

phylo <- function(chr="chr9", start=37800, end=41400, accession=NULL, mutType=NULL, snpSites = NULL) {
  if (exists("snp.lst")){
  }else{
    snp.lst <- read.table("./data/snp.RData.lst", head=T, as.is=T, sep="\t")
  }
  start <- as.numeric(start)
  end <- as.numeric(end)
  reg.gr <- IRanges::IRanges(start, end)
  snp.lst.chr <- snp.lst[snp.lst$chr==chr, ]
  snp.lst.gr <- IRanges::IRanges(start=snp.lst.chr$start, end=snp.lst.chr$end)
  snp.fls <- snp.lst.chr$file[unique(S4Vectors::queryHits(IRanges::findOverlaps(snp.lst.gr, reg.gr)))]
  
  snp.data.lst <- lapply(snp.fls, function(x){
    load(x)
    return(snp.data.inter.Matrix)
  })
  snp.data <- do.call(rbind, snp.data.lst)
  snp.data <- snp.data[order(as.numeric(rownames(snp.data))), ]
  colnames(snp.data) <- acc.info$ID
  
  snpeff.fls <- gsub("snp", "snpeff", snp.fls)
  snpeff.fls.lst <- lapply(snpeff.fls, function(x){
    load(x)
    return(snpeff)
  })
  snpeff <- do.call(rbind, snpeff.fls.lst)
  snpeff <- snpeff[order(as.numeric(snpeff[, 1])), ]
  
  start <- as.numeric(paste0(sprintf("%02d", as.numeric(substr(chr, 4, 4))), sprintf("%09d", start)))
  end <- as.numeric(paste0(sprintf("%02d", as.numeric(substr(chr, 4, 4))), sprintf("%09d", end)))
  
  dat.res <- snp.data[as.numeric(rownames(snp.data))>=start & as.numeric(rownames(snp.data))<=end, , drop=FALSE]
  dat.res <- as.matrix(dat.res)
  
  accession <- gsub(",.+", "", accession)
  accession <- sapply(accession, function(x){
    if (x %in% c("Improved", "Landrace", "Parviglumis")) {
      x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
      return(x.dat)
    } else {
      return(x)
    }
  })
  accession <- unique(unlist(accession))
  
  if (!is.null(accession) && length(accession)>=2) {
    dat.res <- dat.res[, colnames(dat.res) %in% accession, drop=FALSE]
  }
  
  dat.res.row.c <- apply(dat.res, 1, function(x){
    length(unique(x[!is.na(x)]))
  })
  dat.res <- dat.res[dat.res.row.c>1, , drop=FALSE]
  
  if (!is.null(mutType) && length(mutType)>=1 && length(mutType)!=16) {
    snpeff.info <- snpeff[snpeff[, 1] %in% rownames(dat.res),]
    
    snpeff.info[,"eff"][grepl("IT", snpeff.info[,"eff"])] <- "Intergenic"
    snpeff.info[,"eff"][grepl("IR", snpeff.info[,"eff"])] <- "Intron"
    snpeff.info[,"eff"][grepl("IG", snpeff.info[,"eff"])] <- "Start_gained"
    snpeff.info[,"eff"][grepl("IL", snpeff.info[,"eff"])] <- "Start_lost"
    snpeff.info[,"eff"][grepl("SG", snpeff.info[,"eff"])] <- "Stop_gained"
    snpeff.info[,"eff"][grepl("SL", snpeff.info[,"eff"])] <- "Stop_lost"
    snpeff.info[,"eff"][grepl("Up", snpeff.info[,"eff"])] <- "Upstream"
    snpeff.info[,"eff"][grepl("Dn", snpeff.info[,"eff"])] <- "Downstream"
    snpeff.info[,"eff"][grepl("U3", snpeff.info[,"eff"])] <- "three_prime_UTR"
    snpeff.info[,"eff"][grepl("U5", snpeff.info[,"eff"])] <- "five_prime_UTR"
    snpeff.info[,"eff"][grepl("SSA", snpeff.info[,"eff"])] <- "Splice_site_acceptor"
    snpeff.info[,"eff"][grepl("SSD", snpeff.info[,"eff"])] <- "Splice_site_donor"
    snpeff.info[,"eff"][grepl("NSC", snpeff.info[,"eff"])] <- "Non_synonymous_coding"
    snpeff.info[,"eff"][grepl("NSS", snpeff.info[,"eff"])] <- "Non_synonymous_start"
    snpeff.info[,"eff"][grepl("SC", snpeff.info[,"eff"])] <- "Synonymous_coding"
    snpeff.info[,"eff"][grepl("SS", snpeff.info[,"eff"])] <- "Synonymous_stop"
    snpeff.info[,"eff"][grepl("IA", snpeff.info[,"eff"])] <- "Intergenic"
    
    snpeff.info <- snpeff.info[snpeff.info[, "eff"] %in% mutType, , drop=FALSE]
    
    dat.res <- dat.res[rownames(dat.res) %in% snpeff.info[, "id"], , drop=FALSE]
  }
  
  if (!is.null(snpSites) && length(snpSites)>=1) {
    dat.res <- dat.res[rownames(dat.res) %in% snpSites, , drop=FALSE]
  }
  
  dat.res[dat.res == 2] <- 3
  dat.res[dat.res == 1] <- 2
  dat.res[dat.res == 3] <- 1
  
  #### calculate distance matrix
  "%dis%" <- function(x,y){
    return(abs(x-y)/2 + as.numeric((x==1)&(y==1))/2)
  }
  
  dist.mat.nume <- foreach::foreach(x=1:ncol(dat.res),.combine=rbind)%dopar%{colSums(dat.res%dis%dat.res[,x],na.rm=TRUE)}
  dat.res[!is.na(dat.res)] <- 1
  dist.mat.deno <- foreach::foreach(x=1:ncol(dat.res),.combine=rbind)%dopar%{colSums(!is.na(dat.res+dat.res[,x]))}
  
  dist.mat <- dist.mat.nume/dist.mat.deno
  rownames(dist.mat) <- colnames(dist.mat)
  
  ### tree
  dist.mat <- as.dist(dist.mat)
  tre <- ape::nj(dist.mat)
  
  p <- ggtree::ggtree(tre, layout="circular", branch.length="none", size=0.01) + ggplot2::ggtitle("")
  p <- p + ggplot2::theme_void()
  p <- ggtree::gheatmap(p, acc.tree, offset = 1, width=0.1, colnames = FALSE, color=NULL) +
    ggplot2::scale_fill_manual(breaks=c("Improved", "Landrace", "unknown", 
                               "Parviglumis", "Other"), 
                      values=c("blue", "red", "black", 
                               "purple", "gold"))
  figurecp <<- p
  treNwk <<- tre
  return(p)
}

