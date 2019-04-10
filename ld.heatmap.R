
# A function to create the linakge disequilibrium heatmap using SNP data in a specified genomic region.
# Change to the directory of MaizeSNPDB using the setwd function of R.
# Usage: type the next two lines in R Console without the leading #
# source("Global.R")
# ld.plot <- ld.heatmap(chr="chr9", start=37800, end=40400, snp.pos=c(1), gene=FALSE, ld.y=0.64, ld.w=0.80, flip=FALSE, accession=NULL, mutType=NULL, snpSites = NULL)
# Then the linakge disequilibrium heatmap of this region would be displayed in a plotting device.
# For more info, please check the LDheatmap menu of the MaizeSNPDB database.

ld.heatmap <- function(chr="chr9", start=37800, end=40400, snp.pos=c(1), 
                       gene=FALSE, ld.y=0.64, ld.w=0.80, flip=FALSE, accession=NULL, 
                       mutType=NULL, snpSites = NULL, ...){
  start <- as.numeric(start)
  end <- as.numeric(end)
  
  dat <- fetchSnpAllele(chr=chr, start=start, end=end, accession=accession, mutType=mutType)
  
  dat.df <- data.frame(t(dat[[1]]), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(dat.df) <- 1:nrow(dat.df)
  dat.lst <- lapply(dat.df, function(x){genotype(x)})
  dat.snp.mat <- data.frame(dat.lst, check.names = FALSE)
  snp.code.pos <- as.numeric(substr(rownames(dat[[2]]), 3, 11))
  
  if (gene) {
    ll <- LDheatmap(dat.snp.mat, snp.code.pos, 
                    flip=TRUE, title=NULL, ...)
    
    p1 <- geneStru(chr=chr, start=start, end=end)
    
    plot.new()
    llQplot2 <- LDheatmap.addGrob(ll, rectGrob(gp=gpar(col="white")), height=.3)
    pushViewport(viewport(x=0.483, y=ld.y, width=ld.w, height=.1))
    
    grid.draw(ggplotGrob(p1))
  } else {
    if (flip) {
      LDheatmap(dat.snp.mat, snp.code.pos, flip=TRUE, title=NULL, ...)
    } else {
      LDheatmap(dat.snp.mat, snp.code.pos, flip=FALSE, SNP.name = colnames(dat.snp.mat)[snp.pos], title=NULL, ...)
    }
  }
  
}

