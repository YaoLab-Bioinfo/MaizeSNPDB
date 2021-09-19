

anaReg <- function(x=NULL) {
  if (exists("gene.info")){
  }else{
    load("./data/gene.info.RData")
  }
  x <- gsub("\\s+", "", x)
  if (grepl("chr", x)) {
    myChr <- gsub(":.+", "", x)
    myPos <- as.numeric(gsub("\\s","", strsplit(gsub(".+:", "", x),"-")[[1]]))
  } else {
    myChr <- gene.info$chr[gene.info$id==x]
    myPos <- c(gene.info$start[gene.info$id==x], gene.info$end[gene.info$id==x])
  }
  
  chr.size <- c(307041717, 244442276, 235667834, 246994605, 223902240, 
                174033170, 182381542, 181122637, 159769782, 150982314)
  names(chr.size) <- paste0("chr", 1:10)
  
  myPos[1] <- max(1, myPos[1])
  myPos[2] <- min(myPos[2], chr.size[myChr])
  
  return(list(chr=myChr, start=myPos[1], end=myPos[2]))
}

