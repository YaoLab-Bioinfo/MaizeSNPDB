
Args <- commandArgs(TRUE)
in.fl <- as.character(Args[1])

load(in.fl)

snp.data.allele <- apply(snp.data, 1, function(x){
  x <- x[!x %in% c("H", "N")]
  y <- sort(table(x), decreasing=TRUE)
  major <- names(y)[1]
  minor <- names(y)[2]
  return(c(major, minor, "H"))
})

snp.data.allele <- t(snp.data.allele)
colnames(snp.data.allele) <- c("major", "minor", "Het")

snp.data.inter <- apply(snp.data, 1, function(x){
  x.fil <- x[!x %in% c("H", "N")]
  y <- sort(table(x.fil), decreasing=TRUE)
  major <- names(y)[1]
  minor <- names(y)[2]
  z <- x
  z[z==major] <- 0
  z[z==minor] <- 1
  z[z=="H"] <- 2
  z[z=="N"] <- NA
  return(z)
})

snp.data.inter <- t(snp.data.inter)
mode(snp.data.inter) <- "integer"

library(Matrix)
snp.data.inter.Matrix <- as(snp.data.inter, "sparseMatrix")

out.fl <- gsub("snp.RData", "snp.Mat.RData", in.fl)
save(snp.data.allele, snp.data.inter.Matrix, compress="xz", file=out.fl)

