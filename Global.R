
options(warn=-1)

#library(IRanges)
#library(plotly)
#library(LDheatmap)
library(chopsticks)
#library(foreach)
#library(ape)
#library(pegas)
#library(plyr)
#library(dplyr)
#library(tidyr)
library(gridExtra)
#library(ggtree)
#library(grid)
library(snpStats)
#library(genetics)
library(shinycssloaders)
library(shinysky)

source("fetchSnp.R")
source("fetchSnpAllele.R")
source("ld.heatmap.R")
source("phylo.R")
source("nucDiv.R")
source("GBrowser.R")
source("anaReg.R")
source("geneStru.R")
source("snpInfo.R")
source("validReg.R")

acc.info <- read.table("./data/all.acc.txt", head=T, as.is=T, sep="\t", quote="")
#load("./data/gff.AGP.v4.RData")
#snp.lst <- read.table("./data/snp.RData.lst", head=T, as.is=T, sep="\t")
#load("./data/gene.info.RData")

source("chooser.R")
all.acc.cho <- paste(acc.info$ID, acc.info$Species, acc.info$Category, sep=", ")
all.acc.cho <- c("Improved", "Landrace", "Parviglumis", all.acc.cho)

chrInfo <- read.table("./data/chrInfo.txt", head=T, as.is=T, sep="\t")

acc.tree <- read.table("./data/acc.tree.txt", 
                       head=T, as.is=T, sep="\t", row.names = 1)

footerTagList <- list(
  tags$footer(id = "myFooter",
              shiny::includeHTML("www/footer.html")
  )
)

`%>%` <- magrittr::`%>%`
`%dopar%` <- foreach::`%dopar%`
