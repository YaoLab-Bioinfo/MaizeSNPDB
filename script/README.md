
This R script **chara-2-sparse.R** is used to convert a traditional **character genotype matrix** into an **integer sparse genotype matrix**.  
This R script can be used in the system command line with R installed in the system or be used in the R enviroment.  

**test.snp.RData** is the input data to this R script. test.snp.RData is an R data file contained the genotype data of 1210 maize lines across 1000 SNP sites in character matrix. The content of this R data can be viewed in the R enviroment as shown below.  
``` R code
ls()
# [1] "snp.data"
dim(snp.data)
# [1] 1000 1210
snp.data[1:5, 1:5]
#            V5  V6  V7  V8  V9 
# 07180000005 "T" "T" "T" "T" "T"
# 07180000006 "A" "A" "C" "A" "C"
# 07180000011 "H" "A" "C" "C" "C"
# 07180000016 "C" "C" "C" "C" "C"
# 07180000017 "H" "A" "A" "G" "A"
```

We then run this R script using **test.snp.RData** as the input data using the following command in the system command line.  
**R --slave --args test.snp.RData <chara-2-sparse.R**

After running the R script, a new R data file named as **test.snp.Mat.RData** would be generated in the same directory. This R data file contained the genotype data of 1210 maize lines across 1000 SNP sites in integer sparse matrix. The content of this R data can be viewed in the R enviroment as shown below.  
``` R code
ls()
# [1] "snp.data.allele"       "snp.data.inter.Matrix"
head(snp.data.allele)
#             major minor Het
# 07180000005 "T"   "A"   "H"
# 07180000006 "A"   "C"   "H"
# 07180000011 "C"   "A"   "H"
# 07180000016 "C"   "T"   "H"
# 07180000017 "A"   "G"   "H"
# 07180000048 "A"   "T"   "H"
library(Matrix)
snp.data.inter.Matrix[1:5, 1:5]
# 5 x 5 sparse Matrix of class "dgCMatrix"
#             V5 V6 V7 V8 V9
# 07180000005  .  .  .  .  .
# 07180000006  .  .  1  .  1
# 07180000011  2  1  .  .  .
# 07180000016  .  .  .  .  .
# 07180000017  2  .  .  1  .
```

