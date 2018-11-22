
validReg <- function(myPos = list(NULL)) {
  if (myPos$chr %in% paste0("chr", 1:10) && !is.na(myPos$start) && 
      !is.na(myPos$end) && myPos$start>=1 && myPos$end>myPos$start && (myPos$end-myPos$start)<=2e6 ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

