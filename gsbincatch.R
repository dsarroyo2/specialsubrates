gs3 <- function(binned) {
  out <- tryCatch(
    as.DNAbin.list(binned)
    ,
    error=function(cond) {
      return(NA)
    },
    warning=function(cond) {
      return(NULL)
    },
    finally={
      message(paste("Processing..."))
    }
  )
  realign <- tryCatch(
    mafft(out)
    ,
    error=function(cond) {
      return(NA)
    },
    warning=function(cond) {
      return(NULL)
    },
    finally={
      message(paste("Processing..."))
    }
  )
  return(realign)
}

