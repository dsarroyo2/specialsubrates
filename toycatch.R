getseq2 <- function(gene) {
  out <- tryCatch(
  get_seqs(gene, lspecies)
  ,
  error=function(cond) {
    return(NA)
  },
  warning=function(cond) {
    return(NULL)
  },
  finally={
    message(paste("Processed gene:", gene))
  }
  )
  return(out)
}


#gs2 <- lapply(genes[1:3], getseq2)

