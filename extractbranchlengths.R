extract_branchlengths <- function(suminfo_result){
  if(length(suminfo_result) > 2){
    res <- suminfo_result[1:(length(suminfo_result)-2)]
    #THIS IS A HORRIBLE HACK THAT DAVID DID. FIX SUMINFO FXN SO
    #THAT THE NAMES CANT BE NA. SORRY. -DAVID.
    old_names <- names(res)
    names(res) <- ifelse(is.na(old_names), "nope_nope_nope", old_names) 
    names(res) <- sapply(str_split(names(res), "_"), "[[", 3)
    return(res)
  }
  return(list())
}