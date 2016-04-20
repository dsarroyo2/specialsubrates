raxtrees <- function(alignment) {
  if(is.null(alignment)){
    return(NULL)
  }
  p <- sample(1:1000, 1)
  x <- sample(1:1000, 1)
{
    out <- tryCatch(
      raxml(alignment, m = "GTRGAMMA", f = "a", N = 10, p = p, x = x, exec = "/Users/dianaarroyo/Documents/DianasDocuments/ASU/CartwrightLab/HONORSTHESIS/Code/raxmlHPC-AVX-v8/raxml", threads=2)$bestTree
      ,
      error=function(cond) {
        return(NULL)
      }
    )
    return(out)
  }
}
