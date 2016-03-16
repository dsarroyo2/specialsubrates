#David showed me this to avoid trycatch.
is_primate_orth <- function(x){
  if(x$target$species %in% lspecies){
    return( x$method_link_type ==  "ENSEMBL_ORTHOLOGUES")
  }
  FALSE
}

# raw_ali <- mapply(rensembl::homology_symbol, genes, "hsap", format="json", sequence="cdna")
# primate_seqs <- Filter(is_primate_orth, raw_ali[[1]][[1]][[1]])
# 
# seqs <- lapply(primate_seqs, function(x) x$target$align_seq)
# names(seqs) <- sapply(primate_seqs, function(x) x$target$species)
# seqs$homo_sapiens <- primate_seqs[[1]]$source$align_seq
# ALI <- as.DNAbin(lapply(seqs, function(s) tolower(strsplit(s, "")[[1]])))
# write.dna(ALI, "areg_cdnas_aligned.fasta", "fasta")


alignit <- function(gene, write_to_file=FALSE){
    raw_ali <- tryCatch(
      rensembl::homology_symbol("hsap", gene,format="json",  sequence="cdna")
      ,
      error=function(cond) {
        return(NA)
      }
    )
  if(is.na(raw_ali)){
    warning(paste(gene, "has no homologs."))
    return(NA)
  }
  primate_seqs <- Filter(is_primate_orth, raw_ali[[1]][[1]][[1]])
  if (length(primate_seqs)==0){
    warning(paste(gene, "has no primates."))
    return(NULL)
  }

  seqs <- lapply(primate_seqs, function(x) x$target$align_seq)
  names(seqs) <- sapply(primate_seqs, function(x) x$target$species)
  seqs$homo_sapiens <- primate_seqs[[1]]$source$align_seq
  ALI <- as.DNAbin(lapply(seqs, function(s) tolower(strsplit(s, "")[[1]])))
  if (write_to_file){
    write.dna(ALI, paste0(gene, ".fasta"), format="fasta")
  }
  return(ALI)
}

#alignit(genes[1])

#see rprogress3.R to continue
#alignments <- lapply(awesomegenes, alignit)
#really_alignments <- lapply(alignments, function(x) if(is.null(x)) x else mafft(x))
#^Use lapply to use new function to use mafft on alignments except on null.

#Back to rprogress3.R

# gg <- genes[!nulltrees]
# alignments <- lapply(gg, alignit, write_to_file=TRUE)
