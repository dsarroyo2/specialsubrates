#Libraries, working directory, sources, and other functions
library (rensembl)
library (rrelrates) #In order for rrelrates to work (new version), I installed biomaRt and GO.db.

#set your working directory to whatever you need. The example below is how I did it.
setwd("/Users/dianaarroyo/Documents/DianasDocuments/ASU/CartwrightLab/HONORSTHESIS/Code/")

#toycatch, gsbincatch, and avoidcatch were other .R files used.
#They basically have functions that allow for a cleaner version of rprogress4.R.
source("toycatch.R") #toycatch.R contains the function getseqs2
source("gsbincatch.R")
source("avoidcatch.R")
source("raxtreessource.R")
source("extractbranchlengths.R")

coeffv <-function(x){
  (sqrt(x$var)/x$mean)*100
}

#----------------------------------------------------------------
#Genes
genes <- read.table('unique2.txt', stringsAsFactors=FALSE)
genes <- genes[,1]
#awesomegenes.tsv was made in rprogress2.R

#----------------------------------------------------------------
#Realignments
#lspecies can be altered so that it includes whatever species you wish to look at. 
#For this project, a few primates were chosen, so their species names were written into the script.
lspecies <- c('homo_sapiens','pan_troglodytes', 'pan_paniscus', 'gorilla_gorilla', 'gorilla beringei', 'pongo_abelii', 'macaca_mulatta', 'chlorocebus_aethiops', 'sapajus_apella')
#gs <- get_seqs("AGAP2", lspecies) #get seq takes a gene and then a list (of species).

#toycatch.R stuff
gs2 <- lapply(genes, getseq2)
save(gs2, file="gs2.rda")
#load("/Users/dianaarroyo/Documents/DianasDocuments/ASU/CartwrightLab/HONORSTHESIS/Code/gs2.rda")
#Once you have saved the file, you can always load it back instead of doing gs2 again.

#gsbincatch.R stuff
gsbin <- lapply(gs2, gs3)
save(gsbin, file="gsbin.rda")
#load("/Users/dianaarroyo/Documents/DianasDocuments/ASU/CartwrightLab/HONORSTHESIS/Code/gsbin.rda")
#Once you have saved the file, you can always load it back instead of doing gsbin again.

#----------------------------------------------------------------
#Build Trees
trees <-lapply(genes, gettree)
treesnull <- sapply(trees, is.null)
awesomegenes <- genes[!treesnull]
write.table(awesomegenes, file = "awesomegenes.tsv")

#avoidcatch.R stuff
alignments <- lapply(awesomegenes, alignit)
really_alignments <- lapply(alignments, function(x) if(is.null(x)) x else mafft(x))

trees <- lapply(really_alignments, raxtrees)

#----------------------------------------------------------------
#Using filtertrees to get rates & names
nulltrees <- sapply(trees, is.null)
filtertrees <- trees[!nulltrees]

finalrates <-lapply(filtertrees, geneinfo)
names(finalrates) <-awesomegenes[!nulltrees]
sff <- mapply(suminfo, finalrates, filtertrees, SIMPLIFY=FALSE)

#At this point, one can print the new trees.
#Used " plot(filtertrees[[9]]) " for ATP2B2.
#Same thing for FOXF1 with [[43]].
#Same thing for IGFBP5 with [[58]].
#----------------------------------------------------------------
#Coefficient of variation
sapply(sff, coeffv)

#----------------------------------------------------------------
#Printing the table out
nicetableyo <- t(sapply(sff, "[", c("mean", "var")))
coefficient <- sapply(sff, coeffv)
cbind(nicetableyo, coefficient)

branchl <- mapply(head, sff, -2)
exhibitA <- cbind(nicetableyo,coefficient, branchl)
exhibitA
#ExhibitB was made after testing this out again (and replacing 'A' in 'exhibitA' with 'B') just in case on March 8, 2016.

#----------------------------------------------------------------
#Check interesting genes
exhibitA[sapply(exhibitA[, 4], length) > 0, 4]

#DONE
#----------------------------------------------------------------
