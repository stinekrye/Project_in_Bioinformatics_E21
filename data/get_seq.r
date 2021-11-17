#Script: get_seq.r
#Author: Ioannis Kampouris
#Script Task: 1. Assign ASVs names in the sequences and use the ASV names. 
#             2. prepare a phylogenetic tree for downstream analysis.
#Save output for dealing later in downstream analysis. 

library (dplyr)
library (tidyr)
library(stringr)
library(readr)
library(reshape2)
library("tibble")
library(dada2)
library(readxl)

setwd("/home/ioanniskampouris/Desktop/AArhus/combineallseqs")
load("seqtab_taxa.species_overlappingSeqRemoved.Rdata")

## Checking for any inconsistency between the samples.

samples_test = as.data.frame( seqtab.nochim) 
samples_test = rownames_to_column(samples_test, var = "Sample_Name_Seq"  )
samples_test=select(samples_test, Sample_Name_Seq)
map1 <- read_excel("sampleData.xlsx",
                         sheet = "link to R data")

samples_test = full_join(samples_test, map1, by = "Sample_Name_Seq" )


## Relabeling the sequences as ASVs. 
samples_df = as.data.frame( seqtab.nochim)
seq_df = as.data.frame( taxa.species)
seq_df = rownames_to_column(seq_df, var = "Sequence")
sequences = select(seq_df, Sequence)
sequences = rownames_to_column(sequences, var ="ASVs")
sequences$ASVs <- sub("^", "ASV", sequences$ASVs )
#Preparing a file for the tree generation 
seqtreeprep = sequences
seqtreeprep$ASVs <- sub("^", ">", seqtreeprep$ASVs )


write.table(seqtreeprep, file="forfasttree.fasta",  sep = "\n", row.names = FALSE, col.names = FALSE, quote = F)

## Append ASV numbering for each sequence
seq_df= full_join(seq_df, sequences, by="Sequence") %>% column_to_rownames(., var="ASVs") %>% select(., -Sequence)

#### Moving to processing the ASV table.

samples_df= rownames_to_column(samples_df, var = "Sample_Name_Seq")
ASV_df = gather(samples_df, key = "Sequence", value = "Abundance", -Sample_Name_Seq)
ASV_df = full_join(ASV_df, sequences, by="Sequence") %>%  select(., -Sequence)  %>%  dcast(., ASVs~Sample_Name_Seq, value.var = "Abundance", sum)
ASV_df = as.data.frame(ASV_df)

write.table(ASV_df, file="oil_ASVtable.txt", sep = "\t", row.names = F, quote = F)
write.table(seq_df, file="oil_ASVtaxid.txt", sep = "\t", quote = F)


