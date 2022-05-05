#Preprocessing the files

library(shiny)
library(ggplot2)
library(colourpicker)
library(DT)
library('tidyverse')
library(ggbeeswarm)
library("ggplot2")
library("gridExtra") 
library('tidyverse')
library(ggbeeswarm)
library("ggplot2")
library("gridExtra")  
library(DESeq2)
options(shiny.maxRequestSize=30*1024^2) 
#Normalization of the counts
counts<-read_csv("Project/GSE150450_gene_count_matrix.csv")
nonzero_genes <- rowSums(counts[-1])!=0
filtered_counts <- counts[nonzero_genes,]

# DESeq2 requires a counts matrix, column data (sample information), and a formula
# the counts matrix *must be raw counts*
count_mat <- as.matrix(filtered_counts[-1])

row.names(count_mat) <- filtered_counts$gene_id

dds <- DESeqDataSetFromMatrix(
  countData=count_mat,
  colData=tibble(sample_name=colnames(filtered_counts[-1])),
  design=~1 # no formula is needed for normalization, ~1 produces a trivial design matrix
)

# compute normalization factors
dds <- estimateSizeFactors(dds)

# extract the normalized counts
deseq_norm_counts <- as_tibble(counts(dds,normalized=TRUE)) %>%
  mutate(gene=filtered_counts$gene_id) %>%
  relocate(gene)

#norm_counts<-write_csv(deseq_norm_counts, "Project/Normcounts.csv")

#Differential expression using DeSeq2
#reading and filtering out the raw counts which have just zeroes in all rows and converting to matrix
counts<-read.csv("Project/GSE150450_gene_count_matrix.csv")
filtered_counts <- counts[rowSums(counts[-1])>0,]
count_mat <- as.matrix(filtered_counts[-1])
row.names(count_mat) <- filtered_counts$gene

#coldata 
coldata<-read.csv("Project/Coldata.csv")
#design
design <- formula(~ sex+lifestage+treatment)

#diff expression
dds <- DESeqDataSetFromMatrix(
  countData=count_mat,
  colData=coldata,
  design=design
)
dds <- DESeq(dds)
#resultsNames(dds)
#write.csv(results(dds),"Project/differential_expression.csv")

#pre processing for FGSEA

#mapping ids for drosophilla
drosophila <- useMart('ENSEMBL_MART_ENSEMBL', dataset='dmelanogaster_gene_ensembl')
map <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),mart = drosophila)
#write_csv(map,"Project/map_dros.csv")

#gmt pathways
pathways <- gmtPathways("Project/Drosophila_melanogaster_fruit-fly_gmt2.gmt")

y<-read_csv("Project/differential_expression.csv")
gene_ident <- y %>% left_join(map, by=c('Gene' = 'ensembl_gene_id'))

rnks <- gene_ident %>% 
  arrange(desc(log2FoldChange)) %>% 
  dplyr::select(external_gene_name, log2FoldChange) %>% 
  deframe()
fgsea <- fgsea(pathways, rnks, minSize= 15, maxSize= 500) %>% as_tibble()
#write_csv(fgsea, "Project/fgsea.csv")

