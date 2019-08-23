#Biomart annotations
#Goal - input is list with HUGO annotated files, I need  to annotate them to ENSEMBL, and general symbol
library("biomaRt")
setwd("/g/scb2/zaugg/alishaye")
listMarts()
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

###load gene file
remap_names<-read.table("names_remap.tsv")
remap_vector<-as.character(remap_names$V1)

#attributes: is a vector of the output of the query
#filters: is a vector of filters that one wil use as input to the query.
#values: a vector of values for the input 
#mart: is an object of class Mart, which is created by the useMart() function.
length(remap_vector)

biomart_annotated<-getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
      filters = 'hgnc_symbol', 
      values = remap_vector, 
      mart = ensembl)
biomart_annotated[duplicated(biomart_annotated$hgnc_symbol), ]
#10 files that seem to be duplicated
redundant_files<-remap_names[!(remap_names$V1 %in% biomart_annotated$hgnc_symbol), ]
not_redundant_files<-remap_names[(remap_names$V1 %in% biomart_annotated$hgnc_symbol), ]
#manually check these 10 files with the website 

names(biomart_annotated)[2] <- "ENSEMBL"
names(biomart_annotated)[1] <- "HUGO"
write_delim(biomart_annotated,"remap_hgnc_to_ensembl.csv", delim =" ")



#as an input for the GRN thing the table needs to have 3 cols:
#SYMBOL
#ENSEMBLE
#HOCOID -> the file name will be searched by that id!


