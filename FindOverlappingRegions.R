setwd("/g/scb2/zaugg/alishaye")
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(readr)
library(ChIPpeakAnno)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###Prepare 2 gene lists in Granges format

#select significantly DE peaks (441) in the comparison of 3 groups where both genders are included, but only 
#autosoumal chromosomes are considered
#select_005 = order(res$padj)[1:length(which(res$padj < 0.05))]
#df_select_005<-coords[select_005,]
#gr_select_005=makeGRangesFromDataFrame(df_select_005)
#saveRDS(gr_select_005,"441_sign_peaks_granges")

#load significantly DE peaks (441)in the comparison of 3 groups where both genders are included, but only 
#autosoumal chromosomes are considered
both_sexes_sign_peakss<-read_rds("441_sign_peaks_granges")
#load significantly DE peaks (629) from the comparison of 3 groups where only female samples were analysed
fem_sign_genes<-read_rds("629_sign_peaks_granges")


###Find overlapping peaks with Granges function
findOverlaps(both_sexes_sign_peakss,fem_sign_genes,maxgap = 5000)
#Total:17,changeing maxgap or minoverlap doesn't make a difference

#annotate and find overlaps between the gene names
both_genders <- annotatePeak(both_sexes_sign_peakss, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoList_fem <- annotatePeak(fem_sign_genes, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
df_both_genders<-as.data.frame(both_genders)
dfpeakAnnoList_fem<-as.data.frame(peakAnnoList_fem)

###Find overlapping gene names
unique(dfpeakAnnoList_fem$SYMBOL)
intersect(df_both_genders$SYMBOL,dfpeakAnnoList_fem$SYMBOL)
#Total:54. Some peaks are assigned to the same gene (but different region - promoter, exonic ..etc)

