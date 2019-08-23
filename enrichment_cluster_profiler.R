##### Saudat - 2 July 2019 #####
setwd("/g/scb2/zaugg/alishaye/Robjects/females")
library(clusterProfiler)
library(genomation)
library(readr)
library(org.Hs.eg.db)
library(data.table)
library(GenomicRanges)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
OrgDb <- org.Hs.eg.db
library(ChIPseeker)
library(DOSE)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(both_sexes_sign_peaks) <- "UCSC"


#annotate and find overlaps between the gene names
#coords_autosomal<-coords[-grep(coords$CHR,pattern = "chrX"), ]
#coords_autosomal<-coords_autosomal[-grep(coords_autosomal$CHR,pattern = "chrY"), ]
#results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.5, "FDR<0.05", "Not Sig")), row.names=rownames(res))
#granges_res=makeGRangesFromDataFrame(coords_autosomal,keep.extra.columns = T)
#DE_res = order(res$padj)[1:length(which(res$padj < 0.1))]
#DE_res=readRDS("fem_upd_outl_3grps_758_peaks")
#DE_res=coords_autosomal[DE_res,]
#DE_res_granges=makeGRangesFromDataFrame(DE_res, keep.extra.columns = T)

DE_res_granges=readRDS("fem_upd_outl_3grps_758_peaks")
granges_res=readRDS("fem_upd_outl_3grps_ALL_peaks_background")
results=readRDS("fem_upd_outl_3grps_deseq_results_dataframe")
#annot <- annotatePeakInBatch(DE_res_granges, 
#                             AnnotationData=annoData, 
#                             featureType = "TSS",
#                             output="nearestLocation",multiple = FALSE,
#                             PeakLocForDistance = "start")

#annot_res <- annotatePeakInBatch(granges_res, 
 #                                AnnotationData=annoData, 
#                                 featureType = "TSS",
#                                 output="nearestLocation",
#                                 PeakLocForDistance = "start")
annotated_res <- annotatePeak(granges_res, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
annotated_DE_res <- annotatePeak(DE_res_granges, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

df_ann_results<-as.data.frame(annotated_res)
df_ann_DE_results<-as.data.frame(annotated_DE_res)

### gene names
df_ann_results$SYMBOL
df_ann_DE_results$SYMBOL

#somehow you need to add results with p values to the df_ann_results
rownames(df_ann_results)=df_ann_results$id
df_ann_results$id=NULL
head(results)
merged_annotat_results <- merge(results,df_ann_results, by= 0, all=TRUE) #by=0 means it will merge by raw names
#Prepare an input file with ranked genes
geneList <- as.vector(merged_annotat_results$log2FoldChange)
names(geneList) <- as.vector(merged_annotat_results$geneId) #geneId stands for entrezgene ID 
geneList <- sort(geneList, decreasing = T)
gene <- df_ann_DE_results$geneId

###############################################    GO   ###############################################


#GO over-represintation test
ego <- enrichGO(gene          = unique(gene),
                universe      = names(unique(geneList)),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH")
head(ego)
emapplot(ego)
#Gene Set Enrichment Analysis GSEA
ego3 <- gseGO(geneList     = geneList[unique(names(geneList))],
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000)
head(ego3)
emapplot(ego3)

###############################################  KEGG  ###############################################
#kegg over-repres
kk <- enrichKEGG(gene         = unique(gene),
                 universe      = names(unique(geneList)),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
#kegg gsea
kk2 <- gseKEGG(geneList     = unique(geneList),
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
emapplot(kk2)

#disease enrichment
disEnr <- enrichDO(gene    = unique(gene),
              ont           = "DO",
              pAdjustMethod = "BH",
              universe      = names(geneList))
head(disEnr)

