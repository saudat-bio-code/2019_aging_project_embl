############GREAT ANALYSIS############
############Saudat 9 July#############
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rGREAT")
library(DiffBind)
library(rGREAT)
#load and prepare data
both_sexes_sign_peakss<-read_rds("441_sign_peaks_granges")
seqlevelsStyle(both_sexes_sign_peaks) <- "UCSC"
seqlevelsStyle(cons_peaks_granges) <- "UCSC"
granges_res=makeGRangesFromDataFrame(coords,keep.extra.columns = T)

df_ann_results<-as.data.frame(annotated_res)
df_ann_DE_results<-as.data.frame(annotated_DE_res)

ann_gr<-makeGRangesFromDataFrame(df_ann_results)
ans_DE_gr<-makeGRangesFromDataFrame(df_ann_DE_results)

#send the jobs
job = submitGreatJob(ans_DE_gr,bg=ann_gr,species = "hg19")
enrichment_statistics = getEnrichmentTables(job, ontology = c("GO Molecular Function"))
availableOntologies(job)
head(enrichment_statistics[[1]])
par(mfrow = c(1, 2))
res = plotRegionGeneAssociationGraphs(job)
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job, ontology = "GO Molecular Function",termID = "GO:0060589")
res[1:5, ]


pecies = "hg19", rule = "twoClosest", request_interval = 0, adv_twoDistance = 1000)
tb = rGREAT::getEnrichmentTables(job)


