#################Enrichment Heatmap##################
#####################Saudat##########################
library(GenomicRanges)
library(data.table)
library(EnrichedHeatmap)
library(circlize)
library(data.table)
library(ggplot2)
library(GenomicFeatures)
library(annotatr)
library(readr)
library(DiffBind)
ATAC_D = read_rds("ATAC_D")
ATAC_D<-ATAC_D[-c(2,31,5,12), ]
ATAC_D$Condition[ATAC_D$age<55]<-"Young"
ATAC_D$Condition[ATAC_D$age>55]<-"Old"
ATAC_young= ATAC_D[ATAC_D$Condition=="Young", ]
ATAC_old = ATAC_D[ATAC_D$Condition=="Old", ]
dba_young= dba(sampleSheet = ATAC_young)
dba_old= dba(sampleSheet = ATAC_old)
consensus_young = dba.count(DBA = dba_young, summits = 250, score = DBA_SCORE_READS, minOverlap = 2)
consensus_old = dba.count(DBA = dba_old, summits = 250, score = DBA_SCORE_READS, minOverlap = 2)
dds = dba.peakset(consensus_young, bRetrieve = T, DataType = DBA_DATA_FRAME)
coords = as.data.frame(dds[,1:3], stringsAsFactors = F)
gr_young=makeGRangesFromDataFrame(coords)
dds = dba.peakset(consensus_old, bRetrieve = T, DataType = DBA_DATA_FRAME)
coords= as.data.frame(dds[,1:3], stringsAsFactors = F)
gr_old=makeGRangesFromDataFrame(coords)


states_bed <- as.data.frame(read.table("E035_15_coreMarks_dense.bed.bgz",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
states = GRanges(seqnames = states_bed[[1]], 
                 ranges = IRanges(states_bed[[2]] + 1, states_bed[[3]]), 
                 states = states_bed[[4]])
unique(states_bed[[4]])


#Install genomation if the package is not found 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("genomation")

E035<-read_annotations("E035_15_coreMarks_dense.bed.bgz", genome = 'hg19',format ='bed')

###Prepare CHROMHMM annotation file for the CD34 positive mobilized HSC (download in bed format)

bed <- as.data.frame(read.table("E035_15_coreMarks_dense.bed.bgz",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(bed)[names(bed)=="V1"]<-"chr"
names(bed)[names(bed)=="V2"]<-"start"
names(bed)[names(bed)=="V3"]<-"end"
names(bed)[names(bed)=="V4"]<-"id"

extraCols= c(tx_id = "character",
             V5= "numeric",
             V6= "character",
             V7= "numeric",
             V8= "numeric",
             V9= "character")
r_an<-read_annotations(con = "E035_15_coreMarks_dense.bed.bgz", name = 'test',format = 'bed', genome = 'hg19',extraCols= extraCols)
build_annotations(genome = 'hg19', annotations = 'hg19_custom_test')


map = c(
  "1_TssA"      = "TssActive",
  "2_TssAFlnk"  = "TssActive",
  "3_TxFlnk"    = "Transcript",
  "4_Tx"        = "Transcript",
  "5_TxWk"      = "Transcript",
  "6_EnhG"      = "Enhancer",
  "7_Enh"       = "Enhancer",
  "8_ZNF/Rpts"  = "Heterochromatin",
  "9_Het"       = "Heterochromatin",
  "10_TssBiv"   = "TssBivalent",
  "11_BivFlnk"  = "TssBivalent",
  "12_EnhBiv"   = "Enhancer",
  "13_ReprPC"   = "Repressive",
  "14_ReprPCWk" = "Repressive",
  "15_Quies"    = "Quiescent"
)
states$states_simplified = map[states$states]

states_col = c(
  "TssActive"       = "Red",
  "Transcript"      = "Green",
  "Enhancer"        = "Yellow",
  "Heterochromatin" = "PaleTurquoise",
  "TssBivalent"     = "Orange",
  "Repressive"      = "Grey",
  "Quiescent"       = "black"
)
states_name = names(states_col)
n_states = length(states_col)
txdb= makeTxDbFromGFF(file = "gen10.long.gtf",organism = "Homo sapiens")
g = genes(txdb)
tss = promoters(g, upstream = 0, downstream = 1)
# column "states_simplified" is in character mode
mat_states = normalizeToMatrix(states_united, tss, value_column = "states_simplified")
mat_states = normalizeToMatrix(states_young, tss, value_column = "states_simplified")
mat_states = normalizeToMatrix(states_old, tss, value_column = "states_simplified")

mat_states
EnrichedHeatmap(mat_states, name = "states", col = states_col)dev.

###############################################      benign   ###############################################


#sample 1
dm_annotated = annotate_regions(
  regions = gr_young,
  annotations = r_an,
  ignore.strand = TRUE,
  quiet = FALSE)
dm_annotated_form<-as.data.frame(dm_annotated)
states_young = GRanges(seqnames = dm_annotated_form[[1]], 
                        ranges = IRanges(dm_annotated_form[[2]] + 1, dm_annotated_form[[3]]), 
                        states = dm_annotated_form[[12]])
states_young$states_simplified = map[states_young$states]
states_young$states_simplified = factor(states_young$states_simplified, levels = states_name)

#sample 2
dm_annotated = annotate_regions(
  regions = gr_old,
  annotations = r_an,
  ignore.strand = TRUE,
  quiet = FALSE)
dm_annotated_form<-as.data.frame(dm_annotated)
states_old = GRanges(seqnames = dm_annotated_form[[1]], 
                     ranges = IRanges(dm_annotated_form[[2]] + 1, dm_annotated_form[[3]]), 
                     states = dm_annotated_form[[12]])
states_old$states_simplified = map[states_old$states]
states_old$states_simplified = factor(states_old$states_simplified, levels = states_name)

#general
states$states_simplified = map[states$states]
states$states_simplified = factor(states$states_simplified, levels = states_name)

window1 = makeWindows(states, w = 200)
mtch = as.matrix(findOverlaps(window1, states))
window1$states_young = states_young$states_simplified[mtch[, 2]]
window1$states_old = states_old$states_simplified[mtch[, 2]]

########

window1_T = window1[na.omit(window1$states_young == "TssBivalent" | window1$states_old == "TssBivalent")]

transition_mat = table(mcols(window1_T)[ , c("c","states_old")]) * 200
class(transition_mat) = "matrix"
transition_mat = transition_mat[states_name, states_name]
transition_mat

rownames(transition_mat) = paste0("states_old", rownames(transition_mat))
colnames(transition_mat) = paste0("states_old", colnames(transition_mat))
grid.col = c(states_col, states_col)
names(grid.col) = c(rownames(transition_mat), colnames(transition_mat))
chordDiagram(transition_mat, grid.col = grid.col, annotationTrack = c("grid", "axis"),
             directional = TRUE)
circos.clear()
text(0.5, -1, "states_young")
text(0.5, 1, "statates_old")
legend("left", pch = 15, col = states_col, legend = names(states_col))


mat_bivtss = normalizeToMatrix(states[states$states_simplified == "TssBivalent"], tss)
l = rowSums(mat_bivtss[, 40:60]) > 0 # 1kb upstream and downstream
tss_biv = tss[l]
tss_biv


mat_states_ESC = normalizeToMatrix(states, tss_biv, value_column = "states_simplified")
mat_states_lung = normalizeToMatrix(states_lung, tss_biv, value_column = "states_simplified")

mat_meth_ESC = normalizeToMatrix(meth, tss_biv, value_column = "E003", mean_mode = "absolute",
                                 smooth = TRUE)
mat_meth_lung = normalizeToMatrix(meth, tss_biv, value_column = "E096", mean_mode = "absolute",
                                  smooth = TRUE)


mat_meth_diff = mat_meth_ESC - mat_meth_lung
meth_diff_col_fun = colorRamp2(c(-0.25, 0, 0.25), c("#3794bf", "#FFFFFF", "#df8640"))

ht_list = EnrichedHeatmap(mat_states_ESC, name = "states_ESC", col = states_col,
                          top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2), ylim = c(0, 1))),
                          column_title = "States ESC") +
  EnrichedHeatmap(mat_states_lung, name = "states_lung", col = states_col,
                  top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2), ylim = c(0, 1))),
                  show_heatmap_legend = FALSE, column_title = "States lung")

ht_list = ht_list + EnrichedHeatmap(mat_meth_ESC, name = "meth_ESC", col = meth_col_fun,
                                    top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2), ylim = c(0, 1))),
                                    column_title = "Meth ESC") +
  EnrichedHeatmap(mat_meth_lung, name = "meth_lung", col = meth_col_fun,
                  top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2), 
                                                                            ylim = c(0, 1))),
                  show_heatmap_legend = FALSE, column_title = "Meth lung") +
  EnrichedHeatmap(mat_meth_diff, name = "meth_diff", col = meth_diff_col_fun,
                  top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2, 
                                                                                      pos_col = "#df8640", neg_col = "#3794bf"))),
                  column_title = "Meth ESC - lung")


mat_meth_diff_discrete = discretize(mat_meth_diff,
                                    rule = list(
                                      "hypo" = c(-Inf, -0.2),
                                      "hyper" = c(0.2, Inf)
                                    )
)

mat_meth_diff_discrete

ht_list = ht_list + EnrichedHeatmap(mat_meth_diff_discrete, name = "meth_diff_discrete", 
                                    col = c("hyper" = "#df8640", hypo = "#3794bf"),
                                    top_annotation = HeatmapAnnotation(enrich = anno_enriched(gp = gpar(lty = 1:2))))

e = log2(expr[names(tss_biv), c("E003", "E096")] + 1)
ht_list = ht_list + Heatmap(e, name = "expr", 
                            show_row_names = FALSE, width = unit(10, "mm"), cluster_columns = FALSE)

row_order = hclust(dist(cbind(mat_states_ESC[, 40:60], mat_states_lung[, 40:60])))$order
split = ifelse(e[, "E096"] > e[, "E003"], "activation", "repression")
draw(ht_list, row_order = row_order, row_split = split)


