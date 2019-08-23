library("ChIPpeakAnno")
library("GenomicRanges")
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("biomaRt")
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotate genomic intervals in bed format using ChIPpeakAnno
# This script was written for Arabidopsis
# Has to be changed to accomodate any other organism available through Ensembl Biomart

setwd("comparisons/")
list.files()

annotateDiffBindResult <- function(resFile, outFile) 
{
  # resFile <- path to resFile
  # name of the out file that contains annotated results
  
  # Read in results table
  res <- read.csv(resFile, header=T)
  
  # Convert to GRanges
  gr <- makeGRangesFromDataFrame(res, ignore.strand = T, seqnames.field = "chr",
                                 start.field = "start", end.field = "end")
  
  # Give ranges numeric names in order
  names(DE_res_granges) <- c(1:length(DE_res_granges))
  
  # Create GRanges object with annotations from TxDb database
  annoData <- toGRanges(txdb, feature="gene")
  
  # Annotate granges with the nearest TSS
  annot <- annotatePeakInBatch(DE_res_granges, 
                               AnnotationData=annoData, 
                               featureType = "TSS",
                               output="nearestLocation",
                               PeakLocForDistance = "start")
  
  annot_res <- annotatePeakInBatch(granges_res, 
                               AnnotationData=annoData, 
                               featureType = "TSS",
                               output="nearestLocation",
                               PeakLocForDistance = "start")
  
  
  
  # Load mart
  ensembl <- useMart(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     host = "www.ensembl.org",
                     ensemblRedirect = FALSE)
  
  # Add gene information
  annot <- addGeneIDs(annot, mart = ensembl, feature_id_type = "hgnc_symbol",
                      IDs2Add = c("ensembl_gene_id"))
  
  write.table(annot, outFile, sep = "\t", col.names = T, 
              row.names = F, quote = F)
}

annotateDiffBindResult(resFile="DiffBind_result.csv", 
                       outFile="DiffBind_result_annot.txt")