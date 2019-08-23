#change diffTF colours and labels

library("tidyverse")
library("ggplot2")
library("ggrepel")
library("ggpubr")

dir <-'/g/scb2/zaugg/alishaye/'
plotdir <-'/g/scb2/zaugg/alishaye/termporal_r_plots_R/diffTF_plots/'

#HOCOMOCO
#'TF_aml_06_2019/2_groups_1000_bootstraps/output/FINAL_OUTPUT/extension100/'
#'TF_aml_06_2019/Continous_1000_bootstraps/output/FINAL_OUTPUT/extension100/'

#REMAP - SEX
#'TF_aml_06_2019/using_remap/bootstraps/runremap/output/FINAL_OUTPUT/extension100/' #cont boots
#'TF_aml_06_2019/using_remap/bootstraps/mutsepremap/output/FINAL_OUTPUT/extension100/' #sep permutation
#'TF_aml_06_2019/using_remap/permutations/remap_cont/output/FINAL_OUTPUT/extension100/'#permutation

#REMAP - MUTATION
#'TF_aml_06_2019/using_remap/permutations/mutcontremap/output/FINAL_OUTPUT/extension100/' #permutation
#'TF_aml_06_2019/using_remap/bootstraps/cgd_outiers/output/FINAL_OUTPUT/extension100/'






#----the ggplot----
vulcano.plot.class <- function(table, xvector,yvector)
{label <- table %>% filter(log10(-log10(pvalueAdj)) > 1.5)
return(ggplot()+
         geom_rect(mapping=aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf,
                               fill="Young"),alpha=0.3)+
         geom_rect(mapping=aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill="Young"),
                   alpha=0.3)+
         scale_fill_manual(name="TF activity higher in:", values = c("Young"="#ef8a62","Old"="#67a9cf"))+
         geom_point(data=table,
                    mapping = aes(
                      x = weighted_meanDifference,
                      y = log10(-log10(pvalueAdj)),
                      color=TFBS))+
         scale_color_gradient(low = "darkseagreen4", high = "darkolivegreen2")+
         ylim(yvector)+
         xlim(xvector)+
         geom_label_repel(
           data = label,
           mapping = aes(
             x = weighted_meanDifference,
             y = log10(-log10(pvalueAdj)),
             label=TF,
             fill=classification_q0.001_final),
           segment.color="black", segment.alpha = 0.3)+
         scale_fill_manual(name="Classification:", values = c("Young"="#ef8a62","Old"="#67a9cf"))+
         xlab("Weighted Mean Difference")+
         ylab("log10(-log10(Adjusted pvalue))")+
         ggtitle("Weighted Mean Differences of TFs between Young and Old")+
         labs(color='Transcription Factor Binding Sites'))}


vulcano.plot <- function(table, xvector,yvector)
{label <- table %>% filter(log10(-log10(pvalueAdj)) > 1.5)
return(ggplot()+
         geom_rect(mapping=aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf,
                               fill="Young"),alpha=0.3)+
         geom_rect(mapping=aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill="Old"),
                   alpha=0.3)+
         scale_fill_manual(name="TF activity higher in:", values = c("Young"="#ef8a62","Old"="#67a9cf"))+
         geom_point(data=table,
                    mapping = aes(
                      x = weighted_meanDifference,
                      y = log10(-log10(pvalueAdj)),
                      color=TFBS))+
         scale_color_gradient(low = "darkseagreen4", high = "darkolivegreen2")+
         ylim(yvector)+
         xlim(xvector)+
         geom_label_repel(
           data = label,
           mapping = aes(
             x = weighted_meanDifference,
             y = log10(-log10(pvalueAdj)),
             label=TF),
           segment.color="black", segment.alpha = 0.3)+
         xlab("Weighted Mean Difference")+
         ylab("log10(-log10(Adjusted pvalue))")+
         ggtitle("Weighted Mean Differences of TFs between Young and Old")+
         labs(color='Transcription Factor Binding Sites'))}

vulcano.plot.less <- function(table)
{label <- table %>% filter(log10(-log10(pvalueAdj)) > 1.5)
return(ggplot()+
         geom_rect(mapping=aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf,
                               fill="Young"),alpha=0.3)+
         geom_rect(mapping=aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill="Old"),
                   alpha=0.3)+
         scale_fill_manual(name="TF activity higher in:", values = c("Young"="#ef8a62","Old"="#67a9cf"))+
         geom_point(data=table,
                    mapping = aes(
                      x = weighted_meanDifference,
                      y = -log10(pvalueAdj),
                      color=TFBS))+
         scale_color_gradient(low = "darkseagreen4", high = "darkolivegreen2")+
         #ylim(yvector)+
         #xlim(xvector)+
         geom_label_repel(
           data = label,
           mapping = aes(
             x = weighted_meanDifference,
             y = -log10(pvalueAdj),
             label=TF),
           segment.color="black", segment.alpha = 0.3)+
         xlab("Weighted Mean Difference")+
         ylab("-log10(Adjusted pvalue)")+
         ggtitle("Weighted Mean Differences of TFs between Young and Old")+
         labs(color='Transcription Factor Binding Sites')
)}

foldbars <- function(table)
{
  label <- table %>% filter(log10(-log10(pvalueAdj)) > 1.5)
  ggplot(table)+
    geom_bar(mapping=aes(x=reorder(TF, weighted_meanDifference),
                         y= weighted_meanDifference,
                         fill=log10(-log10(pvalueAdj))),
             stat="identity")+
    scale_fill_gradient(high = "darkgoldenrod", low = "wheat")+ 
    geom_label_repel(
      data = label,
      mapping = aes(
        x = TF,
        y = weighted_meanDifference,
        label=TF),
      segment.color="black",
      segment.alpha = 0.3) +
    xlab("Transcription Factor")+
    ylab("Weighted Mean Difference")+
    ggtitle("Significant Transcription Factor's Weighted Mean Differences")+
    labs(fill='log10(-log10(Adjusted pvalue))') +
    theme(axis.text.x = element_text(color="black", 
                                     size=7, angle=90))
}
#----pdf----

file1 <- 'age.summary.tsv.gz'
folder <- 'TF_aml_06_2019/2_groups_1000_bootstraps/output/FINAL_OUTPUT/extension100/'
tsv_scr <- read_tsv(paste0(dir,folder,file1))
filter <- filter(tsv_scr, pvalueAdj < 0.001) 

pdf(paste0(plotdir, "HOCO_BT_2GR_GENDER.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file2 <- 'age.summary.tsv.gz'
folder <- 'TF_aml_06_2019/Continous_1000_bootstraps/output/FINAL_OUTPUT/extension100/'
tsv_scr <- read_tsv(paste0(dir,folder,file2))
filter <- filter(tsv_scr, pvalueAdj < 0.001)

pdf(paste0(plotdir, "HOCO_BT_CON_GENDER.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()


file3 <- 'GetOld.summary.tsv.gz'
folder <- 'TF_aml_06_2019/using_remap/bootstraps/runremap/output/FINAL_OUTPUT/extension100/' #cont boots
tsv_scr <- read_tsv(paste0(dir,folder,file3))
filter <- filter(tsv_scr, pvalueAdj < 0.001)

pdf(paste0(plotdir, "REMAP_BT_CON_GENDER.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()


file4 <- 'GetOld.summary.tsv.gz'
folder <-  'TF_aml_06_2019/using_remap/bootstraps/mutsepremap/output/FINAL_OUTPUT/extension100/' #sep permutation
tsv_scr <- read_tsv(paste0(dir,folder,file4))
filter <- filter(tsv_scr, pvalueAdj < 0.001)
pdf(paste0(plotdir, "REMAP_BT_2GRS_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file5 <- 'PermCont.summary.tsv.gz'
folder <-  'TF_aml_06_2019/using_remap/permutations/remap_cont/output/FINAL_OUTPUT/extension100/'#permutation
tsv_scr <- read_tsv(paste0(dir,folder,file5))
filter <- filter(tsv_scr, pvalueAdj < 0.05)
pdf(paste0(plotdir, "REMAP_PT_CON_GENDER.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file6 <- 'GetOld.summary.tsv.gz'
folder <- 'TF_aml_06_2019/using_remap/permutations/mutcontremap/output/FINAL_OUTPUT/extension100/' #permutation
tsv_scr <- read_tsv(paste0(dir,folder,file6))
filter <- filter(tsv_scr, pvalueAdj < 0.05)

pdf(paste0(plotdir, "REMAP_PT_CON_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file7 <- 'Bootstraps.summary.tsv.gz'
folder <- 'TF_aml_06_2019/using_remap/bootstraps/cgd_outiers/output/FINAL_OUTPUT/extension100/'
tsv_scr<-read_tsv('/g/scb2/zaugg/alishaye/TF_aml_06_2019/using_remap/bootstraps/cgd_outliers/output/FINAL_OUTPUT/extension100/Bootstraps.summary.tsv.gz')
filter <- filter(tsv_scr, pvalueAdj < 0.001)

pdf(paste0(plotdir, "REMAP_BT_up2GR_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file8 <- 'Bootstraps.summary.tsv.gz'
folder <- 'TF_aml_06_2019/HOCO_2grs_new_samples/output/FINAL_OUTPUT/extension100/'
tsv_scr <- read_tsv(paste0(dir,folder,file8))
filter <- filter(tsv_scr, pvalueAdj < 0.001)

pdf(paste0(plotdir, "HOCO_BT_up2GR_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file9 <- 'Permutations.summary.tsv.gz'
folder <- 'TF_aml_06_2019/HOCO_2grs_new_samples/output/FINAL_OUTPUT/extension100/'
tsv_scr <- read_tsv(paste0(dir,folder,file8))
filter <- filter(tsv_scr, pvalueAdj < 0.001)
pdf(paste0(plotdir, "HOCO_PT_up2GR_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

file10 <- 'Permutations.summary.tsv.gz'
folder <- 'TF_aml_06_2019/HOCO_PT_up2grs_new_samples/output/FINAL_OUTPUT/extension100/'
tsv_scr <- read_tsv(paste0(dir,folder,file8))
filter <- filter(tsv_scr, pvalueAdj < 0.001)

pdf(paste0(plotdir, "REMAP_PT_up2GR_MUT.pdf"), width = 12, height = 12 )
#plot all TFS
vulcano.plot(tsv_scr, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(tsv_scr)
#filtered - only 0,001
vulcano.plot(filter, yvector=c(-2.5, 3), xvector=c(-0.3,0.3))
vulcano.plot.less(filter)
foldbars(filter)
dev.off()

