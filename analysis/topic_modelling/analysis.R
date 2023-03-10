library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(plyr)
library(scales)
library(ggpubr)
countdata <- read.table("/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/de_dimension_denoise/new_count_2000/count.matrix",row.names= 1,header=TRUE,check.names=FALSE,comment.char = "",quote="",sep='')
metadata <- read.csv("/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/de_dimension_denoise/new_count_2000/ecDNA_metadata.csv",header = TRUE, sep =',',row.names = 1)

#filter and normalize the matrix (no scaling, after scaling the cluster merge)
countdata_filter <- countdata[which(rowMeans(countdata)>1),]
countdata_filter <- t((t(countdata_filter)/colSums(countdata_filter))*mean(colSums(countdata_filter)))
countdata_filter <- log(countdata_filter+1)
countdata_filter <- countdata_filter[which(rowMeans(countdata_filter) > 0),]
cell_select <- rownames(metadata[which(metadata$Circle_Read_Percent... > 50 & metadata$Mapping_Ratio... > 80),])
countdata_filter <- countdata_filter[,which(colnames(countdata_filter) %in% cell_select)]
metadata <- metadata[cell_select,]

#cistopic
suppressWarnings(library(cisTopic))
cisTopicObject <- createcisTopicObject(countdata_filter, project.name='ecDNA')
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=17, addModels=FALSE)
pdf('model_parameter.pdf')
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()
cisTopicObject <- runUmap(cisTopicObject, target='cell')

# cisTopicObject <- readRDS('cisTopicObject.rds')

pdf('cell_umap_embedding.pdf')
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('Cell_type'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

pdf('cell_umap.pdf')
par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

pdf('cell_topic_heatmap.pdf')
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('Cell_type'))
dev.off()

pdf('cell_drug_topic_heatmap.pdf')
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('Cell_type','Drug'))
dev.off()

#define the topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

#visualize the topic
pdf('topic_parameter.pdf')
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
dev.off()
cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)
pdf('topic_umap.pdf')
par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

#visualize the topic in bigwig
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
getBigwigFiles(cisTopicObject, path='cisTopics_asBW', seqlengths=seqlengths(txdb)) #p

# #define the topic
library(org.Hs.eg.db)
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

pdf('annotateRegions_heatmap.pdf')
par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
dev.off()

saveRDS(cisTopicObject,'cisTopicObject.rds')

# plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)

#annotate the region with rGREAT
cisTopicObject <- GREAT(cisTopicObject, genome='hg38', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
pdf('GO_annotation_pheatmap.pdf')
ontologyDotPlot(cisTopicObject, top=5, topics=c(1:14), var.y='name', order.by='Binom_Adjp_BH')
dev.off()

#enrichment of specific marker
path_to_signatures <- c('/gshare/xielab/chenjx/ecDNA_Project/analysis/TSS_positioning/HeLa/HeLa_loop_middle.bed','/gshare/xielab/chenjx/ecDNA_Project/analysis/TSS_positioning/HeLa/HeLa_hairpin.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/de_dimension_denoise/new_count_2000_including_drug_cells/k562_loop_region.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/de_dimension_denoise/new_count_2000_including_drug_cells/k562_harpin_region.bed','/gshare/xielab/chenjx/ecDNA_Project/resources/cell_line_reference/HeLa/H3K4me3/H3K4me3/H3K4me3/Regions_H3K4me3.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed','/gshare/xielab/chenjx/ecDNA_Project/resources/cell_line_reference/HeLa/H3K9me3/H3K9me3/H3K9me3/Regions_H3K9me3.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed','/gshare/xielab/chenjx/ecDNA_Project/resources/cell_line_reference/HeLa/H3K27ac/H3K27ac/H3K27ac/Regions_H3K27ac.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed','/gshare/xielab/chenjx/ecDNA_Project/resources/cell_line_reference/HeLa/H3K27me3/H3K27me3/H3K27me3/Regions_H3K27me3.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed','/gshare/xielab/chenjx/ecDNA_Project/resources/cell_line_reference/HeLa/H3K36me3/H3K36me3/H3K36me3/Regions_H3K36me3.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed')
labels  <- c('HeLa_loop','HeLa_harpin','K562_loop','K562_harpin','HeLa_H3K4me3','HeLa_H3K9me3','HeLa_H3K27ac','HeLa_H3K27me3','HeLa_H3K36me3')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)

pdf('3D_annotation_heatmap.pdf')
par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject)
dev.off()

#enrichment of specific marker
path_to_signatures <- c('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/K562/H3K9me3/K562_H3K9me3.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/293T/H3K9me3/ENCFF037SXA.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/PC3/H3K9me3/ENCFF022FXG.bed')
labels  <- c('K562_H3K9me3','293T_H3K9me3','PC3_H3K9me3')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)


# Compute cell rankings
pred.matrix <- predictiveDistribution(cisTopicObject)
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

# Plot
pdf('3D_annotation_embedding.pdf')
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('HeLa_loop','HeLa_harpin','K562_loop','K562_harpin'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
dev.off()

pdf('HeLa_histone_embedding.pdf')
par(mfrow=c(2,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('HeLa_H3K4me3','HeLa_H3K9me3','HeLa_H3K27ac','HeLa_H3K27me3','HeLa_H3K36me3'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
dev.off()

pdf('H3K9me3_embedding.pdf')
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('HeLa_H3K9me3','K562_H3K9me3','293T_H3K9me3','PC3_H3K9me3'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
dev.off()

pdf('H3K9me3_heatmap.pdf')
par(mfrow=c(2,2))
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('LineType', 'cellLine'))
dev.off()

#enrichment of specific marker
path_to_signatures <- c('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/PC3/H3K27me3/ENCFF061MHU.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/K562/H3K27me3/ENCFF881ONN.bed')
labels  <- c('PC3_H3K27me3','K562_H3K27me3')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings
pred.matrix <- predictiveDistribution(cisTopicObject)
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

path_to_signatures <- c('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/resources/293T/H3K27me3/H3K27me3/H3K27me3/Regions_H3K27me3.Sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed')
labels  <- c('293T_H3K27me3')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings
pred.matrix <- predictiveDistribution(cisTopicObject)
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)


pdf('H3K27me3_embedding.pdf')
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('HeLa_H3K27me3','K562_H3K27me3','PC3_H3K27me3','293T_H3K27me3'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
dev.off()


cell_data <- cisTopicObject@cell.data
pdf('HeLa_H3K9me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "HeLa_H3K9me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("HeLa H3K9me3 region enrichment")
dev.off()

pdf('K562_H3K9me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "K562_H3K9me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("K562 H3K9me3 region enrichment")
dev.off()

pdf('PC3_H3K9me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "PC3_H3K9me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("PC3 H3K9me3 region enrichment")
dev.off()

pdf('293T_H3K9me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "293T_H3K9me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("293T H3K9me3 region enrichment")
dev.off()

pdf('HeLa_H3K27me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "HeLa_H3K27me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("HeLa H3K27me3 region enrichment")
dev.off()

pdf('K562_H3K27me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "K562_H3K27me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("K562 H3K27me3 region enrichment")
dev.off()

pdf('PC3_H3K27me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "PC3_H3K27me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("PC3 H3K27me3 region enrichment")
dev.off()

pdf('293T_H3K27me3_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "293T_H3K27me3",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("293T H3K27me3 region enrichment")
dev.off()


#RNA expression
path_to_signatures <- c('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/high_expression_gene/Colo320DM_high_gene_region.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/high_expression_gene/PC3_high_gene_region.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/high_expression_gene/HeLa_high_gene_region.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/high_expression_gene/K562_high_gene_region.bed','/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/high_expression_gene/H293T_high_gene_region.bed')
labels  <- c('Colo320DM_high_gene_region','PC3_high_gene_region','HeLa_high_gene_region','K562_high_gene_region','293T_high_gene_region')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings
pred.matrix <- predictiveDistribution(cisTopicObject)
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

pdf('RNA_expression_embedding.pdf')
par(mfrow=c(3,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('Colo320DM_high_gene_region','PC3_high_gene_region','HeLa_high_gene_region','K562_high_gene_region','293T_high_gene_region'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
dev.off()

#RNA expression
path_to_signatures <- c('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/RNA_seq/contrast/res_HeLa_vs_293T.bed')
labels  <- c('HeLa_vs_293T')
cisTopicObject <- getSignaturesRegions(cisTopicObject, path_to_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings
pred.matrix <- predictiveDistribution(cisTopicObject)
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

# pdf('RNA_contrast_embedding.pdf')
# par(mfrow=c(1,2))
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('HeLa_vs_K562','K562_vs_HeLa'), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
# dev.off()

# pdf('HeLa_vs_K562_boxplot.pdf')
# ggboxplot(cell_data, "Cell_type", "HeLa_vs_K562",
#      fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("HeLa_vs_K562_boxplot")
# dev.off()

pdf('HeLa_vs_293T_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "HeLa_vs_293T",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("HeLa_vs_293T_boxplot")
dev.off()


pdf('annotateRegions_heatmap.pdf')
par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject, selected.signatures = c('HeLa_H3K9me3','K562_H3K9me3','293T_H3K9me3','PC3_H3K9me3','HeLa_H3K27me3','K562_H3K27me3','293T_H3K27me3','PC3_H3K27me3'))
dev.off()

pdf('cell_Histone_heatmap.pdf')
par(mfrow=c(1,1))
cell_data <- cisTopicObject@cell.data
cell_data <- cell_data[order(cell_data$Cell_type),]
coldata <- as.data.frame(cell_data$Cell_type)
cell_data <- cell_data[,c('HeLa_H3K9me3','K562_H3K9me3','293T_H3K9me3','PC3_H3K9me3','HeLa_H3K27me3','K562_H3K27me3','293T_H3K27me3','PC3_H3K27me3')]
bk <- c(seq(0,0.38,by=0.0002))
rownames(coldata) <- rownames(cell_data)
pheatmap::pheatmap(cell_data,show_colnames =T,show_rownames = F,
         cluster_cols = T,
         cluster_rows = F,
         annotation_row=coldata,
         breaks = bk,
        #  color =  brewer.pal(9,"YlOrRd"))
         color = c(colorRampPalette(colors = c('#76A199','#EAD68D'))(length(bk)/2),colorRampPalette(colors = c('#EAD68D','#CE6D48'))(length(bk)/2)))
# signaturesHeatmap(cisTopicObject, selected.signatures = c('HeLa_H3K9me3','K562_H3K9me3','293T_H3K9me3','PC3_H3K9me3','HeLa_H3K27me3','K562_H3K27me3','293T_H3K27me3','PC3_H3K27me3'))
dev.off()



cell_data <- cisTopicObject@cell.data
pdf('Colo320DM_high_gene_region_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "Colo320DM_high_gene_region",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("Colo320DM highly expressed gene")
dev.off()

pdf('PC3_high_gene_region_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "PC3_high_gene_region",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("PC3 highly expressed gene")
dev.off()

pdf('HeLa_high_gene_region_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "HeLa_high_gene_region",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("HeLa highly expressed gene")
dev.off()

pdf('K562_high_gene_region_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "K562_high_gene_region",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("K562 highly expressed gene")
dev.off()

pdf('293T_high_gene_region_enrichment_boxplot.pdf')
ggboxplot(cell_data, "Cell_type", "293T_high_gene_region",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("293T highly expressed gene")
dev.off()


countdata <- read.table("/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/de_dimension_denoise/new_count_2000/count.matrix",row.names= 1,header=TRUE,check.names=FALSE,comment.char = "",quote="",sep='')
metadata <- read.csv("/gshare/xielab/chenjx/ecDNA_Project/ecDNA_metadata_with_TE.csv",header = TRUE, sep =',',row.names = 1)
metadata <- metadata[colnames(countdata),]
metadata[,c(13:861)] <- (metadata[,c(13:861)]/colSums(countdata))*mean(colSums(countdata))
metadata$Transposon_enrichment <- unname(rowSums(metadata[,c(13:861)]))

pdf('L1P1_enrichment_boxplot.pdf')
ggboxplot(metadata, "Cell_type", "L1P1",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("L1P1")
dev.off()

pdf('LTR49.int_enrichment_boxplot.pdf')
ggboxplot(metadata, "Cell_type", "LTR49.int",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("LTR49.int")
dev.off()

pdf('Transposon_enrichment_boxplot.pdf')
ggboxplot(metadata, "Cell_type", "Transposon_enrichment",
     fill = "Cell_type", palette = c("#00AFBB", "#E7B800","#FC4E07","#9b3a74","#3cb346")) + xlab("Cell Type") + ylab("Transposon_enrichment")+ylim(0, 100000)
dev.off()

svg('HeLa_H3K27me3_region.svg')
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('HeLa_H3K27me3'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

svg('K562_H3K27me3_region.svg')
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('K562_H3K27me3'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()


svg('ncell.svg')
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('nCells'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

svg('annotation.svg')
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()
