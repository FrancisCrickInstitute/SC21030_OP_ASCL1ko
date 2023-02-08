# Open packages necessary for analysis.
library(tidyverse)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(Seurat)
library(readxl)
library(writexl)
library(GenomeInfoDb)
library(DESeq2)
library(ggrepel)
library(broom)

#Load Seurat objects from same folder as script
temp_seur <- readRDS(file = 'integrated_seurat_object_cell_state_assign_RNA_all.RDS')
#KO-only Seurat object, with label transfer UMAP under 'ref.umap' reduction
ko_seur <- temp_seur[[1]]

#WT + DAPT Seurat object
seur <- readRDS(file = 'integrated_seurat_object_cell_state_assign_RNA.RDS')
#WT + DAPT and ASCL1KO datasets merged and integrated into single Seurat object
temp_seur <- readRDS(file = 'integrated_seurat_object_sampleGroups_celltype.RDS')
double_seur <- temp_seur[[2]]
wt_seur_noDAPT <- temp_seur[[1]]
rm(temp_seur)

#All visualizations were done on RNA assay
DefaultAssay(seur) <- 'RNA'
DefaultAssay(ko_seur) <- 'RNA'
DefaultAssay(wt_seur_noDAPT) <- 'RNA'
DefaultAssay(double_seur) <- 'RNA'

#Renaming ASCL1 Progenitors to Transitional Progenitors as we decided on this name change midway during analysis
seur$cell_classification_res.1 <- str_replace(seur$cell_classification_res.1, 'ASCL1 Progenitors', 'Transitional Progenitors')
seur$group_cell_classification_res.1 <- str_replace(seur$group_cell_classification_res.1, 'ASCL1 Progenitors', 'Transitional Progenitors')
ko_seur$predicted.celltype <- str_replace(ko_seur$predicted.celltype, 'ASCL1 Progenitors', 'Transitional Progenitors')
wt_seur_noDAPT$merged.celltype <- str_replace(wt_seur_noDAPT$merged.celltype, 'ASCL1 Progenitors', 'Transitional Progenitors')
double_seur$merged.celltype <- str_replace(double_seur$merged.celltype, 'ASCL1 Progenitors', 'Transitional Progenitors')

#Genes used for module scores
cycling_prog = c('PAX6','FABP7','VIM','NES','SOX2','MCM2','MKI67','TOP2A','UBE2C','SLC1A3','HES5','HES1')
ascl1_prog = c('ASCL1','SOX4','DLX1','DLL1','RGS16','GADD45G','CDKN1C','SMOC1','DLL3','IGFBP5','HES6')
neurons = c('DCX','ELAVL3','ELAVL4','MAP2','BCL11B','RBFOX3','NEFM','STMN4','STMN2', 'ENSG00000258947')

#Figure 1C, D, E, F, G plots
pdf(file = paste0("Fig 1CDEFG.pdf"), width = 15, height = 8)
Idents(seur) <- 'cell_classification_res.1'
Idents(seur) <- factor(x = Idents(seur), levels = c("15_Cycling Progenitors","5_Cycling Progenitors","13_Cycling Progenitors","1_Cycling Progenitors",
                                                    "19_Cycling Progenitors","10_Cycling Progenitors","17_Cycling Progenitors","18_Cycling Progenitors",
                                                    "3_Cycling Progenitors","12_Cycling Progenitors","4_Cycling Progenitors","21_Cycling Progenitors",
                                                    "0_Transitional Progenitors","2_Transitional Progenitors","7_Transitional Progenitors",
                                                    "11_Neurons","14_Neurons","20_Neurons","9_Neurons","16_Neurons","8_Neurons","6_Neurons"))
plot(DimPlot(seur, group.by = 'integrated_snn_res.1', label = TRUE) + 
       scale_color_discrete(labels = str_sort(levels(seur), numeric = TRUE)) + 
       ylim(-7, 7) + xlim(-10, 10) + scale_y_continuous(breaks = seq(-6,6,3)))
plot(DimPlot(seur, group.by = 'group_cell_classification_res.1', label = TRUE, cols = c('Transitional Progenitors' = '#1f77b4', 'Cycling Progenitors' = '#ff7f0e', 'Neurons' = '#2ca02c')) +
       NoLegend() + ylim(-7, 7) + xlim(-10, 10) + scale_y_continuous(breaks = seq(-6,6,3)))
plot(DimPlot(seur, group.by = 'Phase', label = TRUE, cols = c('G1' = 'gold', 'S' = 'magenta', 'G2M' = 'cyan')) + 
       NoLegend() + ylim(-7, 7) + xlim(-10, 10) + scale_y_continuous(breaks = seq(-6,6,3)))
plot(DotPlot(seur, features = c(cycling_prog, ascl1_prog, neurons)) & 
       scale_color_viridis_c() & RotatedAxis() & xlab('Gene') & ylab('Cluster'))
plot(FeaturePlot(seur, features = 'ASCL1') + scale_color_viridis_c() + ylim(-7, 7) + 
       xlim(-10, 10) + scale_y_continuous(breaks = seq(-6,6,3)))
dev.off()

#Figure S1B, E plots
#Define clusters on the integrated WT+/-DAPT UMAP
wt_seur_noDAPT <- AddModuleScore(wt_seur_noDAPT, features = list(cycling_prog), name = 'Cycling_Progenitors')
wt_seur_noDAPT <- AddModuleScore(wt_seur_noDAPT, features = list(ascl1_prog), name = 'Transitional_Progenitors')
wt_seur_noDAPT <- AddModuleScore(wt_seur_noDAPT, features = list(neurons), name = 'Neurons')
FeaturePlot(wt_seur_noDAPT, features = c('Cycling_Progenitors1', 'Transitional_Progenitors1', 'Neurons1'))
DimPlot(wt_seur_noDAPT, group.by = 'integrated_snn_res.1', label = TRUE) + NoLegend()
wt_seur_noDAPT$group_cell_classification_res.1 <- case_when(as.vector(wt_seur_noDAPT[['integrated_snn_res.1']][[1]]) %in% c('0', '7') ~ 'Transitional Progenitors',
                                                         as.vector(wt_seur_noDAPT[['integrated_snn_res.1']][[1]]) %in% c('2','16','15','19','9','18','13') ~ 'Neurons',
                                                         TRUE ~ 'Cycling Progenitors')
#Define Top 5 DGEs for WT DAPT UMAP clusters
Idents(seur) <- 'cell_classification_res.1'
Idents(seur) <- factor(x = Idents(seur), levels = c("15_Cycling Progenitors","5_Cycling Progenitors","13_Cycling Progenitors","1_Cycling Progenitors",
                                                    "19_Cycling Progenitors","10_Cycling Progenitors","17_Cycling Progenitors","18_Cycling Progenitors",
                                                    "3_Cycling Progenitors","12_Cycling Progenitors","4_Cycling Progenitors","21_Cycling Progenitors",
                                                    "0_Transitional Progenitors","2_Transitional Progenitors","7_Transitional Progenitors",
                                                    "11_Neurons","14_Neurons","20_Neurons","9_Neurons","16_Neurons","8_Neurons","6_Neurons"))
top5dge <- FindAllMarkers(object = seur, only.pos = TRUE, min.pct = 0.25)
top5dge <- top5dge %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf(file = paste0("Fig S1B.pdf"), width = 15, height = 8)
plot(DimPlot(wt_seur_noDAPT, group.by = "group_cell_classification_res.1", label = TRUE,
             label.size = 3, repel = TRUE, split.by = 'treatment',
             cols = c('Transitional Progenitors' = '#1f77b4', 'Cycling Progenitors' = '#ff7f0e', 'Neurons' = '#2ca02c')) +
       NoLegend())
dev.off()
pdf(file = paste0("Fig S1E.pdf"), width = 15, height = 20)
plot(DotPlot(seur, features = unique(top5dge$gene)) + coord_flip() + RotatedAxis() + scale_color_viridis_c() +
       labs(title = 'Top 5 differentially expressed genes by cluster'))
dev.off()
write_xlsx(top5dge, 'top5dge.xlsx')

#Figure 3A, C plots
double_seur$genotype_celltype <- paste(double_seur$genotype, double_seur$merged.celltype)
Idents(double_seur) <- 'genotype_celltype'
levels(double_seur) <- c('ASCL1ko Cycling Progenitors', 'WT Cycling Progenitors',
                    'ASCL1ko Transitional Progenitors', 'WT Transitional Progenitors', 
                    'ASCL1ko Neurons', 'WT Neurons')
pdf(file = paste0("Fig 3AC.pdf"), width = 15, height = 8)
plot(DimPlot(ko_seur, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
             label.size = 3, repel = TRUE, 
             cols = c('Transitional Progenitors' = '#1f77b4', 'Cycling Progenitors' = '#ff7f0e', 'Neurons' = '#2ca02c')) +
       NoLegend() + ggtitle("ASCL1-KO projected UMAP") + ylim(-7, 7) + xlim(-10, 10) + scale_y_continuous(breaks = seq(-6,6,3)))
plot(DotPlot(double_seur, features = c(cycling_prog, ascl1_prog, neurons)) & 
       scale_color_viridis_c() & xlab('Gene') & ylab('Cluster') & coord_flip() & RotatedAxis())
dev.off()

#Figure S2C plot
#Calculate module scores for the 3 celltypes in double_seur, to aid annotation of celltypes
double_seur <- AddModuleScore(double_seur, features = list(cycling_prog), name = 'Cycling_Progenitors')
double_seur <- AddModuleScore(double_seur, features = list(ascl1_prog), name = 'Transitional_Progenitors')
double_seur <- AddModuleScore(double_seur, features = list(neurons), name = 'Neurons')
FeaturePlot(double_seur, features = c('Cycling_Progenitors1', 'Transitional_Progenitors1', 'Neurons1'))
DimPlot(double_seur, group.by = 'integrated_snn_res.1', label = TRUE) + NoLegend()
double_seur$group_cell_classification_res.1 <- case_when(as.vector(double_seur[['integrated_snn_res.1']][[1]]) %in% c('18', '7') ~ 'Transitional Progenitors',
                                                         as.vector(double_seur[['integrated_snn_res.1']][[1]]) %in% c('2') ~ 'Neurons',
                                                         TRUE ~ 'Cycling Progenitors')
pdf(file = paste0("Fig S2C.pdf"), width = 15, height = 8)
plot(DimPlot(double_seur, group.by = "group_cell_classification_res.1", label = TRUE,
             label.size = 3, repel = TRUE, split.by = 'genotype',
             cols = c('Transitional Progenitors' = '#1f77b4', 'Cycling Progenitors' = '#ff7f0e', 'Neurons' = '#2ca02c')) +
       NoLegend())
plot(DimPlot(double_seur, group.by = 'genotype'), split.by = 'group_cell_classification_res.1')
dev.off()

#Figure S3D plot
pdf(file = paste0("Fig S3D.pdf"), width = 15, height = 8)
plot(FeaturePlot(seur, features = 'ACTL6A'))
plot(FeaturePlot(seur, features = 'ACTL6B'))
dev.off()

#Proportion stuff
seur_table <- as_tibble(seur@meta.data) %>% select(orig.ident, group_cell_classification_res.1)
ko_table <- as_tibble(ko_seur@meta.data) %>% select(orig.ident, predicted.celltype)
seur_table %<>% group_by(orig.ident, group_cell_classification_res.1) %>% summarise(n()) %>% add_column(genotype = 'WT')
ko_table %<>% group_by(orig.ident, predicted.celltype) %>% summarise(n()) %>% add_column(genotype = 'ASCL1 KO')
combined_table <- bind_rows(seur_table,ko_table)
write_xlsx(combined_table, 'cluster counts.xlsx')

wt_seur_noDAPT_table <- as_tibble(wt_seur_noDAPT@meta.data) %>% select(orig.ident, group_cell_classification_res.1, treatment)
wt_seur_noDAPT_table %<>% group_by(orig.ident, group_cell_classification_res.1, treatment) %>% summarise(n())
write_xlsx(wt_seur_noDAPT_table, 'integrated_wt_cluster counts.xlsx')

double_table <- as_tibble(double_seur@meta.data) %>% select(orig.ident, group_cell_classification_res.1, genotype)
double_table %<>% group_by(orig.ident, group_cell_classification_res.1, genotype) %>% summarise(n())
write_xlsx(double_table, 'integrated_cluster counts.xlsx')

#Subsetting our Seurat object to obtain pseudobulk count matrices for DESeq2
double_seur$id <- case_when(as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A21' ~ 'WT_1',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A22' ~ 'WT_2',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A23' ~ 'KO_1',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A24' ~ 'KO_2',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A26' ~ 'WT_3',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A27' ~ 'KO_3',
                            TRUE ~ 'Failed')
double_seur$id.celltype <- paste(double_seur$id, double_seur$merged.celltype)
##Warning! Next step may take a while
subsets <- SplitObject(double_seur, split.by = 'id.celltype')
##Summing the RNA counts for each subset - will also take a while
t2 <- NULL
for(i in unique(double_seur$id.celltype)){
  t = apply(subsets[[i]]@assays$RNA@counts, 1, sum)
  t2 = cbind(t2, as.numeric(t))
}
pseudobulk = as_tibble(cbind(names(t), t2))
names(pseudobulk) = c("gene", unique(double_seur$id.celltype))
rm(t2, i, t)
##Making the DESeq2 objects. Will have to round the numbers to obtain integers, as Alevin gives non-integer RNA counts.
genotype <- c('WT', 'WT', 'WT', 'KO', 'KO', 'KO')
CP <- pseudobulk %>% select(`WT_1 Cycling Progenitors`,`WT_2 Cycling Progenitors`,`WT_3 Cycling Progenitors`,
                            `KO_1 Cycling Progenitors`,`KO_2 Cycling Progenitors`,`KO_3 Cycling Progenitors`)
CP_anno <- tibble(genotype)
rownames(CP_anno) <- colnames(CP)
for(i in 1:6){
  CP[[i]] <- as.numeric(CP[[i]])
}
CP <- round(CP)
rownames(CP) <- pseudobulk$gene
CP_bulk <- DESeqDataSetFromMatrix(countData = CP,
                                  colData = CP_anno,
                                  design = ~ genotype)
CP_bulk$genotype <- relevel(CP_bulk$genotype, ref = 'WT')
CP_bulk <- DESeq(CP_bulk)

TP <- pseudobulk %>% select(`WT_1 Transitional Progenitors`,`WT_2 Transitional Progenitors`,`WT_3 Transitional Progenitors`,
                            `KO_1 Transitional Progenitors`,`KO_2 Transitional Progenitors`,`KO_3 Transitional Progenitors`)
TP_anno <- tibble(genotype)
rownames(TP_anno) <- colnames(TP)
for(i in 1:6){
  TP[[i]] <- as.numeric(TP[[i]])
}
TP <- round(TP)
rownames(TP) <- pseudobulk$gene
TP_bulk <- DESeqDataSetFromMatrix(countData = TP,
                                  colData = TP_anno,
                                  design = ~ genotype)
TP_bulk$genotype <- relevel(TP_bulk$genotype, ref = 'WT')
TP_bulk <- DESeq(TP_bulk)

Neu <- pseudobulk %>% select(`WT_1 Neurons`,`WT_2 Neurons`,`WT_3 Neurons`,
                             `KO_1 Neurons`,`KO_2 Neurons`,`KO_3 Neurons`)
Neu_anno <- tibble(genotype)
rownames(Neu_anno) <- colnames(Neu)
for(i in 1:6){
  Neu[[i]] <- as.numeric(Neu[[i]])
}
Neu <- round(Neu)
rownames(Neu) <- pseudobulk$gene
Neu_bulk <- DESeqDataSetFromMatrix(countData = Neu,
                                  colData = Neu_anno,
                                  design = ~ genotype)
Neu_bulk$genotype <- relevel(Neu_bulk$genotype, ref = 'WT')
Neu_bulk <- DESeq(Neu_bulk)

##Saving and plotting the DESeq2 results
write.csv(as.data.frame(results(CP_bulk)), 'KO_vs_WT_CP_DESeq2.csv')
write.csv(as.data.frame(results(TP_bulk)), 'KO_vs_WT_TP_DESeq2.csv')
write.csv(as.data.frame(results(Neu_bulk)), 'KO_vs_WT_Neu_DESeq2.csv')

CP_DGE <- read_csv('KO_vs_WT_CP_DESeq2.csv',) %>% rename(...1 = 'geneID') %>% na.omit()
CP_DGE$diffexpressed <- "NO"
CP_DGE$diffexpressed[CP_DGE$log2FoldChange > 1.5 & CP_DGE$padj < 0.05] <- "UP"
CP_DGE$diffexpressed[CP_DGE$log2FoldChange < -1.5 & CP_DGE$padj < 0.05] <- "DOWN"
pdf(file = paste0("ASCL1KO_vs_WT_CP_DGE.pdf"), width = 15, height = 8)
plot(ggplot(CP_DGE, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) + geom_point(size = 0.5) + 
  geom_hline(yintercept = -log10(0.05), col = 'red') + geom_vline(xintercept = c(-1.5, 1.5), col = 'red') +
  theme_classic() + scale_color_manual(values=c("#4d72c0", "black", "#d8021c")))
dev.off()

TP_DGE <- read_csv('KO_vs_WT_TP_DESeq2.csv',) %>% rename(...1 = 'geneID') %>% na.omit()
TP_DGE$diffexpressed <- "NO"
TP_DGE$diffexpressed[TP_DGE$log2FoldChange > 1.5 & TP_DGE$padj < 0.05] <- "UP"
TP_DGE$diffexpressed[TP_DGE$log2FoldChange < -1.5 & TP_DGE$padj < 0.05] <- "DOWN"
pdf(file = paste0("ASCL1KO_vs_WT_TP_DGE.pdf"), width = 15, height = 8)
ggplot(TP_DGE, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) + geom_point(size = 0.5) + 
  geom_hline(yintercept = -log10(0.05), col = 'red') + geom_vline(xintercept = c(-1.5, 1.5), col = 'red') +
  theme_classic() + scale_color_manual(values=c("#4d72c0", "black", "#d8021c"))
dev.off()

Neu_DGE <- read_csv('KO_vs_WT_Neu_DESeq2.csv',) %>% rename(...1 = 'geneID') %>% na.omit()
Neu_DGE$diffexpressed <- "NO"
Neu_DGE$diffexpressed[Neu_DGE$log2FoldChange > 1.5 & Neu_DGE$padj < 0.05] <- "UP"
Neu_DGE$diffexpressed[Neu_DGE$log2FoldChange < -1.5 & Neu_DGE$padj < 0.05] <- "DOWN"
pdf(file = paste0("ASCL1KO_vs_WT_Neu_DGE.pdf"), width = 15, height = 8)
ggplot(Neu_DGE, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) + geom_point(size = 0.5) + 
  geom_hline(yintercept = -log10(0.05), col = 'red') + geom_vline(xintercept = c(-1.5, 1.5), col = 'red') +
  theme_classic() + scale_color_manual(values=c("#4d72c0", "black", "#d8021c"))
dev.off()

##Correlation plots between scRNAseq DEGs by celltype and bulk RNAseq DEGs
#multiplied the log2fc by -1 due to direction of comparison being WT vs KO instead of KO vs WT
RNA_DGE <- read_delim("NRS_WTvsNRS_ASCL1_KO.deseq2.refseq.results.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(bulk_log2FoldChange = as.numeric(log2FoldChange) * -1) %>%  na.omit() 
CP <- inner_join(CP_DGE, RNA_DGE, by = 'geneID')
TP <- inner_join(TP_DGE, RNA_DGE, by = 'geneID')
Neu <- inner_join(Neu_DGE, RNA_DGE, by = 'geneID')

df <- list(CP=CP, TP=TP, Neu=Neu)
rm(TP, CP, Neu)
for(i in names(df)){
  df[[i]]$diffexpressed <- 'ns'
  df[[i]]$diffexpressed[df[[i]]$padj.x < 0.05 & df[[i]]$padj.y < 0.05] <- 'Both padj<0.05'
  df[[i]]$diffexpressed[df[[i]]$padj.x < 0.05 & df[[i]]$padj.y >= 0.05] <- 'sc padj<0.05'
  df[[i]]$diffexpressed[df[[i]]$padj.x >= 0.05 & df[[i]]$padj.y < 0.05] <- 'bulk padj<0.05'
  }

pdf(file = paste0("ASCL1KO_vs_WT_CP_TP_Neurons_scRNAseq_vs_bulk_correlation.pdf"), width = 15, height = 8)
for(i in names(df)){
  cort <- cor.test(df[[i]]$log2FoldChange.x, df[[i]]$bulk_log2FoldChange)
  plot(ggplot(df[[i]], aes(y = log2FoldChange.x, x = bulk_log2FoldChange, col = diffexpressed)) + geom_point(size = 0.05) +
    geom_hline(yintercept = 0, size = 0.2) + geom_vline(xintercept = 0, size = 0.2) + 
    geom_smooth(method = "lm",
                colour = "black",
                linetype = 2,
                size = 0.2) +
    annotate("text", x = -5.5, y = 5, 
             label = paste("italic(r)(", cort$parameter, ")==", cort$estimate, "~~italic(p)==", cort$p.value), 
             parse = TRUE,
             size = 3) + 
    scale_color_manual(values=c("orange", "#d8021c", "black", 'darkgreen')))
  }
dev.off()

##Re-doing the above, but using pseudobulk of the entire replicate, instead of pseudobulks by celltypes
pseu_DGE <- read_csv('KO_vs_WT_pseudobulk_DESeq2.csv') %>% rename(...1 = 'geneID') %>% na.omit()
#Multiplied the bulk_log2fc as the file's direction of comparison is WT vs KO, and not KO vs WT
RNA_DGE <- read_delim("NRS_WTvsNRS_ASCL1_KO.deseq2.refseq.results.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(bulk_log2FoldChange = as.numeric(log2FoldChange) * -1) %>%  na.omit() 
pseu <- inner_join(pseu_DGE, RNA_DGE, by = 'geneID')
pseu$diffexpressed <- 'ns'
pseu$diffexpressed[pseu$padj.x < 0.05 & abs(pseu$log2FoldChange.x) > log2(1.5)] <- 'sc padj<0.05, log2Fc > log2(1.5)'
pseu$diffexpressed[pseu$padj.y < 0.05 & abs(pseu$bulk_log2FoldChange) > log2(1.5)] <- 'bulk padj<0.05, log2Fc > log2(1.5)'
pseu$diffexpressed[pseu$padj.x < 0.05 & pseu$padj.y < 0.05 & abs(pseu$log2FoldChange.x) > log2(1.5) & abs(pseu$bulk_log2FoldChange) > log2(1.5)] <- 'Both padj<0.05, log2Fc > log2(1.5)'
pdf(file = paste0("ASCL1KO_vs_WT_scRNAseq_vs_bulk_correlation.pdf"), width = 8, height = 8)
cort <- cor.test(pseu$log2FoldChange.x, pseu$bulk_log2FoldChange)
plot(ggplot(pseu, aes(y = log2FoldChange.x, x = bulk_log2FoldChange, col = diffexpressed)) + geom_point(size = 0.05) +
       geom_hline(yintercept = 0, size = 0.2) + geom_vline(xintercept = 0, size = 0.2) + 
       geom_smooth(method = "lm",
                   colour = "black",
                   linetype = 2,
                   size = 0.2) +
       annotate("text", x = -5.5, y = 5, 
                label = paste("italic(r)(", cort$parameter, ")==", cort$estimate, "~~italic(p)==", cort$p.value), 
                parse = TRUE,
                size = 3) + 
       scale_color_manual(values=c("orange", "#d8021c", "black", 'darkgreen')))
dev.off()

### Preparing the WT-only pseudobulk for Cristina to do GO analysis with it as background
temp_seur <- readRDS(file = 'integrated_seurat_object_sampleGroups_celltype.RDS')
double_seur <- temp_seur[[2]]
DefaultAssay(double_seur) <- 'RNA'
double_seur$merged.celltype <- str_replace(double_seur$merged.celltype, 'ASCL1 Progenitors', 'Transitional Progenitors')
double_seur$id <- case_when(as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A21' ~ 'WT_1',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A22' ~ 'WT_2',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A23' ~ 'KO_1',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A24' ~ 'KO_2',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A26' ~ 'WT_3',
                            as.vector(double_seur[['orig.ident']][[1]]) == 'TAN2457A27' ~ 'KO_3',
                            TRUE ~ 'Failed')
subsets <- SplitObject(double_seur, split.by = 'id')
t2 <- NULL
for(i in unique(double_seur$id)){
  t = apply(subsets[[i]]@assays$RNA@counts, 1, sum)
  t2 = cbind(t2, as.numeric(t))
}
pseudobulk = as_tibble(cbind(names(t), t2))
names(pseudobulk) = c("gene", unique(double_seur$id))
rm(t2, i, t)
for(i in 2:7){
  pseudobulk[[i]] <- as.numeric(pseudobulk[[i]]) %>% round()
}
pseudobulk2 <- select(pseudobulk, 'gene','WT_1','WT_2','WT_3') %>% subset(WT_1 > 10|WT_2 > 10|WT_3 > 10)
write_csv(pseudobulk2, 'WT_pseudobulk_for_background_filtered.csv')
