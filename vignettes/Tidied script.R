# Open packages necessary for analysis.
library(tidyverse)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(Seurat)
library(readxl)
library(writexl)

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
