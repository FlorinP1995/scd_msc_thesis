library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpointdensity)
library(scico)
library(scales)
library(cowplot)

combine_datasets <- function(m1) {
  
  m1_seu <- CreateSeuratObject(counts = m1, project = 'ss1')

  m1_seu[['percent.mt']] <- PercentageFeatureSet(m1_seu, pattern = "^MT")
  m1_seu[["percent.rp"]] <- PercentageFeatureSet(m1_seu, pattern = "^RP")
  

  # QUALITY METIRCS
  
  # UMIs per cell
  
  umis_per_cell <- ggplot(m1_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Cell Ranger") + xlab("Number of Molecules Detected") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell <- ggplot(m1_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Alevin-Fry") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell <- ggplot(m1_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Kallisto Bustools") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell <- ggplot(m1_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for STARsolo") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell_combined <- umis_per_cell1 + umis_per_cell2 + umis_per_cell3 + umis_per_cell
  ggsave("results/umis_per_cell.png", umis_per_cell, width = 2560, height = 1369, units = c("px"))
  
  # Genes per cell
  genes_per_cell <- ggplot(m1_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Cell Ranger") + xlab("Number of Genes Detected") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell <- ggplot(m1_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Alevin-Fry") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell <- ggplot(m1_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Kallisto Bustools") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell <- ggplot(m1_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for STARsolo") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell_combined <- genes_per_cell1 + genes_per_cell2 + genes_per_cell3 + genes_per_cell
  ggsave("results/genes_per_cell.png", genes_per_cell, width = 2560, height = 1369, units = c("px"))
  
  # Mitochondrial genes per cell
  mitochondrial_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Cell Ranger") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Alevin-Fry") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Kallisto Bustools") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for STARsolo") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell_combined <- mitochondrial_genes_per_cell1 + mitochondrial_genes_per_cell2 + mitochondrial_genes_per_cell3 + mitochondrial_genes_per_cell 
  ggsave("results/mitochondrial_genes_per_cell1.png", mitochondrial_genes_per_cell, width = 2560, height = 1369, units = c("px"))
  
  # Ribosomal genes per cell
  ribosomal_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for Cell Ranger") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell For Alevin-Fry") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for Kallisto Bustools") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell <- ggplot(m1_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for STARsolo") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell_combined <- ribosomal_genes_per_cell1 + ribosomal_genes_per_cell2 + ribosomal_genes_per_cell3 + ribosomal_genes_per_cell
  ggsave("results/ribosomal_genes_per_cell.png", ribosomal_genes_per_cell, width = 2560, height = 1369, units = c("px"))
  
  # TOTAL UMI COUNTS PER PERCENTAGE MITOCHONDRIAL
  #umi_vs_mito <- ggplot(m1_seu@meta.data, aes(percent.mt, nCount_RNA)) +
  #  geom_pointdensity() +
  #  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  #  labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
  #  ggtitle("Total UMI Counts per Percentage Mitochondrial\nsfor Cell Ranger") +
  #  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  #  scale_y_log10()
  #umi_vs_mito <- ggplot(m1_seu@meta.data, aes(percent.mt, nCount_RNA)) +
  #  geom_pointdensity() +
  #  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  #  labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
  #  ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor Alevin-Fry") +
  #  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  #  scale_y_log10()
  #umi_vs_mito <- ggplot(m1_seu@meta.data, aes(percent.mt, nCount_RNA)) +
  #  geom_pointdensity() +
  #  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  #  labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
  #  ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor Kallisto Bustools") +
  #  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  #  scale_y_log10()
  #umi_vs_mito <- ggplot(m1_seu@meta.data, aes(percent.mt, nCount_RNA)) +
  #  geom_pointdensity() +
  #  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  #  labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
  #  ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor STARsolo") +
  #  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  #  scale_y_log10()
  #umi_vs_mito_combined <- umi_vs_mito1 + umi_vs_mito2 + umi_vs_mito3 + umi_vs_mito
  #ggsave("results/umi_vs_mito_combined.png", umi_vs_mito, width = 2560, height = 1369, units = c("px"))
  
  # PROCESSING
  
  # Normalize the data
  m1_seu <- NormalizeData(m1_seu, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Scale the data
  all.genes <- rownames(m1_seu)
  m1_seu <- ScaleData(m1_seu, features = all.genes)
  
  # IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)
  m1_seu <- FindVariableFeatures(m1_seu, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(m1_seu), 10)
  top2 <- head(VariableFeatures(m1_seu), 2)
  plot1 <- VariableFeaturePlot(m1_seu, log = FALSE)
  LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave("results/top10.jpeg", plot1, width = 2560, height = 1369, units = c("px"))
  
  
  # Perform Linear Dimensional reduction
  m1_seu <- RunPCA(m1_seu, features = VariableFeatures(object = m1_seu), ndims = 10)
  elbow_plot <- ElbowPlot(m1_seu, ndims = 10) + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ggsave("results/elbow_plot.jpeg", elbow_plot, width = 2560, height = 1369, units = c("px"))
  
  print(m1_seu[["pca"]], dims = 1:5, nfeatures = 5)
  
  # Clustering and visualization
  m1_seu <- FindNeighbors(m1_seu, dims = 1:10)
  m1_seu <- FindClusters(m1_seu, dims = 1:10)
  pca_plot <- PCAPlot(m1_seu) + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) 
  ggsave("results/pca_plot.jpeg", pca_plot, width = 2560, height = 1369, units = c("px"))
  
  # tSNE per Library
  m1_seu <- RunTSNE(m1_seu, dims = 1:10, check_duplicates = FALSE)
  tsne_plot <- DimPlot(object = m1_seu, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'indianred1', 'ss2' = 'grey87', 'ss3' = 'grey87', 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Cell Ranger") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne1_plot <- DimPlot(object = m1_seu, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = adjustcolor('royalblue1', alpha.f = 1), 'ss3' = adjustcolor('grey87', alpha.f = 1), 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Alevin-Fry") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne2_plot <- DimPlot(object = m1_seu, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'lightgreen', 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Kallisto Bustools") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne3_plot <- DimPlot(object = m1_seu, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'grey87','ss4' = 'purple')) + ggtitle("tSNE Plot for STARsolo") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne_all <- tsne1_plot + tsne2_plot + tsne3_plot + tsne_plot
  ggsave("results/tsne_all.png", tsne_all, width = 2560, height = 1369, units = c("px"))
  
  # UMAP per Library
  m1_seu <- RunUMAP(m1_seu, dims = 1:10, verbose = FALSE)
  umap_plot <- DimPlot(object = m1_seu, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'indianred1', 'ss2' = 'grey87', 'ss3' = 'grey87', 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Cell Ranger") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss1 <- DimPlot(object = m1_seu, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = adjustcolor('royalblue1', alpha.f = 1), 'ss3' = adjustcolor('grey87', alpha.f = 1), 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Alevin-Fry") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss2 <- DimPlot(object = m1_seu, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'lightgreen', 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Kallisto Bustools") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss3 <- DimPlot(object = m1_seu, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'grey87','ss4' = 'purple'))  + ggtitle("UMAP Plot for STARsolo") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_all <- umap_ss1 + umap_ss2 + umap_ss3 + umap_plot
  ggsave("results/umap_all.png", umap_all, width = 2560, height = 1369, units = c("px"))
  
  # CLUSTERING
  
  clustered_umap <- DimPlot(m1_seu, reduction = "umap")  + ggtitle("Clustered UMAP") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  clustered_tsne <- DimPlot(m1_seu, reduction = "tsne")  + ggtitle("Clustered tSNE") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ggsave("results/clustered_umap.png", clustered_umap, width = 2560, height = 1369, units = c("px"))
  ggsave("results/clustered_tsne.png", clustered_tsne, width = 2560, height = 1369, units = c("px"))
  
  # DE ON CLUSTERS
  
  # find all markers of cluster 2
  #cluster2.markers <- FindMarkers(seu_combined, ident.1 = 2, min.pct = 0.25)
  #head_2 <- head(cluster2.markers, n = 5)
  
  #node <- list()
  
  #for (i in 1:10) {
  
  #  cluster_marker <- FindMarkers(seu_combined, ident.1 = i, min.pct = 0.25)
  #  node[[i]] <- head(cluster_marker, n = 7)
  
  #}
  
  #DefaultAssay(seu_combined) <- "RNA"
  #nk.markers <- FindConservedMarkers(seu_combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
  #head(nk.markers)
  
  m1_seu.markers <- FindAllMarkers(m1_seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  m1_seu.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # Heatmap of top 10 genes per cluster
  
  top10 <- m1_seu.markers %>% group_by(cluster) %>% top_n(n = 2)
  heat_map_all <- DoHeatmap(m1_seu, features = top10$gene) + ggtitle("Heatmap of Top 10 Genes per Cluster") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 5)) #  + NoLegend() 
  ggsave("results/heat_map_all.png", heat_map_all, width = 2560, height = 1369, units = c("px"))
  
  top_gene_clusters <- list()
  
  for (i in 1:(length(top10)-1)) {
    top_gene_clusters[[i]] <- DotPlot(m1_seu, features = top10$gene[(1+10*(i-1)):(10+10*(i-1))]) + RotatedAxis() + coord_flip() + theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5), legend.key.height = unit(0.15, 'cm'), legend.key.width = unit(0.05, 'cm')) #+ NoLegend() 
  }
  
  top_gene_clusters_all <- plot_grid(plotlist = top_gene_clusters) + ggtitle("Average Expression and Percentage Expressed of Top Genes for each Cluster") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) #  + NoLegend() 
  ggsave("results/top_gene_clusters_all.jpeg", top_gene_clusters_all, width = 3000, height = 1369, units = c("px"))
  
  #return(list(umis_per_cell,genes_per_cell,mitochondrial_genes_per_cell,ribosomal_genes_per_cell, umi_vs_mito, elbow_plot, pca_plot, tsne_plot, umap_plot, umap_plot_2))
  return(list(umis_per_cell, genes_per_cell, mitochondrial_genes_per_cell, ribosomal_genes_per_cell, umi_vs_mito, plot1, elbow_plot, pca_plot, clustered_umap, clustered_tsne, heat_map_all, top_gene_clusters_all))
  
}

# Mapping Rate of Kallisto-Bustools

mapping_rate_kb<- function(kb_file) {
  
  # run_info.json
  
  # Load in the file using fromJSON since it's a json file
  library("rjson")
  
  # run_info.json
  features_kb <- fromJSON(file=kb_file)
  
  # Load the p_pseudoaligned 
  features_kb[6]
  
  # In case needed to know, this creates a dataframe of the json file
  features_kb <- as.data.frame(features_kb)
  
  # Load the actual value of p_pseudoaligned
  p_pseudoaligned <- features_kb$p_pseudoaligned
  
  return(p_pseudoaligned)
  
}

# Mapping Rate of Alevin-Fry

mapping_rate_af<- function(af_file) {
  
  # SCD-WP-g001_feature.txt
  
  # Load in the file using read.table since it's a tab delimited file
  features_txt <- read.table(file(af_file), header=TRUE)
  
  # Sum of the Mapped Reads over the sum of the Corrected Reads
  alevin_fry_mapping_rate <- sum(features_txt[1:290,3:3])/sum(features_txt[1:290,2:2])
  
  return(alevin_fry_mapping_rate)
  
}

