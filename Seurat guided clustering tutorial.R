library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#setwd('/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday0')
#setwd('/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday1')
setwd('/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday3')
#setwd('/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday5')

#XtSCI.data <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday0/raw')
#XtSCI.data <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday1/raw')
XtSCI.data <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday3/raw')
#XtSCI.data <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday5/raw')

XtSCI <- CreateSeuratObject(counts = XtSCI.data, project = 'XtSCI', min.cells = 1, min.features = 1)
XtSCI

XtSCI.data

head(XtSCI@meta.data, 5)
vln_raw <- VlnPlot(XtSCI, features = c('nFeature_RNA')) + theme(legend.position = 'none')
#ggsave("vln_raw_day1.png", vln_raw, width = 6, height = 6, dpi = 300)

scatterplot <- FeatureScatter(XtSCI, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scatterplot
ggsave("scatter_raw_day3.png", scatterplot, width = 6, height = 6, dpi = 300)

XtSCIsubset <- subset(XtSCI, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)
vln_subset <- VlnPlot(XtSCIsubset, features = c('nFeature_RNA')) + theme(legend.position = 'none')
#ggsave("vln_subset_day1.png", vln_subset, width = 6, height = 6, dpi = 300)
scatterplot_subset <- FeatureScatter(XtSCIsubset, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scatterplot_subset
ggsave("scatter_subset_day3.png", scatterplot_subset, width = 6, height = 6, dpi = 300)

merged_plot_vln <- wrap_plots(vln_raw, vln_subset, ncol = 2)
merged_plot_scatter <- wrap_plots(scatterplot, scatterplot_subset, ncol = 2)
ggsave('day3scatterrawsubset.png', merged_plot_scatter, width = 6, height = 6, dpi = 300)

XtSCInormalised <- NormalizeData(XtSCIsubset)

XtSCIvariablefeatures <- FindVariableFeatures(XtSCInormalised, selection.method = 'vst', nfeatures = 2000)
top10variable <- head(VariableFeatures(XtSCIvariablefeatures), 10)
top10variable
#write.table(top10variable, file = "XtSCI_day3_top10variable.csv", sep = ",", quote = TRUE, row.names = TRUE, col.names = TRUE)

variablefeatureplot <- VariableFeaturePlot(XtSCIvariablefeatures)
variablefeatureplot <- LabelPoints(plot = variablefeatureplot, points = top10variable, repel = TRUE)
variablefeatureplot
#ggsave("variable_feature_plot_day3.png", variablefeatureplot, width = 6, height = 6, dpi = 300)

all.genes <- rownames(XtSCIvariablefeatures)
XtSCIscaled <- ScaleData(XtSCIvariablefeatures)

XtSCIPCA <- RunPCA(XtSCIscaled)
print(XtSCIPCA[['pca']], dims = 1:10, nfeatures = 5)
pca_variablefeatures <- VizDimLoadings(XtSCIPCA, dims = 1:10, reduction = 'pca')
pcaplot <- DimPlot(XtSCIPCA, reduction = 'pca')
heatmap <- DimHeatmap(XtSCIPCA, dims = 1 : 10, cells = 300, balanced = TRUE)
#ggsave("pca_variablefeatures_day3.png", pca_variablefeatures, width = 6, height = 6, dpi = 300)
#ggsave("pcaplot_day3.png", pcaplot, width = 6, height = 6, dpi = 300)
#ggsave("pcaheatmap_day3.png", heatmap, width = 6, height = 6, dpi = 300)


XtSCIdimensionality <- JackStraw(XtSCIPCA, num.replicate = 100)
XtSCIdimensionality <- ScoreJackStraw(XtSCIdimensionality, dims = 1:10)
jackstraw <- JackStrawPlot(XtSCIdimensionality, dims = 1:10)
elbowplot <- ElbowPlot(XtSCIdimensionality)
#ggsave("jackstraw_day3.png", jackstraw, width = 6, height = 6, dpi = 300)
#ggsave("elbow_day3.png", elbowplot, width = 6, height = 6, dpi = 300)

XtSCIclustered <- FindNeighbors(XtSCIdimensionality, dims = 1:10)
XtSCIclustered <- FindClusters(XtSCIclustered, resolution = 0.5)
head(Idents(XtSCIclustered), 5)

XtSCI_nldr <- RunUMAP(XtSCIclustered, dims = 1:10)
umapplot <- DimPlot(XtSCI_nldr, reduction = 'umap')
#ggsave("umap_day5.png", umapplot, width = 6, height = 6, dpi = 300)

cluster.markers <- FindMarkers(XtSCI_nldr, ident.1 = 0, ident.2 = 4, min.pct = 0.25, logfc.threshold = 1)
head(cluster.markers, n = 5)
#write.table(cluster.markers, file = "XtSCI_day1_0_4.csv", sep = ",", quote = TRUE, row.names = TRUE, col.names = TRUE)

XtSCI.markers <- FindAllMarkers(XtSCI_nldr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)
head(XtSCI.markers, n = 5)
#write.table(XtSCI.markers, file = "XtSCI_day5_allmarkers.csv", sep = ",", quote = TRUE, row.names = FALSE, col.names = TRUE)

XtSCI.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

cluster.markers <- FindMarkers(XtSCI_nldr, ident.1 = 0, logfc.threshold = 0.25, test.use = 'roc', only.pos = TRUE)
VlnPlot(XtSCI_nldr, features = c('snca', 'mpz'))

VlnPlot(XtSCI_nldr, features = c("hbd", "hba4"), slot = "counts", log = TRUE)
FeaturePlot(XtSCI_nldr, features = c("hbd", "hba4", "hba2", "apoc1"))

XtSCI.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
marker_heatmap <- DoHeatmap(XtSCI_nldr, features = top10$gene) + NoLegend()
#ggsave("marker_heatmap_day3.png", marker_heatmap, width = 6, height = 6, dpi = 300)

levels(XtSCI_nldr)
newclusterIDs_day0 <- c('Schwann cells', 'fibroblasts', 'oligodendrocytes', 'neurons', 'radial glia cells', 'astrocytes', 'macrophages', 'endothelial cells', 'erythroid progenitor cells', 'A')
#newclusterIDs_day1 <- c('fibroblasts', 'neurons', 'astrocytes', 'oligodendrocytes', 'schwann cells', 'erythroid progenitor cells')
#newclusterIDs_day3 <- c('macrophages', 'fibroblasts_A', 'astrocytes', 'endothelial cells', 'neurons', 'erythroid progenitor cells', 'fibroblasts_B', 'schwann cells', 'A', 'oligodendrocytes', 'oligodendrocyte precursor cells')
#newclusterIDs_day5 <- c('fibroblasts', 'fibroblasts', 'astrocytes', 'skeletal muscle A', 'macrophages', 'A', 'radial glia cells', 'hepatocytes', 'melanocytes', 'Schwann cells', 'skeletal muscle B', 'neurons')
names(newclusterIDs_day0) <- levels(XtSCI_nldr)
XtSCI_nldr <- RenameIdents(XtSCI_nldr, newclusterIDs_day0)

library(SingleCellExperiment)

# Assuming you have already defined 'XtSCI_nldr', 'neurons', and 'cells.use'

library(SingleCellExperiment)

# Assuming you have already defined 'XtSCI_nldr', 'neurons', and 'cells.use'

expr <- GetAssayData(object = XtSCI_nldr, assay = "RNA", slot = "data")[, cells.use]
expr <- as.matrix(expr)  # Convert to matrix
expr_transposed <- t(expr)  # Transpose the matrix

write.csv(x = expr_transposed, file = "expression_neurons_transposed.csv", quote = FALSE, row.names = FALSE)


