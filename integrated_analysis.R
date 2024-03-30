library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(monocle3)
library(scran)
library(SingleCellExperiment)
library(SeuratWrappers)
library(slingshot)

setwd('/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday0-5')

XtSCI.data0 <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday0/raw')
XtSCI.data1 <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday1/raw')
XtSCI.data3 <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday3/raw')
XtSCI.data5 <- Read10X(data.dir = '/home/ethanlewisbaird/working directory A/spinal cord scRNAseq/XtSCI_raw/XtSCIday5/raw')

XtSCI0 <- CreateSeuratObject(counts = XtSCI.data0, project = 'XtSCI', min.cells = 1, min.features = 1)
XtSCI1 <- CreateSeuratObject(counts = XtSCI.data1, project = 'XtSCI', min.cells = 1, min.features = 1)
XtSCI3 <- CreateSeuratObject(counts = XtSCI.data3, project = 'XtSCI', min.cells = 1, min.features = 1)
XtSCI5 <- CreateSeuratObject(counts = XtSCI.data5, project = 'XtSCI', min.cells = 1, min.features = 1)

XtSCI0 <- AddMetaData(object = XtSCI0, metadata = as.list(rep(0, 727)), col.name = 'timepoint')
XtSCI1 <- AddMetaData(object = XtSCI1, metadata = as.list(rep(1, 373)), col.name = 'timepoint')
XtSCI3 <- AddMetaData(object = XtSCI3, metadata = as.list(rep(3, 532)), col.name = 'timepoint')
XtSCI5 <- AddMetaData(object = XtSCI5, metadata = as.list(rep(5, 1033)), col.name = 'timepoint')

XtSCI_list <- c(XtSCI0, XtSCI1, XtSCI3, XtSCI5)

XtSCI_list <- lapply(X = XtSCI_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = XtSCI_list)

anchors <- FindIntegrationAnchors(object.list = XtSCI_list, anchor.features = features)

XtSCImerged <- IntegrateData(anchorset = anchors)

#count_matrix <- GetAssayData(object = subset, assay = "RNA", slot = "data")
#write.csv(count_matrix, "microglia_count_matrix.csv")

DefaultAssay(XtSCImerged) <- 'integrated'

XtSCIscaled <- ScaleData(XtSCImerged, verbose = FALSE)
XtSCIPCA <- RunPCA(XtSCIscaled, npcs = 10, verbose = FALSE)
XtSCIUMAP <- RunUMAP(XtSCIPCA, reduction = 'pca', dims = 1:10)
XtSCIneighbors <- FindNeighbors(XtSCIUMAP, reduction = 'pca', dims = 1:10)
XtSCIclusters <- FindClusters(XtSCIneighbors, resolution = 0.5)

levels(XtSCIclusters)
newclusterIDs_merged <- c('Microglia', 'Schwann', 'Neural_stem_cells', 'Epithelial_cells', 'Fibroblasts', 'Neurons', 'Microglia', 'Macrophages', 'Oligodendrocytes', 'Microglia', 'Endothelial_cells', 'Erythroid_cells', 'NK_cells', 'Melanocytes', 'Muscle')
names(newclusterIDs_merged) <- levels(XtSCIclusters)
XtSCIclusters <- RenameIdents(XtSCIclusters, newclusterIDs_merged)

merged_count_matrix <- GetAssayData(object = XtSCIclusters, assay = "integrated", slot = "data")
write.csv(merged_count_matrix, "merged_count_matrix.csv")

umapplot <- DimPlot(XtSCIclusters, reduction = 'umap', group.by = 'pseudotime', label = FALSE) + NoLegend()
umapplot
ggsave("umap_merged_split.png", umapplot, width = 10, height = 6, dpi = 300)

FeaturePlot(XtSCIclusters, features = 'pseudotime', label = FALSE)

XtSCI.markers <- FindAllMarkers(XtSCIclusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)
head(XtSCI.markers, n = 5)
write.table(XtSCI.markers, file = "XtSCI_merged_allmarkers.csv", sep = ",", quote = TRUE, row.names = FALSE, col.names = TRUE)

# Feature plot - visualize feature expression in low-dimensional space
features = c('rbp4')
FeaturePlot(subset, features = features, split.by = 'timepoint', label = TRUE)
DotPlot(XtSCIclusters, features = features) + RotatedAxis()
vln_snca <- VlnPlot(XtSCIclusters, features = features, split.by = 'timepoint', log = TRUE)
vln_snca
RidgePlot(XtSCIclusters, features = features, group.by = 'timepoint', log = TRUE)
DoHeatmap(subset(XtSCIclusters, downsample = 100), size = 3)
ggsave("vln_snca.png", vln_snca, width = 10, height = 6, dpi = 300)
#cluster.markers <- FindMarkers(XtSCIclusters, ident.1 = '0', ident.2 = 'Oligodendrocytes', min.pct = 0.25, logfc.threshold = 1)
#head(cluster.markers, n = 5)
#write.table(cluster.markers, file = "XtSCI_merged_0_4.csv", sep = ",", quote = TRUE, row.names = TRUE, col.names = TRUE)

cluster_subset <- c('Microglia', 'Neural_stem_cells', 'Epithelial_cells', 'Fibroblasts', 'Neurons', 'Macrophages', 'Oligodendrocytes', 'Endothelial_cells', 'Erythroid_cells', 'NK_cells', 'Melanocytes', 'Muscle')
timepoint_subset <- 0
#cluster_subset <- c('NK_cells', 'Macrophages')

# Replace "Timepoint1" with the specific timepoint you want to subset for
timepoint_0 <- '0'
timepoint_5 <- '5'

# Subset the Seurat object based on the selected timepoint
subset_seurat_object <- seurat_object[seurat_object$Timepoint == subset_timepoint, ]

Idents(XtSCIclusters) <- XtSCIclusters@meta.data$'timepoint'

subset <- subset(x = XtSCIclusters, idents = cluster_subset)
subset <- subset(subset, idents <- timepoint_subset)
variablefeatures <- FindVariableFeatures(XtSCIclusters, selection.method = 'vst', nfeatures = 2000)
scaled <- ScaleData(variablefeatures, verbose = FALSE)
pca <- RunPCA(scaled, npcs = 10, verbose = FALSE)
umap <- RunUMAP(pca, reduction = 'pca', dims = 1:10)
neighbors <- FindNeighbors(umap, reduction = 'pca', dims = 1:10)
clusters <- FindClusters(neighbors, resolution = 0.05)
umapplot <- DimPlot(clusters, reduction = 'umap', label = TRUE) + NoLegend()
umapplot
#ggsave("neural_cells_umap.png", umapplot, width = 10, height = 6, dpi = 300)
subset.markers <- FindMarkers(clusters, ident.1 = '0', ident.2 = '3', min.pct = 0, logfc.threshold = 1, only.pos = TRUE)
head(subset.markers, n = 5)
write.table(subset.markers, file = "markers_0_5.csv", sep = ",", quote = TRUE, row.names = TRUE, col.names = TRUE)

levels(clusters)
newIDs_subset <- c('Neural_stem_cells', 'Neurons')
names(newIDs_subset) <- levels(clusters)
clusters <- RenameIdents(clusters, newIDs_subset)

# XtSCIclusters$celltype_timepoint <- paste(Idents(XtSCIclusters), XtSCIclusters$timepoint, sep = '_')
# Idents(XtSCIclusters) <- 'celltype_timepoint'
# 
# NK_cells_zero_five <- FindMarkers(XtSCIclusters, ident.1 = 'Neural_stem_cells_0', ident.2 = 'Neural_stem_cells_1', logfc.threshold = 1)
# write.table(NK_cells_zero_five, file = "Neural_stem_cells_zero_one_markers.csv", sep = ",", quote = TRUE, row.names = TRUE, col.names = TRUE)

expressiontable <- AverageExpression(XtSCIclusters)
head(expressiontable)
RNAexpression <- expressiontable$RNA
integratedexpression <- expressiontable$integrated

write.table(RNAexpression, file = "rnaaverage_expression.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(integratedexpression, file = "integaverage_expression.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

markers_to_plot <- c('tmem119', 'rbp4', 'col1a1', 'mpz', 'snca', 'pmp22', 'krt18.1', 'fabp7', 'tecta', 'csta', 'epcam', 'cdh6', 'col3a1', 'serpinf1', 'elavl4', 'stmn2', 'gap43', 'cd74', 'lgals1.1', 'ctsb', 'cldn11', 'ablim3', 'plp1', 'hbd', 'hba4', 'hbg2', 'gzmh', 'ctsh', 'ccl11', 'tyrp1', 'pmel', 'mlana', 'mylpf', 'actc1', 'myl1')
DotPlot(XtSCIclusters, features = rev(markers_to_plot), cols = c('blue', 'red', 'green', 'yellow'), dot.scale = 8, split.by = 'timepoint') + RotatedAxis()

selected_timepoints <- c(0, 5)
XtSCIclusters_subset <- subset(XtSCIclusters, timepoint %in% selected_timepoints)
DotPlot(XtSCIclusters_subset, features = rev(markers_to_plot), cols = c('blue', 'red'), dot.scale = 8, split.by = 'timepoint') + RotatedAxis()

timepoints <- c(0, 1, 3, 5)
cluster_assignments <- Idents(XtSCIclusters)

cluster_abundance <- data.frame(timepoint = character(),
                                cluster = character(),
                                relative_abundance = numeric(),
                                stringsAsFactors = FALSE)

for (time in timepoints) {
  total_cells_at_timepoint <- sum(XtSCIclusters$timepoint == time)
  for (cluster in unique(cluster_assignments)) {
    abundance <- sum(cluster_assignments == cluster & XtSCIclusters$timepoint == time)
    relative_abundance <- abundance / total_cells_at_timepoint
    cluster_abundance <- rbind(cluster_abundance, data.frame(timepoint = time,
                                                             cluster = cluster,
                                                             relative_abundance = relative_abundance))
  }
}

cluster_abundance$cluster <- factor(cluster_abundance$cluster, levels = unique(cluster_abundance$cluster))

relative_abundance <- ggplot(cluster_abundance, aes(x = factor(timepoint), y = relative_abundance, fill = cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Timepoint", y = "Relative Abundance", title = "Cluster Relative Abundance Over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Cell Type")
relative_abundance
ggsave("relative_abundance.png", relative_abundance, width = 10, height = 6, dpi = 300)

cluster_ids <- Seurat::Idents(object = clusters)
head(cluster_ids)

cds <- as.cell_data_set(XtSCIclusters)

seurat_clusters <- XtSCIclusters$seurat_clusters
colData(cds)$seurat_clusters <- seurat_clusters

cds <- clusterCells(cds, assay.type = 'logcounts')
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
cds_plot <- plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds_plot

# Perform differential gene expression analysis

#Extract data, phenotype data, and feature data from the SeuratObject
timepoint_subset <- c('0', '5')
subset <- subset(x = XtSCIclusters, idents = timepoint_subset)
data <- as(as.matrix(XtSCIclusters@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = XtSCIclusters@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~timepoint", cores = 8)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Hours")

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))
plot_cell_trajectory(cds, color_by = "Pseudotime")


XtSCIclusters$pseudotime <- cds$Pseudotime
FeaturePlot(data, features = "pseudotime")

XtSCIclusters_cds <- exportCDS(cds, export_to = 'Seurat', export_all = T)

ggsave("cds_plot.png", cds_plot, width = 8, height = 6, dpi = 300)


get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
# cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
# 
# plot_cells(cds,
#            color_cells_by = "pseudotime",
#            group_cells_by = "cluster",
#            label_cell_groups = FALSE,
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            label_roots = FALSE,
#            trajectory_graph_color = "grey60",
#            cell_size = 1.5)

cds <- order_cells(cds)
nscs_neurons <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
ggsave("nscs_neurons.png", nscs_neurons, width = 8, height = 6, dpi = 300)

cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(clusters[["RNA"]])


expression_matrix <- counts(cds)  # Get the expression matrix
cell_metadata <- colData(cds)      # Get the cell metadata
seurat_obj <- CreateSeuratObject(counts = expression_matrix, meta.data = cell_metadata)

# Assuming you have a SingleCellExperiment object named 'cds' and you want to add metadata
metadata_values <- clusters[["RNA"]][, "features"]  # Replace with your metadata values
colData(cds)$gene_short_name <- metadata_values

cds <- AddMetaData(cds, metadata = clusters[["RNA"]][, "gene_short_name"], col.name = "gene_short_name")

gene_short_names <- clusters[["RNA"]][, "gene_short_name"]
cds$RNA$gene_short_name <- gene_short_names

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)
head(cds_graph_test_results)

XtSCIclusters <- AddMetaData(
  object = XtSCIclusters,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "npcs_neurons"
)

nscs_neurons2 <- FeaturePlot(XtSCIclusters, c("npcs_neurons"), pt.size = 0.1) & scale_color_viridis_c()
ggsave("nscs_neurons2.png", nscs_neurons2, width = 8, height = 6, dpi = 300)
