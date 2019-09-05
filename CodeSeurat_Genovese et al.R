#---------------------------- load libraries -------------------------------
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

#---------------------------- load data ----------------------------------------
## Need to open the outs filtered gene matrices from 10X 
CTRL.data <- Read10X(data.dir = "~/IBDM/Single cell RNA-seq/poxpros/OutsFiles//outs/filtered_gene_bc_matrices/Drosophila_melanogaster.BDGP6")

#---------------------------- Setup the Seurat Object --------------------------------
dense.size <- object.size(x = as.matrix(x = CTRL.data))
dense.size

sparse.size <- object.size(x = CTRL.data)
sparse.size

dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
CTRL <- CreateSeuratObject(counts = CTRL.data, min.cells = 0, min.features = 200, project = "Control")
CTRL


#-------------- QC and selecting cells for further analysis---------------------
CTRL[["percent.mt"]] <- PercentageFeatureSet(CTRL, pattern = "^mt:")
CTRL.QC1 <- VlnPlot(CTRL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("1-Seurat_First_Steps/CTRL_VlnPlot_MitoGenes.pdf", width = 8.5, height = 6)
CTRL.QC1
dev.off()
CTRL.QC1


plot1 <- FeatureScatter(CTRL, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTRL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CTRL.QC2 <- CombinePlots(plots = list(plot1, plot2))
pdf("1-Seurat_First_Steps/CTRL_GenePlot.pdf", width = 11.4, height = 8)
CTRL.QC2
dev.off()
CTRL.QC2

CTRL <- subset(CTRL, subset = nFeature_RNA > 200 & nFeature_RNA < 3700 & percent.mt < 5)


#---------------------Normalizing the data-----------------------------
CTRL <- NormalizeData(CTRL, normalization.method = "LogNormalize", scale.factor = 10000)


#------------------------- Detection of variable genes across the single cells --------------------------
CTRL <- FindVariableFeatures(CTRL, selection.method = "vst", nfeatures = 2000)
CTRL.top10var <- head(VariableFeatures(CTRL), 10)
CTRL.top10var

CTRL.plot1 <- VariableFeaturePlot(CTRL)
CTRL.plot2 <- LabelPoints(plot = CTRL.plot1, points = CTRL.top10var, repel = TRUE)
CTRL.VarGenes <- CombinePlots(plots = list(CTRL.plot1, CTRL.plot2))
pdf("1-Seurat_First_Steps/CTRL_VariableGenes.pdf", width = 11.4, height = 8)
CTRL.VarGenes
dev.off()
CTRL.VarGenes

#---------------------Remove cell-cycling genes from the analysis-------------------------
## Need to open the files containing the list of cell cycling genes
s.genes <- readLines(con = "~/IBDM/Single cell RNA-seq/analysis with R/For Genovese et al/Seurat UMAP/cell cycle genes to regress/S_Phase_Genes_Droso.csv")
g2m.genes <- readLines(con = "~/IBDM/Single cell RNA-seq/analysis with R/For Genovese et al/Seurat UMAP/cell cycle genes to regress/G2-M_Phase_Genes_Droso.csv")

CTRL <- CellCycleScoring(CTRL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = CTRL@meta.data)


#----------------- Scaling the data and removing unwanted sources of variation ---------------------------
all.genes <- rownames(CTRL)
CTRL <- ScaleData(object = CTRL, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mt","S.Score", "G2M.Score"), display.progress = FALSE)


#-------------------- Perform linear dimensional reduction ---------------------
CTRL <- RunPCA(CTRL, features = VariableFeatures(object = CTRL), nfeatures.print = 10, ndims.print = 1:5)

VizDimLoadings(CTRL, dims = 1:2, reduction = "pca")
CTRL.CellCyclingRegression <- DimPlot(CTRL, reduction = "pca")
CTRL.CellCyclingRegression


#------------------PCheatmap 12 first PCs and custom color palette-----------------------
myPalette <- CustomPalette(low = "#195a8c", high = "#e0891f", mid = "#e9eae3")
myPalette
PCHeatmap(object = poxpros, pc.use = 1:12, cells.use = 500, col = myPalette(),  do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


myPalette <- CustomPalette(low = "#195a8c", high = "#e0891f", mid = "#e9eae3")
CTRL.PCA <- DimHeatmap(CTRL, dims = 1:12, cells = 500, balanced = TRUE, ncol = 3)

CTRL.PCA.1 <- DimHeatmap(CTRL, dims = 1, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.2 <- DimHeatmap(CTRL, dims = 2, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.3 <- DimHeatmap(CTRL, dims = 3, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.4 <- DimHeatmap(CTRL, dims = 4, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.5 <- DimHeatmap(CTRL, dims = 5, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.6 <- DimHeatmap(CTRL, dims = 6, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.7 <- DimHeatmap(CTRL, dims = 7, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.8 <- DimHeatmap(CTRL, dims = 8, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.9 <- DimHeatmap(CTRL, dims = 9, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.10 <- DimHeatmap(CTRL, dims = 10, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.11 <- DimHeatmap(CTRL, dims = 11, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
CTRL.PCA.12 <- DimHeatmap(CTRL, dims = 12, cells = 500, balanced = TRUE, fast = FALSE) + scale_fill_gradientn(colors = myPalette)
p1 <- plot(CTRL.PCA.1, col = myPalette())
p2 <- plot(CTRL.PCA.2, col = myPalette())
p3 <- plot(CTRL.PCA.3, col = myPalette())
p4 <- plot(CTRL.PCA.4, col = myPalette())
p5 <- plot(CTRL.PCA.5, col = myPalette())
p6 <- plot(CTRL.PCA.6, col = myPalette())
p7 <- plot(CTRL.PCA.7, col = myPalette())
p8 <- plot(CTRL.PCA.8, col = myPalette())
p9 <- plot(CTRL.PCA.9, col = myPalette())
p10 <- plot(CTRL.PCA.10, col = myPalette())
p11 <- plot(CTRL.PCA.11, col = myPalette())
p12 <- plot(CTRL.PCA.12, col = myPalette())
p1 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 1")
p2 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 2")
p3 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 3")
p4 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 4")
p5 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 5")
p6 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 6")
p7 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 7")
p8 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 8")
p9 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 9")
p10 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 10")
p11 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 11")
p12 + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 12")

p12 + theme(legend.text = element_text(hjust = 1), axis.text.y.left = element_text(hjust = 0)) + labs(title = "PC 12")


#-----------------Determine statistically significant principal components----------------------
CTRL <- JackStraw(CTRL, dims = 50, num.replicate = 100, verbose =TRUE)
CTRL <- ScoreJackStraw(CTRL, dims = 1:45)
JackStrawPlot(CTRL, dims = 1:45)
ElbowPlot(CTRL, ndims = 45)


#-----------Cluster the cells ---------------
CTRL <- FindNeighbors(CTRL, dims = 1:20)
CTRL <- FindClusters(CTRL, resolution = 0.4)
head(Idents(CTRL), 5)

CTRL <- RunUMAP(CTRL, dims = 1:20)
CTRL.UMAP <-DimPlot(CTRL, reduction = 'umap', pt.size = 0.4, cols = c('0' = '#c384d1', '1' = '#6d95bf', '2' = '#f74848', '3' = '#48cef7', '4' = 'blue', '5' = 'green', '6' = 'brown' ))
CTRL.UMAP


#----------------------------------Save UMAP and computation---------------------------------
saveRDS(CTRL, file = "3-Seurat_Clusterisation_on_UMAP/CTRL_Analysis.rds")


#-------------------Finding differentially expressed genes (cluster biomarkers)--------------------
# find markers for every cluster compared to all remaining cells, report only the positive ones
CTRL.markers <- FindAllMarkers(CTRL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTRL.top50 <- CTRL.markers %>% group_by(cluster) %>% top_n(n = Inf, wt = avg_logFC)
write.csv(CTRL.top50, file = "3-Seurat_Clusterisation_on_UMAP/CTRL_AllDiffExpGenes_per_cluster.csv", row.names = FALSE)

#heatmap differentially expressed genes for all clusters
CTRL.top10 <- CTRL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(CTRL, features = CTRL.top10$gene) + NoLegend() 
PlotCluster <- DoHeatmap(CTRL, features = CTRL.top10$gene) + scale_fill_gradientn(colors = myPalette)
PlotCluster + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0))

#heatmap differentially expressed genes for all clusters
CTRL.top15 <- CTRL.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
DoHeatmap(CTRL, features = CTRL.top15$gene) + NoLegend() 
PlotCluster <- DoHeatmap(CTRL, features = CTRL.top15$gene) + scale_fill_gradientn(colors = myPalette)
PlotCluster + theme(legend.position = "none", axis.text.y.left = element_text(hjust = 0))

# find markers for cluster 2 compared to all remaining cells
CTRL.cluster2.markers <- FindMarkers(CTRL, ident.1 = 2, only.pos = TRUE, min.pct = 0.25)
head(CTRL.cluster2.markers, n = 20)
View(CTRL.cluster2.markers)
write.csv(CTRL.cluster2.markers, file = "3-Seurat_Clusterisation_on_UMAP/CTRL_TopMarkers_of_Cluster3.csv")

CTRL.cluster2.markers <- FindMarkers(CTRL, ident.1 = 2, ident.2 = c(0, 1, 3, 4), min.pct = 0.25)
head(CTRL.cluster2.markers, n = 100)
write.csv(CTRL.cluster2.markers, file = "~/IBDM/Single cell RNA-seq/analysis with R/For Genovese et al/Seurat UMAP/CTRL_TopMarkers_of_Cluster2vsCluster0,1,3,4.csv")




#---------------Save MetaData----------------------------
CTRLmetaData <- CTRL@meta.data
write.csv(CTRLmetaData, file = "5-Metadata/CTRL_MetaData.csv")


#------------------------------number of cells per cluster--------------------
CTRL.NumbCellsCluster <- table(CTRLmetaData$seurat_clusters)
CTRL.NumbCellsCluster
write.csv(CTRL.NumbCellsCluster, file = "3-Seurat_Clusterisation_on_UMAP/CTRL_Number_of_Cells_per_Cluster.csv", row.names = FALSE)


#---------------- Distribution of genes on UMAP----------------
FeaturePlot(CTRL, features = c("Imp"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("mira"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Oatp74D"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Eip93F"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("plum"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("stg"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Syp"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG15646"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Chd64"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("chinmo"), cols = c("#dbd5d2","#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("lin-28"), cols = c("#dbd5d2","#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG10939"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("lncRNA:CR43334"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG5953"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("SP1173"), cols = c("#dbd5d2", "#b32c18"), pt.size = 0.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

FeaturePlot(CTRL, features = c("dpn"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("grh"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("VAChT"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("VGlut"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Gad1"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("repo"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Gadd45"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Irbp18"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Xrp1"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("GstE6"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

FeaturePlot(CTRL, features = c("abd-A"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Abd-B"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Rcd-1r"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("kek1"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("bi"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Optix"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("glec"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Hsp68"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Hsp70Bc"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Pvf3"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("robo2"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("cas"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("svp"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("D"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("trx"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("E(z)"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("esc"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Su(z)12"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Psc"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("E23"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Ubx"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("wg"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("en"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("in"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("ci"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("tsh"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("bru3"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Thor"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("hig"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("nAChRalpha5"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("nAChRalpha3"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG2852"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG43088"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CG4250"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Set1"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("trr"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Mnn1"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("CBP"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Ten-a"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Antp"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("Kr"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("nub"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("hb"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
FeaturePlot(CTRL, features = c("elav"), cols = c("#dbd5d2", "#b32c18"), pt.size = 2.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())


FeaturePlot(CTRL, features = c("PCNA", "Mcm2", "stg", "Cdk1", "CycB"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Imp", "lin-28", "abd-A", "Abd-B", "Ubx", "Antp"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Imp", "lin-28", "grim", "rpr", "hid", "Diap1"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Imp", "lin-28", "Psc", "Su(z)12", "esc", "Pc"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Imp", "lin-28", "cas", "D", "svp", "Kr"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("mira", "dpn", "grh", "nab", "wor", "ase"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Gadd45", "Irbp18", "Xrp1", "GstE6"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("VAChT", "VGlut", "Gad1", "repo", "elav"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("trx", "E(z)", "Set1", "trr", "CBP", "brm", "Mnn1"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("Imp", "lin-28", "chinmo", "E23", "Eip93F", "Syp"), cols = c("#dbd5d2", "#e0681e"))
FeaturePlot(CTRL, features = c("pros", "Snr1", "aPKC", "numb"), cols = c("#dbd5d2", "#e0681e"))


#--------------Expression of two genes on the same UMAP-------------------------------
FeaturePlot(CTRL, features = c("E23", "lin-28"), blend = TRUE, combine = FALSE, blend.threshold = 0.5, pt.size = 2)


#----------Label clusters on UMAP-----------------
new.cluster.ids <- c(" ", " ", " ", "Imp+", "Neurons", "Stressed cells")
names(new.cluster.ids) <- levels(CTRL)
CTRL <- RenameIdents(CTRL, new.cluster.ids)
DimPlot(CTRL, reduction = "umap", label = TRUE) + NoLegend()

#----------------------------------Save UMAP and computation---------------------------------
saveRDS(CTRL, file = "3-Seurat_Clusterisation_on_UMAP/CTRL_Analysis_withLabels.rds")

#-----------number of single, double or more positive cells----------------------------
E23.cutoff <- 0
lin28.cutoff <- 0
cas.cutoff <- 0
length(which(FetchData(CTRL, vars='lin-28', slot = "data" ) > lin28.cutoff))
length(which(FetchData(CTRL, vars='E23') > E23.cutoff))
length(which(FetchData(CTRL, vars='cas') > cas.cutoff))
length(which(FetchData(CTRL, vars='lin-28') > lin28.cutoff & FetchData(CTRL, vars='E23') > E23.cutoff & FetchData(CTRL, vars='cas') > cas.cutoff))
length(which(FetchData(CTRL, vars='lin-28') > lin28.cutoff & FetchData(CTRL, vars='E23') > E23.cutoff & FetchData(CTRL, vars='cas') == cas.cutoff))
length(which(FetchData(CTRL, vars='lin-28') > lin28.cutoff & FetchData(CTRL, vars='cas') > cas.cutoff & FetchData(CTRL, vars='E23') == E23.cutoff))
length(which(FetchData(CTRL, vars='cas') > cas.cutoff & FetchData(CTRL, vars='E23') > E23.cutoff & FetchData(CTRL, vars='lin-28') == lin28.cutoff))
