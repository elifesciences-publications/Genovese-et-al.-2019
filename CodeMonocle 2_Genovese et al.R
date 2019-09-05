# ----------------------------------- load library -------------------------
library(cellrangerRkit)
library(monocle)
library(ggplot2)
library(cowplot)
library(dplyr)
library(hexbin)
library(colormap)
library(Matrix)

# ----------------------------------- Data loading -------------------------
matrix_dir = "~/IBDM/Single cell RNA-seq/poxpros/OutsFiles/outs/filtered_gene_bc_matrices/Drosophila_melanogaster.BDGP6/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "genes.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1



my_dir <- "~/IBDM/Single cell RNA-seq/poxpros/OutsFiles/"
gbm <- load_cellranger_matrix(my_dir)


# -----------------------------------Cell Ranger data conversion -------------------------
dim(exprs(gbm))
dim(pData(gbm))
head(fData(gbm))
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = fData(gbm)),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
#warning

my_feat <- fData(gbm)
names(my_feat) <- c('id', 'gene_short_name')
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = my_feat),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


# ----------------------------------- Normalization & variance estimation -------------------------
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)


# ----------------------------------- UMI & gene expression -------------------------
my_cds <- detectGenes(my_cds, min_expr = 0.1)
head(fData(my_cds))
summary(fData(my_cds)$num_cells_expressed)
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))


# ----------------------------------- Pseudotime Computation -------------------------
head(pData(my_cds))
ggplot(pData(my_cds), aes(num_genes_expressed, UMI)) + geom_point()
disp_table <- dispersionTable(my_cds)
head(disp_table)
table(disp_table$mean_expression>=0.1)


#Ordering cells using known marker genes

Imp_id <- row.names(subset(fData(my_cds), gene_short_name == "Imp"))
Eip93F_id <- row.names(subset(fData(my_cds), gene_short_name == "Eip93F"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Imp+tNBs", classify_func=function(x) {x[Imp_id,] >= 1})
cth <- addCellType(cth, "E93+tNBs", classify_func=function(x) {x[Eip93F_id,] >=1})
cth <- addCellType(cth, "Reserve cell", classify_func=function(x) {x[Eip93F_id,] == 0 & x[Imp_id,] == 0})
my_cds <- classifyCells(my_cds, cth)

#Now we select the set of genes that co-vary (in either direction) with these two "bellweather" genes:
my_cds_expressed_genes <-  row.names(subset(fData(my_cds),
                                            num_cells_expressed >= 10))
marker_diff <- markerDiffTable(my_cds[my_cds_expressed_genes,],
                               cth,
                               cores = 1)

#semisup_clustering_genes <-
#row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <-
  row.names(marker_diff)[order(marker_diff$qval)][1:100]



my_cds <- setOrderingFilter(my_cds, semisup_clustering_genes)

#plot_ordering_genes(my_cds)
my_cds <- reduceDimension(my_cds, max_components = 2,
                          method = 'DDRTree', norm_method = 'log')
my_cds <- orderCells(my_cds)

plot_cell_trajectory(my_cds, color_by = "CellType") +
  theme(legend.position = "right")

# to reverse trajectory
  my_cds <- orderCells(my_cds, reverse = TRUE)
  
  
# ------------------------------ GRAPHIQUES -----------------------------------------
  
  #gene curves along pseudotime
  to_be_tested <- row.names(subset(fData(my_cds),
                                   gene_short_name %in% c("Imp", "Eip93F", "Oatp74D", "Pgi", "Hex-A", "Pfk", "Ald", "Tpi", "Gapdh1", "Gapdh2", "Pgk", "CG9961", "Pglym78", "Eno", "PyK", "ImpL3" )))
  cds_subset <- my_cds[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  plot_genes_in_pseudotime(cds_subset, color_by = "State") +   theme(legend.position="none")
  
  
  #for plot_genes_branched_pseudotime
  #use: plot_genes_branched_pseudotime(cds_subset, branch_point = 4, color_by = "State") +   theme(legend.position="none")
  
  
  
  #FBgn0264490          Eip93F  0.000000e+00  0.000000e+00
  #FBgn0039431            plum  2.164841e-48  2.783367e-48
  #FBgn0003525             stg  3.925092e-33  4.415728e-33
  #FBgn0035499           Chd64 5.993983e-146 8.990974e-146
  #FBgn0036732         Oatp74D  0.000000e+00  0.000000e+00
  #FBgn0010620         CG10939  0.000000e+00  0.000000e+00
  #FBgn0285926             Imp  0.000000e+00  0.000000e+00
  #FBgn0030665         CG15646 5.460538e-281 9.828969e-281
  #FBgn0086758          chinmo  2.790154e-05  2.790154e-05#
  
  #FBgn0000064             Ald 8.545422e-157 5.127253e-156
  #FBgn0014869         Pglym78 2.916335e-150 1.166534e-149
  #FBgn0086355             Tpi  2.984726e-42  5.969453e-42
  #FBgn0001258           ImpL3  2.414792e-29  3.622188e-29
  #FBgn0001091          Gapdh1  2.141179e-76  5.138829e-76
  #FBgn0003074             Pgi  1.439659e-01  1.727591e-01
  #FBgn0003071             Pfk  1.346085e-21  1.794780e-21
  #FBgn0001186           Hex-A  1.774405e-01  1.935715e-01
  #FBgn0001092          Gapdh2  0.000000e+00  0.000000e+00
  #FBgn0000579             Eno 1.133611e-127 3.400832e-127
  #FBgn0031451          CG9961  7.760592e-01  7.760592e-01
  #FBgn0250906             Pgk  6.090158e-33  1.044027e-32
  
  #differential gene expression by states with violin plots
  to_be_tested <- row.names(subset(fData(my_cds),
                                   gene_short_name %in% c("Imp")))
  cds_subset <- my_cds[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(State)")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  plot_genes_violin(cds_subset, grouping = "State", color_by = "State")
  
  
  #plot_cell_trajectory --> color by state
  plot_cell_trajectory(my_cds, color_by = "State", show_branch_points = FALSE) + theme(legend.position = "right")
  #plot_cell_trajectory --> color for gene expression with gradient 
  plot_cell_trajectory(my_cds, markers = c("PCNA"), show_tree = FALSE, use_color_gradient = TRUE, show_branch_points = FALSE, show_backbone = FALSE, cell_size = 1.5) + scale_colour_gradient2(low = "#f9f7f4", high = "#e53b06", mid = "#c4c4c4") + theme(legend.position = "right")
  
  
  
  
  # cluster the top 200 genes that vary as a function of pseudotime on a heatmap
  my_pseudotime_de <- differentialGeneTest(my_cds,
                                           fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                           cores = 8)
  my_pseudotime_de %>% arrange(qval) %>% head(200) %>% select(id) -> gene_to_cluster
  gene_to_cluster <- gene_to_cluster$id
  
  my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds[gene_to_cluster,],
                                                   num_clusters = 3,
                                                   cores = 8,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)
  
  # to collect gene names for each clusters
  t <- plot_pseudotime_heatmap(my_cds[gene_to_cluster,], return_heatmap=T)
  t <- as.data.frame(cutree(t$tree_row, k=3))
  colnames(t) <- "Cluster"
  t$Gene <- rownames(t)
  dplyr::filter
  filter(t)
  
  #Analyzing Branches in Single-Cell Trajectories (BEAM)
  BEAM_res <- BEAM(my_cds, branch_point = 4, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  my_branched_heatmap <- plot_genes_branched_heatmap(my_cds[row.names(subset(BEAM_res,
                                                                             qval < 1e-20)),],
                                                     branch_point = 4,
                                                     num_clusters = 4,
                                                     cores = 1,
                                                     use_gene_short_name = T,
                                                     show_rownames = T,
                                                     return_heatmap = TRUE)
  print(my_branched_heatmap$annotation_row)
  