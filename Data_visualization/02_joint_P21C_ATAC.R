
rm(list=ls())

library(paletteer)#
col <- paletteer_d("ggthemes::Classic_20")[c(20,6,9,19,3,18,8,7,5,1,13,2,17,15,10,11,12,16)]
col
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

rm(list=ls())

getwd()

setwd("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/ATAC/ArchR")
rm(list=ls())

library(SeuratObject)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(Matrix)

sampleNames <- 'P21C_ATAC-DBiT'
P21_C <- readRDS("/seuratobject/01_joint_P21C_RNA.rds")
P21_C

DefaultAssay(P21_C) <- 'ATAC'
P21_C
P21_C <- FindTopFeatures(P21_C, min.cutoff = 10)
P21_C <- RunTFIDF(P21_C)
P21_C <- RunSVD(P21_C)
P21_C <- FindNeighbors(
  object = P21_C,
  reduction = 'lsi',
  dims = 2:20
)
P21_C <- FindClusters(
  object = P21_C,
  algorithm = 3,
  resolution = 0.6,
  verbose = FALSE
)
P21_C <- RunUMAP(P21_C, reduction = 'lsi', dims = 2:20, reduction.name = 'umap.cut')

p2 <- DimPlot(P21_C, reduction = 'umap.cut', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")
p2



#Export the ATAC UMAP coordinates
#write.csv(P21_C@reductions$umap.cut@cell.embeddings,
#          "/data/P21_C_all_ATAC_UMAP_coordinates.csv", row.names = T, col.names = T, quote=F)

#Find different peaks
# diff.LR_peaks = FindAllMarkers(P21_C, assay = 'ATAC',test.use = "LR")
# write.csv(diff.LR_peaks, 
#           "/data/P21C_all_diff_peaks_LR.csv", row.names = F)


DefaultAssay(P21_C) <- "peaks"
P21_C <- FindTopFeatures(P21_C, min.cutoff = 5)
P21_C <- RunTFIDF(P21_C)
P21_C <- RunSVD(P21_C)


# build a joint neighbor graph using both assays
P21_C <- FindMultiModalNeighbors(
  object = P21_C,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:10, 2:20),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)


P21_C <- RunUMAP(
  object = P21_C,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE,
  reduction.name = "wnn.umap", 
  reduction.key = "wnnUMAP_"
)


P21_C<- FindClusters(P21_C, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)



######################################################################################
getwd()
setwd("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/ATAC/outs")
meta.data <- P21_C@meta.data
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names
spatial_archr_data <- read.csv(file = "/data/ArchR/P21_brain_ATAC-DBiT_gene_score_matrix.csv")
row.names(spatial_archr_data) <- spatial_archr_data$X
spatial_archr_data <- spatial_archr_data[,-1]

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = spatial_archr_data, assay = assay, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object


table(P21_C$wsnn_res.0.8)
table(P21_C$SCT_snn_res.0.8)
table(P21_C$ATAC_snn_res.0.6)

p6 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'wsnn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
p6$layers[[1]]$aes_params <- c(p6$layers[[1]]$aes_params, shape=22)
p6

ggsave("/integration/plot/Spatial-wsnn_res.0.8-UMP.png", plot = p6, width = 9, height = 9)
p7 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'SCT_snn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
p7$layers[[1]]$aes_params <- c(p7$layers[[1]]$aes_params, shape=22)
p7

ggsave("/integration/plot/Spatial-RNA_res.0.8-UMP.png", plot = p7, width = 9, height = 9)

p8 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'ATAC_snn_res.0.6', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
p8$layers[[1]]$aes_params <- c(p8$layers[[1]]$aes_params, shape=22)
p8

ggsave("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/integration/plot/ATAC_spatial-UMP.png", plot = p8, width = 9, height = 9)


p2 <- p7+p8+p6
p2
ggsave("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/integration/plot/wnn-3-spatial-UMP.png", plot = p2, width = 15, height = 9)





p1 <- DimPlot(P21_C, label = TRUE,cols = col,group.by = "SCT_snn_res.0.8",pt.size = 1.8) + NoLegend() + ggtitle("RNA UMAP")
p2 <- DimPlot(P21_C, label = TRUE,cols = col,group.by = "ATAC_snn_res.0.6",pt.size = 1.8) + NoLegend() + ggtitle("ATAC UMAP")
p3 <- DimPlot(P21_C, label = TRUE,cols = col,group.by = "wsnn_res.0.8",pt.size = 1.8) + NoLegend() + ggtitle("WNN UMAP")

p <- p1+p2+p3
p


ggsave("/integration/plot/3-UMP-wnn.png", plot = p, width = 15, height = 9)

write.csv(P21_C@reductions$umap.rna@cell.embeddings, 
          "/integration/data/P21C_RNA_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)
write.csv(P21_C@reductions$umap.cut@cell.embeddings, 
          "/integration/data/P21C_ATAC_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)
write.csv(P21_C@reductions$wnn.umap@cell.embeddings, 
          "/integration/data/P21C_wnn_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)



