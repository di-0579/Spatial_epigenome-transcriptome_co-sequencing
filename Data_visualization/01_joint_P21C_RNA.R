rm(list=ls())


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

getwd()

setwd("./P21C_20um/")
rm(list=ls())

library(SeuratObject)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(Matrix)

sampleNames <- 'P21C_ATAC-DBiT'


# load the RNA and ATAC data

data <- read.table(file = "P21_ATAC-DBiT_matrix.tsv",  header = TRUE, row.names = 1, as.is = TRUE)
data <- as.data.frame(t(data))

data$V1 <- paste0(row.names(data))
data$V1 <- paste0(data$V1, "-1")
data$V1
row.names(data) <- data$V1
data <- data[,!(colnames(data) %in% c("V1"))]

RNA_counts <- as.data.frame(t(data))

P21_C <- CreateSeuratObject(counts = RNA_counts, assay = "RNA")

meat.data <- P21_C@meta.data

#############################################

data1 <- read.table("position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(data1[1,])
x = x[-1]
data1 = data1[-1]
data1 = as.data.frame(t(data1))

my_data <- read.table(file = 'spatial_barcodes.txt', 
                      sep = '\t', stringsAsFactors=FALSE)
my_data$V4 <- paste0(my_data$V2,"x",my_data$V3)
row.names(my_data) <- my_data$V4
my_data_in_tissue <- my_data[data1$V1, ]
my_data_in_tissue$V1 <- paste0(my_data_in_tissue$V1, "-1")
bc_in_tissue <- my_data_in_tissue$V1 
###############
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

P21_C_counts <- Read10X_h5("./P21C_20um/ATAC/outs/raw_peak_bc_matrix.h5")
P21_C_counts <- P21_C_counts[,bc_in_tissue]

P21_C[["ATAC"]] <- CreateChromatinAssay(
  counts = P21_C_counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = "./P21C_20um/ATAC/outs/fragments.tsv.gz",
  min.cells = 0
)

Annotation(P21_C[["ATAC"]]) <- annotations

###############################
P21_C
DefaultAssay(P21_C) <- "ATAC"

P21_C <- NucleosomeSignal(P21_C)
P21_C <- TSSEnrichment(P21_C)



P21_C$blacklist_fraction <- FractionCountsInRegion(
  object = P21_C,
  assay = 'ATAC',
  regions = blacklist_mm10
)


# call peaks using MACS2
peaks <- CallPeaks(P21_C, macs2.path = "./conda_envs/ArchR/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(P21_C),
  features = peaks,
  cells = colnames(P21_C)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
P21_C[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "./P21C_20um/ATAC/outs/fragments.tsv.gz",
  genome = "mm10"
)
Annotation(P21_C[["peaks"]]) <- annotations


#########################################

DefaultAssay(P21_C) <- "RNA"


P21_C<- FindVariableFeatures(P21_C, nfeatures = 3000)
P21_C <- SCTransform(P21_C)
P21_C <- RunPCA(P21_C, assay = "SCT", verbose = FALSE)

P21_C <- RunUMAP(P21_C, dims = 1:10, reduction.name = "umap.rna")
P21_C <- FindNeighbors(P21_C, dims = 1:10)
P21_C <- FindClusters(P21_C, resolution = 0.8, algorithm = 3)


table(P21_C$seurat_clusters)

p1 <- DimPlot(P21_C, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")
p1 

#ggsave("jihe/RNA_UMP-3.png", plot = p1, width = 9, height = 9)
# Export the RNA UMAP coordinates
#write.csv(P21_C@reductions$umap.rna@cell.embeddings, 
#          "./data/P21C_all_RNA_UMAP_coordinates.csv", row.names = T, col.names = T, quote=F)


#########################################################
meta.data <- P21_C@meta.data
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names
spatial_data <- read.table(file = "./P21_ATAC-DBiT_matrix.tsv",  header = TRUE, row.names = 1, as.is = TRUE)


data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = spatial_data, assay = assay, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object


library(paletteer)#
col <- paletteer_d("ggthemes::Classic_20")[c(20,6,9,19,3,18,8,7,5,1,13,2,17,15,10,11,12,16)]
col

p1 <- DimPlot(P21_C, label = TRUE,cols = col,group.by = "SCT_snn_res.0.8",pt.size = 1.8) + NoLegend() + ggtitle("RNA UMAP")
p1

p6 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'SCT_snn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
p6$layers[[1]]$aes_params <- c(p6$layers[[1]]$aes_params, shape=22)
p6



