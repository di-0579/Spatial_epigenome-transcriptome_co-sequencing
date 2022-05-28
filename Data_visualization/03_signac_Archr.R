rm(list=ls())

library(Seurat)
library(ArchR)
library(grid)

#set.seed(1)

getwd()
setwd("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/ATAC/ArchR/")


data1 <- read.table("position_ATAC_3.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(data1[1,])
x = x[-1]
data1 = data1[-1]
data1 = as.data.frame(t(data1))



my_data <- read.table(file = '/gpfs/ycga/project/fan/dz286/useful/spatial_barcodes.txt', 
                      sep = '\t', stringsAsFactors=FALSE)

my_data$V4 <- paste0(my_data$V2,"x",my_data$V3)

row.names(my_data) <- my_data$V4

my_data_in_tissue <- my_data[data1$V1, ]

my_data_in_tissue$V1 <- paste0(my_data_in_tissue$V1, "-1")
my_data_in_tissue$V1 <- paste0( "P21_ATAC-DBiT#", my_data_in_tissue$V1)

bc_in_tissue <- my_data_in_tissue$V1 


#########################
#library(parallel)
threads = 8
addArchRThreads(threads = threads)

library("BSgenome.Mmusculus.UCSC.mm10")
addArchRGenome("mm10")

inputFiles <- "fragments.tsv.gz"
inputFiles
sampleNames <- 'P21_ATAC-DBiT'

## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj




proj_in_tissue <- proj[bc_in_tissue, ]

proj_in_tissue


getAvailableMatrices(proj_in_tissue)

head(proj_in_tissue$cellNames)
quantile(proj_in_tissue$TSSEnrichment)

meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))


##################################################################

signacdata <- read.table("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/ATAC/outs/data/Signac_meta_data.txt",sep = '\t', stringsAsFactors=FALSE)
signacdata$seurat_clusters

seurat_clusters <- signacdata[, "seurat_clusters", drop=FALSE]
seurat_clusters$seurat_clusters <- paste0('C', seurat_clusters$seurat_clusters)

row.names(seurat_clusters) <- paste0("P21_ATAC-DBiT#", row.names(seurat_clusters))

seurat_clusters <- seurat_clusters[proj_in_tissue$cellNames, , drop=FALSE]
all(proj_in_tissue$cellNames == row.names(seurat_clusters))
proj_in_tissue$seurat_clusters <- seurat_clusters$seurat_clusters

getCellColData(ArchRProj = proj_in_tissue)

###################################################################


proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 20000, 
    n.start = 10
  ), 
  varFeatures = 30000, 
  dimsToUse = 2:30,
  force = TRUE
)

table(proj_in_tissue$seurat_clusters)


proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p3 <- plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "seurat_clusters", embedding = "UMAP", size = 1)
p3

#ggsave("./plot/seurat_Clusters_ArchR-UMAP.png", plot = p3, width = 9, height = 10)



proj_in_tissue <- addImputeWeights(proj_in_tissue)

## Identify the marker genes for each cluster 
markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "seurat_clusters",
  testMethod = "wilcoxon"
)
###################################
#all gene
# dd <- getMarkers(markersGS, cutOff = "FDR <= 1")
# dd
# markerGenes_dd <- list()
# for (i in seq_len(length(dd))) {
#   markerGenes_dd <- c(markerGenes_dd, dd[[i]]$name)
# }
# 
# markerGenes_dd <- unlist(markerGenes_dd)
# markerGenes_dd <- unique(markerGenes_dd)
# all_gene_score <- getGeneScore_ArchR(ArchRProj = proj_in_tissue, name = markerGenes_dd, imputeWeights = getImputeWeights(proj_in_tissue))
# 
# write.csv(all_gene_score, "data/all_gene_score_matrix.csv", row.names = T, col.names = T, quote=F)
# 

##########################################


markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.1 & Log2FC >=   0.5")
markerList_neg <- getMarkers(markersGS, cutOff = "FDR <= 0.1 & Log2FC >=  -0.5")

df <-as.data.frame(markerList_pos@listData$C0) 




markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}

markerGenes <- unlist(markerGenes)


# for (i in seq_len(length(markerList_pos))) {
#   write.table(markerList_pos[[i]], file=paste0('./data/signac_markers_list/', sampleNames, '_C+', i, '+markers.txt'),
#               quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
# }

#################

markerGenes <- list()
for (i in seq_len(length(markerList_neg))) {
  markerGenes <- c(markerGenes, markerList_neg[[i]]$name)
}

markerGenes <- unlist(markerGenes)


for (i in seq_len(length(markerList_neg))) {
  write.table(markerList_neg[[i]], file=paste0('./markers_list_8/', sampleNames, '_C-', i, '-markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}



#######################################################################################

markerGenes  <- c("Ano1","Emx2","Khdrbs3","Mobp","Pde10a","Prdm12")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >=   0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

reorder_row = c("C1","C2","C4","C3","C8","C5","C6","C7","C0") 
ComplexHeatmap::draw(heatmapGS,row_order=reorder_row,heatmap_legend_side = "bot", annotation_legend_side = "bot")
###########################################################################################


p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue, 
  groupBy = "seurat_clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)


grid::grid.newpage()
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)

ArchRBrowser(proj_in_tissue)


feature_track <- 'Prdm12'
p_track <- plotBrowserTrack(
  ArchRProj = proj_in_tissue,
  plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
  sizes = c(6, 6, 2),
  groupBy = "seurat_clusters", 
  useGroups = 'C4',
  geneSymbol = feature_track,
  upstream = 50000*2,
  downstream = 50000*2,
  scCellsMax = 300,
  scTileSize = 0.2,
  pal = "#D51F26"
)

grid::grid.newpage()
#png(filename = 'Sptb_upReg_track.png', width = 1200, height = 1200, res = 300)
grid::grid.draw(p_track$Prdm12)
#dev.off()
#
plotPDF(plotList = p_track,
        name = paste0(sampleNames, "_liver_brain_Prdm12_Tracks.pdf"),
        ArchRProj = proj_in_tissue,
        addDOC = FALSE, width = 10, height = 5)

feature_track <- 'Prdm12'
p_track <- plotBrowserTrack(
  ArchRProj = proj_in_tissue,
  plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
  sizes = c(6, 6, 2),
  groupBy = "seurat_clusters", 
  useGroups = c('C1', 'C4'),
  geneSymbol = feature_track,
  upstream = 50000*2,
  downstream = 50000*2,
  scCellsMax = 153,
  scTileSize = 0.2,
  pal = c("#208A42", "#D51F26")
)



###########################################################################################
###########################################################################################
## Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "seurat_clusters")

#pathToMacs2 <- findMacs2()
pathToMacs2 <- "/gpfs/ycga/project/fan/dz286/conda_envs/ArchR/bin/macs2"

#library(presto)
proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "seurat_clusters", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

getPeakSet(proj_in_tissue)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getAvailableMatrices(proj_in_tissue)

#devtools::install_github("GreenleafLab/chromVARmotifs")
proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)



proj_in_tissue@peakAnnotation
#######################
#if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
# proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
#}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "PeakMatrix", 
  groupBy = "seurat_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# 
# pma <- markerPlot(seMarker = markersPeaks, name = "C9", cutOff = "FDR <= 0.05 & Log2FC >= 0.1", plotAs = "MA")
# pma
# pv <- markerPlot(seMarker = markersPeaks, name = "C11", cutOff = "FDR <= 0.05 & Log2FC >= 0.1", plotAs = "Volcano")
# pv
# 
# plotPDF(pma, pv, name = "Erythroid-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj_in_tissue, addDOC = FALSE)
# 


enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj_in_tissue,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_in_tissue, addDOC = FALSE)


####################################
## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj_in_tissue, name = "MotifMatrix", plot = TRUE)
motifdata <- plotVarDev$data
write.csv(motifdata, "data/all-motifdata.csv", row.names = T, col.names = T, quote=F)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_in_tissue, addDOC = FALSE)


markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "seurat_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)
markersMotifs 
markersMotifs@metadata

markersMotifs 



motifPositions <- getPositions(proj_in_tissue)
motifPositions



#saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-ATAC0525_clusters", load = FALSE)




markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR <= 0.2 & MeanDiff >= 0.1")

# for (i in seq_len(length(markerMotifsList))) {
#   write.table(markerMotifsList[[i]], file=paste0('./data/motifs_list/', sampleNames, '_C', i, '_motifs.txt'),
#               quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
# }




motifs <- list()
for (i in seq_len(length(markerMotifsList))) {
  motifs <- c(motifs, markerMotifsList[[i]]$name)
}
motifs <- unlist(motifs)
motifs <- paste0('z:', motifs)

# p <- plotEmbedding(
#   ArchRProj = projCUTA, 
#   colorBy = "MotifMatrix", 
#   name = motif_feature, 
#   embedding = "UMAP",
#   imputeWeights = getImputeWeights(projCUTA)
# )
# p

## Spatial plots
library(ggplot2)
library(patchwork)
library(dplyr)

source("/gpfs/ycga/project/fan/dz286/useful/code/getDeviation_ArchR.R")
source("/gpfs/ycga/project/fan/dz286/useful/code/SpatialPlot_new.R")

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

proj_in_tissue <- addImputeWeights(proj_in_tissue)
dev_score <- getDeviation_ArchR(ArchRProj = proj_in_tissue, name = motifs, imputeWeights = getImputeWeights(proj_in_tissue))
dev_score[is.na(dev_score)] <- 0 #min(dev_score, na.rm = TRUE)


# spatial_archr_data <- read.csv(file = "/gpfs/ycga/scratch60/fan/dz286/brainF/integration/data/all_brainF_gene_score_matrix.csv")
# row.names(spatial_archr_data) <- spatial_archr_data$X
# spatial_archr_data <- spatial_archr_data[,-1]
# 

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = dev_score, assay = assay, meta.data = meta.data)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

motifs_cluster <- as.data.frame(markerMotifsList$C6)
motifs <- motifs_cluster$name[1:5]
motifs


motifs <- c("CIC-749",  "NRL-123" , "ARID3C-10",  "EMX2-525",   "FOSB-121" ,
            "ID2-35","KLF16-205","NR4A1-671","POU3F3-620","IRF3-630","SNAI3-268",
            "SOX10-751","SOX12-768","BCL11B-825")
motif_feature <- motifs[14]
motif_feature <- ("Sox10-735")




# p <- SpatialPlot_new(spatial.obj, features = motif_feature, pt.size.factor = 4, image.alpha = 0, stroke = 0) + 
#   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p <- SpatialPlot(spatial.obj, features = motif_feature, pt.size.factor = 1, 
                     image.alpha = 1, stroke = 0, alpha = c(1, 1),  min.cutoff = "q5", max.cutoff = "q95") + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('./motif_ATAC/',motif_feature, '_spatial.png'), width = 1200, height = 1200, res = 300)
p
dev.off()


#saveRDS(spatial.obj, file = "motif_spatial_ATAC.rds")

#spatial.obj<-readRDS("/gpfs/ycga/scratch60/fan/dz286/brainF/ATAC/outs/motif_spatial_ATAC.rds")

df <- as.data.frame(spatial.obj@assays$Spatial@counts)

write.table(df, file = 'brain_E_2-motif_spatial__matrix.tsv', sep = '\t',col.names=TRUE, row.names = TRUE,quote = FALSE)





























motifPositions <- getPositions(proj_in_tissue)
motifPositions



saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-signac_clusters", load = FALSE)



proj_in_tissue <- readRDS("/gpfs/ycga/scratch60/fan/dz286/P21C_20um/ATAC/ArchR/Save-signac_clusters/Save-ArchR-Project.rds")
proj_in_tissue
#######################################################
# proj_in_tissue <- addCoAccessibility(
#   ArchRProj = proj_in_tissue,
#   reducedDims = "IterativeLSI"
# )
# cA <- getCoAccessibility(
#   ArchRProj = proj_in_tissue,
#   corCutOff = 0.5,
#   resolution = 1,
#   returnLoops = FALSE
# )
# cA
# ###############################################


P21C_DBiT_data <- read.table(file = "/gpfs/ycga/scratch60/fan/dz286/P21C_20um/DBiT/output/P21_ATAC-DBiT_matrix.tsv",  header = TRUE, row.names = 1, as.is = TRUE)
P21C_DBiT <- CreateSeuratObject(counts = P21C_DBiT_data, project = 'P21C_DBiT', assay = "RNA")
P21C_DBiT<- FindVariableFeatures(P21C_DBiT, nfeatures = 3000)
P21C_DBiT <- SCTransform(P21C_DBiT)
P21C_DBiT <- RunPCA(P21C_DBiT, assay = "SCT", verbose = FALSE)

P21C_DBiT <- RunUMAP(P21C_DBiT, dims = 1:10, reduction.name = "umap.rna")
P21C_DBiT<- FindNeighbors(P21C_DBiT, dims = 1:10)
P21C_DBiT <- FindClusters(P21C_DBiT, resolution = 0.8)


table(P21C_DBiT$seurat_clusters)






proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = P21C_DBiT,
  addToArrow = TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

proj_in_tissue <- addPeak2GeneLinks(
  ArchRProj = proj_in_tissue,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj_in_tissue,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p2g






p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue, 
  groupBy = "seurat_clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj_in_tissue)
)
grid::grid.newpage()

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)


p <- plotPeak2GeneHeatmap(ArchRProj = proj_in_tissue, groupBy = "seurat_clusters")
p
