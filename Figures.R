library("Seurat")
library("ggplot2")
library("plyr")
library("dplyr")
library("easyGgplot2")
library("ComplexHeatmap")
library("circlize")
library("pcaExplorer")
library("topGO")
library("data.table")
library("tools") 
library("SCENIC")
library("AUCell")
library("SeuratWrappers")
library("monocle3")
library("scales")
library("ggrepel")
library("cowplot")
library("glmGamPoi")
library("patchwork")
library("ggthemes")
library("paletteer")
library("rstatix")
library("ggpubr")
library("gridExtra")
library("SummarizedExperiment")
library("MetaNeighbor")
library("plotly")
library("shiny")
library("fmsb")
library("RColorBrewer")
library("raster")
library("iTALK")
library("viridis")
library("ggeasy")
library("cowplot")
library("patchwork")
library("scater")
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("Matrix")
library("DESeq2")
library("rcartocolor")
library("stringr")
library("AnnotationDbi")
library("factoextra")
library("plot3D")
library("UpSetR")
library("viridisLite")
library("gplots")
library("SCP")

meta.data.list <- list("orig.treatment" = setNames(c("Naive", "LPS", "aGr1.LPS"), 
                                                   c("#4285F4", "#DB4437", "#F4B400")),
                       "orig.ident" = setNames(c("mLu.Naive1", "mLu.Naive2", "mLu.LPS1", "mLu.LPS2", "mLu.aGr1.LPS1", "mLu.aGr1.LPS2"), 
                                               c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#B85A0D", "#FFECB3")), 
                       "Cell.Types" = setNames(c("IMs", "ncIM", "recMac", "ncMono", "DC1", "CD301b-.DC2", "CD301b+.DC2", "infDC2", "migDC", "Cycling"), 
                                               c("#FB8072FF", "#FDB462FF", "#BC80BDFF", "#FCCDE5FF", "#BEBADAFF", "#B3DE69FF", "#CCEBC5FF", "#80B1D3FF", "#FFED6FFF", "#B3B3B3FF")),
                       "IMs.Subtypes" = setNames(c("CD206hi.IMs", "CD206lo.IMs", "ncIM", "recMac", "ncMono", "DC1", "CD301b-.DC2", "CD301b+.DC2", "infDC2", "migDC", "Cycling"), 
                                                 c("#FB8072FF", "#FEE8E5", "#FDB462FF", "#BC80BDFF", "#FCCDE5FF", "#BEBADAFF", "#B3DE69FF", "#CCEBC5FF", "#80B1D3FF", "#FFED6FFF", "#B3B3B3FF")),
                       "IMs.Chemokine.Subtypes" = setNames(paste0("IMck", 0:9), 
                                                           c("#999999", "#FFB6DBFF", "#D8152F", "#B66DFFFF", "#B6DBFFFF", "#006DDBFF", "#009999", "#008600", "#FF9932", "#920000FF")),
                       "rec.Macs.Chemokine.Subtypes" = setNames(paste0("rec.Mac.ck", 0:4), 
                                                                c("#999999", "#E41A1CFF", "#377EB8FF", "#4DAF4AFF", "#FF7F00FF"))) 

Chemokine.Ligands <- unique(AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0008009", columns = c('SYMBOL'), keytype = "GOALL")$SYMBOL)
IMs.Chemokine.List <- list("a" = c("Ccl12", "Ccl7", "Ccl2", "Cxcl14", "Ccl5", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl8", "Ccl6", "Ccl9", "Cxcl10", "Cxcl9", "Cxcl13", "Ccl24"), "b" = c("Pf4", "Cxcl16", "Cklf")) 
IMs.Chemokine.List.Extend <- list("a" = c("Ccl12", "Ccl7", "Ccl2", "Cxcl14", "Ccl5", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl8", "Ccl6", "Ccl9", "Cxcl10", "Cxcl9", "Cxcl13", "Ccl24"), "b" = c("Pf4", "Cxcl16", "Cklf"), "c" = setdiff(Chemokine.Ligands, unlist(IMs.Chemokine.List))) 

Chemokine.Ligands.Hs <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0008009", columns = c('SYMBOL'), keytype = "GOALL")$SYMBOL)
IMs.Chemokine.List.Hs <- list("a" = c("CCL12", "CCL7", "CCL2", "CXCL14", "CCL5", "CCL3", "CCL4", "CXCL1", "CXCL2", "CXCL3", "CCL8", "CCL6", "CCL9", "CXCL10", "CXCL9", "CXCL13", "CCL24"), "b" = c("PF4", "CXCL16", "CKLF"))
IMs.Chemokine.List.Extend.Hs <- list("a" = c("CCL12", "CCL7", "CCL2", "CXCL14", "CCL5", "CCL3", "CCL4", "CXCL1", "CXCL2", "CXCL3", "CCL8", "CCL6", "CCL9", "CXCL10", "CXCL9", "CXCL13", "CCL24"), "b" = c("PF4", "CXCL16", "CKLF"), "c" = setdiff(Chemokine.Ligands, unlist(IMs.Chemokine.List))) 

##Data.Processing##mLu.combined########################################################################################################################################################################################################################

##mLu.combined
mLu.Naive1.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/Naive-A/outs/filtered_feature_bc_matrix")
mLu.Naive2.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/Naive-B/outs/filtered_feature_bc_matrix")
mLu.aGr1.LPS1.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/aGr1-LPS-A/outs/filtered_feature_bc_matrix")
mLu.aGr1.LPS2.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/aGr1-LPS-B/outs/filtered_feature_bc_matrix")
mLu.LPS1.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/LPS-A/outs/filtered_feature_bc_matrix")
mLu.LPS2.data <- Read10X(data.dir = "/Volumes/Jakubzick/Mouse_10x_3-23-20/analysis/LPS-B/outs/filtered_feature_bc_matrix")
mLu.data <- list(mLu.Naive1.data, mLu.Naive2.data, mLu.aGr1.LPS1.data, mLu.aGr1.LPS2.data, mLu.LPS1.data, mLu.LPS2.data)
mLu.list.vector <- c("mLu.Naive1", "mLu.Naive2", "mLu.aGr1.LPS1", "mLu.aGr1.LPS2", "mLu.LPS1", "mLu.LPS2")
percent.mt.max.y <- 20
nFeature_RNA.min <- c(1400, 1200, 1000, 1600, 2000, 1500) # based on distribution
nFeature_RNA.max <- c(3700, 4200, 3500, 4600, 6000, 5500) # based on distribution
percent.mt.max <- c(3, 5.5, 3, 3, 3.5, 4.5) # based on distribution
mLu.list <- list()
for (i in 1:length(mLu.list.vector)){
  mLu.list[[i]] <- CreateSeuratObject(counts = mLu.data[[i]], project = mLu.list.vector[i], min.cells = 3, min.features = 200)
  mLu.list[[i]]$"orig.ident" <- mLu.list.vector[i]
  mLu.list[[i]]$"orig.treatment" <- substr(mLu.list.vector, 5, nchar(mLu.list.vector))[i]
  mLu.list[[i]]$"orig.tissue" <- "mLu"
  mLu.list[[i]]$"orig.species" <- "mouse"
  mLu.list[[i]]$"percent.mt" <- PercentageFeatureSet(mLu.list[[i]], pattern = "^mt-")
  mLu.list[[i]] <- subset(mLu.list[[i]], subset = nFeature_RNA > nFeature_RNA.min[i] & nFeature_RNA < nFeature_RNA.max[i] & percent.mt < percent.mt.max[i]) ##Cell Filter Adjustment##
  mLu.list[[i]] <- RenameCells(mLu.list[[i]], add.cell.id = mLu.list.vector[i])
  mLu.list[[i]] <- mLu.list[[i]] %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters()
}
mLu.combined <- merge(mLu.combined, y = c(mLu.list[[2]], mLu.list[[3]], mLu.list[[4]], mLu.list[[5]], mLu.list[[6]]), merge.data = TRUE)
mLu.combined <- mLu.combined %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
#Cell.Types
#Resolution of 0.2 and 0.4#cluster 10 lowest nCount_RNA and nFeature_RNA, not doublet, decided to remove cause ribosomal RNA enriched in DEGs
resolution <- c(0.2, 0.4)
mLu.combined <- FindClusters(mLu.combined, resolution = resolution)
mLu.combined$Cell.Types <- factor(mapvalues(mLu.combined$SCT_snn_res.0.2, c(0, 1, 3, 4, 6, 10,  2,  8,  9,  5,  7), c(rep("IMs", 6), "recMac", "ncMono", "DC1", "DC2", "Cycling")), levels = c("IMs", "recMac", "ncMono", "DC1", "DC2", "infDC2", "Cycling"))
mLu.combined$Cell.Types[mLu.combined$SCT_snn_res.0.4 == 11 & mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 0] <- "infDC2"
mLu.combined <- subset(mLu.combined, subset = SCT_snn_res.0.2 == "10", invert = TRUE)
mLu.combined <- mLu.combined %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters(resolution = 4)
Idents(mLu.combined) <- "Cell.Types"
mLu.combined <- FindSubCluster(mLu.combined, graph.name = "SCT_snn", cluster = "DC2", resolution = 0.1, subcluster.name = "Cell.Types")
mLu.combined$Cell.Types <- factor(mapvalues(mLu.combined$Cell.Types, from = c("DC2_0", "DC2_1"), to = c("DC2", "migDC")), levels = c("IMs", "recMac", "ncMono", "DC1", "DC2", "infDC2", "migDC", "Cycling"))
Idents(mLu.combined) <- "Cell.Types"
mLu.combined <- FindSubCluster(mLu.combined, graph.name = "SCT_snn", cluster = "IMs", resolution = 1, subcluster.name = "Cell.Types")
mLu.combined$Cell.Types <- factor(mapvalues(mLu.combined$Cell.Types, from = c(paste0("IMs_", 0:21), "DC2"), c(rep("IMs", 17), "ncIM", "CD301b+.DC2", rep("IMs", 3), "CD301b−.DC2")), levels = c("IMs", "ncIM", "recMac", "ncMono", "DC1", "CD301b−.DC2", "CD301b+.DC2", "infDC2", "migDC", "Cycling"))
mLu.combined$orig.ident <- factor(mLu.combined$orig.ident, levels = paste0("mLu.", rep(c("Naive", "LPS", "aGr1.LPS"), each = 2), 1:2))
mLu.combined$orig.treatment <- factor(gsub("mLu.", ", ", mLu.combined$orig.treatment), levels = c("Naive", "LPS", "aGr1.LPS"))
#IMs.Chemokine.Subtypes
mLu.combined <- FindClusters(mLu.combined, resolution = 4)
mLu.combined.IMs <- subset(mLu.combined, subset = Cell.Types == "IMs")
mLu.combined.IMs <- mLu.combined.IMs %>% FindNeighbors() %>% RunUMAP()
mLu.combined.IMs$IMs.Chemokine.Subtypes <- factor(mapvalues(mLu.combined$SCT_snn_res.4, c(19, 11, 30,   49, 0, 42,   1,   53,   31, 34,   40, 50,  15, 18, 21,   20,   14, 32, 44, 45, 48,   37, 39,  2, 3, 4, 5, 6, 8, 17, 22, 23, 25, 26, 27, 28, 29, 33, 35, 36, 38, 43), c(rep("IMck1", 3), rep("IMck2", 3), "IMck3", "IMck4", rep("IMck5", 2), rep("IMck6", 5), "IMck7", rep("IMck8", 5), rep("IMck9", 2), rep("IMck0", 19))), levels = paste0("IMck", 0:9))
mLu.combined$IMs.Chemokine.Subtypes <- factor("Others", levels = c("Others", paste0("IMck", 0:9)))
for (i in 1:length(paste0("IMck", 0:9))){
  mLu.combined$IMs.Chemokine.Subtypes[colnames(mLu.combined.IMs)[mLu.combined.IMs$IMs.Chemokine.Subtypes == paste0("IMck", 0:9)[i]]] <- paste0("IMck", 0:9)[i]
}
#IMs.Subtypes
CD206.p <- colnames(mLu.combined.IMs)[mLu.combined.IMs[["SCT"]]@data["Mrc1", ] > 2]
mLu.combined$IMs.Subtypes <- factor(mapvalues(mLu.combined$Cell.Types, from = "IMs", to = "CD206lo.IMs"), levels = c("CD206hi.IMs", "CD206lo.IMs", "ncIM", "recMac", "ncMono", "DC1", "CD301b-.DC2", "CD301b+.DC2", "infDC2", "migDC", "Cycling"))
mLu.combined$IMs.Subtypes[CD206.p] <- "CD206hi.IMs"

##Data.Processing##mLu.combined.IMs##RNA.velocity##############################################################################################################################################################################

read.loom.matrices <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}

# https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/Velocyto_fixed
# https://github.com/zhanghao-njmu/SCP
sample_id <- paste0(rep(c("aGr1-LPS-", "LPS-"), each = 2), rep(c("A", "B"), 2))
seurat_sample_id <- paste0("mLu.", rep(c("aGr1.LPS", "LPS"), each = 2), rep(1:2, 2))
reformatLoomColnames <- function(loomObject, assay){
  paste0(stri_replace_all_regex(colnames(loomObject[[assay]]), pattern = c("x", paste0(sample_id, ":")), replacement = c("", paste0(seurat_sample_id, "_")), vectorize = FALSE), "-1")
}
loom_data <- lapply(sample_id, function(sample_id){
  loom_object <- read.loom.matrices(paste0("/Volumes/mLu/Analysis/scRNA_seq/Datasets/velocity/", sample_id, ".loom"))
  colnames(loom_object$spliced) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "spliced")
  colnames(loom_object$unspliced) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "unspliced")
  colnames(loom_object$ambiguous) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "ambiguous")
  loom_object
})
names(loom_data) <- sample_id

loom_aggregate <- lapply(c("spliced", "unspliced"), function(assay){
  do.call("cbind", lapply(loom_data,"[[", assay))
})
names(loom_aggregate) <- c("spliced", "unspliced")

mLu.combined.velocity <- CreateSeuratObject(counts = loom_aggregate$spliced, project = "Velocity", assay = "spliced", min.cells = 0, min.features = 0, names.field = 2, names.delim = "\\-")
mLu.combined.velocity[["unspliced"]] <- CreateAssayObject(counts = loom_aggregate$unspliced)
mLu.combined.velocity <- mLu.combined.velocity[rownames(mLu.combined), colnames(mLu.combined)]
mLu.combined.velocity[["RNA"]] <- CreateAssayObject(counts = mLu.combined[["RNA"]]@counts[rownames(mLu.combined.velocity), intersect(colnames(mLu.combined.velocity), colnames(mLu.combined))])
mLu.combined.velocity[["SCT"]] <- CreateAssayObject(counts = mLu.combined[["SCT"]]@counts[rownames(mLu.combined.velocity), intersect(colnames(mLu.combined.velocity), colnames(mLu.combined))])
mLu.combined.velocity$IMs.Chemokine.Subtypes <- as.vector(mLu.combined$IMs.Chemokine.Subtypes[colnames(mLu.combined.velocity)])
mLu.combined.velocity <- subset(mLu.combined.velocity, subset = IMs.Chemokine.Subtypes != "Others")
DefaultAssay(mLu.combined.velocity) <- "SCT"
mLu.combined.velocity <- mLu.combined.velocity %>% RunPCA() %>% FindNeighbors() %>% RunUMAP()

##Data.Processing##mHr.combined################################################################################################################################################################################

##mHr.combined
mHr.1.data <- Read10X(data.dir = "/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/mHt/mHr.1")
mHr.2.data <- Read10X(data.dir = "/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/mHt/mHr.2")
mHr.data <- list(mHr.1.data, mHr.2.data)
mHr.list.vector <- c("mHr.1", "mHr.2")
nFeature_RNA.min <- c(900, 900)
nFeature_RNA.max <- c(4500, 4500)
percent.mt.max <- c(8, 8)
mHr.list <- list()
for (i in 1:length(mHr.list.vector)){
  mHr.list[[i]] <- CreateSeuratObject(counts = mHr.data[[i]][[1]], project = mHr.list.vector[i], min.cells = 3, min.features = 200)
  mHr.list[[i]]$"orig.ident" <- mHr.list.vector[i]
  mHr.list[[i]]$"orig.treatment" <- "Naive"
  mHr.list[[i]]$"orig.tissue" <- "mHr"
  mHr.list[[i]]$"orig.species" <- "mouse"
  mHr.list[[i]]$"percent.mt" <- PercentageFeatureSet(mHr.list[[i]], pattern = "^mt-")
  mHr.list[[i]] <- subset(mHr.list[[i]], subset = nFeature_RNA > nFeature_RNA.min[i] & nFeature_RNA < nFeature_RNA.max[i] & percent.mt < percent.mt.max[i])
  mHr.list[[i]] <- RenameCells(mHr.list[[i]], add.cell.id = mHr.list.vector[i])
  mHr.list[[i]] <- mHr.list[[i]] %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters()
}
mHr.features <- SelectIntegrationFeatures(object.list = mHr.list, nfeatures = 10000)
mHr.features.list <- PrepSCTIntegration(object.list = mHr.list, anchor.features = mHr.features)
mHr.anchors <- FindIntegrationAnchors(object.list = mHr.features.list, normalization.method = "SCT", anchor.features = mHr.features)
mHr.combined <- IntegrateData(anchorset = mHr.anchors, normalization.method = "SCT", new.assay.name = "ITG")
mHr.combined <- mHr.combined %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
##Cell.Types
mHr.combined <- mHr.combined %>% FindClusters(graph.name = "ITG_snn", resolution = 0.2)
#0, 1, 5, 10 are IMs, based on c("Cd163", "Cx3cr1", "H2-Aa", "Timd4", "Ccr2", "Folr2", "Mrc1", "Itgax", "C1qc", "Vcan", "Dpp4", "Pf4", "Cxcl16", "C1qa")
mHr.combined$Cell.Types <- "Others"
mHr.combined$Cell.Types[mHr.combined$ITG_snn_res.0.2 %in% c(0, 1, 5, 10)] <- "IMs"
#mHr.combined.IMs
mHr.combined.IMs <- subset(mHr.combined, subset = Cell.Types == "IMs")
DefaultAssay(mHr.combined.IMs) <- "SCT"
mHr.combined.IMs <- mHr.combined.IMs %>% FindVariableFeatures() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
mHr.combined.IMs$Chemokine.Cell.Types <- factor(mapvalues(mHr.combined.IMs$SCT_snn_res.0.8, from = c(14, 9, 16, 20, 11, 17, 15, 7, 21), to = paste0("IMck", 1:9)), levels = c(paste0("IMck", 0:9), 0:22))
mHr.combined.IMs$Chemokine.Cell.Types <- factor(mapvalues(mHr.combined.IMs$Chemokine.Cell.Types, from = 0:22, to = paste0("IMck", rep(0, 23))), levels = paste0("IMck", 0:9))

##Figure.1##################################################################################################################################################################################

#Figure.1A.DimPlot
#Figure.1B.DimPlot
mLu.combined.orig.treatment.list <- append(append(list(setNames(meta.data.list$orig.treatment, rep("#E5E5E5", 3))), list(meta.data.list$orig.treatment)), split(meta.data.list$orig.treatment, names(meta.data.list$orig.treatment)))
mLu.combined.UAMP_1.length <- max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]) - min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])
mLu.combined.UAMP_2.length <- max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]) - min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])
for ( i in 1:length(mLu.combined.orig.treatment.list)){
  Figure.1.DimPlot <- DimPlot(mLu.combined, group.by = "orig.treatment", cells = colnames(mLu.combined)[mLu.combined$orig.treatment %in% mLu.combined.orig.treatment.list[[i]]], cols = names(mLu.combined.orig.treatment.list[[i]]), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.1.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.1", c("A", paste0("B.", 1:4))[i], ".DimPlot.pdf")
  pdf(Figure.1.DimPlot.directory, width = 18, height = 15)
  print(Figure.1.DimPlot)
  dev.off() 
}
#Figure.1C.FeaturePlot
Figure.1.FeaturePlot.Features <- c("Dpp4", "C5ar1", "Ly6c2")
for (i in 1:length(Figure.1.FeaturePlot.Features)){
  Figure.1.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.1.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.1.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.1C.", Figure.1.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.1.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.1.FeaturePlot)
  dev.off() 
}
#Figure.1D.DimPlot
Figure.1.DimPlot <- DimPlot(mLu.combined, group.by = "Cell.Types", cols = names(meta.data.list$Cell.Types), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
Figure.1.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.1D", ".DimPlot.pdf")
pdf(Figure.1.DimPlot.directory, width = 18, height = 15)
print(Figure.1.DimPlot)
dev.off() 
#Figure.1E.DotPlot
#Figure.1.Curated.Genes.DotPlot
Figure.1.Curated.Genes.DotPlot.Features <- list("IMs" = c("Lyve1", "Folr2", "C1qa", "Adgre1"), "ncIM" = c("Ctsz", "Ctsd", "Esd", "Mmp12"), "recMacs" = c("Cd14", "Fn1", "Ccr2", "Ly6c2"), "ncMono" = c("Plac8", "Cx3cr1", "Spn", "Itgal"), "DC1" = c("Xcr1", "Clec9a", "Itgae", "Irf8"), "CD301b-.DC2" = c("S100a4", "S100a6", "H2-DMb2", "Irf4"), "CD301b+.DC2" = c("Mgl2", "Cd209a", "H2-Oa", "Siglecg"), "infDC2" = c("Ifit1", "Phf11d", "Ifit3", "Ifi205"), "migDC" = c("Ccl17", "Ccl22", "Ccr7", "Ccl5"), "Cycling" = c("Mki67", "Stmn1", "Top2a", "Cdk1"))
Figure.1.Curated.Genes.DotPlot <- DotPlot(mLu.combined, group.by = "Cell.Types", assay = "SCT", features = Figure.1.Curated.Genes.DotPlot.Features, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
Figure.1.Curated.Genes.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.1E.Curated.Genes.DotPlot.pdf")
pdf(Figure.1.Curated.Genes.DotPlot.directory, width = 14, height = 3.7)
print(Figure.1.Curated.Genes.DotPlot)
dev.off()
#Figure.1F.DotPlot
#Figure.1.DEGs.DotPlot
mLu.combined <- PrepSCTFindMarkers(mLu.combined, assay = "SCT", verbose = TRUE)
mLu.combined.Cell.Types.DEGs <- FindAllMarkers(mLu.combined, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
mLu.combined.Cell.Types.DEGs.top4 <- mLu.combined.Cell.Types.DEGs %>% group_by(cluster) %>% slice_head(n = 4)
Figure.1.DEGs.DotPlot.Features <- setNames(mLu.combined.Cell.Types.DEGs.top4$gene, mLu.combined.Cell.Types.DEGs.top4$cluster)
Figure.1.DEGs.DotPlot.Features <- Figure.1.DEGs.DotPlot.Features[!duplicated(Figure.1.DEGs.DotPlot.Features)]
Figure.1.DEGs.DotPlot.Features <- lapply(split(Figure.1.DEGs.DotPlot.Features, names(Figure.1.DEGs.DotPlot.Features)), unname)
Figure.1.DEGs.DotPlot.Features <- Figure.1.DEGs.DotPlot.Features[meta.data.list$Cell.Types]
# Figure.1.DEGs.DotPlot.Features <- list("IMs" = c("Pf4", "Apoe", "Ccl12", "C1qc"), "ncIM" = c("Slc7a11", "Lgals3", "AA467197", "Cd9"), "recMacs" = c("Fn1", "Il1b", "Plac8", "Ccr2"), "ncMono" = c("Ifitm6", "Itgb7", "Adgre4", "Atp1a3"), "DC1" = c("Cd24a", "Itgae", "Xcr1", "Rnase6"), "CD301b-.DC2" = c("Ccl17", "S100a6", "S100a4", "Napsa"), "CD301b+.DC2" = c("Cd209a", "Ear2", "Clec4b1", "Bhlhe40"), "infDC2" = c("Ifit3", "Rsad2", "Ifit1", "Ifit2"), "migDC" = c("Ccl22", "Ccr7", "Il4i1", "Traf1"), "Cycling" = c("Stmn1", "2810417H13Rik", "Birc5", "Top2a"))
Figure.1.DEGs.DotPlot <- DotPlot(mLu.combined, group.by = "Cell.Types", assay = "SCT", features = Figure.1.DEGs.DotPlot.Features, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
Figure.1.DEGs.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.1F.DEGs.DotPlot.New.pdf")
pdf(Figure.1.DEGs.DotPlot.directory, width = 14, height = 4)
print(Figure.1.DEGs.DotPlot)
dev.off() 

##Figure.2##################################################################################################################################################################################

#Figure.2A.FeaturePlot
#Figure.2B.FeaturePlot
Figure.2.FeaturePlot.Features <- c("Cd163", "Cx3cr1", "H2-Aa", "Timd4", "Ccr2", "Folr2", "Mrc1", "Itgax", "C1qc")
for (i in 1:length(Figure.2.FeaturePlot.Features)){
  Figure.2.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.2.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.2.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.2", ifelse(Figure.2.FeaturePlot.Features[i] == "C1qc", "B", "A"), ".", Figure.2.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.2.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.2.FeaturePlot)
  dev.off() 
}
#Figure.2.Features
mLu.combined.IMs <- subset(mLu.combined, subset = IMs.Subtypes %in% c("CD206hi.IMs", "CD206lo.IMs"))
mLu.combined.IMs$Orig.Treatment <- factor(mapvalues(mLu.combined.IMs$orig.treatment, "aGr1.LPS", "LPS"), levels = c("Naive", "LPS"))
mLu.combined.IMs <- PrepSCTFindMarkers(mLu.combined.IMs)
Idents(mLu.combined.IMs) <- mLu.combined.IMs$IMs.Subtypes
mLu.combined.IMs.Cell.Types.DEGs <- FindAllMarkers(mLu.combined.IMs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mLu.combined.IMs.Cell.Types.DEGs.top20 <- mLu.combined.IMs.Cell.Types.DEGs %>% group_by(cluster) %>% slice_head(n = 20)
Figure.2.DEGs.DotPlot.Features <- setNames(mLu.combined.IMs.Cell.Types.DEGs.top20$gene, mLu.combined.IMs.Cell.Types.DEGs.top20$cluster)
Figure.2.DEGs.DotPlot.Features <- Figure.2.DEGs.DotPlot.Features[!duplicated(Figure.2.DEGs.DotPlot.Features)]
Figure.2.DEGs.DotPlot.Features <- lapply(split(Figure.2.DEGs.DotPlot.Features, names(Figure.2.DEGs.DotPlot.Features)), unname)
Figure.2.DEGs.DotPlot.Features <- Figure.2.DEGs.DotPlot.Features[c("CD206hi.IMs", "CD206lo.IMs")]
Figure.2.DEGs.DotPlot.Features <- list("CD206hi.IMs" = c("Sepp1", "F13a1", "Apoe", "Cbr2", "Mrc1", "Folr2", "Lyve1", "Cd163", "Gas6", "Ltc4s", "Ly6e", "Ccl6", "Ccl8", "Cfh", "Trf", "Fcrls", "Fcgrt", "Klf2", "Clec10a", "Stab1"), "CD206lo.IMs" = c("Ccl5", "H2-Eb1", "Il1b", "H2-Ab1", "Lgals3", "Cd74", "Cd9", "Mmp12", "Lpl", "Cd72", "Cd52", "Coro1a", "Bcl2a1a", "Zmynd15", "Sod2", "Clec4n", "Bcl2a1b", "Lsp1", "Acp5", "Capg"))
#Figure.2D.DEGs.Heatmap
Figure.2.DEGs.Heatmap.data <- mLu.combined.IMs[["SCT"]]@scale.data[unlist(Figure.2.DEGs.DotPlot.Features), ]
Figure.2.DEGs.Heatmap.data.cell.names <- c(sample(colnames(mLu.combined.IMs)[mLu.combined.IMs$IMs.Subtypes == "CD206hi.IMs"], 2000), sample(colnames(mLu.combined.IMs)[mLu.combined.IMs$IMs.Subtypes == "CD206lo.IMs"], 2000))
Figure.2.DEGs.Heatmap.data <- Figure.2.DEGs.Heatmap.data[, Figure.2.DEGs.Heatmap.data.cell.names]
colors <- c("#F96650", "#FEE8E5")
column_split <- factor(c(rep("CD206hi.IMs", 2000), rep("CD206lo.IMs", 2000)), levels = c(meta.data.list$IMs.Subtypes[1:2]))
row_split <- factor(c(rep("CD206hi.IMs", 20), rep("CD206lo.IMs", 20)), levels = c(meta.data.list$IMs.Subtypes[1:2]))
top_annotation <- columnAnnotation(IMs.Cell.Types = c(rep("CD206hi.IMs", 2000), rep("CD206lo.IMs", 2000)), col = list(IMs.Cell.Types = setNames(colors, meta.data.list$IMs.Subtypes[1:2])), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation <- rowAnnotation(IMs.Cell.Types = c(rep("CD206hi.IMs", 20), rep("CD206lo.IMs", 20)), col = list(IMs.Cell.Types = setNames(colors, meta.data.list$IMs.Subtypes[1:2])), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.2.DEGs.Heatmap <- Heatmap(Figure.2.DEGs.Heatmap.data, use_raster = TRUE, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = FALSE, column_split = column_split, column_gap = unit(2, "mm"), column_title = NULL, row_split = row_split, row_gap = unit(2, "mm"), row_title = NULL, cluster_column_slices = FALSE, 
                                 row_names_gp = gpar(fontsize = 12, fontface = "italic"), show_heatmap_legend = FALSE, width = unit(3.75, "in"), height = unit(8, "in"), top_annotation = top_annotation, right_annotation = right_annotation)
Figure.2.DEGs.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.2D.IMs.top20.DEGs.Heatmap.pdf")
pdf(Figure.2.DEGs.Heatmap.directory, width = 5.25, height = 9)
print(Figure.2.DEGs.Heatmap)
dev.off()

#Figure.2E.GO.Heatmap
mLu.combined.IMs.bg <- rownames(mLu.combined.IMs)[rowSums(mLu.combined.IMs) > 0]
mLu.combined.IMs.topGO.list <- vector()
for (i in c(meta.data.list$IMs.Subtypes[1:2])){
  mLu.combined.IMs.DEG <- mLu.combined.IMs.Cell.Types.DEGs[mLu.combined.IMs.Cell.Types.DEGs$cluster == i, "gene"]
  mLu.combined.IMs.topGO.xls.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Datasets/topGO/Cell.Types/mLu.combined.IMs_", i, ".topGO.xls")
  mLu.combined.IMs.topGO <- topGOtable(mLu.combined.IMs.DEG, mLu.combined.IMs.bg, mapping = "org.Mm.eg.db", ontology = "BP", writeOutput = TRUE, topTablerows = 500, outputFile = mLu.combined.IMs.topGO.xls.directory)
  mLu.combined.IMs.topGO.list <- c(mLu.combined.IMs.topGO.list, mLu.combined.IMs.topGO$GO.ID[1:20])
}
names(mLu.combined.IMs.topGO.list) <- rep(meta.data.list$IMs.Subtypes[1:2], each = 20)
mLu.combined.IMs.topGO.list <- mLu.combined.IMs.topGO.list[!duplicated(mLu.combined.IMs.topGO.list)]
mLu.combined.IMs.topGO.targetGO.Heatmap.table <- data.frame(matrix(ncol = 0, nrow = length(mLu.combined.IMs.topGO.list)))
for (i in c(meta.data.list$IMs.Subtypes[1:2])){
  mLu.combined.IMs.DEG <- mLu.combined.IMs.Cell.Types.DEGs[mLu.combined.IMs.Cell.Types.DEGs$cluster == i, "gene"]
  mLu.combined.IMs.topGO.targetGO.xls.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Datasets/topGO/Cell.Types/mLu.combined.IMs_", i, ".targetGO.xls")
  mLu.combined.IMs.topGO.targetGO <- topGOtable.targetGO(mLu.combined.IMs.DEG, mLu.combined.IMs.bg, mapping = "org.Mm.eg.db", ontology = "BP", TargetGO = mLu.combined.IMs.topGO.list, writeOutput = TRUE, outputFile = mLu.combined.IMs.topGO.targetGO.xls.directory)
  setnafill(mLu.combined.IMs.topGO.targetGO, "nocb", cols = "p.value_elim")
  mLu.combined.IMs.topGO.targetGO.List <- mLu.combined.IMs.topGO.targetGO[match(mLu.combined.IMs.topGO.list, mLu.combined.IMs.topGO.targetGO$GO.ID), ]
  mLu.combined.IMs.topGO.targetGO.Heatmap.table <- cbind(mLu.combined.IMs.topGO.targetGO.Heatmap.table, mLu.combined.IMs.topGO.targetGO.List$p.value_elim)
}
rownames(mLu.combined.IMs.topGO.targetGO.Heatmap.table) <- toTitleCase(mLu.combined.IMs.topGO.targetGO.List$Term)
colnames(mLu.combined.IMs.topGO.targetGO.Heatmap.table) <- meta.data.list$IMs.Subtypes[1:2]
mLu.combined.IMs.topGO.targetGO.Heatmap.table <- data.matrix(-log10(mLu.combined.IMs.topGO.targetGO.Heatmap.table), rownames.force = TRUE)
colors <- c("#F96650", "#FEE8E5")
rownames(mLu.combined.IMs.topGO.targetGO.Heatmap.table) <- str_to_sentence(rownames(mLu.combined.IMs.topGO.targetGO.Heatmap.table))
fontface <- rep("plain", nrow(mLu.combined.IMs.topGO.targetGO.Heatmap.table))
fontface[grep("chemotaxis", rownames(mLu.combined.IMs.topGO.targetGO.Heatmap.table))] <- "bold"
col_fun = colorRamp2(range(mLu.combined.IMs.topGO.targetGO.Heatmap.table), c("white", "red"))
column_split <- factor(c(rep("CD206hi.IMs", 1), rep("CD206lo.IMs", 1)), levels = c(meta.data.list$IMs.Subtypes[1:2]))
row_split <- factor(c(rep("CD206hi.IMs", 20), rep("CD206lo.IMs", 18)), levels = c(meta.data.list$IMs.Subtypes[1:2]))
top_annotation <- columnAnnotation(IMs.Cell.Types = c("CD206hi.IMs", "CD206lo.IMs"), col = list(IMs.Cell.Types = setNames(colors, meta.data.list$IMs.Subtypes[1:2])), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation <- rowAnnotation(IMs.Cell.Types = c(rep("CD206hi.IMs", 20), rep("CD206lo.IMs", 18)), col = list(IMs.Cell.Types = setNames(colors, meta.data.list$IMs.Subtypes[1:2])), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.2E.GO.Heatmap <- Heatmap(mLu.combined.IMs.topGO.targetGO.Heatmap.table, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, column_title = NULL, show_heatmap_legend = FALSE, heatmap_width = unit(4.5, "in"), heatmap_height = unit(8, "in"), col = col_fun, column_split = column_split, column_gap = unit(2, "mm"), right_annotation = right_annotation, top_annotation = top_annotation, 
                                row_names_gp = gpar(fontsize = 12, fontface = fontface), column_names_rot = 0, column_names_centered = TRUE, column_names_gp = gpar(fontsize = 14))
Figure.2E.GO.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.2E.IMs.top20.topGO.Heatmap.pdf")
pdf(Figure.2E.GO.Heatmap.directory, width = 10, height = 9)
print(Figure.2E.GO.Heatmap)
dev.off()

##Figure.3##Figure.S4################################################################################################################################################################################

#Figure.3A.DimPlot
Figure.3.DimPlot <- DimPlot(mLu.combined, group.by = "IMs.Chemokine.Subtypes", cols = c("#E5E5E5", names(meta.data.list$IMs.Chemokine.Subtypes)), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
Figure.3.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3A", ".DimPlot.pdf")
pdf(Figure.3.DimPlot.directory, width = 18, height = 15)
print(Figure.3.DimPlot)
dev.off()
#Figure.3A.Re-cluster.DimPlot
Figure.1.DimPlot <- DimPlot(mLu.combined.IMs, group.by = "IMs.Chemokine.Subtypes", cols = names(meta.data.list$IMs.Chemokine.Subtypes), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])*1.01) + ylim(range(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])*1.01) + coord_cartesian(expand = FALSE)
Figure.1.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3A.Re-cluster", ".DimPlot.pdf")
pdf(Figure.1.DimPlot.directory, width = 18, height = 15)
print(Figure.1.DimPlot)
dev.off() 
#Figure.3B.DotPlot
Chemokine.DotPlot <- DotPlot(mLu.combined, idents = paste0("IMck", 0:10), group.by = "IMs.Chemokine.Subtypes", assay = "SCT", features = IMs.Chemokine.List, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.3B.mLu.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 5.6, height = 3.3)
print(Chemokine.DotPlot)
dev.off() 
#Figure.3C.FeaturePlot
#Figure.S4E.FeaturePlot
FeaturePlot.Features.Vector <- c("Ccl12", "Ccl7", "Ccl2", "Cxcl14", "Ccl5", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl8", "Ccl6", "Ccl9", "Cxcl10", "Cxcl9", "Cxcl13", "Ccl24", "Pf4","Cxcl16","Cklf")
for (i in 1:length(FeaturePlot.Features.Vector)){
  Figure.1.FeaturePlot <- FeaturePlot(mLu.combined.IMs, features = FeaturePlot.Features.Vector[i], raster = FALSE, order = TRUE, pt.size = 4) + scale_colour_gradientn(colours = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000")) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])*1.01) + ylim(range(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])*1.01) + coord_cartesian(expand = FALSE)
  Figure.1.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.", ifelse(FeaturePlot.Features.Vector[i] %in% c("Pf4", "Cxcl16", "Cklf"), "S4E", "3C"), ".", FeaturePlot.Features.Vector[i], ".FeaturePlot.pdf")
  pdf(Figure.1.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.1.FeaturePlot)
  dev.off() 
}
#Figure.3D.All.PCAPlot
#Figure.3D.IMs.PCAPlot
mLu.combined$IMs.Chemokine.Subtypes.Cell.Types <- factor(mLu.combined$Cell.Types, levels = c(paste0("IMck", 0:9), levels(mLu.combined$Cell.Types)[-(1:2)]))
for (i in paste0("IMck", 0:9)){
  mLu.combined$IMs.Chemokine.Subtypes.Cell.Types[mLu.combined$IMs.Chemokine.Subtypes == i] <- i
}
mLu.combined.AverageExpression <- AverageExpression(mLu.combined, assays = "SCT", slot = "data", group.by = "IMs.Chemokine.Subtypes.Cell.Types")
PCA.Result <- prcomp(t(mLu.combined.AverageExpression$SCT), scale = FALSE)
PCAPlot.Data <- as.data.frame(PCA.Result$x)

PCA.Percentage <- c(PCA.Result$sdev[[1]]^2/sum(PCA.Result$sdev^2), PCA.Result$sdev[[2]]^2/sum(PCA.Result$sdev^2), PCA.Result$sdev[[3]]^2/sum(PCA.Result$sdev^2))
PCA.Percentage <- percent(PCA.Percentage, accuracy = 0.01)
#Figure.3D.All.PCAPlot
Figure.3.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3D", ".All", ".PCAPlot.pdf")
pdf(Figure.3.FeaturePlot.directory, width = 7.5, height = 7.5)
print(text3D(PCAPlot.Data$PC1, PCAPlot.Data$PC2, PCAPlot.Data$PC3, xlim = range(PCAPlot.Data$PC1) + c(-10, 10), ylim = range(PCAPlot.Data$PC2) + c(-10, 10), zlim = range(PCAPlot.Data$PC3) + c(-10, 10),
             labels = rownames(PCAPlot.Data), col = rep("black", 18), theta = 30, phi = 15,
             xlab = paste0("PC1 (", PCA.Percentage[1], ")"), ylab = paste0("PC2 (", PCA.Percentage[2], ")"), zlab = paste0("PC3 (", PCA.Percentage[3], ")"), 
             main = "", cex = 1.25, cex.lab = 1, bty = "b2", ticktype = "simple", d = 2, adj = 0.5, font = 1))
print(scatter3D(PCAPlot.Data$PC1[1:10], PCAPlot.Data$PC2[1:10], PCAPlot.Data$PC3[1:10] - 18, colvar = 1:10, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)]))[1:10], colkey = FALSE, 
                type = "p", cex = 1.5, pch = 19, add = TRUE))
print(scatter3D(PCAPlot.Data$PC1[11:18], PCAPlot.Data$PC2[11:18], PCAPlot.Data$PC3[11:18] - 18, colvar = 1:8, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)]))[11:18], colkey = FALSE, 
                type = "p", cex = 1.5, pch = 18, add = TRUE))
print(scatter3D(PCAPlot.Data$PC1, PCAPlot.Data$PC2, PCAPlot.Data$PC3 - 18, colvar = 1:18, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)])), colkey = FALSE, 
                type = "h", pch = 1, add = T))
dev.off() 
#Figure.3D.IMs.PCAPlot
PCAPlot.Data.IMs <- PCAPlot.Data[1:10, ]
Figure.3.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3D", ".IMs", ".PCAPlot.pdf")
pdf(Figure.3.FeaturePlot.directory, width = 7.5, height = 7.5)
print(text3D(PCAPlot.Data$PC1, PCAPlot.Data$PC2, PCAPlot.Data$PC3, xlim = range(PCAPlot.Data$PC1) + c(-10, 10), ylim = range(PCAPlot.Data$PC2) + c(-10, 10), zlim = range(PCAPlot.Data$PC3) + c(-10, 10),
             labels = rownames(PCAPlot.Data), col = rep("black", 18), theta = 30, phi = 15,
             xlab = paste0("PC1 (", PCA.Percentage[1], ")"), ylab = paste0("PC2 (", PCA.Percentage[2], ")"), zlab = paste0("PC3 (", PCA.Percentage[3], ")"), 
             main = "", cex = 1.25, cex.lab = 1, bty = "b2", ticktype = "simple", d = 2, adj = 0.5, font = 1))
print(scatter3D(PCAPlot.Data$PC1[1:10], PCAPlot.Data$PC2[1:10], PCAPlot.Data$PC3[1:10] - 18, colvar = 1:10, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)]))[1:10], colkey = FALSE, 
                type = "p", cex = 1.5, pch = 19, add = TRUE))
print(scatter3D(PCAPlot.Data$PC1[11:18], PCAPlot.Data$PC2[11:18], PCAPlot.Data$PC3[11:18] - 18, colvar = 1:8, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)]))[11:18], colkey = FALSE, 
                type = "p", cex = 1.5, pch = 18, add = TRUE))
print(scatter3D(PCAPlot.Data$PC1, PCAPlot.Data$PC2, PCAPlot.Data$PC3 - 18, colvar = 1:18, col = names(c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)])), colkey = FALSE, 
                type = "h", pch = 1, add = T))
dev.off() 
#Figure.3E.Heatmap
Gene.Matrix_by.Cell.Types <- mLu.combined[["SCT"]]@data
Gene.Matrix_by.Cell.Types <- Gene.Matrix_by.Cell.Types[, !is.na(mLu.combined$IMs.Chemokine.Subtypes.Cell.Types)]
colData_by.Cell.Types <- mLu.combined@meta.data[, c("IMs.Chemokine.Subtypes.Cell.Types", "orig.ident")]
colData_by.Cell.Types <- tibble::rownames_to_column(colData_by.Cell.Types, "Cell_ID")
colData_by.Cell.Types <- colData_by.Cell.Types[!is.na(mLu.combined$IMs.Chemokine.Subtypes.Cell.Types), ]
SE.Combined_by.Cell.Types <- SummarizedExperiment(assays = SimpleList(gene_matrix = Gene.Matrix_by.Cell.Types), colData = colData_by.Cell.Types)
Var.Genes_by.Cell.Types <- variableGenes(dat = SE.Combined_by.Cell.Types, exp_labels = SE.Combined_by.Cell.Types$IMs.Chemokine.Subtypes.Cell.Types)
Cell.Types.NV_by.Cell.Types <- MetaNeighborUS(var_genes = Var.Genes_by.Cell.Types, dat = SE.Combined_by.Cell.Types, study_id = rep("", sum(!is.na(mLu.combined$IMs.Chemokine.Subtypes.Cell.Types))), cell_type = SE.Combined_by.Cell.Types$IMs.Chemokine.Subtypes.Cell.Types)
rownames(Cell.Types.NV_by.Cell.Types) <- gsub("\\|", "", rownames(Cell.Types.NV_by.Cell.Types))
colnames(Cell.Types.NV_by.Cell.Types) <- gsub("\\|", "", colnames(Cell.Types.NV_by.Cell.Types))
Figure.3.Heatmap.Cell.Types <- c(meta.data.list$IMs.Chemokine.Subtypes, meta.data.list$Cell.Types[-(1:2)])
column_names_gp.col  <- names(Figure.3.Heatmap.Cell.Types)[match(rownames(Cell.Types.NV_by.Cell.Types), Figure.3.Heatmap.Cell.Types)]
left_annotation.Colors <- setNames(column_names_gp.col, column_names_gp.col)
left_annotation <- rowAnnotation(Figure.3.Heatmap.Cell.Types = left_annotation.Colors,
                                 col = list(Figure.3.Heatmap.Cell.Types = left_annotation.Colors), simple_anno_size = unit(0.25, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.3.Heatmap <- Heatmap(Cell.Types.NV_by.Cell.Types, left_annotation = left_annotation, rect_gp = gpar(type = "none"), column_names_gp = gpar(col = column_names_gp.col, fontsize = 15, fontface = "bold"), column_names_rot = 45, show_column_dend = FALSE, show_row_names = FALSE, show_heatmap_legend = FALSE, width = unit(7, "in"), height = unit(7, "in"),
                            cell_fun = function(j, i, x, y, w, h, fill) {if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))}})
Figure.3.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3E", ".Heatmap.pdf")
pdf(Figure.3.Heatmap.directory, width = 8, height = 8.5)
print(Figure.3.Heatmap)
dev.off() 

#Figure.3G.UpSetPlot
Idents(mLu.combined.IMs) <- "IMs.Chemokine.Subtypes"
mLu.IMs.Chemokine.Subtypes.DEGs <- FindAllMarkers(mLu.combined.IMs, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
UpSetPlot.Data <- list()
for (i in 1:length(paste0("IMck", 0:9))){
  UpSetPlot.Data[[i]] <- mLu.IMs.Chemokine.Subtypes.DEGs$gene[mLu.IMs.Chemokine.Subtypes.DEGs$cluster == paste0("IMck", 0:9)[i] & mLu.IMs.Chemokine.Subtypes.DEGs$p_val_adj < 0.05]
}
names(UpSetPlot.Data) <- paste0("IMck", 0:9)
UpSetPlot.Data <- make_comb_mat(UpSetPlot.Data, mode = "intersect")
UpSetPlot.Data <- UpSetPlot.Data[order(comb_size(UpSetPlot.Data), decreasing = TRUE)[1:40]]
comb_name_order <- character(0)
for (i in 1:10) {
  vec <- rep("0", 10)
  vec[i] <- "1"
  comb_name_order <- c(comb_name_order, paste(vec, collapse = ""))
}
comb_col <- ifelse(comb_name(UpSetPlot.Data) %in% comb_name_order, comb_name(UpSetPlot.Data), "grey60")
comb_col <- mapvalues(comb_col, from = comb_name_order, to = names(meta.data.list$IMs.Chemokine.Subtypes))
Figure.3.UpSetPlot <- UpSet(UpSetPlot.Data, comb_order = order(comb_size(UpSetPlot.Data), decreasing = TRUE), pt_size = unit(4, "mm"), row_names_gp = gpar(fontsize = 10), heatmap_width = unit(8, "in"), heatmap_height = unit(3, "in"), comb_col = comb_col,
                            right_annotation = upset_right_annotation(UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "top", gp = gpar(fontsize = 8)), width = unit(2, "cm"), gp = gpar(fill = names(meta.data.list$IMs.Chemokine.Subtypes))),
                            top_annotation = upset_top_annotation(UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "left", gp = gpar(fontsize = 8)), height = unit(3.75, "cm"), anno_size = unit(2, "in"), gp = gpar(fill = gsub("grey60", "grey90", comb_col))))
Figure.3.UpSetPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3G", ".UpSetPlot.pdf")
pdf(Figure.3.UpSetPlot.directory, width = 8.5, height = 3.25)
print(Figure.3.UpSetPlot)
dev.off()

#Figure.3F.Heatmap
Heatmap.Features <- c()
for (i in 1:length(paste0("IMck", 0:9))){
  No.Chemokine.Features <- setdiff(extract_comb(UpSetPlot.Data, comb_name_order[i]), Chemokine.Ligands)
  Heatmap.Features <- c(Heatmap.Features, mLu.IMs.Chemokine.Subtypes.DEGs$gene[mLu.IMs.Chemokine.Subtypes.DEGs$cluster == paste0("IMck", 0:9)[i] & mLu.IMs.Chemokine.Subtypes.DEGs$gene %in% No.Chemokine.Features][1:5])
}
Heatmap.Data <- AverageExpression(mLu.combined.IMs, assays = "SCT", slot = "data", group.by = "IMs.Chemokine.Subtypes")
Heatmap.Data <- Heatmap.Data[["SCT"]][Heatmap.Features, ]
Heatmap.Data <- t(scale(t(Heatmap.Data)))
top_annotation <- columnAnnotation(IMs.Cell.Types = as.vector(meta.data.list$IMs.Chemokine.Subtypes), col = list(IMs.Cell.Types = setNames(colors, meta.data.list$IMs.Chemokine.Subtypes)), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation.Colors <- setNames(names(meta.data.list$IMs.Chemokine.Subtypes), as.vector(meta.data.list$IMs.Chemokine.Subtypes))
right_annotation <- rowAnnotation(IMs.Chemokine.Subtypes = rep(as.vector((meta.data.list$IMs.Chemokine.Subtypes)), each = 5),
                                  col = list(IMs.Chemokine.Subtypes = right_annotation.Colors), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.3.Heatmap <- Heatmap(Heatmap.Data, show_heatmap_legend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, top_annotation = top_annotation, right_annotation = right_annotation, row_names_gp = gpar(fontface = "italic"))
Figure.3.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.3F", ".Heatmap.pdf")
pdf(Figure.3.Heatmap.directory, width = 6, height = 8)
print(Figure.3.Heatmap)
dev.off() 

##Figure.4##################################################################################################################################################################################

# Figure.4B.scVelo.DimPlot
mLu.combined.velocity <- RunSCVELO(
  srt = mLu.combined.velocity, group_by = "IMs.Chemokine.Subtypes", 
  linear_reduction = "pca", nonlinear_reduction = "umap", save = TRUE, palcolor = names(meta.data.list$IMs.Chemokine.Subtypes)
)
rangediff <- function(x) {diff(range(x))}
range1.01 <- function(x) {c(min(x) - rangediff(x)*0.01, max(x) + rangediff(x)*0.01)}
Figure.4.DimPlot <- DimPlot(mLu.combined.velocity, reduction = "umap", group.by = "IMs.Chemokine.Subtypes", cols = alpha(names(meta.data.list$IMs.Chemokine.Subtypes), 0.75), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = 4) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(mLu.combined.velocity@reductions$umap@cell.embeddings[, "UMAP_1"])) + ylim(range1.01(mLu.combined.velocity@reductions$umap@cell.embeddings[, "UMAP_2"])) + coord_cartesian(expand = FALSE) & theme(panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA))
Figure.4.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.4B.DimPlot.png")
png(Figure.4.DimPlot.directory, width = 1800, height = 1500, bg = "transparent")
print(Figure.4.DimPlot)
dev.off() 
#Figure.4B.PAGA.DimPlot
PAGAPlot_build <- ggplot_build(PAGAPlot(srt = mLu.combined.velocity, reduction = "umap", label = TRUE, label_insitu = TRUE, label_repel = TRUE) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
                                 xlim(range1.01(mLu.combined.velocity@reductions$umap@cell.embeddings[, "UMAP_1"])) + ylim(range1.01(mLu.combined.velocity@reductions$umap@cell.embeddings[, "UMAP_2"])) + coord_cartesian(expand = FALSE) & theme(panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA)))
dot.x <- PAGAPlot_build$data[[5]]$x
PAGAPlot_build$data[[5]] <- NULL
PAGAPlot_build$data[[4]] <- NULL
PAGAPlot_build$data[[2]] <- NULL
size <- log(table(mLu.combined.velocity$IMs.Chemokine.Subtypes))/3
PAGAPlot_build$data[[2]]$size <- as.numeric(size)^4
PAGAPlot_build$data[[2]]$size[10] <- 15
PAGAPlot_build$data[[2]]$size[5] <- 20
PAGAPlot_build$data[[2]]$fill <- PAGAPlot_build$data[[2]]$colour
PAGAPlot_build$data[[2]]$shape <- 21
PAGAPlot_build$data[[2]]$colour <- "grey"
linewidth_new <- seq(from = 5, to = 0.5, length.out = nrow(PAGAPlot_build$data[[1]]))
dot.cluster <- list(c(2, 7), c(6, 9), c(6, 8), c(0, 8), c(0, 1), c(0, 5), c(0, 3), c(0, 2),    c(1, 7), c(1, 4), c(2, 4), c(4, 6), c(8, 9), c(1, 2), c(5, 6), c(1, 8), c(0, 9), c(1, 5))
for (i in 1:length(dot.cluster)){
  line.order <- which((PAGAPlot_build$data[[1]]$x == dot.x[dot.cluster[[i]][1] + 1] & PAGAPlot_build$data[[1]]$xend == dot.x[dot.cluster[[i]][2] + 1]) | (PAGAPlot_build$data[[1]]$x == dot.x[dot.cluster[[i]][2] + 1] & PAGAPlot_build$data[[1]]$xend == dot.x[dot.cluster[[i]][1] + 1]))
  PAGAPlot_build$data[[1]]$linewidth_new[line.order] <- linewidth_new[i]
  if (i %in% 1:8){
    PAGAPlot_build$data[[1]]$colour[line.order] <- "black"
  } else {
    PAGAPlot_build$data[[1]]$linetype[line.order] <- 2
  }
}
PAGAPlot_build$data[[1]]$alpha <- 1
PAGAPlot_build$data[[2]]$fill <- names(meta.data.list$IMs.Chemokine.Subtypes)
PAGAPlot_build$data[[2]]$stroke <- 5
PAGAPlot_build$data[[1]]$colour <- gsub("grey40", "grey70", PAGAPlot_build$data[[1]]$colour)
other.line.order <- sample(which(PAGAPlot_build$data[[1]]$colour == "grey90" & PAGAPlot_build$data[[1]]$linetype == 1))
PAGAPlot_build$data[[1]]$linewidth_new[other.line.order] <- linewidth_new[-(1:length(dot.cluster))]
PAGAPlot_build$data[[1]]$linetype[other.line.order] <- 2
PAGAPlot_build$data[[1]] <- PAGAPlot_build$data[[1]][c(which(PAGAPlot_build$data[[1]]$colour != "black"), which(PAGAPlot_build$data[[1]]$colour == "black")), ]
Figure.4.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.4B.PAGA.DimPlot.png")
png(Figure.4.DimPlot.directory, width = 1800, height = 1500, bg = "transparent")
print(grid::grid.draw(ggplot_gtable(PAGAPlot_build)))
dev.off() 
#Figure.4C.ScatterPlot
DefaultAssay(mLu.combined.velocity) <- "spliced"
mLu.combined.velocity <- NormalizeData(mLu.combined.velocity)
DefaultAssay(mLu.combined.velocity) <- "unspliced"
mLu.combined.velocity <- NormalizeData(mLu.combined.velocity)
spliced.data <- t(as.matrix(mLu.combined.velocity[["spliced"]]@data))
colnames(spliced.data) <- paste0(colnames(spliced.data), ".spliced")
unspliced.data <- t(as.matrix(mLu.combined.velocity[["unspliced"]]@data))
colnames(unspliced.data) <- paste0(colnames(unspliced.data), ".unspliced")
mLu.combined.velocity.df <- cbind(as.data.frame(cbind(spliced.data, unspliced.data)), "Cell.Types" = mLu.combined.velocity$IMs.Chemokine.Subtypes)
write.table(t(mLu.combined.velocity.df), "/Volumes/mLu/Analysis/scRNA_seq/Datasets/velocity/mLu.combined.velocity.exprMat.txt", row.names = TRUE, quote = FALSE, sep = "\t")
#Magic in python
#https://magic.readthedocs.io/en/stable/
Figure.4.ScatterPlot.Data <- read.csv("/Volumes/mLu/Analysis/scRNA_seq/Datasets/velocity/mLu.Magic.emt_magic.csv", row.names = 1)
Figure.4.ScatterPlot.Data <- cbind(Figure.4.ScatterPlot.Data, "Cell.Types" = mLu.combined.velocity$IMs.Chemokine.Subtypes)

Figure.4.ScatterPlot.Features <- c("Ccl12", "Ccl5", "Ccl4", "Cxcl2", "Ccl8", "Ccl6", "Cxcl10", "Cxcl13", "Ccl24")
Figure.4.ScatterPlot.Features.Order <- c()
for (i in 1:length(Figure.4.ScatterPlot.Features)) {
  Figure.4.ScatterPlot.Features.Order <- c(Figure.4.ScatterPlot.Features.Order, grep(Figure.4.ScatterPlot.Features[i], colnames(Figure.4.ScatterPlot.Data)))
}
for (i in 1:length(Figure.4.ScatterPlot.Features)) {
  Figure.4.ScatterPlot <- ggplot(Figure.4.ScatterPlot.Data, aes(x = Figure.4.ScatterPlot.Data[, Figure.4.ScatterPlot.Features.Order[i*2-1]], y = Figure.4.ScatterPlot.Data[, Figure.4.ScatterPlot.Features.Order[i*2]], color = Cell.Types)) + geom_point(size = 3) + scale_color_manual(values = names(meta.data.list$IMs.Chemokine.Subtypes)) + theme_light() + 
    theme(axis.text.x = element_text(size = 30, color = "black"), axis.text.y = element_text(size = 30, color = "black"), axis.title = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(0.15, "cm"), axis.ticks = element_line(linewidth = 1, colour = "black"), legend.position = "none")
  Figure.4.ScatterPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Revision.Figure.4.", Figure.4.ScatterPlot.Features[i], ".ScatterPlot", ".pdf")
  pdf(Figure.4.ScatterPlot.directory, width = 8, height = 8)
  print(Figure.4.ScatterPlot)
  dev.off() 
}

##Figure.5##Figure.S6##Figure.S7##############################################################################################################################################################################

##SENIC Regulon Analysis
mLu.combined.IMs <- subset(mLu.combined, subset = IMs.Chemokine.Subtypes %in% paste0("IMck", 0:9))
mLu.combined.IMs.exprMat <- GetAssayData(object = mLu.combined.IMs, assay= "SCT", slot = "counts")
mLu.combined.IMs.exprMat <- as(Class = "matrix", object = mLu.combined.IMs.exprMat)
write.table(mLu.combined.IMs.exprMat, "/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/mLu.combined.IMs.exprMat.txt", row.names = TRUE, sep = "\t")
#pyscenic in python
#https://pyscenic.readthedocs.io/en/latest/
#Python3 analysis under pyscenic environment using Jupyter Notebook
#https://mp.weixin.qq.com/s/gTOKdqawzqBPNokc-aEomQ
#https://www.jianshu.com/p/eccfe2d1b2c7
#https://pyscenic.readthedocs.io/en/latest/tutorial.html
mLu.combined.IMs.regulon.AUC <- importAUCfromText(file.path("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/",  "mLu.combined.IMs.10k.10k.auc.tsv"))
# mLu.combined.IMs.regulon.AUC <- importAUCfromText(file.path("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/",  "mLu.combined.IMs.500.100.auc.tsv"))
mLu.combined.IMs.regulon.AUC <- mLu.combined.IMs.regulon.AUC[onlyNonDuplicatedExtended(rownames(mLu.combined.IMs.regulon.AUC)),]
cellInfo <- data.frame(IMs.Chemokine.Subtypes = factor(mLu.combined.IMs$IMs.Chemokine.Subtypes, levels = paste0("IMck", 0:9)))
mLu.combined.IMs.regulonActivity <- sapply(split(rownames(cellInfo), cellInfo$IMs.Chemokine.Subtypes),
                                           function(cells) rowMeans(getAUC(mLu.combined.IMs.regulon.AUC)[,cells]))
mLu.combined.IMs.regulonActivity_Scaled <- t(scale(t(mLu.combined.IMs.regulonActivity), center = T, scale=T))
mLu.combined.IMs.regulonActivity_Scaled[is.na(mLu.combined.IMs.regulonActivity_Scaled)] = 0
rownames(mLu.combined.IMs.regulonActivity_Scaled) <- gsub("\\s*\\([^\\)]+\\)", "", rownames(mLu.combined.IMs.regulonActivity_Scaled))
Top10.Regulons <- vector()
for (i in 1:ncol(mLu.combined.IMs.regulonActivity_Scaled)){
  Top10.Regulons <- c(Top10.Regulons, names(sort(mLu.combined.IMs.regulonActivity_Scaled[, i], decreasing = TRUE))[1:10])
}
names(Top10.Regulons) <- rep(colnames(mLu.combined.IMs.regulonActivity_Scaled), each = 10)
Top10.Regulons <- Top10.Regulons[!duplicated(Top10.Regulons)]
mLu.combined.IMs.regulon <- read.csv("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/mLu.combined.IMs.10k.10k.regulons.csv", header = FALSE)
# mLu.combined.IMs.regulon <- read.csv("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/mLu.combined.IMs.500.100.regulons.csv", header = FALSE)
mLu.combined.IMs.regulon.Vector <- vector()
mLu.combined.IMs.regulon.List <- list()
for (i in 1:nrow(mLu.combined.IMs.regulon)){
  mLu.combined.IMs.regulon.Vector <- gsub("\\[|\\]", "", mLu.combined.IMs.regulon[i, 1]) # remove square brackets
  mLu.combined.IMs.regulon.List[[i]] <- strsplit(mLu.combined.IMs.regulon.Vector, ",\\s*")[[1]] # split by comma and trim whitespace
}
names(mLu.combined.IMs.regulon.List) <- rownames(mLu.combined.IMs.regulon.AUC)
max_length <- max(lengths(mLu.combined.IMs.regulon.List)) # fill the uneven column with NA
mLu.combined.IMs.regulon.List.csv <- lapply(mLu.combined.IMs.regulon.List, function(x) {
  if(length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))
  } else {
    x
  }
})
# Do the same for mLu.combined.IMs.regulon.Chemokine.List
Chemokine.Ligands <- c("Ccl12", "Ccl7", "Ccl2", "Cxcl14", "Ccl5", "Ccl3", "Ccl4", "Cxcl2", "Cxcl3", "Ccl8", "Ccl6", "Ccl9", "Cxcl10", "Cxcl9", "Cxcl13", "Ccl24", "Pf4", "Cxcl16", "Cklf") 
mLu.combined.IMs.regulon.Chemokine.List <- lapply(mLu.combined.IMs.regulon.List, function(x) {
  matches <- which(x %in% Chemokine.Ligands)
  if (length(matches) > 0) {
    return(setNames(x[matches], names(x)))
  } else {
    return(NULL)
  }
})
mLu.combined.IMs.regulon.Chemokine.List <- mLu.combined.IMs.regulon.Chemokine.List[!sapply(mLu.combined.IMs.regulon.Chemokine.List, is.null)]
names(mLu.combined.IMs.regulon.Chemokine.List) <- gsub("\\s*\\([^\\)]+\\)", "", names(mLu.combined.IMs.regulon.Chemokine.List))
# col_fun = colorRamp2(c(min(mLu.combined.IMs.regulonActivity_Scaled[Top10.Regulons, ]), 0, max(mLu.combined.IMs.regulonActivity_Scaled[Top10.Regulons, ])), c("blue", "white", "red"))
row.IMs.Chemokine.Subtypes <- vector()
for (i in 1:ncol(mLu.combined.IMs.regulonActivity_Scaled)){
  row.IMs.Chemokine.Subtypes <- c(row.IMs.Chemokine.Subtypes, rep(colnames(mLu.combined.IMs.regulonActivity_Scaled)[i], table(names(Top10.Regulons))[colnames(mLu.combined.IMs.regulonActivity_Scaled)][i]))
}
mLu.combined.IMs.AverageExpression <- AverageExpression(mLu.combined.IMs, assays = "SCT", slot = "data", group.by = "IMs.Chemokine.Subtypes")
Figure.Regulons.Heatmap.data <- list("Heatmap.data" = mLu.combined.IMs.regulonActivity_Scaled, "Top10.Regulons" = Top10.Regulons, "row.IMs.Chemokine.Subtypes" = row.IMs.Chemokine.Subtypes, "mLu.combined.IMs.regulon.Chemokine.List" = mLu.combined.IMs.regulon.Chemokine.List)
mLu.combined.IMs_Regulon.Heatmap.data.CK.table <- lapply(Figure.Regulons.Heatmap.data$mLu.combined.IMs.regulon.Chemokine.List, function(x) {
  if (is.character(x)) {
    toString(x)
  } else {
    x}})

# Figure.S6A_10k.10k_top10_Regulons.Heatmap.pdf
# Figure.5A_500.100_top10_Regulons.Heatmap.pdf
col_fun = colorRamp2(c(-2.4, 0, 2.8), c("blue", "white", "red"))
Heatmap.anno_mark.at <- which(Figure.Regulons.Heatmap.data[["Top10.Regulons"]] %in% names(Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]]))
Heatmap.anno_mark.labels <- Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]][Figure.Regulons.Heatmap.data[["Top10.Regulons"]][Heatmap.anno_mark.at]]
Heatmap.anno_mark.labels <- as.vector(unlist(lapply(Heatmap.anno_mark.labels, function(x) paste(x, collapse = "\n"))))
Heatmap.btm <- 1 - (Heatmap.anno_mark.at / length(Figure.Regulons.Heatmap.data[["Top10.Regulons"]]))
Heatmap.top <- Heatmap.btm + (1/length(Figure.Regulons.Heatmap.data[["Top10.Regulons"]]))
Heatmap.box_col <- names(meta.data.list$IMs.Chemokine.Subtypes[match(Figure.Regulons.Heatmap.data[["row.IMs.Chemokine.Subtypes"]][Heatmap.anno_mark.at], meta.data.list$IMs.Chemokine.Subtypes)])
Heatmap.rownames.col <- rep("grey50", length(Figure.Regulons.Heatmap.data[["Top10.Regulons"]]))
Heatmap.rownames.col[Heatmap.anno_mark.at] <- "black"
Heatmap.rownames.fontface <- rep("italic", length(Figure.Regulons.Heatmap.data[["Top10.Regulons"]]))
Heatmap.rownames.fontface[Heatmap.anno_mark.at] <- "bold.italic"
Heatmap.ca <- columnAnnotation(IMs.Cell.Types = as.vector(meta.data.list$IMs.Chemokine.Subtypes), col = list(IMs.Cell.Types = setNames(names(meta.data.list$IMs.Chemokine.Subtypes), meta.data.list$IMs.Chemokine.Subtypes)), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Heatmap.ra.Colors <- setNames(names(meta.data.list$IMs.Chemokine.Subtypes), as.vector(meta.data.list$IMs.Chemokine.Subtypes))
#Figure.5A_500.100_top10_Regulons.Heatmap
#Figure.S6A_10k.10k_top10_Regulons.Heatmap
Figure.Regulons.Heatmap.ra <- rowAnnotation(IMs.Chemokine.Subtypes = Figure.Regulons.Heatmap.data[["row.IMs.Chemokine.Subtypes"]], Chemokine = anno_mark(at = Heatmap.anno_mark.at, labels = Heatmap.anno_mark.labels, labels_gp = gpar(fontsize = 12, fontface = "bold.italic"), padding = unit(5, "mm"), link_width = unit(15, "mm")), 
                                            col = list(IMs.Chemokine.Subtypes = Heatmap.ra.Colors), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.Regulons.Heatmap <- Heatmap(Figure.Regulons.Heatmap.data[["Heatmap.data"]][Figure.Regulons.Heatmap.data[["Top10.Regulons"]], ], name = "Regulon", cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, heatmap_width = unit(7, "in"), heatmap_height = unit(12, "in"), col = col_fun, right_annotation = Figure.Regulons.Heatmap.ra, top_annotation = Heatmap.ca, 
                                   row_names_gp = gpar(col = Heatmap.rownames.col, fontsize = 10, fontface = Heatmap.rownames.fontface), show_column_names = FALSE)
Figure.Regulons.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S6A_10k.10k_top10_Regulons.Heatmap.pdf")
# Figure.Regulons.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.5A_500.100_top10_Regulons.Heatmap.pdf")
pdf(Figure.Regulons.Heatmap.directory, width = 8, height = 13)
print(Figure.Regulons.Heatmap)
print(decorate_heatmap_body("Regulon", { for (i in 1:length(Heatmap.anno_mark.at)) {
  grid.lines(c(0, 1), c(Heatmap.top[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(0, 1), c(Heatmap.btm[i],Heatmap.btm[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(0, 0), c(Heatmap.btm[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(1, 1), c(Heatmap.btm[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
}}))
dev.off()

#Figure.5B_500.100_top10_AverageExpression.Heatmap
#Figure.S6B_10k.10k_top10_AverageExpression.Heatmap
mLu.combined.IMs.AverageExpression <- AverageExpression(mLu.combined.IMs, assays = "SCT", slot = "data", group.by = "IMs.Chemokine.Subtypes")
mLu.combined.IMs.AverageExpression.Heatmap.data <- mLu.combined.IMs.AverageExpression[["SCT"]][Figure.Regulons.Heatmap.data[["Top10.Regulons"]], ]
mLu.combined.IMs.AverageExpression.Heatmap.data <- t(scale(t(mLu.combined.IMs.AverageExpression.Heatmap.data)))
Figure.AverageExpression.Heatmap <- Heatmap(mLu.combined.IMs.AverageExpression.Heatmap.data, name = "AverageExpression", cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, heatmap_width = unit(7, "in"), heatmap_height = unit(12, "in"), left_annotation = Figure.Regulons.Heatmap.ra, top_annotation = Heatmap.ca, 
                                            row_names_gp = gpar(col = Heatmap.rownames.col, fontsize = 10, fontface = Heatmap.rownames.fontface), row_names_side = "left", show_column_names = FALSE)
Figure.AverageExpression.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S6B_10k.10k_top10_AverageExpression.Heatmap.pdf")
# Figure.AverageExpression.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.5B_500.100_top10_AverageExpression.Heatmap.pdf")
pdf(Figure.AverageExpression.Heatmap.directory, width = 8, height = 13)
print(Figure.AverageExpression.Heatmap)
print(decorate_heatmap_body("AverageExpression", { for (i in 1:length(Heatmap.anno_mark.at)) {
  grid.lines(c(0, 1), c(Heatmap.top[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(0, 1), c(Heatmap.btm[i],Heatmap.btm[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(0, 0), c(Heatmap.btm[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
  grid.lines(c(1, 1), c(Heatmap.btm[i],Heatmap.top[i]), gp = gpar(lty = 1, lwd = 2, col = Heatmap.box_col[i]))
}}))
dev.off()

#Figure.S7A_500.100_CK_Regulons.Heatmap
#Figure.S7B_10k.10k_CK_Regulons.Heatmap
mLu.combined.IMs.CK.Regulon.Heatmap.data <- Figure.Regulons.Heatmap.data[["Heatmap.data"]][names(Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]]), ]
Figure.Regulons.CK.Heatmap.anno_mark.at <- 1:length(Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]])
Figure.Regulons.CK.Heatmap.anno_mark.labels <- as.vector(unlist(lapply(Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]], function(x) paste(x, collapse = "\n"))))
# Figure.Regulons.CK.Heatmap.anno_mark.labels <- gsub("\n", " ", Figure.Regulons.CK.Heatmap.anno_mark.labels) # Figure.S7A_500.100_CK_Regulons.Heatmap
# Figure.Regulons.CK.Heatmap.anno_mark.labels <- as.vector(sapply(Figure.Regulons.CK.Heatmap.anno_mark.labels, FUN = function(str) {gsub(paste0("([^ ]+( +[^ ]+){",3-1,"}) +"), "\\1\n", str)})) # Figure.S7A_500.100_CK_Regulons.Heatmap
Figure.Regulons.CK.Heatmap.ra <- rowAnnotation(Chemokine = anno_mark(at = Figure.Regulons.CK.Heatmap.anno_mark.at, labels = Figure.Regulons.CK.Heatmap.anno_mark.labels, labels_gp = gpar(fontsize = 12, fontface = "bold.italic"), padding = unit(4, "mm"), link_width = unit(15, "mm")), show_legend = FALSE, show_annotation_name = FALSE)
# Figure.Regulons.CK.Heatmap.ra <- rowAnnotation(Chemokine = anno_mark(at = Figure.Regulons.CK.Heatmap.anno_mark.at, labels = Figure.Regulons.CK.Heatmap.anno_mark.labels, labels_gp = gpar(fontsize = 12, fontface = "bold.italic"), padding = unit(3, "mm"), link_width = unit(15, "mm")), show_legend = FALSE, show_annotation_name = FALSE)
Figure.Regulons.CK.Heatmap <- Heatmap(mLu.combined.IMs.CK.Regulon.Heatmap.data, cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, heatmap_width = unit(7, "in"), heatmap_height = unit(12, "in"), col = col_fun, right_annotation = Figure.Regulons.CK.Heatmap.ra, top_annotation = Heatmap.ca, 
                                      row_names_gp = gpar(col = "black", fontsize = 12, fontface = "bold.italic"), show_column_names = FALSE)
Figure.Regulons.CK.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S7B_10k.10k_CK_Regulons.Heatmap.pdf")
# Figure.Regulons.CK.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S7A_500.100_CK_Regulons.Heatmap.pdf")
pdf(Figure.Regulons.CK.Heatmap.directory, width = 8, height = 13)
print(Figure.Regulons.CK.Heatmap)
dev.off()

#Figure.S7A_500.100_CK_AverageExpression.Heatmap
#Figure.S7B_10k.10k_CK_AverageExpression.Heatmap
mLu.combined.IMs.CK.AverageExpression.Heatmap.data <- mLu.combined.IMs.AverageExpression[["SCT"]][names(Figure.Regulons.Heatmap.data[["mLu.combined.IMs.regulon.Chemokine.List"]]), ]
mLu.combined.IMs.CK.AverageExpression.Heatmap.data <- t(scale(t(mLu.combined.IMs.CK.AverageExpression.Heatmap.data)))
Figure.AverageExpression.CK.Heatmap <- Heatmap(mLu.combined.IMs.CK.AverageExpression.Heatmap.data, name = "AverageExpression", cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, heatmap_width = unit(7, "in"), heatmap_height = unit(12, "in"), left_annotation = Figure.Regulons.CK.Heatmap.ra, top_annotation = Heatmap.ca, 
                                               row_names_gp = gpar(col = "black", fontsize = 12, fontface = "bold.italic"), row_names_side = "left", show_column_names = FALSE)
Figure.AverageExpression.CK.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S7B_10k.10k_CK_AverageExpression.Heatmap.pdf")
# Figure.AverageExpression.CK.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S7A_500.100_CK_AverageExpression.Heatmap.pdf")
pdf(Figure.AverageExpression.CK.Heatmap.directory, width = 8, height = 13)
print(Figure.AverageExpression.CK.Heatmap)
dev.off()

##Figure.6##################################################################################################################################################################################

#Figure.6A.FeaturePlot
Figure.6.FeaturePlot.Features <- c("Lyve1", "Pdpn", "Nes", "Pf4")
for (i in 1:length(Figure.6.FeaturePlot.Features)){
  Figure.6.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.6.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.6.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.6A", ".", Figure.6.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.6.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.6.FeaturePlot)
  dev.off() 
}
#Figure.6C.Heatmap
GSE127267.raw.data <- readRDS("/Volumes/Public/ImmGen/GSE127267.SummarizedExperiment.rds")
GSE127267.target.raw.data <- GSE127267.raw.data@assays@data@listData$logcounts[c("Lyve1", "Pdpn", "Nes", "Pf4"), ]
GSE127267.target.average.data <- sapply(names(table(colData(GSE127267.raw.data)$Population)), function(x) apply(as.matrix(GSE127267.target.raw.data[, grep(x, colnames(GSE127267.target.raw.data))]), 1, mean))
GSE127267.target.average.data.colnames <- colnames(GSE127267.target.average.data)[endsWith(colnames(GSE127267.target.average.data), "Lu") | (!grepl("^MF.|MG.|BEC.|Ep.|FRC.|IAP.|LEC.|LTHSC.|STHSC.|MMP.", colnames(GSE127267.target.average.data)))] # no stromal and stem cells "BEC", "Ep", "FRC", "IAP", "LEC", "LTHSC", "STHSC", "MMP"
GSE127267.target.average.data.colnames <- GSE127267.target.average.data.colnames[order(factor(str_extract(GSE127267.target.average.data.colnames, pattern = "^[a-zA-Z]+"), levels = c("MF", "Mo", "DC", "MC", "GN", "Eo", "Ba", "EMP", "ILC", "NK", "NKT", "T", "preT", "Tgd", "Treg", "B", "proB")))] # no MG also no MF in other tissues than lung
Figure.6.Heatmap.data <- GSE127267.target.average.data[, GSE127267.target.average.data.colnames]
Figure.6.Heatmap.column_split <- factor(c(rep("Macrophages", 15), rep("Monocytes", 22), rep("Dendritic Cells", 43), rep("Other Myeloid Cells", 23), rep("Innate Lymphocytes", 36), rep("T Cells", 65), rep("B Cells", 39)), levels = c("Macrophages", "Monocytes", "Dendritic Cells", "Other Myeloid Cells", "Innate Lymphocytes", "T Cells", "B Cells"))
Figure.6.Heatmap.row_split <- factor(c(rep("Others", 3), "Pf4"), levels = c("Others", "Pf4"))
Figure.6.Heatmap.Annotation <- columnAnnotation(Cell.Types = anno_block(gp = gpar(col = hue_pal()(7), fill = hue_pal()(7))))
Figure.6.Heatmap <- Heatmap(Figure.6.Heatmap.data, cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = FALSE, column_split = Figure.6.Heatmap.column_split, row_split = Figure.6.Heatmap.row_split, column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), row_title = NULL, row_names_gp = gpar(fontsize = 18), column_title_gp = gpar(fontsize = 18), column_title_rot = 90, show_heatmap_legend = FALSE, width = unit(8, "in"), height = unit(4, "in"), top_annotation = Figure.6.Heatmap.Annotation)
Figure.6.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.6C.Heatmap.pdf")
pdf(Figure.6.Heatmap.directory, width = 12, height = 13)
print(Figure.6.Heatmap)
dev.off()
#Figure.6C.Pf4.Heatmap
Figure.6.Pf4.Heatmap.data <- t(as.matrix(GSE127267.target.average.data["Pf4", GSE127267.target.average.data.colnames[grep("^MF.|Mo.|DC.", GSE127267.target.average.data.colnames)]]))
rownames(Figure.6.Pf4.Heatmap.data) <- "Pf4"
Figure.6.Pf4.Heatmap.column_split <- factor(c(rep("Macrophages", 15), rep("Monocytes", 22), rep("Dendritic Cells", 43)), levels = c("Macrophages", "Monocytes", "Dendritic Cells"))
Figure.6.Pf4.Heatmap.Annotation <- columnAnnotation(Cell.Types = anno_block(gp = gpar(col = hue_pal()(7)[1:3], fill = hue_pal()(7)[1:3])))
Figure.6.Pf4.Heatmap.column.names.color <- c(rep("black", 10), rep("grey50", 5), rep("grey50", 65))
Figure.6.Pf4.Heatmap.column.names.fontface <- c(rep("bold", 10), rep("plain", 5), rep("plain", 65))
Figure.6.Pf4.Heatmap <- Heatmap(Figure.6.Pf4.Heatmap.data, cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = TRUE, column_split = Figure.6.Pf4.Heatmap.column_split, column_gap = unit(2, "mm"), row_title = NULL, row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(col = Figure.6.Pf4.Heatmap.column.names.color, fontface = Figure.6.Pf4.Heatmap.column.names.fontface, fontsize = 7), column_names_side = "bottom", column_names_rot = 60, column_title_gp = gpar(fontsize = 14), show_heatmap_legend = FALSE, width = unit(8, "in"), height = unit(1, "in"), top_annotation = Figure.6.Pf4.Heatmap.Annotation)
Figure.6.Pf4.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.6C.Pf4.Heatmap.pdf")
pdf(Figure.6.Pf4.Heatmap.directory, width = 10, height = 4)
print(Figure.6.Pf4.Heatmap)
dev.off()

##Figure.S1##################################################################################################################################################################################

#Figure.S1.Features
mLu.combined.Mac <- subset(mLu.combined, subset = IMs.Subtypes %in% c("CD206hi.IMs", "CD206lo.IMs", "ncIM", "recMac"))
mLu.combined.Mac$Orig.Treatment <- factor(mapvalues(mLu.combined.Mac$orig.treatment, "aGr1.LPS", "LPS"), levels = c("Naive", "LPS"))
mLu.combined.Mac <- PrepSCTFindMarkers(mLu.combined.Mac)
Idents(mLu.combined.Mac) <- mLu.combined.Mac$IMs.Subtypes
mLu.combined.Mac.Cell.Types.DEGs <- FindAllMarkers(mLu.combined.Mac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # for Figure.S1E # repeat for other comparison in Figure.S1
mLu.combined.Mac.Cell.Types.DEGs.top10 <- mLu.combined.Mac.Cell.Types.DEGs %>% group_by(cluster) %>% slice_head(n = 10)
Figure.S1.DEGs.DotPlot.Features <- setNames(mLu.combined.Mac.Cell.Types.DEGs.top10$gene, mLu.combined.Mac.Cell.Types.DEGs.top10$cluster)
Figure.S1.DEGs.DotPlot.Features <- Figure.S1.DEGs.DotPlot.Features[!duplicated(Figure.S1.DEGs.DotPlot.Features)]
Figure.S1.DEGs.DotPlot.Features <- lapply(split(Figure.S1.DEGs.DotPlot.Features, names(Figure.S1.DEGs.DotPlot.Features)), unname)
Figure.S1.DEGs.DotPlot.Features <- Figure.S1.DEGs.DotPlot.Features[c("CD206hi.IMs", "CD206lo.IMs", "ncIM", "recMac")]
Figure.S1.DEGs.DotPlot.Features <- list("CD206hi.IMs" = c("Apoe", "Sepp1", "Cbr2", "F13a1", "Mrc1", "Folr2", "Lyve1", "Cd163", "Ccl8", "Gas6"), "CD206lo.IMs" = c("Ccl5", "Ccl12", "Ccl2", "Ccl7", "Cd72", "Zmynd15", "Lpl", "Ch25h", "Lmna", "Coro1a"), "ncIM" = c("Ear2", "AA467197", "Cxcl3", "Slc7a11", "Lgals3", "Cd9", "Mmp12", "Psap", "Slc7a2", "Plxnd1"), "recMac" = c("Fn1", "Plac8", "S100a6", "Il1b", "S100a4", "Spp1", "Ccr2", "Ly6c2", "Arg1"))
Figure.S1.DEGs.DotPlot.Features.subset <- list("Naive" = list("CD206hi.IMs" = c("Folr2", "F13a1", "Cbr2", "Gas6", "Pf4", "Cd163", "Mrc1", "Ly6e", "Ninj1", "Ltc4s"), "CD206lo.IMs" = c("H2-Eb1", "H2-Ab1", "Cd74", "H2-Aa", "H2-DMb1", "Mmp12", "Ccr2", "Lsp1", "Cd9", "H2-DMa"), "ncIM" = c("Gpnmb", "Cd200r4", "Ahnak2", "Trib3", "Ctsk", "Rgcc", "Cd300lf", "Arhgap24", "Vegfa", "Cd2"), "recMac" = c("Upb1", "Hr", "F10", "Cd226", "Ltb4r1", "Fn1", "Ear2", "Cd24a", "Bhlhe40", "Atp1a3")), 
                                               "LPS" = list("CD206hi.IMs" = c("Cbr2", "F13a1", "Ccl8", "Lyve1", "Mrc1", "Sepp1", "Folr2", "Cd163", "Apoe", "Fcgrt"), "CD206lo.IMs" = c("Ccl5", "Ccl12", "Ccl2", "Ccl7", "Cd72", "Zmynd15", "Ch25h", "Cd81", "Ms4a7", "Lpl"), "ncIM" = c("Ear2", "Cxcl3", "Slc7a11", "AA467197", "Plxnd1", "Lgals3", "Mmp12", "Cd9", "Itgax", "Mcemp1"), "recMac" = c("Fn1", "Plac8", "S100a6", "S100a4", "Il1b", "Spp1", "Ccr2", "H2-Ab1", "Ly6c2")))
#Figure.S1A.DEGs.DotPlot
#Figure.S1B.DEGs.DotPlot
mLu.combined.Mac.subset <- list()
SCTModel.list.list <- list(c(3:6), c(1:2))
for (i in 1:length(c("Naive", "LPS"))){
  mLu.combined.Mac.subset[[i]] <- subset(mLu.combined.Mac, subset = Orig.Treatment == c("Naive", "LPS")[i])
  mLu.combined.Mac.subset[[i]][["SCT"]]@SCTModel.list[SCTModel.list.list[[i]]] <- NULL
  mLu.combined.Mac.subset[[i]] <- PrepSCTFindMarkers(mLu.combined.Mac.subset[[i]])
  Figure.S1.DEGs.DotPlot <- DotPlot(mLu.combined.Mac.subset[[i]], group.by = "IMs.Subtypes", assay = "SCT", features = Figure.S1.DEGs.DotPlot.Features.subset[[i]], dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
  Figure.S1.DEGs.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S1", c("A", "B")[i], ".Mac.", c("Naive", "LPS")[i], ".top10.DEGs.DotPlot.pdf")
  pdf(Figure.S1.DEGs.DotPlot.directory, width = 14, height = 2.6)
  print(Figure.S1.DEGs.DotPlot)
  dev.off() 
}
#Figure.S1E.DEGs.Heatmap
scale.data <- mLu.combined.Mac[["SCT"]]@scale.data[intersect(unlist(Figure.S1.DEGs.DotPlot.Features), rownames(mLu.combined.Mac[["SCT"]]@scale.data)), ]
data <- mLu.combined.Mac[["SCT"]]@data[setdiff(unlist(Figure.S1.DEGs.DotPlot.Features), rownames(mLu.combined.Mac[["SCT"]]@scale.data)), ]
Figure.S1.DEGs.Heatmap.data <- rbind(scale.data, data)[unlist(Figure.S1.DEGs.DotPlot.Features), ]
Figure.S1.DEGs.Heatmap.data.cell.names <- c(sample(colnames(mLu.combined.Mac)[mLu.combined.Mac$IMs.Subtypes == "CD206hi.IMs"], 2000), sample(colnames(mLu.combined.Mac)[mLu.combined.Mac$IMs.Subtypes == "CD206lo.IMs"], 2000), colnames(mLu.combined.Mac)[mLu.combined.Mac$IMs.Subtypes == "ncIM"], sample(colnames(mLu.combined.Mac)[mLu.combined.Mac$IMs.Subtypes == "recMac"], 2000))
Figure.S1.DEGs.Heatmap.data <- Figure.S1.DEGs.Heatmap.data[, Figure.S1.DEGs.Heatmap.data.cell.names]
Figure.S1.DEGs.Heatmap.column_split <- factor(c(rep("CD206hi.IMs", 2000), rep("CD206lo.IMs", 2000), rep("ncIM", 532), rep("recMac", 2000)), levels = c(meta.data.list$IMs.Subtypes[1:4]))
Figure.S1.DEGs.Heatmap.row_split <- factor(c(rep("CD206hi.IMs", 10), rep("CD206lo.IMs", 10), rep("ncIM", 10), rep("recMac", 9)), levels = c(meta.data.list$IMs.Subtypes[1:4]))
colors <- names(meta.data.list$IMs.Subtypes[1:4])
cells <- levels(Figure.S1.DEGs.Heatmap.column_split)
Figure.S1.DEGs.Heatmap.column.Annotation <- columnAnnotation(IMs.Cell.Types = c(rep("CD206hi.IMs", 2000), rep("CD206lo.IMs", 2000), rep("ncIM", 532), rep("recMac", 2000)), col = list(IMs.Cell.Types = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S1.DEGs.Heatmap <- Heatmap(Figure.S1.DEGs.Heatmap.data, use_raster = TRUE, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = FALSE, column_split = Figure.S1.DEGs.Heatmap.column_split, column_gap = unit(2, "mm"), column_title = NULL, row_split = Figure.S1.DEGs.Heatmap.row_split, row_gap = unit(2, "mm"), row_title = NULL, cluster_column_slices = FALSE, row_names_gp = gpar(fontsize = 20, fontface = "italic"), show_heatmap_legend = FALSE, width = unit(10, "in"), height = unit(12, "in"), top_annotation = Figure.S1.DEGs.Heatmap.column.Annotation)
Figure.S1.DEGs.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S1E.Mac.top10.DEGs.Heatmap.pdf")
pdf(Figure.S1.DEGs.Heatmap.directory, width = 12, height = 13)
print(Figure.S1.DEGs.Heatmap)
dev.off()
#Figure.S1C.DEGs.Heatmap
#Figure.S1D.DEGs.Heatmap
Figure.S1.DEGs.Heatmap.data <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/FindAllMarkers/Figure.S1.Mac.Naive.LPS.top10.DEGs.Heatmap.data.rds")
colors <- names(meta.data.list$IMs.Subtypes[1:4])
for (i in 1:length(c("Naive", "LPS"))){
  cells <- names(table(Figure.S1.DEGs.Heatmap.data[["Figure.S1.DEGs.Heatmap.column.split"]][[i]]))
  top_annotation <- columnAnnotation(IMs.Cell.Types = Figure.S1.DEGs.Heatmap.data[["Figure.S1.DEGs.Heatmap.column.split"]][[i]], col = list(IMs.Cell.Types = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
  Figure.S1.DEGs.Heatmap <- Heatmap(as.matrix(Figure.S1.DEGs.Heatmap.data[["Figure.S1.DEGs.Heatmap.data"]][[i]]), use_raster = TRUE, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = FALSE, column_split = Figure.S1.DEGs.Heatmap.data[["Figure.S1.DEGs.Heatmap.column.split"]][[i]], column_gap = unit(2, "mm"), column_title = NULL, row_split = Figure.S1.DEGs.Heatmap.data[["Figure.S1.DEGs.Heatmap.row.split"]][[i]], row_gap = unit(2, "mm"), row_title = NULL, cluster_column_slices = FALSE, row_names_gp = gpar(fontsize = 20, fontface = "italic"), show_heatmap_legend = FALSE, width = unit(10, "in"), height = unit(12, "in"), top_annotation = top_annotation)
  Figure.S1.DEGs.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S1", c("C", "D")[i], ".Mac.", c("Naive", "LPS")[i], ".top10.DEGs.Heatmap.pdf")
  pdf(Figure.S1.DEGs.Heatmap.directory, width = 12, height = 13)
  print(Figure.S1.DEGs.Heatmap)
  dev.off()
}
#Figure.S1F.topGO.Heatmap
mLu.combined.Mac.bg <- rownames(mLu.combined.Mac)[rowSums(mLu.combined.Mac) > 0]
mLu.combined.Mac.topGO.list <- vector()
for (i in c(meta.data.list$IMs.Subtypes[1:4])){
  mLu.combined.Mac.DEG <- mLu.combined.Mac.Cell.Types.DEGs[mLu.combined.Mac.Cell.Types.DEGs$cluster == i, "gene"]
  mLu.combined.Mac.topGO.xls.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Datasets/topGO/Cell.Types/mLu.combined.Mac_", i, ".topGO.xls")
  mLu.combined.Mac.topGO <- topGOtable(mLu.combined.Mac.DEG, mLu.combined.Mac.bg, mapping = "org.Mm.eg.db", ontology = "BP", writeOutput = TRUE, topTablerows = 500, outputFile = mLu.combined.Mac.topGO.xls.directory)
  mLu.combined.Mac.topGO.list <- c(mLu.combined.Mac.topGO.list, mLu.combined.Mac.topGO$GO.ID[1:10])
}
names(mLu.combined.Mac.topGO.list) <- rep(meta.data.list$IMs.Subtypes[1:4], each = 10)
mLu.combined.Mac.topGO.list <- mLu.combined.Mac.topGO.list[!duplicated(mLu.combined.Mac.topGO.list)]
mLu.combined.Mac.topGO.targetGO.Heatmap.table <- data.frame(matrix(ncol = 0, nrow = length(mLu.combined.Mac.topGO.list)))
for (i in c(meta.data.list$IMs.Subtypes[1:4])){
  mLu.combined.Mac.DEG <- mLu.combined.Mac.Cell.Types.DEGs[mLu.combined.Mac.Cell.Types.DEGs$cluster == i, "gene"]
  mLu.combined.Mac.topGO.targetGO.xls.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Datasets/topGO/Cell.Types/mLu.combined.Mac_", i, ".targetGO.xls")
  mLu.combined.Mac.topGO.targetGO <- topGOtable.targetGO(mLu.combined.Mac.DEG, mLu.combined.Mac.bg, mapping = "org.Mm.eg.db", ontology = "BP", TargetGO = mLu.combined.Mac.topGO.list, writeOutput = TRUE, outputFile = mLu.combined.Mac.topGO.targetGO.xls.directory)
  setnafill(mLu.combined.Mac.topGO.targetGO, "nocb", cols = "p.value_elim")
  mLu.combined.Mac.topGO.targetGO.List <- mLu.combined.Mac.topGO.targetGO[match(mLu.combined.Mac.topGO.list, mLu.combined.Mac.topGO.targetGO$GO.ID), ]
  mLu.combined.Mac.topGO.targetGO.Heatmap.table <- cbind(mLu.combined.Mac.topGO.targetGO.Heatmap.table, mLu.combined.Mac.topGO.targetGO.List$p.value_elim)
}
rownames(mLu.combined.Mac.topGO.targetGO.Heatmap.table) <- str_to_sentence(rownames(mLu.combined.Mac.topGO.targetGO.Heatmap.table))
colnames(mLu.combined.Mac.topGO.targetGO.Heatmap.table) <- meta.data.list$IMs.Subtypes[1:4]
mLu.combined.Mac.topGO.targetGO.Heatmap.table <- data.matrix(-log10(mLu.combined.Mac.topGO.targetGO.Heatmap.table), rownames.force = TRUE)
col_fun = colorRamp2(range(mLu.combined.Mac.topGO.targetGO.Heatmap.table), c("white", "red"))
mLu.combined.Mac.topGO.list <- c(rep("CD206hi.IMs", 10), rep("CD206lo.IMs", 10), rep("ncIM", 8), rep("recMac", 8))
colors <- names(meta.data.list$IMs.Subtypes[1:4])
cells <- c("CD206hi.IMs", "CD206lo.IMs", "ncIM", "recMac")
column_annotation <- columnAnnotation(IMs.Cell.Types = c(rep("CD206hi.IMs", 1), rep("CD206lo.IMs", 1), rep("ncIM", 1), rep("recMac", 1)), col = list(IMs.Cell.Types = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation <- rowAnnotation(Mac.Cell.Types = mLu.combined.Mac.topGO.list, col = list(Mac.Cell.Types = setNames(colors, meta.data.list$IMs.Subtypes[1:4])), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S1.Mac.top10.topGO.Heatmap <- Heatmap(mLu.combined.Mac.topGO.targetGO.Heatmap.table, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 20), show_heatmap_legend = FALSE, heatmap_width = unit(10, "in"), heatmap_height = unit(12, "in"), col = col_fun, right_annotation = right_annotation, top_annotation = column_annotation)
Figure.S1.Mac.top10.topGO.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S1F.Mac.top10.topGO.Heatmap.pdf")
pdf(Figure.S1.Mac.top10.topGO.Heatmap.directory, width = 26, height = 13)
print(Figure.S1.Mac.top10.topGO.Heatmap)
dev.off()

##Figure.S2##################################################################################################################################################################################

#Figure.S2A.DotPlot
Growth.Factor.Activity <- sort(unique(AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0008083", columns = c('SYMBOL'), keytype = "GOALL")$SYMBOL))
Growth.Factor.Activity <- setdiff(intersect(Growth.Factor.Activity, rownames(mLu.combined[["RNA"]]@data)), c("Bmp7", "Btc", "Csf2", "Epgn", "Fgf12", "Fgf17", "Gdf10", "Mia", "Ntf5"))
Figure.S2.DotPlot <- DotPlot(mLu.combined, group.by = "Cell.Types", assay = "SCT", features = Growth.Factor.Activity, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
Figure.S2.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S2A.DotPlot.pdf")
pdf(Figure.S2.DotPlot.directory, width = 22.5, height = 3.5)
print(Figure.S2.DotPlot)
dev.off() 
#Figure.S2B.FeaturePlot
Figure.S2.FeaturePlot.Features <- c("Igf1", "Osm", "Tgfb1", "Vegfa")
for (i in 1:length(Figure.S2.FeaturePlot.Features)){
  Figure.S2.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.S2.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.S2.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S2B", ".", Figure.S2.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.S2.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.S2.FeaturePlot)
  dev.off() 
}
#Figure.S2B.VlnPlot # Set min due to long tail
Figure.S2.VlnPlot.Idents <- rep(c("IMs", "ncIM", "recMac", "ncMono"), 4)
Figure.S2.VlnPlot.Features <- rep(c("Igf1", "Osm", "Tgfb1", "Vegfa"), each = 4)
y.max.increment <- c(rep(0.5, 8), rep(0, 4), rep(0.01, 4))
for (i in 1:16){
  Figure.S2.VlnPlot.build <- ggplot_build(VlnPlot(mLu.combined, idents = Figure.S2.VlnPlot.Idents[i], features = Figure.S2.VlnPlot.Features[i], assay = "SCT", slot = c(rep("scale.data", 8), rep("data", 8), rep("scale.data", 4))[i], log = TRUE))
  Figure.S2.VlnPlot.range <- c(min(Figure.S2.VlnPlot.build$data[[2]][, "y"][!is.na(Figure.S2.VlnPlot.build$data[[2]][, "y"])]), max(Figure.S2.VlnPlot.build$data[[2]][, "y"][!is.na(Figure.S2.VlnPlot.build$data[[2]][, "y"])]))
  Figure.S2.VlnPlot <- VlnPlot(mLu.combined, idents = Figure.S2.VlnPlot.Idents[i], features = Figure.S2.VlnPlot.Features[i], assay = "SCT", slot = c(rep("scale.data", 8), rep("data", 8), rep("scale.data", 4))[i], log = TRUE, cols = names(meta.data.list$orig.treatment), pt.size = 0, group.by = "orig.treatment", y.max = (Figure.S2.VlnPlot.range[2] + y.max.increment[i])) + geom_violin(trim = FALSE, scale = "width") + geom_jitter(width = 0.4, size = 1) + 
    theme(plot.title = element_blank(), legend.position = "none", axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) + scale_y_log10(limits = 10^(Figure.S2.VlnPlot.range + y.max.increment[i]), expand = c(0, 0)) +
    annotate(x = 0.42, xend = 0.42, y = 10^quantile(Figure.S2.VlnPlot.range + y.max.increment[i], 0.2), yend = 10^quantile(Figure.S2.VlnPlot.range + y.max.increment[i], 0.8), colour = "black", lwd = 2.5, geom = "segment", arrow = arrow()) + 
    annotate("text", x = 0.27, y = 10^quantile(Figure.S2.VlnPlot.range + y.max.increment[i], 0.5), label = "Expression", angle = 90, size = 15) +
    scale_x_discrete(expand = expansion(mult = c(0.4, 0)))  
  Figure.S2.VlnPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S2B.", Figure.S2.VlnPlot.Idents[i], ".", Figure.S2.VlnPlot.Features[i], ".VlnPlot.pdf")
  pdf(Figure.S2.VlnPlot.directory, width = 18, height = 10)
  print(Figure.S2.VlnPlot)
  dev.off() 
}
#Figure.S2C.FeaturePlot
Figure.S2.FeaturePlot.Features <- c("Arg1", "Nos2")
for (i in 1:length(Figure.S2.FeaturePlot.Features)){
  Figure.S2.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.S2.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.S2.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S2C", ".", Figure.S2.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.S2.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.S2.FeaturePlot)
  dev.off() 
}
#Figure.S2D.Heatmap
M2.vs.M1.DEGs <- read.csv("/Volumes/mLu/Analysis/scRNA_seq/Datasets/FindAllMarkers/M2.vs.M1.DEGs.csv") # PMID: 31178859
M2.vs.M1.DEGs <- M2.vs.M1.DEGs[M2.vs.M1.DEGs$Gene.symbol != " ", ]
M2.vs.M1.DEGs <- M2.vs.M1.DEGs[M2.vs.M1.DEGs$Gene.symbol %in% rownames(mLu.combined[["SCT"]]@scale.data), ]
Figure.S2.Heatmap.data <- mLu.combined[["SCT"]]@scale.data[c(head(M2.vs.M1.DEGs$Gene.symbol, 20), rev(tail(M2.vs.M1.DEGs$Gene.symbol, 20))), ]
Figure.S2.Heatmap.data.cell.names <- c(sample(colnames(mLu.combined)[mLu.combined$Cell.Types == "IMs"], 2000), colnames(mLu.combined)[mLu.combined$Cell.Types == "ncIM"], sample(colnames(mLu.combined)[mLu.combined$Cell.Types == "recMac"], 2000), colnames(mLu.combined)[mLu.combined$Cell.Types == "ncMono"])
Figure.S2.Heatmap.data <- Figure.S2.Heatmap.data[, Figure.S2.Heatmap.data.cell.names]
Figure.S2.Heatmap.data <- Figure.S2.Heatmap.data[c(21:40, 1:20), ]
Figure.S2.Heatmap.column_split <- factor(c(rep("IMs", 2000), rep("ncIM", 532), rep("recMac", 2000), rep("ncMono", 224)), levels = c(meta.data.list$Cell.Types[1:4]))
Figure.S2.Heatmap.row_split <- factor(c(rep("M1 Signatures", 20), rep("M2 Signatures", 20)), levels = c("M1 Signatures", "M2 Signatures"))
Figure.S2.Heatmap.column.Annotation <- columnAnnotation(Cell.Types = c(rep("IMs", 2000), rep("ncIM", 532), rep("recMac", 2000), rep("ncMono", 224)), col = list(Cell.Types = setNames(names(meta.data.list$Cell.Types[1:4]), meta.data.list$Cell.Types[1:4])), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S2.Heatmap.row.Annotation <- rowAnnotation(Polarization = c(rep("M1 Signatures", 20), rep("M2 Signatures", 20)), col = list(Polarization = setNames(c("grey70", "black"), c("M1 Signatures", "M2 Signatures"))), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S2.Heatmap <- Heatmap(Figure.S2.Heatmap.data, use_raster = TRUE, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = FALSE, column_split = Figure.S2.Heatmap.column_split, column_gap = unit(2, "mm"), column_title = NULL, row_split = Figure.S2.Heatmap.row_split, row_gap = unit(2, "mm"), row_title = NULL, cluster_column_slices = TRUE, row_names_gp = gpar(fontsize = 20, fontface = "italic"), show_heatmap_legend = FALSE, width = unit(12, "in"), height = unit(12, "in"), top_annotation = Figure.S2.Heatmap.column.Annotation, left_annotation = Figure.S2.Heatmap.row.Annotation)
Figure.S2.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S2D.Heatmap.pdf")
pdf(Figure.S2.Heatmap.directory, width = 14, height = 13)
print(Figure.S2.Heatmap)
dev.off()

##Figure.S3##########################################################################################################################################################################################################################

#Figure.S3A.FeaturePlot
#Figure.S3C.FeaturePlot
for (j in c("Nos2", "Arg1")){
  for (i in 1:length(meta.data.list$orig.treatment)){
    Figure.S3.FeaturePlot <- FeaturePlot(mLu.combined, features = j, cells = colnames(mLu.combined)[mLu.combined$orig.treatment == meta.data.list$orig.treatment[i]], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 4) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
      xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
    Figure.S3.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S3", ifelse(j == "Nos2", "A", "C"), ".", j, ".FeaturePlot_", meta.data.list$orig.treatment[i], ".pdf")
    pdf(Figure.S3.FeaturePlot.directory, width = 18, height = 15)
    print(Figure.S3.FeaturePlot)
    dev.off() 
  }
}
#Figure.S3B.FeaturePlot
#Figure.S3D.FeaturePlot
Figure.S3.FeaturePlot.Features <- c("Chil3", "Retnla", "Clec10a", "Mgl2", "Cdh1", "Selp", "Pafah1b1", "Pafah1b2", "Psap", "Pmp22", "Folr2", "Trem2")
Figure.S3.FeaturePlot.Features <- c("Il31ra", "Tbx21", "Slc6a9", "Itk", "Gbp5", "Chac1", "Bspry", "Sspn", "Trim2", "Clu", "Ablim1", "Ubd")
for (i in 1:length(Figure.S3.FeaturePlot.Features)){
  Figure.S3.FeaturePlot <- FeaturePlot(mLu.combined, features = Figure.S3.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 4) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.UAMP_1.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.UAMP_2.length), max(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.S3.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S3", ".", Figure.S3.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.S3.FeaturePlot.directory, width = 18, height = 15)
  print(Figure.S3.FeaturePlot)
  dev.off() 
}

##Figure.S4##################################################################################################################################################################################

#Figure.S4A.DimPlot
Figure.S4.DimPlot <- DimPlot(mLu.combined, group.by = "SCT_snn_res.4", raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"])*1.01) + ylim(range(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"])*1.01) + coord_cartesian(expand = FALSE)
Figure.S4.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S4A", ".DimPlot.pdf")
pdf(Figure.S4.DimPlot.directory, width = 18, height = 15)
print(Figure.S4.DimPlot)
dev.off() 
#Figure.S4B.Heatmap
Idents(mLu.combined) <- "SCT_snn_res.4"
mLu.combined.SCT_snn_res.4.DEGs <- FindAllMarkers(mLu.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mLu.combined.SCT_snn_res.4.DEGs.top10 <- mLu.combined.SCT_snn_res.4.DEGs %>% group_by(cluster) %>% slice_head(n = 10)
mLu.combined.SCT_snn_res.4.DEGs.top10 <- setNames(mLu.combined.SCT_snn_res.4.DEGs.top10$gene, mLu.combined.SCT_snn_res.4.DEGs.top10$cluster)
mLu.combined.SCT_snn_res.4.DEGs.top10 <- mLu.combined.SCT_snn_res.4.DEGs.top10[!duplicated(mLu.combined.SCT_snn_res.4.DEGs.top10)]
Heatmap.data <- AverageExpression(mLu.combined, assays = "SCT", slot = "data", group.by = "SCT_snn_res.4")
Heatmap.data <- Heatmap.data[["SCT"]][mLu.combined.SCT_snn_res.4.DEGs.top10, ]
Heatmap.data <- t(scale(t(Heatmap.data)))
Figure.S4.df <- table(mLu.combined$Cell.Types, mLu.combined$SCT_snn_res.4)
Cell.Types <- vector()
for (i in 0:57){
  Cell.Types <- c(Cell.Types, which.max(Figure.S4.df[, i+1]))
}
Figure.S4.Heatmap.data <- list("Heatmap.data" = Heatmap.data, "Cell.Types" = Cell.Types, "mLu.combined.SCT_snn_res.4.DEGs.top10" = mLu.combined.SCT_snn_res.4.DEGs.top10)
Heatmap.Annotation.Cell.Types.Colors <- meta.data.list$Cell.Types[Figure.S4.Heatmap.data[["Cell.Types"]]]
Figure.S4.Heatmap.columnAnnotation <- columnAnnotation(Cell.Types = as.vector(Heatmap.Annotation.Cell.Types.Colors), col = list(Cell.Types = setNames(names(meta.data.list$Cell.Types), meta.data.list$Cell.Types)), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Heatmap.Annotation.SCT_snn_res.4.Colors <- setNames(0:57, names(meta.data.list$Cell.Types)[as.vector(Figure.S4.Heatmap.data[["Cell.Types"]])])
Figure.S4.anno_mark.at <- which(rownames(Figure.S4.Heatmap.data[["Heatmap.data"]]) %in% Chemokine.Ligands)
Figure.S4.ha <- rowAnnotation(SCT_snn_res.4 = names(Figure.S4.Heatmap.data[["mLu.combined.SCT_snn_res.4.DEGs.top10"]]), Chemokine = anno_mark(at = Figure.S4.anno_mark.at, labels = rownames(Figure.S4.Heatmap.data[["Heatmap.data"]])[Figure.S4.anno_mark.at], labels_gp = gpar(fontsize = 15, fontface = "italic"), padding = unit(2, "mm")), 
                              col = list(SCT_snn_res.4 = setNames(names(Heatmap.Annotation.SCT_snn_res.4.Colors), Heatmap.Annotation.SCT_snn_res.4.Colors)), simple_anno_size = unit(0.2, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S4.Heatmap <- Heatmap(Figure.S4.Heatmap.data[["Heatmap.data"]], cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 12), column_names_side = "bottom", column_names_rot = 90, show_heatmap_legend = FALSE, width = unit(8, "in"), height = unit(6, "in"), top_annotation = Figure.S4.Heatmap.columnAnnotation, right_annotation = Figure.S4.ha)
Figure.S4.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S4B.Heatmap.pdf")
pdf(Figure.S4.Heatmap.directory, width = 11, height = 7)
print(Figure.S4.Heatmap)
dev.off()
#Figure.S4C.FeaturePlot
mLu.combined.IMs.UAMP_1.length <- max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]) - min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])
mLu.combined.IMs.UAMP_2.length <- max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]) - min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])
mLu.combined.IMs.UAMP_2.length.times <- (mLu.combined.UAMP_1.length/mLu.combined.UAMP_2.length)*mLu.combined.IMs.UAMP_2.length/mLu.combined.IMs.UAMP_1.length
Figure.S4.FeaturePlot.Features <- c("Ccl24", "Cxcl13", "Cxcl9", "Cxcl10", "Ccl9", "Ccl6", "Ccl8", "Cxcl3", "Cxcl2", "Cxcl1", "Ccl4", "Ccl3", "Ccl5", "Cxcl14", "Ccl2", "Ccl7", "Ccl12", "Pf4", "Cxcl16", "Cklf")
for (i in 1:length(Figure.S4.FeaturePlot.Features)){
  Figure.S4.FeaturePlot <- FeaturePlot(mLu.combined.IMs, features = Figure.S4.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 4) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.IMs.UAMP_1.length), max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.IMs.UAMP_2.length), max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.S4.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S4C.", Figure.S4.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.S4.FeaturePlot.directory, width = 18, height = 15*mLu.combined.IMs.UAMP_2.length.times)
  print(Figure.S4.FeaturePlot)
  dev.off() 
}

##Figure.S5##################################################################################################################################################################################

#Figure.S5A.Trajectory.DimPlot
mLu.combined.CD206hiIMs <- subset(mLu.combined.IMs, subset = IMs.Subtypes == "CD206hi.IMs")
mLu.combined.CD206loIMs <- subset(mLu.combined.IMs, subset = IMs.Subtypes == "CD206lo.IMs")
mLu.combined.CD206hiIMs.cds <- as.cell_data_set(mLu.combined.CD206hiIMs)
mLu.combined.CD206hiIMs.cds <- cluster_cells(mLu.combined.CD206hiIMs.cds, resolution = 5.2e-4)
mLu.combined.CD206hiIMs.cds <- learn_graph(mLu.combined.CD206hiIMs.cds, close_loop = FALSE, use_partition = FALSE)
# root_pr_nodes <- c("Y_140") for mLu.combined.CD206hiIMs.cds
# root_pr_nodes <- c("Y_8") for mLu.combined.CD206loIMs.cds
mLu.combined.CD206hiIMs.cds <- order_cells(mLu.combined.CD206hiIMs.cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)
Other.Data <- mLu.combined.IMs@reductions$umap@cell.embeddings[mLu.combined.IMs$IMs.Subtypes != "CD206hi.IMs", ]
Background <- ggplot() + geom_point(data = data.frame(Other.Data), aes(x = UMAP_1, y = UMAP_2), colour = "#E5E5E5", size = 2.5, stroke = 0) + scale_shape(solid = T)
Figure.S5A <- plot_cells(mLu.combined.CD206hiIMs.cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE, cell_size = 2.5, trajectory_graph_segment_size = 4, trajectory_graph_color = "black", graph_label_size = 5, alpha = 0.7, cell_stroke = 0) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(c(min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mLu.combined.IMs.UAMP_1.length), max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]) + (0.005*mLu.combined.IMs.UAMP_1.length))) + ylim(c(min(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mLu.combined.IMs.UAMP_2.length), max(mLu.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]) + (0.005*mLu.combined.IMs.UAMP_2.length))) + coord_cartesian(expand = FALSE)
Figure.S5A.Build <- ggplot_build(Figure.S5A)
Figure.S5A.Build[[1]][[1]] <- rbind(ggplot_build(Background)[[1]][[1]], Figure.S5A.Build[[1]][[1]])
Figure.S5A.Build[[1]][[4]]$label <- NULL
Figure.S5A.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S5A.", "mLu.combined.CD206hiIMs.Monocle3_", "Trajectory.DimPlot", ".pdf")
Figure.S5A.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S5C.", "mLu.combined.CD206loIMs.Monocle3_", "Trajectory.DimPlot", ".pdf")
pdf(Figure.S5A.directory, width = 18*mLu.combined.IMs.UAMP_1.length/mLu.combined.UAMP_1.length, height = 15*mLu.combined.IMs.UAMP_2.length/mLu.combined.UAMP_2.length)
print(grid::grid.draw(ggplot_gtable(Figure.S5A.Build)))
dev.off() 

#Figure.S5B.Pseudotime.Heatmap
mLu.combined.IMs$IMs.Chemokine.Subtypes <- mLu.combined@meta.data[colnames(mLu.combined.IMs), "IMs.Chemokine.Subtypes"]
mLu.combined.CD206hiIMs$IMs.Chemokine.Subtypes <- mLu.combined@meta.data[colnames(mLu.combined.CD206hiIMs), "IMs.Chemokine.Subtypes"]
mLu.combined.CD206loIMs$IMs.Chemokine.Subtypes <- mLu.combined@meta.data[colnames(mLu.combined.CD206loIMs), "IMs.Chemokine.Subtypes"]
mLu.combined.CD206hiIMs.Pseudotime.Gene.List <- graph_test(mLu.combined.CD206hiIMs.cds, neighbor_graph = "principal_graph", cores = 1)
Pseudotime.data <- list()
Pseudotime.data[["Heatmap.Annotation.Pseudotime"]] <- pseudotime(mLu.combined.CD206hiIMs.cds)[colnames(mLu.combined.CD206hiIMs)]
mLu.combined.CD206hiIMs$Pseudotime.Part <- "End"
mLu.combined.CD206hiIMs$Pseudotime.Part[names(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]])[1:(length(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]])/2)]] <- "Start"
Pseudotime.AverageExpression <- AverageExpression(mLu.combined.CD206hiIMs, assays = "SCT", slot = "scale.data", group.by = "Pseudotime.Part")
mLu.combined.CD206hiIMs.Pseudotime.Gene.List <- mLu.combined.CD206hiIMs.Pseudotime.Gene.List[rownames(Pseudotime.AverageExpression$SCT), ]
mLu.combined.CD206hiIMs.Pseudotime.Gene.List <- mLu.combined.CD206hiIMs.Pseudotime.Gene.List[(Pseudotime.AverageExpression$SCT[, "End"] > Pseudotime.AverageExpression$SCT[, "Start"]), ]
Pseudotime.Gene.Top120 <- rownames(mLu.combined.CD206hiIMs.Pseudotime.Gene.List)[order(mLu.combined.CD206hiIMs.Pseudotime.Gene.List$morans_I, decreasing = TRUE)][1:120]
Pseudotime.data[["Heatmap"]] <- mLu.combined.CD206hiIMs[["SCT"]]@scale.data[intersect(Pseudotime.Gene.Top120, rownames(mLu.combined.CD206hiIMs[["SCT"]]@scale.data))[1:100], ]
Pseudotime.data[["Heatmap.Annotation.IMs.Chemokine.Subtypes"]] <- mLu.combined.CD206hiIMs$IMs.Chemokine.Subtypes
Pseudotime.col_fun <- colorRamp2(seq(from = min(Pseudotime.data[["Heatmap"]]), to = 15, len = 9), brewer.pal(9, "Purples"))
Heatmap.Annotation.Pseudotime.Color <- colorRamp2(seq(from = min(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]]), to = max(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]]), len = 5), c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))
Heatmap.Annotation.Pseudotime.Order <- names(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]][order(Pseudotime.data[["Heatmap.Annotation.Pseudotime"]])])
Figure.S5.Heatmap.Annotation.Pseudotime.Cell.Types = HeatmapAnnotation(Pseudotime = Pseudotime.data[["Heatmap.Annotation.Pseudotime"]], IMs.Chemokine.Subtypes = Pseudotime.data[["Heatmap.Annotation.IMs.Chemokine.Subtypes"]], 
                                                                       col = list(Pseudotime = Heatmap.Annotation.Pseudotime.Color, IMs.Chemokine.Subtypes = setNames(names(meta.data.list$IMs.Chemokine.Subtypes), meta.data.list$IMs.Chemokine.Subtypes)), show_legend = FALSE, annotation_label = c("Pseudotime", "Chemokine-expressing IMs"), annotation_name_gp = gpar(fontsize = 16))
# Figure.S5.Heatmap.Annotation.Pseudotime.Cell.Types = HeatmapAnnotation(Pseudotime = Pseudotime.data[["Heatmap.Annotation.Pseudotime"]], IMs.Chemokine.Subtypes = Pseudotime.data[["Heatmap.Annotation.IMs.Chemokine.Subtypes"]], IMs.Subtypes = Pseudotime.data[["Heatmap.Annotation.IMs.Subtypes"]],  # mLu.combined.IMs
#                                                                                     col = list(Pseudotime = Heatmap.Annotation.Pseudotime.Color, IMs.Chemokine.Subtypes = setNames(names(meta.data.list$IMs.Chemokine.Subtypes), meta.data.list$IMs.Chemokine.Subtypes), IMs.Subtypes = setNames(names(meta.data.list$IMs.Subtypes[1:2]), meta.data.list$IMs.Subtypes[1:2])), show_legend = FALSE, annotation_label = c("Pseudotime", "Chemokine-expressing IMs", "IM subtypes"), annotation_name_gp = gpar(fontsize = 16))
Figure.S5.Heatmap.Annotation.Pseudotime.Gene = rowAnnotation(Pseudotime.Gene = anno_mark(at = 1:20, labels = rownames(Pseudotime.data[["Heatmap"]])[1:20], labels_gp = gpar(fontsize = 20), padding = unit(2, "mm")), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S5.Pseudotime.Heatmap <- Heatmap(pmin(Pseudotime.data[["Heatmap"]], 15), column_order = Heatmap.Annotation.Pseudotime.Order, clustering_distance_rows = "spearman", show_column_names = FALSE, show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE, width = unit(8, "in"), height = unit(8, "in"), col = Pseudotime.col_fun, top_annotation = Figure.S5.Heatmap.Annotation.Pseudotime.Cell.Types, right_annotation = Figure.S5.Heatmap.Annotation.Pseudotime.Gene)
Figure.S5.Pseudotime.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S5B.", "mLu.combined.CD206hiIMs.Monocle3_Pseudotime", ".Heatmap.pdf")
Figure.S5.Pseudotime.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S5D.", "mLu.combined.CD206loIMs.Monocle3_Pseudotime", ".Heatmap.pdf")
pdf(Figure.S5.Pseudotime.Heatmap.directory, width = 12, height = 9)
print(Figure.S5.Pseudotime.Heatmap)
dev.off()

##Figure.S8##################################################################################################################################################################################

#Figure.S8A.mLu.IM.DotPlot
Chemokine.DotPlot <- DotPlot(mLu.combined, idents = paste0("IMck", 0:10), group.by = "IMs.Chemokine.Subtypes", assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8A.mLu.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 3.8)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8B.mTumor.IM.DotPlot
Idents(mTumor.combined.IMs) <- mTumor.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(mTumor.combined.IMs, idents = paste0("IMck", 0:10), assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8B.mTumor.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 3.8)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8C.mPL.IM.DotPlot
Idents(mPL.combined.IMs) <- mPL.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(mPL.combined.IMs, assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8C.mPL.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 3.36)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8D.mSk.IM.DotPlot
Idents(mSk.combined.IMs) <- mSk.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(mSk.combined.IMs, assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8D.mSk.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 3.36)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8E.mBAL.IM.DotPlot
Idents(mBAL.combined.IMs) <- mBAL.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(mBAL.combined.IMs, idents = paste0("IMck", 0:7), assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8E.mBAL.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 2.7)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8F.hSk.IM.DotPlot
Idents(hSk.combined.IMs) <- hSk.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(hSk.combined.IMs, assay = "SCT", features = IMs.Chemokine.List.Extend.Hs, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8F.hSk.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 13, height = 2.74)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S8G.hBAL.IM.DotPlot
Idents(hBAL.combined.IMs) <- hBAL.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(hBAL.combined.IMs, assay = "SCT", features = IMs.Chemokine.List.Extend.Hs, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S8G.hBAL.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 13, height = 2.74)
print(Chemokine.DotPlot)
dev.off() 

##Figure.S9##################################################################################################################################################################################

#Figure.S9A.mHr.IM.DotPlot
Idents(mHr.combined.IMs) <- mHr.combined.IMs$Chemokine.Cell.Types
Chemokine.DotPlot <- DotPlot(mHr.combined.IMs, assay = "SCT", features = IMs.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S9A.mHr.IM.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.55, height = 4.02)
print(Chemokine.DotPlot)
dev.off()
#Figure.S9B.DimPlot
mHr.combined.IMs.UAMP_1.length <- max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]) - min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])
mHr.combined.IMs.UAMP_2.length <- max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]) - min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])
Figure.S9.DimPlot <- DimPlot(mHr.combined.IMs, group.by = "Chemokine.Cell.Types", cols = c(names(meta.data.list$IMs.Chemokine.Subtypes)[c(1,2,2,3,4,4,5,6,7,8,10)]), raster = FALSE, shuffle = TRUE, seed = 1, pt.size = c(0.5, 2, 2, 0.5)[2]) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(c(min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mHr.combined.IMs.UAMP_1.length), max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mHr.combined.IMs.UAMP_2.length), max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
Figure.S9.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S9B", ".DimPlot.pdf")
pdf(Figure.S9.DimPlot.directory, width = 15, height = 15*0.92)
print(Figure.S9.DimPlot)
dev.off() 
#Figure.S9C.FeaturePlot
#Figure.S9D.FeaturePlot
Figure.S9.FeaturePlot.Features <- IMs.Chemokine.List
for (i in 1:length(Figure.S9.FeaturePlot.Features)){
  Figure.S9.FeaturePlot <- FeaturePlot(mHr.combined.IMs, features = Figure.S9.FeaturePlot.Features[i], cols = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000"), raster = FALSE, order = TRUE, pt.size = 3) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(c(min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"])-(0.075*mHr.combined.IMs.UAMP_1.length), max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_1"]))) + ylim(c(min(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"])-(0.075*mHr.combined.IMs.UAMP_2.length), max(mHr.combined.IMs@reductions$umap@cell.embeddings[, "UMAP_2"]))) + coord_cartesian(expand = FALSE)
  Figure.S9.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S9", ifelse(Figure.S9.FeaturePlot.Features[i] %in% c("Pf4", "Cxcl16", "Cklf"), "C", "D"), ".", Figure.S9.FeaturePlot.Features[i], ".FeaturePlot.pdf")
  pdf(Figure.S9.FeaturePlot.directory, width = 15, height = 15*0.92)
  print(Figure.S9.FeaturePlot)
  dev.off() 
}

#Figure.S9E.Heatmap
mLu.Heatmap.Data <- Heatmap.Data # from Figure.3F
Heatmap.Features <- rownames(mLu.Heatmap.Data)
# mTumor.combined.IMs
Idents(mTumor.combined.IMs) <- "Chemokine.Cell.Types"
mTumor.Heatmap.Data <- AverageExpression(mTumor.combined.IMs, assays = "SCT", slot = "data", group.by = "Chemokine.Cell.Types")
mTumor.Heatmap.Data <- mTumor.Heatmap.Data[["SCT"]][Heatmap.Features, ]
mTumor.Heatmap.Data <- t(scale(t(mTumor.Heatmap.Data)))
# mHr.combined.IMs
Idents(mHr.combined.IMs) <- "Chemokine.Cell.Types"
mHr.Heatmap.Data <- AverageExpression(mHr.combined.IMs, assays = "SCT", slot = "data", group.by = "Chemokine.Cell.Types")
mHr.Heatmap.Data <- mHr.Heatmap.Data[["SCT"]][Heatmap.Features, ]
mHr.Heatmap.Data <- t(scale(t(mHr.Heatmap.Data)))
# mSk.combined.IMs
Idents(mSk.combined.IMs) <- "Chemokine.Cell.Types"
mSk.Heatmap.Data <- AverageExpression(mSk.combined.IMs, assays = "SCT", slot = "data", group.by = "Chemokine.Cell.Types")
mSk.Heatmap.Data <- mSk.Heatmap.Data[Heatmap.Features, ]
mSk.Heatmap.Data <- t(scale(t(mSk.Heatmap.Data)))
# hSk.combined.IMs
Heatmap.Features <- mapvalues(toupper(Heatmap.Features), from = toupper(c("Irg1", "Sepp1", "Fcna", "Cd209f", "Cd209g", "Cd209d")), to = c("ACOD1", "SELENOP", "FCN1", "CD209", "CD209", "CD209")) # Human homologs
hSk.Heatmap.Data <- AverageExpression(hSk.combined.IMs, assays = "SCT", slot = "data", group.by = "Chemokine.Cell.Types")
Idents(hSk.combined.IMs) <- "Chemokine.Cell.Types"
hSk.Heatmap.Data <- hSk.Heatmap.Data[Heatmap.Features, ]
hSk.Heatmap.Data <- t(scale(t(hSk.Heatmap.Data)))
# Heatmap.Data
Heatmap.Data.Vector <- c("mTumor", "mHr", "mSk", "hSk")
Heatmap.Data <- list(mTumor.Heatmap.Data, mHr.Heatmap.Data, mSk.Heatmap.Data, hSk.Heatmap.Data)
for (i in 1:length(Heatmap.Data.Vector)) {
  colnames(Heatmap.Data[[i]]) <- paste0(Heatmap.Data.Vector[i], ".", colnames(Heatmap.Data[[i]]))
}
for (i in 1:length(Heatmap.Data.Vector)) {
  rownames(Heatmap.Data[[i]]) <- paste0(rownames(mLu.Heatmap.Data), " (", rownames(Heatmap.Data[[4]]), ")")
}
Heatmap.Data <- do.call(cbind, Heatmap.Data)
Heatmap.Data <- Heatmap.Data[, paste0(Heatmap.Data.Vector, ".IMck", rep(0:9, each = 4))]
Heatmap.Data[grep("CD209", rownames(Heatmap.Data)), grep("hSk", colnames(Heatmap.Data))] <- rep(Heatmap.Data[grep("Cd209d", rownames(Heatmap.Data)), grep("hSk", colnames(Heatmap.Data))], each = 3) # Same for all CD209 in hSk
top_annotation <- HeatmapAnnotation(IMs.Chemokine.Subtypes = anno_text(rep("", each = 40), location = 0.5, just = "center", rot = 0, height = unit(4, "mm"),
                                                                       gp = gpar(fill = rep(names(meta.data.list$IMs.Chemokine.Subtypes), each = 4), col = "white", fontsize = 10, fontface = "bold")),
                                    Dataset.Subtypes = anno_text(rep("", each = 40), location = 0.5, just = "center", rot = 0, height = unit(2, "mm"),
                                                                 gp = gpar(fill = rep(paste0("grey", round(seq(from = 90, to = 0, length.out = 4))), 4), col = "white", fontsize = 10, fontface = "bold")))
right_annotation.Colors <- setNames(names(meta.data.list$IMs.Chemokine.Subtypes), as.vector(meta.data.list$IMs.Chemokine.Subtypes))
right_annotation <- rowAnnotation(IMs.Chemokine.Subtypes = rep(as.vector((meta.data.list$IMs.Chemokine.Subtypes)), each = 5),
                                  col = list(IMs.Chemokine.Subtypes = right_annotation.Colors), simple_anno_size = unit(2, "mm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.S9E.Heatmap <- Heatmap(Heatmap.Data, show_heatmap_legend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = "italic"), top_annotation = top_annotation, right_annotation = right_annotation)
Figure.S9E.Heatmap.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S9E", ".Heatmap.pdf")
pdf(Figure.S9E.Heatmap.directory, width = 7.5, height = 8)
print(Figure.S9E.Heatmap)
dev.off() 

#Figure.S9F.BarGraph
# DEG.Vector
Figure.S9F.BarGraph.Data <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/Figure.Datasets/Figure.S9F.BarGraph.Data.rds")
Figure.S9F.BarGraph.Data <- melt(Figure.S9F.BarGraph.Data)
colnames(Figure.S9F.BarGraph.Data) <- c("Overlap", "IMck", "Number")
Figure.S9F.BarGraph.Data$IMck <- factor(Figure.S9F.BarGraph.Data$IMck, levels = paste0("IMck", 9:1))
Figure.S9F.BarGraph <- ggplot(data = Figure.S9F.BarGraph.Data, aes(x = Number, y = IMck, fill = Overlap, width = 0.55)) + geom_bar(stat = "identity") + scale_fill_manual(values = paste0("grey", round(seq(from = 0, to = 90, length.out = 5)))) + 
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 15), axis.ticks.y = element_blank(), axis.line.x = element_line(linewidth = 0.5), axis.ticks.x = element_line(color = "black", size = 1), axis.ticks.length = unit(.25, "cm"), 
        panel.background = element_blank(), panel.grid.major.x = element_line(color = "grey80", size = 1, linetype = 2), panel.grid.minor.x = element_line(color = "grey80", size = 0.5, linetype = 2), panel.grid.major.y = element_blank(), panel.ontop = FALSE, legend.position = "none") + scale_x_reverse(position = "top", expand = c(0, 0))
Figure.S9F.BarGraph.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S9F", ".BarGraph.pdf")
pdf(Figure.S9F.BarGraph.directory, width = 3.5, height = 8)
print(Figure.S9F.BarGraph)
dev.off() 

# Figure.9G.DonutChart
Figure.9.DonutChart.Data <- cbind(as.data.frame(Figure.S9F.BarGraph.Data), "Overlap" = rownames(Figure.S9F.BarGraph.Data))
NaturalLevel <- function(x) factor(x, levels = unique(x))
Figure.9.DonutChart.Data$Overlap <- NaturalLevel(Figure.9.DonutChart.Data$Overlap)
reticulate::py_run_string("import sys")
for (i in 1:(length(colnames(Figure.9.DonutChart.Data))-1)){
  Figure.9.DonutChart <- Figure.9.DonutChart.Data %>% plot_ly(labels = ~Overlap, values = ~Figure.9.DonutChart.Data[[colnames(Figure.9.DonutChart.Data)[i]]], marker = list(colors = col2hex(paste0("grey", round(seq(from = 0, to = 90, length.out = 5))))), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 100))
  Figure.9.DonutChart <- Figure.9.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  Figure.9.DonutChart.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S9G", ".DonutChart.", colnames(Figure.9.DonutChart.Data)[i], ".pdf")
  save_image(Figure.9.DonutChart, Figure.9.DonutChart.directory, scale = 5)
}

##Figure.S10##################################################################################################################################################################################

#Figure.S10A.mPL.Macs.LPMs.DotPlot
Idents(mPL.combined.LPMs) <- mPL.combined.LPMs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("Cxcl1", "Cxcl2", "Cxcl13", "Cxcl14", "Ccl2", "Ccl3", "Cxcl12", "Cxcl16", "Ccl24", "Ccl6", "Ccl9", "Pf4", "Cklf", "Ccl25")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(mPL.combined.LPMs, assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10A.mPL.Macs.LPMs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 8.15, height = 4.7)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S10B.mPL.Macs.SPMs.DotPlot
Idents(mPL.combined.SPMs) <- mPL.combined.SPMs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("Ccl25", "Cxcl16", "Cxcl10", "Ccl2", "Ccl7", "Ccl27a", "Ccl5", "Ccl4", "Xcl1", "Cxcl1", "Cxcl2", "Cxcl12", "Cxcl13", "Pf4", "Ccl24", "Ccl6", "Ccl9", "Cklf", "Cxcl14")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(mPL.combined.SPMs, idents = paste0("SPMck", 0:9), assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10B.mPL.Macs.SPMs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 8.05, height = 3.82)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S10C.mSk.Macs.DotPlot
Idents(mSk.combined.LCs) <- mSk.combined.LCs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("Ccl7", "2610528A11Rik", "Pf4", "Ccl24", "Ccl4", "Ccl2", "Ccl25", "Ccl6", "Ccl9", "Ccl3", "Ccl27a", "Cxcl10", "Ccl22", "Ccl8", "Ccl5", "Cxcl3", "Cxcl2", "Cxcl1", "Cklf", "Cxcl16")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(mSk.combined.LCs, assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10C.mSk.Macs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 7.88, height = 3.6)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S10D.mBAL.Macs.DotPlot
Idents(mBAL.combined.AMs) <- mBAL.combined.AMs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("Ccl12", "Cxcl10", "Cxcl9", "Cxcl3", "Cxcl2", "Cxcl1", "Ccl3", "Ccl9", "Cxcl14", "Ccl4", "Pf4", "Cxcl16", "Ccl25", "Ccl6", "Cklf", "Cx3cl1")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(mBAL.combined.AMs, idents = paste0("AMck", 0:9), assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10D.mBAL.Macs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 7.94, height = 3.82)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S10E.hSk.Macs.DotPlot
Idents(hSk.combined.LCs) <- hSk.combined.LCs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("CXCL14", "CCL22", "CCL20", "CCL4", "CCL3L1", "CCL3", "CCL18", "CCL13", "CCL2", "CXCL12", "CCL14", "CCL8", "CXCL5", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CKLF", "CXCL16")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands.Hs, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(hSk.combined.LCs, assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10E.hSk.Macs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 10.56, height = 3.17)
print(Chemokine.DotPlot)
dev.off() 
#Figure.S10F.hBAL.Macs.DotPlot
Idents(hBAL.combined.AMs) <- hBAL.combined.AMs$Chemokine.Cell.Types
Mac.Chemokine.List <- c("CCL2", "CCL5", "CCL15", "CXCL9", "CXCL10", "CXCL11", "CCL24", "PPBP", "CXCL5", "CXCL8", "CCL3L1", "CCL3", "CCL4", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CCL23", "CCL18", "CXCL16", "CKLF")
Mac.Chemokine.List.Extend <- list("a" = Mac.Chemokine.List, "b" = setdiff(Chemokine.Ligands.Hs, Mac.Chemokine.List))
Chemokine.DotPlot <- DotPlot(hBAL.combined.AMs, assay = "SCT", features = Mac.Chemokine.List.Extend, dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 15), strip.text = element_text(size = 13), legend.position = "none")
Chemokine.DotPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Chemokine.DotPlot/Figure.S10F.hBAL.Macs.DotPlot.pdf")
pdf(Chemokine.DotPlot.directory, width = 11.15, height = 4.27)
print(Chemokine.DotPlot)
dev.off() 

##Figure.S13##################################################################################################################################################################################

GSE147668.df <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/GSE147668/GSE147668.df.rds")
CRA004586.df <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/CRA004586/CRA004586.df.rds")
GSE149563.df <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/GSE149563/GSE149563.df.rds")
E.MTAB.10026.df <- readRDS("/Volumes/mLu/Analysis/scRNA_seq/Datasets/scRNA_seq/E.MTAB.10026/E.MTAB.10026.df.rds")
Figure.S13.Data <- list(GSE147668.df, CRA004586.df, GSE149563.df, E.MTAB.10026.df)
Figure.S13.Study <- c("GSE147668", "CRA004586", "GSE149563", "E.MTAB.10026")
Figure.S13.Background <- c("#FAFAFA", "white", "#FAFAFA", "white")
Figure.S13.DimPlot.Color <- list(c("#FFC72F", "#1A00F8", "#0955F9", "#CD00F9", "#FF00C2", "#6200F8", "#FF5C21", "#00FF31", "#00FFC6", "#00FF66", "#FF011C", "#FF005D", "#BFFF36", "#00C4FB", "#41FF32"), 
                                 c("#4700F8", "#DE00F9", "#00D5FC", "#00D5FC", "#FF00D2", "#00FF8D", "#FF011C", "#79FF33", "#8E00F9", "#1A00F8", "#FF0045", "#00FF50", "#FF0087", "#00FF32", "#FFD832", "#00FFD6", "#1337F8", "#FF8826", "#D1FF37", "#00FF31", "#0084FA"), 
                                 c("#99C7DD", "#5B9DC5", "#2275AD", "#3E8D9B", "#7DBF8C", "#87CD71", "#4BAE4B", "#2E9636", "#90925D", "#F58F8B", "#F55B5F", "#E92931", "#ED4034", "#F28253", "#FFAC5C", "#FF8E37", "#FF7726", "#E58D60", "#CAA3BC", "#9E7AB5", "#754798", "#845F8D", "#C9B991", "#F3E987", "#CF9956", "#AC4F2B"), 
                                 c("#87D1D4", "#A398D5", "#89D79B", "#D8E19A", "#D76F42", "#FFAD88", "#FFA2BF", "#EC7CDC", "#264785", "#CC494C", "#00A64E", "#3FB8BD", "#428FD4", "#F49195", "#C9D899", "#C069DD", "#F1282B", "#006193"))
# Figure.S13C.DimPlot
for (i in 1:length(Figure.S13.Data)) {
  Figure.S13.DimPlot <- ggplot(Figure.S13.Data[[i]], aes(x = Reduction_1, y = Reduction_2, color = Cell.Types)) + geom_point(size = c(4, 2, 2, 1)[i]) + scale_color_manual(values = alpha(Figure.S13.DimPlot.Color[[i]], c(0.75, 0.75, 1, 0.75)[i])) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range1.01(Figure.S13.Data[[i]]$Reduction_1)) + ylim(range1.01(Figure.S13.Data[[i]]$Reduction_2)) + coord_cartesian(expand = FALSE) + theme(panel.background = element_rect(fill = Figure.S13.Background[i], color = Figure.S13.Background[i]))
  Figure.S13.DimPlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Figure.S13C.", Figure.S13.Study[i], ".DimPlot.png")
  png(Figure.S13.DimPlot.directory, width = 1800, height = 1500)
  print(Figure.S13.DimPlot)
  dev.off() 
}
# Figure.S13C.FeaturePlot
for (i in 1:length(Figure.S13.Data)) {
  if (i %in% 1:3) {FeaturePlot.Features.Vector <- c("Pf4", "Cx3cr1", "C5ar1", "C1qc")} else {FeaturePlot.Features.Vector <- c("PF4", "CX3CR1", "C5AR1", "C1QC")}
  for (j in 1:length(FeaturePlot.Features.Vector)){
    Figure.S13.Data[[i]] <- Figure.S13.Data[[i]][order(Figure.S13.Data[[i]][, FeaturePlot.Features.Vector[j]]), ]
    Figure.S13.FeaturePlot <- ggplot(Figure.S13.Data[[i]], aes(x = Reduction_1, y = Reduction_2, color = Figure.S13.Data[[i]][, FeaturePlot.Features.Vector[j]])) + geom_point(size = c(5, 4, 3, 1.5)[i]) + scale_colour_gradientn(colours = c("#E5E5E5", "#F6D746", "#E55C30", "#84206B", "#140B34", "#000000")) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
      xlim(range1.01(Figure.S13.Data[[i]]$Reduction_1)) + ylim(range1.01(Figure.S13.Data[[i]]$Reduction_2)) + coord_cartesian(expand = FALSE) + theme(panel.background = element_rect(fill = Figure.S13.Background[i], color = Figure.S13.Background[i]))
    Figure.S13.FeaturePlot.directory <- paste0("/Volumes/mLu/Analysis/scRNA_seq/Figures/Revision.Figure.S13C.", Figure.S13.Study[i], ".", FeaturePlot.Features.Vector[j], ".FeaturePlot.png")
    png(Figure.S13.FeaturePlot.directory, width = 1800, height = 1500)
    print(Figure.S13.FeaturePlot)
    dev.off() 
  }
}


