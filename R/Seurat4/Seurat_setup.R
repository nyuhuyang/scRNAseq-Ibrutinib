########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","harmony","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples$patient %<>% as.character()
df_samples$`cite-seq` %<>% as.character()
df_samples %>% kable() %>% kableExtra::kable_styling()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/OSU_20210821.rds")
meta.data = object@meta.data
for(i in 1:length(df_samples$sample.id)){
    cells <- meta.data$orig.ident %in% df_samples$sample.id[i]
    print(df_samples$sample.id[i])
    print(table(cells))
    meta.data[cells,"sample"] = df_samples$sample[i]
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"timepoint"] = df_samples$timepoint[i]
    meta.data[cells,"treatment"] = df_samples$treatment[i]
    meta.data[cells,"response"] = df_samples$response[i]
    meta.data[cells,"disease"] = df_samples$disease[i]
    meta.data[cells,"cite-seq"] = df_samples$`cite-seq`[i]
    #meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    #meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    #meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$sample.id)
meta.data$sample %<>% factor(levels = df_samples$sample)
meta.data$sample %<>% factor(levels = df_samples$sample)

table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% pblapply(SCTransform,method = "glmGamPoi",vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(object.list = object_list)

options(future.globals.maxSize= object.size(object_list)*1.5)
object_list %<>% pblapply(FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
})
anchors <- FindIntegrationAnchors(object.list = object_list,reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:50)
remove(object_list);GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors,normalization.method = "SCT", dims = 1:50)
remove(anchors);GC()
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/OSU_20210821.rds")
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
    
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard workflow for visualization and clustering
# Run the standard cca workflow for umap & tsne visualization
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
npcs = 100

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object[["cca.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "ccaUMAP_", assay = DefaultAssay(object))
colnames(object[["cca.umap"]]@cell.embeddings) %<>% paste0("cca-",.)

object[["cca.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "cca-SNE_", assay = DefaultAssay(object))
colnames(object[["cca.tsne"]]@cell.embeddings) %<>% paste0("cca-",.)

saveRDS(object, file = "data/OSU_20210821.rds")


#======1.7 UMAP from raw pca =========================
format(object.size(object),unit = "GB")
DefaultAssay(object) = "RNA"
object[['SCT']] = NULL
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

saveRDS(object, file = "data/OSU_20210821.rds")

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"

npcs = 100
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
colnames(object[["harmony.umap"]]@cell.embeddings) %<>% paste0("harmony-",.)

object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                                 key = "harmonytSNE_", assay = DefaultAssay(object))
colnames(object[["harmony.tsne"]]@cell.embeddings) %<>% paste0("harmony-",.)


npcs = 100
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
object %<>% FindClusters(resolution = 0.8)
object %<>% FindClusters(resolution = 0.02)
UMAPPlot.1(object,group.by = "SCT_snn_res.0.02",do.print = T,label = T,label.repel = T)


saveRDS(object, file = "data/OSU_20210821.rds")


#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/OSU_SCT_20210821.rds")


