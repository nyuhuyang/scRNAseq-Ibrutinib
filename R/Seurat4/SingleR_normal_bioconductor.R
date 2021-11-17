#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(TCGAbiolinks)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============
object = readRDS("data/OSU_SCT_20210821.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ======= load azimuth PBMC data ==============================
path = "../seurat_resources/azimuth/PBMC/"
counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
table(rownames(meta.data) == colnames(counts))
PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                             colData=DataFrame(meta.data))
rm(counts);GC()

# ====== conbime data =============

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(PBMC)
))
length(common)
table(PBMC$celltype.l3)
system.time(trained <- trainSingleR(ref = PBMC[common,],
                                    labels=PBMC$celltype.l3))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/OSU_SCT_20211117_singleR_pred.rds")