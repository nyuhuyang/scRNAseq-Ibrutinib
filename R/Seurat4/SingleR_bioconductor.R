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

# ====== load reference =============
blue_encode <- BlueprintEncodeData()
#remove = grepl("CD4|CD8|Tregs|B-cells|Monocytes",blue_encode$label.fine)
#blue_encode = blue_encode[,!remove]

#immue_exp <- DatabaseImmuneCellExpressionData()

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(blue_encode)
))
length(common)
table(blue_encode$label.fine)
system.time(trained <- trainSingleR(ref = blue_encode[common,],
                                    labels=blue_encode$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/OSU_SCT_20210821_singleR_pred.rds")