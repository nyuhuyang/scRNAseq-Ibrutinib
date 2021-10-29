library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
object = readRDS("data/OSU_SCT_20210821.rds")
meta_data$`cite-seq` = meta.data$`cite-seq`
meta_data %<>% cbind(object[["umap"]]@cell.embeddings)
MDSCs <- meta_data$UMAP_1 < 0 &  meta_data$UMAP_2 < -6
meta_data[MDSCs,"label.fine"] = "MDSCs"
meta_data[MDSCs,"cell.types"] = "MDSCs"
meta_data[MDSCs,"cell.types.colors"] = "#1F78B4"

object %<>% FindClusters(resolution = 0.02)
Idents(object) = "SCT_snn_res.0.02"
object %<>% RenameIdents('0' = 'FTH1 low',
                         '1' = 'FTH1 low',
                         '2' = 'FTH1 high',
                         '3' = 'FTH1 low',
                         '4' = 'FTH1 high',
                         '5' = 'FTH1 high',
                         '6' = 'FTH1 low',
                         '7' = 'FTH1 high',
                         '8' = 'FTH1 low',
                         '9' = 'FTH1 low',
                         '10' = 'FTH1 high',
                         '11' = 'FTH1 high',
                         '12' = 'FTH1 low',
                         '13' = 'FTH1 low',
                         '14' = 'FTH1 high',
                         '15' = 'FTH1 high',
                         '16' = 'FTH1 high'
                         )
object[["FTH1_lvl"]] = Idents(object)
object$disease %<>% gsub("melanoma","Melanoma",.)
meta_data = object@meta.data
saveRDS(meta_data,"output/20211022/meta_data.rds")
