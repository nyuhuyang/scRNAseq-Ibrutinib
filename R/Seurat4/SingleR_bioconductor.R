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

# 2. download TCGA BLCA data ==============================
query.exp <- GDCquery(project = "TCGA-BLCA",
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results", 
                      #sample.type = "Primary solid Tumor",
                      legacy = TRUE)

GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
exp
exp = exp[, exp$definition == "Primary solid Tumor"]
assays(exp)[[1]][1:4,1:4]
genes = gsub("\\|.*","",rownames(assays(exp)$raw_count))
raw_count <- assays(exp)$raw_count[!duplicated(genes),]
rownames(raw_count) = gsub("\\|.*","",rownames(raw_count))

libsizes <- colSums(raw_count)
size.factors <- libsizes/mean(libsizes)
log_count <- log2(t(t(raw_count)/size.factors) + 1)

meta.data = data.frame("label.main"= rep("BLCA",ncol(raw_count)),
                       "label.fine" = paste("BLCA",exp$tumor_stage),
                       "label.ont" = colnames(exp))
BLCA <- SummarizedExperiment(list(logcounts=log_count),
                                  colData=DataFrame(meta.data))

# ====== load reference =============
blue_encode <- BlueprintEncodeData()
#remove = grepl("CD4|CD8|Tregs|B-cells|Monocytes",blue_encode$label.fine)
#blue_encode = blue_encode[,!remove]

#immue_exp <- DatabaseImmuneCellExpressionData()

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(BLCA),
                                 rownames(blue_encode)
))
length(common)
combine_ref = do.call("cbind", list(blue_encode[common,],
                                    BLCA[common,]))
table(combine_ref$label.fine)
system.time(trained <- trainSingleR(ref = combine_ref,
                                    labels=combine_ref$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/OSU_SCT_20210821_singleR_pred.rds")