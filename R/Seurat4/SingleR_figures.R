# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/OSU_SCT_20210821_singleR_pred.rds")
pred1 = readRDS("output/OSU_SCT_20211117_singleR_pred.rds")
object = readRDS("data/OSU_SCT_20210821.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       "celltype.l3" = pred1$pruned.labels,
                       row.names = rownames(pred1))
table(is.na(singlerDF$label.fine))
table(is.na(singlerDF$celltype.l3))

singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"
singlerDF$celltype.l3[is.na(singlerDF$celltype.l3)]= "unknown"

##############################
# adjust cell label
##############################
# combine cell types
singlerDF$cell.types = singlerDF$label.fine
grep("B-cells",singlerDF$cell.types, value = T) %>% table
singlerDF[grep("B-cells",singlerDF$cell.types),"cell.types"] ="B-cells"
singlerDF[grep("CD4+",singlerDF$cell.types),"cell.types"] ="T-cells:CD4+"
singlerDF[grep("CD8+",singlerDF$cell.types),"cell.types"] ="T-cells:CD8+"
singlerDF$cell.types %<>% gsub("Tregs","T-cells:regs",.)
singlerDF$cell.types %<>% gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",.)
singlerDF$cell.types %<>% gsub("DC|Macrophages|Macrophages M1","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Eosinophils|Megakaryocytes","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Astrocytes|Adipocytes|Fibroblasts|mv Endothelial cells|Keratinocytes|Epithelial cells","Nonhematopoietic cells",.)

# ==========================
path = "../seurat_resources/azimuth/PBMC/"
meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
meta.data = meta.data[!duplicated(meta.data$celltype.l3),]
singlerDF$celltype.l2 = plyr::mapvalues(singlerDF$celltype.l3,from  = meta.data$celltype.l3,
                                        to = meta.data$celltype.l2)
singlerDF$celltype.l1 = plyr::mapvalues(singlerDF$celltype.l3,from  = meta.data$celltype.l3,
                                        to = meta.data$celltype.l1)
##############################
# process color scheme
##############################
singlerDF %<>% cbind(object[["umap"]]@cell.embeddings)
MDSCs <- singlerDF$UMAP_1 < 0 &  singlerDF$UMAP_2 < -6
singlerDF[MDSCs,"label.fine"] = "MDSCs"
singlerDF[MDSCs,"cell.types"] = "MDSCs"

singlerDF$cell.types.colors = singlerDF$cell.types
singlerDF$cell.types.colors %<>% plyr::mapvalues(from = c("B-cells","Erythrocytes","HSC/progenitors",
                                                          "MDSCs", "MCL","Monocytes","NK cells",
                                                          "Nonhematopoietic cells","other Myeloid cells",
                                                          "Plasma cells","T-cells:CD4+","T-cells:CD8+",
                                                          "T-cells:regs","unknown"),
                                                 to = c("#E6AB02", "#ff0000", "#6A3D9A",
                                                        "#1F78B4", "#2055da", "#ADDFEE","#A65628",
                                                        "#B3B3B3","#FDDAEC",
                                                        "#1B9E77","#B3DE69","#F0027F",
                                                        "#7570B3","#F2F2F2"))
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.fine","cell.types", "cell.types.colors")])
object@meta.data %<>% cbind(singlerDF[,c("celltype.l1","celltype.l2", "celltype.l3")])


lapply(c(TSNEPlot.1,UMAPPlot.1), function(fun)
    fun(object = object, label = T, label.repel = T,group.by = "cell.types",
        no.legend = T,
        pt.size = 0.1,label.size = 3,
        do.print = T,do.return = F,
        title ="labeling by blue_encode and TCGA-BLCA RNA-seq"))

saveRDS(object, file = "data/OSU_SCT_20211123.rds")


# by barplot
cell_Freq <- table(object$label.fine) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 30 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()