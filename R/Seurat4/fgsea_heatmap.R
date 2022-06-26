########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(progress)
library(tidyr)
library(kableExtra)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
#=========================

hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.4.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

#===============
object = readRDS("data/OSU_SCT_20210821.rds")
meta_data = readRDS("output/20220309/meta_data.rds")
if(all(colnames(object) == rownames(meta_data))) print("!") 
object %<>% subset(subset = Doublets %in% "Singlet")
object$sample %<>% factor(levels = df_samples$sample)
object$FTH1_lvl = plyr::mapvalues(object$FTH1_lvl,
                                  from = c("FTH1 high","FTH1 low"),
                                  to = c("subCluster1","subCluster2"))
#object$cell.types %<>% gsub("Monocytes|MDSCs","Monocytes+MDSCs",.)
object$cell.types %<>% gsub(":",",",.)
object %<>% subset(subset = cell.types %in% c("B-cells","Monocytes","MDSCs","NK cells","T-cells,CD4+",
                                              "T-cells,CD8+"))
object$cell.type_subCluster = paste0(object$cell.types,"_", object$FTH1_lvl)

object$sample_cell.type_subCluster = paste0(object$sample,"_",object$response,"_", object$cell.type_subCluster)

cell_type_list = list("B-cells",c("Monocytes","MDSCs"),"NK cells","T-cells,CD4+",
                      "T-cells,CD8+")
Fgsea_list <- list()
for(cell_type in cell_type_list){
    sub_object = subset(object, subset = cell.types %in% cell_type)
    sub_EXP = AverageExpression(sub_object,assays = "SCT", group.by = "sample_cell.type_subCluster")
    sub_EXP = sub_EXP$SCT
    df = table(sub_object$sample_cell.type_subCluster) %>% as.data.frame()
    df$sample = gsub("_PR_.*|_SD_.*|_PD_.*","",df$Var1)
    df_list <- split(df, f = df$sample)
    for(i in seq_along(df_list)) {
        if(nrow(df_list[[i]]) > 1) {
            df_list[[i]] = df_list[[i]][order(df_list[[i]]$Freq,decreasing = T)[1],]
        }
    }
    select_sample = bind_rows(df_list)
    sub_EXP = sub_EXP[,select_sample$Var1]
    stats = avg_log2FC(sub_EXP) %>% as.data.frame() %>%
        tibble::rownames_to_column("gene") %>% pivot_longer(-gene,
                                                               names_to = "cluster",
                                                               values_to = "avg_log2FC")
    Fgsea_list[[i]] <- FgseaDotPlot(stats=stats, pathways=hallmark,Rowv = T,Colv = T,
                                    size = " -log10(padj)",
                                    #order.xaxis = groups[[i]],
                              title = cell_type,
                              plot.title = element_text(hjust = 0.5,size = 15),
                              axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
                              axis.text.y = element_text(size = 10),
                              text = element_text(family = "Avenir"),
                              height=10,
                              width = 12,
                              file.name = paste0("Dotplot_",cell_type,"_0.25_0.5"),
                              do.return = T)
}
openxlsx::write.xlsx(Fgsea_list,
                     file =  paste0(path,"gsea_summary.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding")

avg_log2FC <- function(exp){
    log2_exp = log2(exp+1)
    row_Sums <- rowSums(log2_exp)
    Ncol = ncol(log2_exp)
    log2_fc = apply(log2_exp, 2, function(col) col- (row_Sums - col)/(Ncol-1))
    return(log2_fc)
}

