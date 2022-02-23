########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(progress)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("shinyApp/PBMC/util.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

csv_list <- list.files(path = "output/20211101",
                        pattern = "Ibrutinib vs Baseline.*csv",full.names = T)
deg_total <- pbapply::pblapply(csv_list, function(csv){
    tmp <- read.csv(csv,row.names = 1)
    tmp$gene =rownames(tmp)
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    tmp
}) %>% bind_rows()


deg_total$cell.type %<>% plyr::mapvalues(from = c("B-cells",
                                            "MDSCs",
                                            "Monocytes",
                                            "NK cells",
                                            "T-cells:CD4+",
                                            "T-cells:CD8+",
                                            "T-cells:regs"),
                                   to = c("B cells",
                                          "MDSCs",
                                          "monocytes",
                                          "NK cells",
                                          "CD8T",
                                          "CD4T",
                                          "Treg"))

deg_total$cell.type %<>% factor(levels = rev(c("B cells","MDSCs","monocytes","NK cells","CD4T","CD8T","Treg")))

deg_total %<>% filter(p_val < 0.05)
table(abs(deg_total$avg_log2FC) > 0.1)
df_logFC = as.data.frame.matrix(table(deg_total$cell.type,deg_total$avg_log2FC > 0.1))
colnames(df_logFC) = c("log2FC < -0.1", "log2FC > 0.1")
df_logFC$`log2FC < -0.1` = -(df_logFC$`log2FC < -0.1`)
df_logFC %<>% tibble::rownames_to_column("cell.types")
df_logFC %<>% pivot_longer(!cell.types, names_to = "change", values_to = "gene number")
df_logFC = df_logFC[order(df_logFC$change),]

#======================
csv_list <- list.files(path = "output/20211123",
                       pattern = "FC0.25_cell.types_.*csv",full.names = T)
deg_total <- pbapply::pblapply(csv_list, function(csv){
    tmp <- read.csv(csv,row.names = 1)
    tmp$gene =rownames(tmp)
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    tmp
}) %>% bind_rows()

openxlsx::write.xlsx(deg_total, file = paste0(path,"DEGs_cell.types.xlsx"),
                     colNames = TRUE,rowNames = FALSE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
#======================
csv_list <- list.files(path = "output/20211123",
                       pattern = "_FC0.25_cluster_.*csv",full.names = T)
deg_total <- pbapply::pblapply(csv_list, function(csv){
    tmp <- read.csv(csv,row.names = 1)
    tmp$gene =rownames(tmp)
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    tmp
}) %>% bind_rows()

openxlsx::write.xlsx(deg_total, file = paste0(path,"DEGs_Cluster_resolution.0.8.xlsx"),
                     colNames = TRUE,rowNames = FALSE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



jpeg(paste0(path,"gene_number_change.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(df_logFC, x="cell.types", y= "gene number",fill="change",
          palette = c("#A6CEE3","#B15928"),lab.size = 20)+
    facet_wrap(~ change, scales = "free_x") + 
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.spacing.x = unit(0, "mm"))+
    sctheme(base_size = sList[inpfsz], Xang = 90, XjusH = 1)

dev.off()

jpeg(paste0(path,"gene_number_change1.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(df_logFC, x="cell.types", y= "gene number",fill="change",
          palette = c("#A6CEE3"),lab.size = 20)+
    facet_wrap(~ change, scales = "free_x") + 
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),limits = c(-600,0)) +
    theme(panel.spacing.x = unit(0, "mm"))+
    sctheme(base_size = sList[inpfsz], Xang = 90, XjusH = 1)

dev.off()

jpeg(paste0(path,"gene_number_change2.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(df_logFC, x="cell.types", y= "gene number",fill="change",
          palette = c("#B15928"),lab.size = 20)+
    facet_wrap(~ change, scales = "free_x") + 
    coord_flip() +
    scale_y_continuous(limits = c(0,600),breaks = c(1:6*100)) +
    theme(panel.spacing.x = unit(0, "mm"))+
    sctheme(base_size = sList[inpfsz], Xang = 90, XjusH = 1)

dev.off()


# Plot theme
sList = c(12,18,24,30)
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){
    oupTheme = theme(
        text =                         element_text(size = base_size, family = "Helvetica"),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line =     element_line(colour = "black"),
        axis.ticks =    element_line(colour = "black", size = base_size / 20),
        axis.title =    element_text(face = "bold"),
        axis.text =     element_text(size = base_size),
        axis.text.x = element_text(angle = Xang, hjust = XjusH),
        legend.position = "bottom",
        legend.text = element_text(size = base_size),
        legend.key = element_rect(colour = NA, fill = NA)
    )
    if(!XYval){
        oupTheme = oupTheme + theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    return(oupTheme)
}
