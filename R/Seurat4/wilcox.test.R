library(dplyr)
library(tidyr)
library(magrittr)
library(ggpubr)
# cell percentage =======
cell_number <- read.csv("output/partition_Cell_counts_BB.csv")
total_cell_number= rowSums(cell_number[,3:11])
cell_number  = cell_number[,c("patient","treatment","MDSC")]
cell_number$MDSC = cell_number$MDSC/total_cell_number*100
cell_number = pivot_wider(cell_number,names_from = "treatment", values_from = "MDSC")
cell_number$total_patent_cell_number = rowSums(cell_number[,c("baseline","ibrutinib")])
cell_number = filter(cell_number, total_patent_cell_number > 10)

wilcox.test(x = cell_number$baseline,
            y = cell_number$ibrutinib,paired = TRUE)


# cell number =======
cell_number <- read.csv("output/partition_Cell_counts_BB.csv")
cell_number  = cell_number[,c("patient","treatment","MDSC","Mono")]
cell_number = pivot_wider(cell_number,names_from = "treatment", values_from = c("MDSC","Mono"))
cell_number %<>% tibble::column_to_rownames("patient")
cell_number$total_patent_cell_number = rowSums(cell_number)
cell_number = filter(cell_number, total_patent_cell_number > 100)

wilcox.test(x = cell_number$baseline,
            y = cell_number$ibrutinib,paired = FALSE)

colSums(cell_number[,2:ncol(cell_number)],na.rm = T)
log2(972/4430) = -2.188278
log2(19264/10380) = 0.8921009

# cell log =======
cell_number <- read.csv("output/partition_Cell_counts_BB.csv")
cell_number  = cell_number[,-grep("doublet",colnames(cell_number))]
cell_number = pivot_wider(cell_number,names_from = "treatment", values_from = colnames(cell_number)[3:11])
cell_number %<>% tibble::column_to_rownames("patient")
Colsum = colSums(cell_number,na.rm = T)
names(Colsum) = colnames(cell_number)
df = data.frame("cell.type" = c("MDSC","Mono"),
                "log2FC" = c(log2(Colsum["MDSC_ibrutinib"]/Colsum["MDSC_baseline"]),
                             log2(Colsum["Mono_ibrutinib"]/Colsum["Mono_baseline"])))
ggbarplot(df, x = "cell.type", y = "log2FC",fill = c("#00AFBB", "#FC4E07"),
          ylab = "log2(C1D+1/C1D-7")+
    font("xy.text", size = 20)+
    font("xlab", size = 24)+
    font("ylab", size = 24)



# ======== subcluster =========
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()


df_samples <- df_samples %>%
    reorder_levels(response, order = c("PR", "SD", "PD"))
df_samples %>% 
    group_by(response) %>%
    get_summary_stats(cluster, type = "common")




res.kruskal <- df_samples %>% kruskal_test(response ~ cluster)
res.kruskal

jpeg(paste0(path,"Subcluster.jpeg"), units="in", width=8, height=6,res=600)
ggboxplot(df_samples, x = "response", y = "cluster",add = "dotplot",
          color = "timepoint") + scale_y_continuous(limits = c(1, 2),
                                                    breaks = c(1, 2))
dev.off()
