########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","eulerr","openxlsx","kableExtra",
                   "magrittr","readxl"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# =======================================
df_samples <- readxl::read_excel("doc/20211020_scRNAseq_info.xlsx")

xlsx_names = list.files("output/20211027",pattern = "Ibrutinib vs Baseline.*xlsx",full.names = T)

DEG <- pbapply::pblapply(xlsx_names, function(xlsx){
        tmp <- readxl::read_excel(xlsx)
        type = sub(".*_","",xlsx) %>% sub("\\.xlsx","",.)
        tmp$type = type
        tmp
}) %>% bind_rows()

table(DEG$p_val_adj <0.05)
table(abs(DEG$avg_log2FC) >0.5 )
table(abs(DEG$avg_log2FC) >0.25 )
table(abs(DEG$avg_log2FC) >0.1 )

DEG %<>% as.data.frame()

#============================
for(fc_value in c(0.1,0.25)){
        pos_genes <- eulerr(DEG, 
                            group.by = "type", shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                            cut_off = "avg_log2FC", cut_off_value = fc_value,do.print = T, return.raw = T,
                            save.path = path, file.name =paste0("Venn_treatment_log2fc_",fc_value,".jpeg"))
        euler_df <- eulerr::euler(pos_genes,shape = "circle")
        pos_genes_list <- as.list(euler_df$original.values)
        names(pos_genes_list) %<>% paste(":",pos_genes_list)
        id <- eulerr:::bit_indexr(8)
        for (i in nrow(id):1) {
                pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                        setdiff(Reduce(union, pos_genes[!id[i,]]))
        }
        pos_genes_df <- list2df(pos_genes_list)
        pos_genes_df = pos_genes_df[,sapply(pos_genes_list,length) != 0]
        write.xlsx(pos_genes_df, asTable = F,
                   file = paste0(path,"Venn_treatment_log2fc_",fc_value,".xlsx"),
                   borders = "surrounding")
        print(fc_value)
}


cell_types = unique(DEG$type)
for(cell in c("TASC","BC","IC", "S1")){
        eulerr(DEG, key = grep(cell,unique(DEG$type),value = T),
               group.by = "type", shape =  "circle",
               cut_off = "avg_log2FC", cut_off_value = 0.5,do.print = T,
               save.path = path, file.name =paste0("Venn_vivo_vs_vitro_",cell,"_log2fc_1.jpeg"))
}

DEG$abs_log2FC = abs(DEG$avg_log2FC)
DEG %>% group_by(type) %>% top_n(8, abs_log2FC) %>% kable() %>% kableExtra::kable_styling()
deg <- DEG %>% group_by(type) %>% top_n(8, abs_log2FC)
deg <- deg[!duplicated(deg$gene),]
write.csv(deg, paste0(path,"Venn_treatment_top42.csv"))
