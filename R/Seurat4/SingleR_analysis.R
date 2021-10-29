########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# ====== load Seurat =============
load(file = "data/Lorenzo-LS6_20210408_SCT.Rda")
# comapre cell types difference between WT and KO

# by UMAP
UMAPPlot.1(object, group.by="label.fine",split.by = "conditions",
           pt.size = 0.1,label = T,
           label.repel = T,alpha = 0.8,
           no.legend = T,label.size = 3, repel = T, title = "Cell types differece between WT and KO",
           do.print = T, do.return = F)

Idents(object) = "label.fine"
# by t test
cell_Freq <- table(object$label.fine, object$orig.ident) %>% as.data.frame.matrix()

df = data.frame("orig.ident" = c("LS2","LS3","LS4","LS5","LS6","LS7"),
               "conditions" = c("KO", "WT", "KO", "WT", "WT", "KO"))
df$conditions %<>% factor(levels = c("WT","KO"))
df = df[order(df$conditions),]
cell_Freq = cell_Freq[,df$orig.ident]
cell_Freq[,"ttest"] = NA
cell_Freq[,"wilcox"] = NA
cell_Freq %<>% as.matrix()
n = ncol(cell_Freq)/2


for(i in 1:nrow(cell_Freq)) {
    ttest = t.test(x = cell_Freq[i,1:n], y = cell_Freq[i,(n+1):(2*n)],paired = F)
    wilcox = wilcox.test(x = cell_Freq[i,1:n],y = cell_Freq[i,(n+1):(2*n)],paired = F)
    
    cell_Freq[i,"ttest"] = ttest$p.value
    cell_Freq[i,"wilcox"] = wilcox$p.value
    svMisc::progress(i/nrow(cell_Freq)*100)
}
cell_Freq %<>% as.data.frame()
cell_Freq[,"ttest.adj"] = p.adjust(cell_Freq[,"ttest"], method = "BH")
cell_Freq[,"wilcox.adj"] = p.adjust(cell_Freq[,"wilcox"], method = "BH")

write.csv(cell_Freq,paste0(path,"compare_cell_types.csv"))


# by fisher exact test
cell_Freq_split <- table(object$label.fine, object$conditions) %>% as.data.frame.matrix()

cell_Freq_split$p_value = NULL
ColSum <- colSums(cell_Freq_split)

for(i in 1:nrow(cell_Freq_split)){
    conting <- rbind(cell_Freq_split[i,1:2],ColSum-cell_Freq_split[i,1:2])
    FISH <- fisher.test(conting,conf.int = T)
    cell_Freq_split[i,"p_value"] = FISH$p.value
}

cell_Freq_split$p_val_adj = p.adjust(p = cell_Freq_split$p_value, method = "bonferroni", 
                                     n = nrow(cell_Freq_split))
write.csv(cell_Freq_split,paste0(path,"cell_Freq_split.csv"))

Rshiny_path <- "Rshiny/Lorenzo-LS-9829/"
samples <- c("All_samples","WT", "KO","LS2","LS3","LS4","LS5","LS6","LS7")

PrepareShiny(object, samples, Rshiny_path, split.by = c("conditions","orig.ident"),
             reduction = "umap",assay = "SCT")