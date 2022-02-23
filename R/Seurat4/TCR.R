########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#=========================== after scRepertoire ===================
meta_data = readRDS("output/20211029/meta_data.rds")
# common tcr
common_tcr_list <- split(meta_data[,c("patient","CTaa")],f = meta_data$patient)
patients = levels(meta_data$patient)
common_tcr <- c()
for(i in 1:(length(common_tcr_list)-1)){
    for(k in (i+1):length(patients)){
        new_common_tcr = intersect(common_tcr_list[[i]]$CTaa,
                                   common_tcr_list[[k]]$CTaa)
        print(paste( patients[i],"and",  patients[k], "share",length(new_common_tcr),"tcr"))
        common_tcr = c(common_tcr, new_common_tcr)
        
    }
}
rm(common_tcr_list);GC()
common_tcr = unique(common_tcr)
common_tcr = common_tcr[!is.na(common_tcr)]
# Persistent tcr

persistent_tcr_list <- split(meta_data[,c("treatment","CTaa")],f = meta_data$treatment)
timepoints = unique(meta_data$treatment)
persistent_tcr = intersect(persistent_tcr_list[["Baseline"]]$CTaa,
                           persistent_tcr_list[["Ibrutinib"]]$CTaa)
rm(persistent_tcr_list);GC()
persistent_tcr = unique(persistent_tcr);length(persistent_tcr)
persistent_tcr = persistent_tcr[!is.na(persistent_tcr)]
meta_data$persistent_tcr = meta_data$CTaa %in% persistent_tcr

tcr_counter = as.data.frame(table(meta_data$CTaa))
tcr_counter$Var1 %<>% as.character()
meta_data$enriched_tcr = meta_data$CTaa %in% tcr_counter[tcr_counter$Freq > 1,"Var1"] & !is.na(meta_data$CTaa)
meta_data$singlets_tcr = meta_data$CTaa %in% tcr_counter[tcr_counter$Freq == 1,"Var1"]  & !is.na(meta_data$CTaa)

meta_data$tcr_clonetype = meta_data$CTaa
meta_data$tcr_clonetype[meta_data$singlets_tcr] = "singlets"
meta_data$tcr_clonetype[meta_data$enriched_tcr] = "enriched"
meta_data$tcr_clonetype[meta_data$persistent_tcr] = "persistent"
meta_data$tcr_clonetype[meta_data$common_tcr] = "common"
meta_data$tcr_clonetype[meta_data$common_tcr & meta_data$persistent_tcr] = "persistent & common"
meta_data$persistent_tcr = NULL;
meta_data$enriched_tcr = NULL;
meta_data$singlets_tcr = NULL;

meta_data = meta_data[colnames(object),]
meta_data$tcr_clonetype %<>% factor(levels = c("singlets","persistent","enriched","common","persistent & common"))
saveRDS(meta_data,paste0(path,"meta_data.rds")) #20220207
write.csv(meta_data[c("orig.ident","patient","timepoint","treatment","response","cell.types",
                      "frequency","proportion","cdr3s_aa","tcr_clonetype")],file = paste0(path,"meta_data.csv"))




########################################################################
#
#  1 read file
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples$patient %<>% as.character()
df_samples$`cite-seq` %<>% as.character()

#======1.2 load  Seurat =========================
object = readRDS("data/OSU_SCT_20210821.rds")

meta.data = object@meta.data
meta.data$barcode = rownames(meta.data)
meta.data_list <- list()
# for each sample
for(i in 1:length(df_samples$sample.id)){
    s <- df_samples$sample.id[i]
    cells <- meta.data$orig.ident %in% s
    print(s)
    print(table(cells))
    tcr_folder = paste0("data/counts/",s)
    tcr  <- add_tcr_clonotype(tcr_folder)
    tcr$barcode %<>% paste0(s,"-",.)
    
    sub_meta.data = meta.data[cells,]
    sub_meta.data %<>% left_join(tcr,by = "barcode")
    rownames(sub_meta.data) = sub_meta.data$barcode
    meta.data_list[[i]] = sub_meta.data
}
meta_data = bind_rows(meta.data_list)
rownames(meta_data) = meta_data$barcode
meta_data = meta_data[colnames(object),]
table(rownames(meta_data) == rownames(meta.data))
meta_data[is.na(meta_data$cdr3s_aa),"cdr3s_aa"] = "none"
meta_data[is.na(meta_data)] = 0

# for each patient
patients = unique(df_samples$patient)
meta.data_list <- list()
common_tcr_list <- list()
for(i in 1:length(patients)){
    s <- patients[i]
    cells <- meta_data$patient %in% s
    print(s)
    print(table(cells))
    sub_meta.data = meta_data[cells,]
    common_tcr_list[[i]] = unique(sub_meta.data$cdr3s_aa)
    sub_meta.dat_list <- split(sub_meta.data,f = sub_meta.data$orig.ident)
    if(length(sub_meta.dat_list) >1) {
        shared_tcr = intersect(sub_meta.dat_list[[1]]$cdr3s_aa,
                               sub_meta.dat_list[[2]]$cdr3s_aa)
        shared_tcr = shared_tcr[!shared_tcr %in% "none"]
        sub_meta.data$shared_tcr = sub_meta.data$cdr3s_aa %in% shared_tcr
    } else sub_meta.data$shared_tcr = FALSE
    meta.data_list[[i]] = sub_meta.data
}
meta_data = bind_rows(meta.data_list)


add_tcr_clonotype <- function(tcr_folder){
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep="/"))
    
    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr$barcode <- gsub("-1", "", tcr$barcode)
    tcr <- tcr[!duplicated(tcr$barcode), ]
    
    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    
    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep="/"))
    
    # Slap the AA sequences onto our original table by clonotype_id.
    tcr %<>% left_join(clono[, c("clonotype_id","frequency","proportion","cdr3s_aa")],by="clonotype_id")
    tcr$clonotype_id = NULL
    return(tcr)
}
