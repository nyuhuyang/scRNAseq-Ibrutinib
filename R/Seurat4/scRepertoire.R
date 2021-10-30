########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("pbapply","magrittr","scRepertoire"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

########################################################################
#
#  1 read file
#
# ######################################################################
#======1.1 reac TCR file =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples$patient %<>% as.character()
df_samples$`cite-seq` %<>% as.character()


# for each sample
contig_list <- lapply(df_samples$sample, function(s){
    s_id = df_samples[df_samples$sample %in% s,"sample.id"]
    tcr_folder = paste0("data/counts/",s_id,"/")
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep="/"))
    tcr$barcode <- gsub("-1", "", tcr$barcode)
    tcr
})

combined <- combineTCR(contig_list,sample=df_samples$sample, cells = "T-AB")
combined <- lapply(combined, function(obj){
    s = unique(obj$sample)
    obj$barcode = gsub(paste0(s,"_"),paste0(s,"-"), obj$barcode)
    obj
})

saveRDS(combined, "output/TCR.rds")

#======1.2 load  Seurat =========================
object = readRDS("data/OSU_SCT_20210821.rds")
object$barcode = gsub(".*-","",colnames(object))
object$barcode %<>% paste0(object$sample,"-",.)
object %<>% RenameCells(new.names = object$barcode)

object <- combineExpression(combined, object, cloneCall="gene+nt")
object@meta.data[is.na(object$Frequency),"Frequency"] = 0

saveRDS(object,"data/OSU_SCT_20210821.rds")

saveRDS(object@meta.data,"output/20211029/meta_data.rds")

meta_data = readRDS("output/20211029/meta_data.rds")
combined <- split(meta_data,f = meta_data$treatment)

jpeg(paste0(path,"clonalDiversity_treatment.jpeg"), units="in", width=10, height=7,res=600)
clonalDiversity(combined, cloneCall = "gene", group = "treatment")
dev.off()

meta_data = readRDS("output/20211029/meta_data.rds")
combined <- split(meta_data,f = meta_data$sample)


jpeg(paste0(path,"quantContig.jpeg"), units="in", width=10, height=7,res=600)
quantContig(combined, cloneCall="gene+nt", scale = TRUE)+
    theme(axis.text.x =    element_text(angle = as.numeric(45),
                                        hjust = switch (as.character(45),
                                                        "0" = 0.5,
                                                        "30" = 1,
                                                        "45" = 1,
                                                        "90" = 1
                                        ),
                                        vjust = switch (as.character(45),
                                                        "0" = 0,
                                                        "30" = 1,
                                                        "45" = 1,
                                                        "90" = 0.5))
          )
dev.off()

jpeg(paste0(path,"compareClonotypes_1.jpeg"), units="in", width=10, height=7,res=600)
compareClonotypes(combined, numbers = 10, samples = c("1_1", "1_2"), 
                  cloneCall="aa", graph = "alluvial")
dev.off()

jpeg(paste0(path,"compareClonotypes_2.jpeg"), units="in", width=10, height=7,res=600)
compareClonotypes(combined, numbers = 10, samples = c("2_1", "2_2"), 
                  cloneCall="aa", graph = "alluvial")
dev.off()

jpeg(paste0(path,"compareClonotypes_3.jpeg"), units="in", width=10, height=7,res=600)
compareClonotypes(combined, numbers = 10, samples = c("4_1", "4_2"), 
                  cloneCall="aa", graph = "alluvial")
dev.off()

jpeg(paste0(path,"compareClonotypes_5.jpeg"), units="in", width=10, height=7,res=600)
compareClonotypes(combined, numbers = 10, samples = c("5_1", "5_2"), 
                  cloneCall="aa", graph = "alluvial")
dev.off()

sub_combined <- subsetContig(combined, name = "sample", variables = c("1_1", "1_2","2_1","2_2","4_1","4_2"))

jpeg(paste0(path,"abundanceContig.jpeg"), units="in", width=10, height=7,res=600)
abundanceContig(combined, cloneCall = "gene", scale = F)+ 
    scale_y_continuous(trans='log10')+
    scale_color_manual(values=color_generator("colorBlind",30))
dev.off()

jpeg(paste0(path,"scatterClonotype_1.jpeg"), units="in", width=10, height=7,res=600)
scatterClonotype(combined, x.axis = "1_1", y.axis = "1_2")
dev.off()

jpeg(paste0(path,"scatterClonotype_2.jpeg"), units="in", width=10, height=7,res=600)
scatterClonotype(combined, x.axis = "2_1", y.axis = "2_2")
dev.off()

jpeg(paste0(path,"scatterClonotype_3.jpeg"), units="in", width=10, height=7,res=600)
scatterClonotype(combined, x.axis = "4_1", y.axis = "4_2")
dev.off()

jpeg(paste0(path,"scatterClonotype_5.jpeg"), units="in", width=10, height=7,res=600)
scatterClonotype(combined, x.axis = "5_1", y.axis = "5_2")
dev.off()

jpeg(paste0(path,"clonalHomeostasis_1.jpeg"), units="in", width=10, height=7,res=600)
clonalHomeostasis(sub_combined, cloneCall = "gene")
dev.off()

jpeg(paste0(path,"clonalDiversity_1.jpeg"), units="in", width=10, height=7,res=600)
clonalDiversity(combined, cloneCall = "gene", group = "treatment")
dev.off()



write.csv(object@meta.data, file = paste0(path,"TCR_meta_data.csv"))
