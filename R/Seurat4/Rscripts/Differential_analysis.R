########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
#======1.2 load  Seurat =========================
object = readRDS("data/OSU_SCT_20210821.rds")
# Need 32GB
#DefaultAssay(object) = "SCT"
#Idents(object) = "Doublets"
#object <- subset(object, idents = "Singlet")

opts = data.frame(cell.types = c(rep("B-cells",3),
                            rep("MDSCs",3),
                            rep("Monocytes",3),
                            rep("NK cells",3),
                            rep("T-cells:CD4+",3),
                            rep("T-cells:CD8+",3),
                            rep("T-cells:regs",3)),
                  response = c(rep(c("PD","PR","SD"),times = 7))
)

print(opt <- opts[args,])
object %<>% subset(subset = cell.types %in% opt$cell.types
                   & response %in% opt$response)
Idents(object) = "treatment"
system.time(markers <- FindMarkers_UMI(object, 
                                       ident.1 = "Ibrutinib",
                                       group.by = "treatment",
                                       logfc.threshold = 0.1, 
                                    only.pos = F,
                                       latent.vars = "nFeature_SCT",
                                       test.use = "MAST"))

markers$cell.type = opt$cell.types
markers$response = opt$response
markers$cluster = "Ibrutinib vs Baseline"

arg = args
if(args < 10) arg = paste0("0", args)
write.csv(markers,paste0(path,arg,"_FC0.1_",opt$cell.types,"_",opt$response,".csv"))