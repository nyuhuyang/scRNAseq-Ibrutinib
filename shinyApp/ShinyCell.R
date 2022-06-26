library(Seurat)
library(ShinyCell)
library(magrittr)
library(SeuratData)
library(SeuratDisk)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

object = readRDS("data/OSU_SCT_20210821.rds")


#=============
pd <- reticulate::import("pandas")
antigen = pd$read_hdf("output/antigen_mcpas_all.h5", key = "scirpy")
antigen %<>% as.matrix()
antigen[is.na(antigen)] = "NA"
antigen %<>% as.data.frame()

antigen[,"barcode"] = rownames(antigen)
object@meta.data$barcode = colnames(object)

meta.data = object@meta.data
meta.data %<>% full_join(antigen,by = "barcode")
rownames(meta.data) = meta.data$barcode
meta.data = meta.data[,c("Category","antigen.species",
             "antigen.protein","antigen.epitope")]
table(meta.data$Category)
colnames(meta.data) %<>% paste0("mcpas.all.",.)
if(all(colnames(object) == rownames(meta.data))) object@meta.data %<>% cbind(meta.data)
#==============
antigen = pd$read_hdf("output/antigen_mcpas_Cancer.h5", key = "scirpy")
antigen %<>% as.matrix()
antigen[is.na(antigen)] = "unknown"
antigen %<>% as.data.frame()

antigen[,"barcode"] = rownames(antigen)
object@meta.data$barcode = colnames(object)

meta.data = object@meta.data
meta.data %<>% full_join(antigen,by = "barcode")
rownames(meta.data) = meta.data$barcode
meta.data = meta.data[,c("Category","antigen.species",
                         "antigen.protein","antigen.epitope")]
table(meta.data$Category)
colnames(meta.data) %<>% paste0("mcpas.cancer.",.)
if(all(colnames(object) == rownames(meta.data))) object@meta.data %<>% cbind(meta.data)
object$mcpas.cancer.Category.value = plyr::mapvalues(object$mcpas.cancer.Category,
                                                     from = c("Cancer","unknown"),
                                                     to  = c(1,0))
object$mcpas.cancer.Category.value %<>% as.integer()
# ========
# manually annotation based on mcpas
meta_data = readRDS("output/20220309/meta_data.rds")

meta_data[meta_data$Pathology == "none","Pathology"] = "unknown" 
meta_data[meta_data$Category == "none","Category"] = "unknown" 

if(all(colnames(object) == rownames(meta_data))) {
    object$tcr_clonetype = meta_data$tcr_clonetype
    object$Pathology = meta_data$Pathology
    object$Category = meta_data$Category
}
#s.genes <- cc.genes$s.genes
#s.genes %<>% gsub("MLF1IP","CENPU",.)
#g2m.genes <- cc.genes$g2m.genes
#g2m.genes %<>% plyr::mapvalues(from = c("FAM64A", "HN1"),
#                               to = c("PIMREG","JPT1"))
#object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#colnames(object@meta.data)[grep("Phase",colnames(object@meta.data))]="cell.cycle.phase"
#colnames(object@meta.data)[grep("Frequency",colnames(object@meta.data))]="tcr.frequency"
#object$response %<>% factor(levels = c("PR","SD","PD"))
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
#df_samples %<>% filter(paired %in% "T" )
#object %<>% subset(subset = patient %in% c("8","9"), invert = T)

#object$sample %<>% paste0("PT",.)
#object$sample %<>% gsub("PT","-pre",.)
#object$sample %<>% gsub("_2","-post",.)
object$sample %<>% factor(levels = df_samples$sample)
object %<>% subset(subset = Doublets %in% "Singlet")

colnames(object@meta.data) %<>% sub("^sample$","Sample",.)
object$Sample %<>% droplevels()
object$patient %<>% droplevels()
object$combine_cell.types = gsub("T-cells.*","T-cells",object$cell.types)
#object$combine_response = gsub("PR|SD","PR+SD",object$response)
#object$combine_response %<>% factor(levels = c("PR+SD","PD"))
#object$cell.types_response = paste0(object$cell.types,"_",object$response)
object$acid = plyr::mapvalues(object$FTH1_lvl,
                                  from = c("FTH1 high","FTH1 low"),
                                  to = c("acid","normal"))
object$acid_cell.types = paste0(object$acid," ",object$cell.types)
object$treatment_response = paste0(object$treatment,"_",object$response)
object$treatment_patient = paste0(object$treatment,"_",object$patient)


meta.to.include =c("cell.types","Sample","acid","acid_cell.types","tcr_clonetype",
                   "mcpas.cancer.Category","mcpas.cancer.Category.value",
                   "mcpas.cancer.antigen.species",
                   "mcpas.cancer.antigen.protein","mcpas.cancer.antigen.epitope",
                   "mcpas.all.Category","mcpas.all.antigen.species",
                   "Pathology","Category",
                   "patient","timepoint","treatment","response","disease",
                   "label.fine","orig.ident","treatment_response","treatment_patient",
                   "celltype.l1","celltype.l2","celltype.l3",
                   "combine_cell.types","tcr.frequency","cloneType",
                   "S.Score","G2M.Score","cell.cycle.phase","SCT_snn_res.0.8",
                   "nCount_SCT","nFeature_SCT","percent.mt",
                   "Mean.Reads.per.Cell",
                   "Number.of.Reads","Sequencing.Saturation",
                   "CTgene","CTnt","CTaa","CTstrict"
                   )
table(meta.to.include %in% colnames(object@meta.data))
meta.to.include[!meta.to.include %in% colnames(object@meta.data)]
scConf = createConfig(object, meta.to.include = meta.to.include, maxLevels = length(unique(object$CTstrict))+1)

for(ui in c("CTgene","CTnt","CTaa","CTstrict")){
    scConf[UI == ui]$fID = scConf[UI == ui]$fUI = scConf[UI == ui]$fCL =  NA
    scConf[UI == ui]$grp =  FALSE
    scConf[UI == ui]$fRow =  NA
    
    
}

makeShinyApp(object, scConf, gex.assay = "SCT",gene.mapping = TRUE,
             gex.slot = "data",default.gene1 = "FTH1",default.gene2 = "PTPRC",
             default.multigene = c("CD19","CD3D","CD4","CD8A","CD14","FCGR3A","FCGR1A",
                                   "VEGFA","FTL","S100A6", "LYZ",
                                   "GNLY","KLRC1","NCAM1"),
            default.dimred = c("UMAP_1","UMAP_2"),shiny.dir = "shinyApp/Ibrutinib/",
             shiny.title = "PBMC from Ibrutinib plus Nivolumab on a phase 1")


max_exp = qlcMatrix::rowMax(object@assays[["SCT"]]@data)
max_exp_df = data.frame("val"= as.vector(max_exp),row.names = rownames(object))
saveRDS(max_exp_df,"shinyApp/Ibrutinib/sc1maxlvl.rds")


meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell.types),]
meta.data = meta.data[order(meta.data$cell.types),]
sc1conf = readRDS("shinyApp/Ibrutinib/sc1conf.rds")
sc1conf$fCL[1] = paste(meta.data$cell.types.colors,collapse = "|")
sc1conf[UI == "mcpas.cancer.Category"]$fCL = "#570b90|#bcbcbc"
#sc1conf$fCL[grep("Doublets",sc1conf$ID)] = "red|orange|black"
saveRDS(sc1conf,"shinyApp/Ibrutinib/sc1conf.rds")


sc1def  = readRDS("shinyApp/Ibrutinib/sc1def.rds")
sc1def$grp1 = "cell.types"
sc1def$grp2 = "treatment"
saveRDS(sc1def,"shinyApp/Ibrutinib/sc1def.rds")


format(object.size(object),unit = "GB")
format(object.size(object@assays$SCT),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
DefaultAssay(object) = "SCT"
object[["SCT"]]@counts = matrix(0,0,0)
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
format(object.size(object@meta.data),unit = "GB")
object@meta.data = object[["orig.ident"]]


file.remove("shinyApp/Ibrutinib/sc1csr_gexpr.h5ad")
SaveH5Seurat(object, filename = "shinyApp/Ibrutinib/sc1csr_gexpr.h5Seurat")
Convert("shinyApp/Ibrutinib/sc1csr_gexpr.h5Seurat", dest = "h5ad")
file.remove("shinyApp/Ibrutinib/sc1csr_gexpr.h5Seurat")



#=============================
object = readRDS("data/OSU_normal_20211108.rds")
meta.data = object@meta.data
cancer = meta.data[meta.data$orig.ident != "SeuratProject",]
meta.data[meta.data$orig.ident != "SeuratProject","celltype.l3"] = cancer$label.fine
meta.data[meta.data$orig.ident != "SeuratProject","celltype.l1"] = cancer$cell.types
meta.data[meta.data$orig.ident %in% "SeuratProject","patient"] = "Normal"
meta.data[meta.data$orig.ident %in% "SeuratProject","treatment"] = "Normal"
meta.data[meta.data$orig.ident %in% "SeuratProject","response"] = "Normal"
meta.data[meta.data$orig.ident %in% "SeuratProject","FTH1_lvl"] = "Normal"
object@meta.data = meta.data
object$treatment %<>% factor(levels = c("Normal", "Baseline","Ibrutinib"))

meta.to.include =c("celltype.l1","celltype.l3",
                   "patient","timepoint","treatment","response")
table(meta.to.include %in% colnames(object@meta.data))
meta.to.include[!meta.to.include %in% colnames(object@meta.data)]
scConf = createConfig(object, meta.to.include = meta.to.include, maxLevels = length(unique(object$celltype.l3))+1)

makeShinyApp(object, scConf, gex.assay = "SCT",gene.mapping = TRUE,
             gex.slot = "data",default.gene1 = "FTH1",default.gene2 = "FTL",
             default.multigene = c("CD19","CD3D","CD4","CD8A","CD14","FCGR3A","FCGR1A",
                                   "FTH1","FTL","S100A6", "LYZ",
                                   "GNLY","KLRC1","NCAM1"),
             default.dimred = c("UMAP_1","UMAP_2"),shiny.dir = "shinyApp/Ibrutinib_normal/",
             shiny.title = "PBMC from Ibrutinib plus Nivolumab on a phase 1 + normal PBMC")


max_exp = qlcMatrix::rowMax(object@assays[["SCT"]]@data)
max_exp_df = data.frame("val"= as.vector(max_exp),row.names = rownames(object@assays[["SCT"]]@data))
saveRDS(max_exp_df,"shinyApp/Ibrutinib_normal/sc1maxlvl.rds")


meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell.types),]
meta.data = meta.data[order(meta.data$cell.types),]
sc1conf = readRDS("shinyApp/Ibrutinib_normal/sc1conf.rds")
sc1conf$fCL[1] = paste(meta.data$cell.types.colors,collapse = "|")
#sc1conf$fCL[grep("Doublets",sc1conf$ID)] = "red|orange|black"
saveRDS(sc1conf,"shinyApp/Ibrutinib_normal/sc1conf.rds")


sc1def  = readRDS("shinyApp/Ibrutinib_normal/sc1def.rds")
sc1def$grp1 = "celltype.l1"
sc1def$grp2 = "treatment"
saveRDS(sc1def,"shinyApp/Ibrutinib_normal/sc1def.rds")


format(object.size(object),unit = "GB")
format(object.size(object@assays$SCT),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
object[["RNA"]]@counts = matrix(0,0,0)
object[["RNA"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
format(object.size(object@reductions),unit = "GB")
object@meta.data = object[["orig.ident"]]


file.remove("shinyApp/Ibrutinib_normal/sc1csr_gexpr.h5ad")
SaveH5Seurat(object, filename = "shinyApp/Ibrutinib_normal/sc1csr_gexpr.h5Seurat")
Convert("shinyApp/Ibrutinib_normal/sc1csr_gexpr.h5Seurat", dest = "h5ad")
file.remove("shinyApp/Ibrutinib_normal/sc1csr_gexpr.h5Seurat")


