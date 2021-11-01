########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
require(fields)
require(parallel)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  2. DoubletFinder 
# 
# ######################################################################

# samples

object = readRDS("data/OSU_SCT_20210821.rds")
object_list <- SplitObject(object,split.by = "orig.ident")
remove(object);GC()

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
Sys.setenv("OMP_NUM_THREADS" = 16)
npcs <- 100
sweep.res_list <- pbapply::pblapply(object_list,function(obj){
    paramSweep_v4(obj, PCs = 1:npcs, sct = T)
})

save(sweep.res_list,file = "output/OSU_SCT_20211030_sweep.res_list.Rda")
load("output/Lung_30_20211030_sweep.res_list.Rda")
sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)
# find histgram local maximam
find.localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }
    which(x == max(x[y]))
}

(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))
maximal_pk

# http://rstudio-pubs-static.s3.amazonaws.com/329613_f53e84d1a18840d5a1df55efb90739d9.html
qplot_2axis <- function(data,x = "pK", y1 = "MeanBC", y2 = "BCmetric"){
    if(class(data[,x]) == "factor") data[,x] <- as.numeric(as.character(data[,x]))
    data_y1 <- data[,y1]
    data_y2 <- data[,y2]
    a <- range(data_y1)
    b <- range(data_y2)
    scale_factor <- diff(a)/diff(b)
    data_y2 <- ((data_y2 - b[1]) * scale_factor) + a[1]
    trans <- ~ ((. - a[1]) / scale_factor) + b[1]
    
    g <- ggplot(data = data, aes_string(x = x, y = y1))+
        geom_line()+geom_point()+
        geom_point(aes(y = data_y2),colour = "blue")+
        geom_line(aes(y = data_y2),colour = "blue")+
        scale_y_continuous(name = y1,
                           sec.axis = sec_axis(trans=trans, name=y2))+
        theme(axis.text.y.right = element_text(color = "blue"))
    
    g
    
}
#qplot_2axis(data = bcmvn_list[[2]])

Multiplet_Rate <- function(object, numBatches = 1, num10xRuns = 1){
    
    numCellsRecovered = 1.0 * ncol(object)
    m = 4.597701e-06
    r = 0.5714286
    
    numCellsLoaded = numCellsRecovered / r
    multipletRate = m * numCellsLoaded / num10xRuns
    
    singletRate = 1.0 - multipletRate;
    numSinglet = singletRate * numCellsRecovered
    numMultiplet = numCellsRecovered - numSinglet
    numIdentMultiplet = numMultiplet * (numBatches - 1) / numBatches
    numNonIdentMultiplet = numMultiplet - numIdentMultiplet
    numCells = numSinglet + numNonIdentMultiplet
    
    return(numNonIdentMultiplet/numCells)
}
Multiplet_Rate(object_list[[1]])
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
for(i in 1:length(object_list)){
    print(paste("processing",unique(object_list[[i]]$orig.ident)))
    homotypic.prop <- modelHomotypic(object_list[[i]]@meta.data$cell.types)
    nExp_poi <- round(Multiplet_Rate(object_list[[i]])*length(colnames(object_list[[i]])))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    object_list[[i]] <- doubletFinder_v4(object_list[[i]], PCs = 1:50, 
                                         pN = 0.25, pK = maximal_pk[i], 
                                         nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
    object_list[[i]] <- doubletFinder_v4(object_list[[i]], PCs = 1:50,
                                         pN = 0.25, pK = maximal_pk[i],
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = grep("pANN",colnames(object_list[[i]]@meta.data),value = T),
                                         sct = TRUE)
    colName = colnames(object_list[[i]]@meta.data)
    colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                    "High_confident_doublets")
    colnames(object_list[[i]]@meta.data) = colName
    Progress(i,length(object_list))
}

for(i in 1:length(object_list)){
    object_list[[i]]@meta.data$row.names = rownames(object_list[[i]]@meta.data)
}
meta.data_list <- lapply(object_list, function(x) {
    temp <- x@meta.data
    temp$row.names = rownames(temp)
    return(temp)
    })
meta.data = bind_rows(meta.data_list)
rownames(meta.data) = meta.data$row.names

object = readRDS("data/OSU_SCT_20210821.rds")
meta.data = meta.data[rownames(object@meta.data),]
meta.data$doublets = gsub("Doublet","Doublet-Low Confidence",meta.data$Low_confident_doublets)
meta.data[meta.data$High_confident_doublets %in% "Doublet","doublets"] = "Doublet-High Confidence"
meta.data = cbind(object@meta.data,meta.data$doublets)
colnames(meta.data)[ncol(meta.data)] = "Doublets"
table(meta.data$Doublets)


object@meta.data = meta.data
saveRDS(object, file = paste0("data/OSU_SCT_20210821.rds"))

object = readRDS(file = "data/Lung_SCT_30_20210831.rds")

UMAPPlot.1(object, group.by = "Doublets",cols = c("red","orange","black"), 
           title = "Singlets and possible Doublets", do.print = T,
           do.return = F,pt.size = 0.3)

# 1. How many doublets compared to the total cell number? 
cell.number = table(object$Doublets) 
percentage = prop.table(cell.number)*100
cbind(cell.number,percentage) %>% kable() %>% kable_styling()

cell.number = table(object$orig.ident,object$Doublets) 
percentage = prop.table(cell.number,margin = 1)*100
cbind(cell.number,percentage) %>% kable() %>% kable_styling()

cell.number = table(object$SCT_snn_res.0.06,object$Doublets) 
percentage = prop.table(cell.number,margin = 1)*100
cbind(cell.number,percentage) %>% kable() %>% kable_styling()

cell.number = table(object$conditions,object$Doublets) 
percentage = prop.table(cell.number,1)*100
cbind(cell.number,percentage) %>% kable() %>% kable_styling()
