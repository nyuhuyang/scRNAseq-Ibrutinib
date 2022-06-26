# Library
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(magrittr)
# load data
df_samples <- readxl::read_excel("doc/20210927_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples = df_samples[df_samples$paired == "T",]
data = df_samples[,c("patient","sample","treatment","disease","response","cluster","cluster_change",
                     "tcrs_change","batch")]
data= data[!duplicated(data$patient),]
chisq_test(data$cluster_change, y = data$batch,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$treatment,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$disease,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$response,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$tcrs_change,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$cluster,simulate.p.value = T)

chisq_test(data$cluster, y = data$batch,simulate.p.value = T)
chisq_test(data$cluster, y = data$treatment,simulate.p.value = T)
chisq_test(data$cluster, y = data$disease,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$response,simulate.p.value = T)
chisq_test(data$cluster_change, y = data$tcrs_change,simulate.p.value = T)
value_data = apply(data,2,function(x) as.integer(as.factor(as.character(x)))) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("sample")
    #pivot_longer(!patient, names_to = "factors",values_to = "value")

res.chisq = pairwise_chisq_test_against_p(as.matrix(value_data))
heatplot(value_data[,-c(1,3)], scale = "none",dualScale=FALSE)
# Color Brewer palette
library(viridis)
ggplot(value_data, aes(factors, value, fill= value)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    theme_ipsum()


if(is.null(order.yaxis)){
    if(isTRUE(Rowv)) {
        hcr <- hclust(as.dist(1-cor(t(mtx_fgseaRes), method="spearman")),
                      method="ward.D2")
        ddr <- as.dendrogram(hcr)
        rowInd <- order.dendrogram(ddr)
        order.yaxis = rownames(mtx_fgseaRes)[rowInd]
    } else {
        order.yaxis = rownames(mtx_fgseaRes)
        order.yaxis = order.yaxis[order.yaxis %in% df_fgseaRes[,"pathway"]]
    }
}
order.yaxis = order.yaxis[order.yaxis %in% rownames(mtx_fgseaRes)]
df_fgseaRes %<>% filter(pathway %in% order.yaxis)
df_fgseaRes[,"pathway"] %<>% factor(levels = rev(order.yaxis))

if(isTRUE(Colv)) {
    hcc <- hclust(as.dist(1-cor(mtx_fgseaRes, method="spearman")),
                  method="ward.D2")
    ddc <- as.dendrogram(hcc)
    colInd <- order.dendrogram(ddc)
    order.xaxis = colnames(mtx_fgseaRes)[colInd]
}