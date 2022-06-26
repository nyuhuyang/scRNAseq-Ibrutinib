# prepare Colour palette
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(scales)

### Useful stuff
# Colour palette

cList = list(c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59",
                "#EF6548","#D7301F","#B30000","#7F0000"),
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF",
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)],
             c("blue","cyan3","green","greenyellow","yellow","orange","chocolate1","red","darkred")[c(1,1:9,9)],
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C",
               "#2C728E","#3B528B","#472D7B","#440154"),
             c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
             c("#7FC97F","#BEAED4","#FDC086","#386CB0","#F0027F",
               "#BF5B17","#666666","#1B9E77","#7570B3","#66A61E",
               "#E6AB02","#A6761D","#A6CEE3","#B2DF8A","#FB9A99",
               "#E31A1C","#FF7F00","#6A3D9A","#8DA0CB",
               "#4DAF4A","#984EA3","#c6c386","#999999","#66C2A5",
               "#FC8D62","#A6D854","#FFD92F","#BEBADA",
               "#FB8072","#80B1D3","#FDB462","#BC80BD","#B3B3B3",
               "#33A02C","#B3DE69","#4038b0","#ee7576","#e94749","#E78AC3","#ff0000",
               "#A65628","#d80172","#F781BF","#D95F02","#E7298A",
               "#1F78B4","#FDBF6F","#CAB2D6","#B15928","#FBB4AE",
               "#B3CDE3",
               '#0173b2','#de8f05','#029e73','#d55e00','#cc78bc','#ca9161','#fbafe4','#949494','#ece133','#56b4e9', # seaborn.color_palette colorblind
               "#00AFBB", "#E7B800", "#FC4E07",
               "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"))

names(cList) = c("default","Blue-White-Red","Blue-Yellow-Red","Yellow-Green-Purple","colorBlind","Singler")
custom.pal.info =data.frame("maxcolors" = sapply(cList,length),
                            "category" = c("seq","div","div","div","qual","qual"))
rownames(custom.pal.info) = names(cList)
custom.pal.info$palette = NA
custom.pal.info$package = "custom"


ggsci.pal.info = t(data.frame(
    c("NPG", 10, "qual", "nrc"),
    c("AAAS", 10, "qual", "default"),
    c("NEJM", 8, "qual", "default"),
    c("Lancet", 9, "qual", "lanonc"),
    c("JAMA", 7, "qual", "default"),
    c("JCO", 10, "qual", "default"),
    c("UCSCGB", 26, "qual", "default"),
    c("D3.10", 10, "qual", "category10"),
    c("D3.20", 20, "qual", "category20"),
    c("D3.20b", 20, "qual", "category20b"),
    c("D3.20c", 20, "qual", "category20c"),
    c("LocusZoom", 7, "qual", "default"),
    c("IGV", 51, "qual", "default"),
    c("IGV.alternating", 2, "div", "alternating"),
    c("UChicago", 9, "qual", "default"),
    c("Uchicago.light", 9, "qual", "light"),
    c("Uchicago.dark", 9, "qual", "dark"),
    c("StarTrek", 7, "qual", "uniform"),
    c("Tron", 7, "qual", "legacy"),
    c("Futurama", 12, "qual", "planetexpress"),
    c("RickandMorty", 12, "qual", "schwifty"),
    c("Simpsons", 16, "qual", "springfield"),
    c("GSEA", 12, "div", "default"),
    c("Material.red", 10, "seq", "red"),
    c("Material.pink", 10, "seq", "pink"),
    c("Material.purple", 10, "seq", "purple"),
    c("Material.indigo", 10, "seq", "indigo"),
    c("Material.blue", 10, "seq", "blue"),
    c("Material.light-blue", 10, "seq", "light-blue"),
    c("Material.cyan", 10, "seq", "cyan"),
    c("Material.teal", 10, "seq", "teal"),
    c("Material.green", 10, "seq", "green"),
    c("Material.light-green", 10, "seq", "light-green"),
    c("Material.lime", 10, "seq", "lime"),
    c("Material.yellow", 10, "seq", "yellow"),
    c("Material.amber", 10, "seq", "amber"),
    c("Material.orange", 10, "seq", "orange"),
    c("Material.deep-orange", 10, "seq", "deep-orange"),
    c("Material.brown", 10, "seq", "brown"),
    c("Material.grey", 10, "seq", "grey"),
    c("Material.blue-grey", 10, "seq", "blue-grey"))
)
ggsci.pal.info %<>% as.data.frame()
colnames(ggsci.pal.info) = c("Name","maxcolors","category","palette")
rownames(ggsci.pal.info) = ggsci.pal.info[,"Name"]
ggsci.pal.info[,"maxcolors"] %<>% as.integer()
ggsci.pal.info = ggsci.pal.info[, -1]
ggsci.pal.info$package = "ggsci"

brewer.pal.info$package = "RColorBrewer"
brewer.pal.info$palette = NA
brewer.pal.info$colorblind = NULL
brewer.pal.info = brewer.pal.info[,colnames(ggsci.pal.info)]
pal.info = bind_rows(list(custom.pal.info, ggsci.pal.info,brewer.pal.info))
