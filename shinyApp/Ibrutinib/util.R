library(shiny)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
library(ggsci)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(digest)
library(shinyhelper)
library(openxlsx)
library(ggplotify)
library(grid)
library(tidyr)
library(ggridges)
library(ggpubr)
library(vegan)
library(stringr)
library(rstatix)
library(ggpmisc)
library(forcats)
library(igraph)
library(networkD3)
#library(ggcorrplot)
#library(plotly)
#library(anndata)

source("util_palette.R")
set.seed(101)
color_generator <- function(palette.name, n=NULL, alpha = 1){
    if(palette.name == "white") return(rep("white",max(n,5)))
    N = as.integer(pal.info[palette.name,"maxcolors"])
    if(is.null(n)) n = N
    #enough.color <- as.character(pal.info[palette.name,"maxcolors"] >= n)
    switch (pal.info[palette.name,"package"],
            "custom" = colorRampPalette(cList[[palette.name]])(n),
            "RColorBrewer" = colorRampPalette(brewer.pal(n = pal.info[palette.name,"maxcolors"],
                                                         palette.name))(n),
            "ggsci" = {
                ggsci_f <- get(paste0("pal_",tolower(sub("\\..*","",palette.name))))
                colorRampPalette(ggsci_f(palette = pal.info[palette.name,"palette"],alpha = alpha)(pal.info[palette.name,"maxcolors"]))(n)
            }
    )
}

# Panel sizes
pList = c("400px", "600px", "800px","1000px")
names(pList) = c("Small", "Medium", "Large", "Extra Large")
pList2 = c("500px", "700px", "900px")
names(pList2) = c("Small", "Medium", "Large")
sList = c(12,18,24,30)
names(sList) = c("Extra Small","Small", "Medium", "Large")
lList = c(4,5,6,7)
names(lList) = c("Extra Small","Small", "Medium", "Large")

doFactoring <- function(ggData,inpConf, col = "X", inp, inpcols) {
    if(class(ggData[[col]]) == "factor") {
        gglvl = levels(ggData[[col]])[levels(ggData[[col]]) %in% unique(ggData[[col]])]
        ggData[[col]] %<>% droplevels()
        ggData[[col]] %<>% factor(levels = gglvl)
    } else {
        gglvl = sort(unique(ggData[[col]]))
    }

    ggCol = strsplit(inpConf[UI == inp]$fCL, "\\|")[[1]]
    ggCol = switch(inpcols,
                   "default" = ggCol,
                   color_generator(inpcols,length(ggCol)))
    names(ggCol) = strsplit(inpConf[UI == inp]$fID, "\\|")[[1]]

    res = list(ggData, gglvl, ggCol)
    names(res) = c("ggData", "gglvl", "ggCol")
    return(res)
}

# Function to extract legend
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

# Plot theme
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){
    oupTheme = theme(
        text = element_text(size = base_size, family = "Avenir"),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line =     element_line(colour = "black"),
        axis.ticks =    element_line(colour = "black", size = base_size / 20),
        axis.title =    element_text(colour = "black"),
        axis.text =     element_text(size = base_size),
        axis.text.x = element_text(angle = Xang, hjust = XjusH),
        legend.position = "bottom",
        legend.text = element_text(size = base_size),
        legend.key = element_rect(colour = NA, fill = NA)
    )
    if(!XYval){
        oupTheme = oupTheme + theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    return(oupTheme)
}

### Common plotting functions
# Plot cell information on dimred
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsplt,
                     inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                     inpsiz, inpcols, inpord, inpfsz, inpasp, inptxt, inplab,inpleg = FALSE,
                     inplegpos = "bottom", inparrange = "auto"){
    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
                         inpConf[UI == inp1]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID,
                         inpConf[UI %in% inpsplt]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "Y", "group", "sub1","sub2","sub3")[1:ncol(ggData)]
    if(inpsplt %in% inpConf$ID) colnames(ggData)[ncol(ggData)] = "split"
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
    #bgCells = FALSE
    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)

    if(inpord == "Max-1st"){
        ggData = ggData[order(group)]
    } else if(inpord == "Min-1st"){
        ggData = ggData[order(-group)]
    } else if(inpord == "Random"){
        ggData = ggData[sample(nrow(ggData))]
    }
    temp <- doFactoring(ggData,inpConf, col = "group", inp1, inpcols)
    gglvl = temp$gglvl; ggCol = temp$ggCol[gglvl]
    ggData = temp$ggData

    # Actual ggplot
    ggOut = ggplot(ggData, aes(X, Y, color = group))
    #if(bgCells){
    #    ggOut = ggOut +
    #        geom_point(data = ggData, color = "snow2", size = inpsiz, shape = 16)
    #}
    ggOut = ggOut +
        geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) +
        sctheme(base_size = sList[inpfsz], XYval = inptxt)
    if(inpsplt %in% inpConf$ID){
        ggOut = ggOut + switch(inparrange,
                               "auto" = facet_wrap(facets = . ~ split, nrow = NULL,shrink = FALSE),
                               "1row" = facet_wrap(facets = . ~ split, nrow = 1,shrink = FALSE),
                               "1column" = facet_wrap(facets = split ~., ncol = 1,shrink = FALSE)
                               )
    }
    if(is.na(inpConf[UI == inp1]$fCL)) # continous variable
    {
        ggOut = ggOut + scale_color_gradientn("", colours = color_generator(inpcols,length(gglvl))) +
            guides(color = guide_colorbar(barwidth = switch(inplegpos,
                                                            "top" =15,
                                                            "right"=1.5,
                                                            "bottom" =15)))
    } else { # categorical variable
        sListX = min(nchar(paste0(levels(ggData$group), collapse = "")), 200)
        sListX = 0.75 * (sList - (1.5 * floor(sListX/50)))
        ggOut = ggOut + scale_color_manual("", values = switch(inpcols,
                                                               "default" =ggCol,
                                                               color_generator(inpcols,length(gglvl))))
        ggOut = ggOut + guides(color = guide_legend(override.aes = list(size = 5),
                                                    nrow = inpConf[UI == inp1]$fRow)) +
            theme(legend.text = element_text(size = sListX[inpfsz]))
        if(inplab != "No labels"){
            ggData = ggData[, .(X = mean(X), Y = mean(Y)), by = "group"]
            lListX = min(nchar(paste0(ggData$group, collapse = "")), 200)
            lListX = lList - (0.25 * floor(lListX/50))

            ggOut = ggOut + switch(inplab,
                                   "No labels" = NULL,
                                   "black text" = geom_text_repel(data = ggData, aes(X, Y, label = group),
                                                                  color = "grey10",
                                                                  box.padding = unit(0.5, "lines"),
                                                                  point.padding = unit(0.8, "lines"),
                                                                  size = lListX[inpfsz], seed = 42),
                                   "black labels" = geom_label_repel(data = ggData, aes(X, Y, label = group),
                                                                     color = "grey10",
                                                                     box.padding = unit(0.5, "lines"),
                                                                     point.padding = unit(0.8, "lines"),
                                                                     size = lListX[inpfsz], seed = 42),
                                   "color labels" = geom_label_repel(data = ggData, aes(X, Y, label = group),
                                                                     box.padding = unit(0.5, "lines"),
                                                                     point.padding = unit(0.8, "lines"),
                                                                     size = lListX[inpfsz], seed = 42)
            )
        }
    }

    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }
    ggOut = ggOut + theme(legend.position = switch (as.character(inpleg),
                                                    "FALSE" = "none",
                                                    "TRUE" = inplegpos))
    return(ggOut)
}

scDRcellnum <- function(inpConf, inpMeta, inp1, inpsplt,
                        inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2){
    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID,
                         inpConf[UI %in% inpsplt]$ID),
                     with = FALSE]
    colnames(ggData) = c("group", "sub1","sub2","sub3","split")[1:ncol(ggData)]
    if(inpsplt %in% inpConf$ID) colnames(ggData)[ncol(ggData)] = "split"
    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)

    # Actual data.table
    if(inpsplt %in% inpConf$ID){
        ggData = ggData[, .(nCells = .N), by = c("group","split")]
        ggData = dcast(ggData, group ~ split, value.var = "nCells")
    } else ggData = ggData[, .(nCells = .N), by = "group"]
    ggData = ggData[order(group)]

    return(ggData)
}

# Plot gene expression on dimred
scDRgeneData <- function(inpConf, inpMeta, inpMax, inpdrX, inpdrY, inp1,inpsplt,inpGrp,
                     inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                     inpH5, inpGene){
    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID,
                         inpConf[UI == inpdrY]$ID,
                         inpConf[UI == inpGrp]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID,
                         inpConf[UI %in% inpsplt]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "Y", "grpBy", "sub1","sub2","sub3","split")[1:ncol(ggData)]
    if(inpsplt %in% inpConf$ID) colnames(ggData)[ncol(ggData)] = "split"
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=)))
    ggData[val < 0]$val = 0
    h5file$close_all()

    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)

    return(ggData)
}

scDRgene <- function(ggData,inpConf, inpMax, inpdrX, inpdrY, inp1,inpsplt,
                     inpsiz, inpmax, inpcols,inpbground, inpord, inpfsz, inpasp,
                     inptxt, inpleg = FALSE, inplegpos = "bottom",inparrange = "auto"){

    # record coordinates
    max_x = max(ggData$X)
    min_x = min(ggData$X)
    max_y = max(ggData$Y)
    min_y = min(ggData$Y)
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
    #bgCells = FALSE

    if(inpmax) {
        ggData = rbindlist(list(ggData, head(ggData,1)))
        ggData[nrow(ggData),"X"] = max_x*20
        ggData[nrow(ggData),"Y"] = max_y*20
        ggData[nrow(ggData),"val"] = inpMax[inp1,"val"]
    }

    if(inpord == "Max-1st"){
        ggData = ggData[order(val)]
    } else if(inpord == "Min-1st"){
        ggData = ggData[order(-val)]
    } else if(inpord == "Random"){
        ggData = ggData[sample(nrow(ggData))]
    }

    # Actual ggplot
    ggOut = ggplot(ggData, aes(X, Y, color = val))
    #if(bgCells){
    #   ggOut = ggOut +
    #        geom_point(data = ggData, color = "snow", size = inpsiz, shape = 16)
    #}
    if(inpsplt %in% inpConf$ID){
        ggOut = ggOut + switch(inparrange,
                               "auto" = facet_wrap(facets = . ~ split, nrow = NULL,shrink = FALSE),
                               "1row" = facet_wrap(facets = . ~ split, nrow = 1,shrink = FALSE),
                               "1column" = facet_wrap(facets = split ~., ncol = 1,shrink = FALSE)
        )
    }
    ggOut = ggOut +
        geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) +
        sctheme(base_size = sList[inpfsz], XYval = inptxt) +
        theme(legend.position = switch (as.character(inpleg),
                                                        "FALSE" = "none",
                                                        "TRUE" = inplegpos))+
        scale_color_gradientn(inp1, colours = c(inpbground,color_generator(inpcols))) +
        guides(color = guide_colorbar(barwidth = switch (as.character(inplegpos),
                                                         "top" = sList[inpfsz]/1.5,
                                                         "right" = sList[inpfsz]/15,
                                                         "bottom" = sList[inpfsz]/1.5),
                                      barheight = switch (as.character(inplegpos),
                                                          "top" = sList[inpfsz]/15,
                                                          "right" = sList[inpfsz]/1.5,
                                                          "bottom" = sList[inpfsz]/15)
                                      ))

    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat,xlim = c(min_x, max_x),ylim = c(min_y,max_y))
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed(xlim = c(min_x, max_x),ylim = c(min_y,max_y))
    }

    return(ggOut)
}


DensityPlot <- function(sub_ggData,inpConf, inp1, inpGrp,
                        inpsiz, inpcol2,inpfsz, inpasp, inpyrev,
                        inptxt, density_type = "density_ridges"){

    # Do factoring
    ggCol = strsplit(inpConf[UI == inpGrp]$fCL, "\\|")[[1]]
    names(ggCol) = levels(sub_ggData$grpBy)
    gglvl = levels(sub_ggData$grpBy)[levels(sub_ggData$grpBy) %in% unique(sub_ggData$grpBy)]
    gglvl = switch(inpyrev,"As-it-is" = gglvl, "Reverse" = rev(gglvl))
    sub_ggData$grpBy = factor(sub_ggData$grpBy, levels = gglvl)

    # record coordinates
    max_x = max(sub_ggData$val)
    min_x = min(sub_ggData$val)
    max_y = length(gglvl)
    rat = (max(sub_ggData$X) - min(sub_ggData$X)) / max_y

    # Main plot
    colnames(sub_ggData) %<>% sub("val", inp1,.)
    colnames(sub_ggData) %<>% sub("grpBy", inpGrp,.)

    geom_x <- switch (density_type,
                      "density" =  list(geom_density(data = sub_ggData, aes_string(x = inp1, fill = inpGrp),
                                                     alpha = 0.5, size = 0.5)),
                      "density_ridges"= list(ggridges::geom_density_ridges(data = sub_ggData,
                                                                           aes_string(x = inp1, y = inpGrp,fill = inpGrp),
                                                                           alpha = 0.5, size = 0.5),
                                             scale_y_discrete(expand = expansion(mult = c(0.01, 0.5))))
    )
    ggOut <- do.call(what = '+', args = list(ggplot(sub_ggData, aes_string(x = inp1)), geom_x))+
        sctheme(base_size = sList[inpfsz], XYval = inptxt) +
        scale_fill_manual("", values = switch(inpcol2,
                                              "default" =ggCol,
                                              color_generator(inpcol2,length(gglvl))))
    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat/4,xlim = c(-0.2, max_x),ylim = c(0,max_y))
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed(xlim = c(-0.2, max_x),ylim = c(0,max_y))
    }

    return(ggOut)
}


scDRexNum <- function(ggData, inp1){
    if(nrow(ggData) == 0) return(data.table(`expression > 0`="none", nCells=0, percent=0))
    # Actual data.table
    ggData$express = "none"
    ggData[val > 0]$express = inp1
    ggData$express = factor(ggData$express, levels = unique(c(inp1, "none")))
    ggData = ggData[, .(nCells = .N), by = "express"]
    ggData$percent = 100 * ggData$nCells / sum(ggData$nCells)
    ggData = ggData[order(express)]
    colnames(ggData)[1] = "expression > 0"
    return(ggData)
}

# Plot gene coexpression on dimred
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){
    oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22
    oup = oup / (xy*xy)
    return(oup)
}


scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpGrp,
                     inpsub1_1, inpsub1_2, inpsub2_1, inpsub2_2,
                     inpH5, inpGene){
    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID,
                         inpConf[UI == inpdrY]$ID,
                         inpConf[UI == inpGrp]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "Y", "grpBy","sub1","sub2")
    #ad <- read_h5ad(inpad)
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=)))
    ggData[val1 < 0]$val1 = 0
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val2 < 0]$val2 = 0
    h5file$close_all()
    ggData = Subset(ggData, inpsub1_2, inpsub2_2)

    return(ggData)
}

scDRcoexPlot <- function(ggData, inpdrX, inpdrY,inp1, inp2,
                         inpsiz, inpcols, inpord, inpfsz, inpasp, inptxt, plot_brush){
    bgCells = TRUE
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

    # Generate coex color palette
    cInp = strsplit(inpcols, "; ")[[1]]
    if(cInp[1] == "Red (Gene1)"){
        c10 = c(255,0,0)
    } else if(cInp[1] == "Orange (Gene1)"){
        c10 = c(255,140,0)
    } else {
        c10 = c(0,255,0)
    }
    if(cInp[2] == "Green (Gene2)"){
        c01 = c(0,255,0)
    } else {
        c01 = c(0,0,255)
    }
    c00 = c(217,217,217) ; c11 = c10 + c01
    nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2
    gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1)))
    gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid
    gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid
    gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
    gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
    gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
    gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
    gg = gg[, c("v1", "v2", "cMix")]

    # Map colours
    ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1))
    ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2))
    ggData[is.na(ggData)] = 0
    ggData$v0 = ggData$v1 + ggData$v2
    ggData = gg[ggData, on = c("v1", "v2")]
    if(inpord == "Max-1st"){
        ggData = ggData[order(v0)]
    } else if(inpord == "Min-1st"){
        ggData = ggData[order(-v0)]
    } else if(inpord == "Random"){
        ggData = ggData[sample(nrow(ggData))]
    }

    # Actual ggplot
    ggOut = ggplot(ggData, aes(X, Y))
    if(bgCells){
        ggOut = ggOut +
            geom_point(data = ggData, color = "snow2", size = inpsiz, shape = 16)
    }
    ggOut = ggOut +
        geom_point(size = inpsiz, shape = 16, color = ggData$cMix) +
        xlab(inpdrX) + ylab(inpdrY) +
        sctheme(base_size = sList[inpfsz], XYval = inptxt) #+
    guides(color = guide_colorbar(barwidth = 15))

    if(is.null(plot_brush) & file.exists("tempData/plot_brush_tmp.rds")) {
        plot_brush = readRDS(file = "tempData/plot_brush_tmp.rds")
    }
    if(!is.null(plot_brush)) {
        ggOut = ggOut + rectangle(plot_brush$xmin,plot_brush$xmax, plot_brush$ymin, plot_brush$ymax,
                                  colour = "grey")
    }

    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }
    return(ggOut)
}



scDRcoexPlot3 <- function(ggData, inpdrX, inpdrY,inp1, inp2,
                          inpsiz, inpcols, inpord, inpfsz, inpasp, inptxt, plot_brush){
    bgCells = TRUE
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

    # Generate coex color palette
    cInp = strsplit(inpcols, "; ")[[1]] %>% gsub(" .*","",.) %>% tolower()
    gg = colorRamp2D(col1 = cInp[1], col2 = cInp[2],col3 = cInp[3], col0 = "snow2",nTot= 100)
    # Map colours
    nTot = 100
    ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1))
    ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2))
    ggData[is.na(ggData)] = 0
    ggData$v0 = ggData$v1 + ggData$v2
    ggData = gg[ggData, on = c("v1", "v2")]
    if(inpord == "Max-1st"){
        ggData = ggData[order(v0)]
    } else if(inpord == "Min-1st"){
        ggData = ggData[order(-v0)]
    } else if(inpord == "Random"){
        ggData = ggData[sample(nrow(ggData))]
    }

    # Actual ggplot
    ggOut = ggplot(ggData, aes(X, Y))
    if(bgCells){
        ggOut = ggOut +
            geom_point(data = ggData, color = "snow2", size = inpsiz, shape = 16)
    }
    ggOut = ggOut +
        geom_point(size = inpsiz, shape = 16, color = ggData$cMix) +
        xlab(inpdrX) + ylab(inpdrY) +
        sctheme(base_size = sList[inpfsz], XYval = inptxt)

    if(is.null(plot_brush) & file.exists("tempData/plot_brush_tmp.rds")) {
        plot_brush = readRDS(file = "tempData/plot_brush_tmp.rds")
    }
    if(!is.null(plot_brush)) {
        ggOut = ggOut + rectangle(plot_brush$xmin,plot_brush$xmax, plot_brush$ymin, plot_brush$ymax,
                                  colour = "grey")
    }

    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }
    return(ggOut)
}

# draw rectangle by four corner coordinates
#' @param x_left,x_right,y_bottom,y_top, coordinate values anti-clockwise from left bottom
#' @param ... argments pass to geom_segment
#' @example rectangle(-4.6, -4.1, -4.7, -5.5,colour = "blue")
rectangle <- function(x_left = -Inf, x_right = Inf, y_bottom = -Inf, y_top = Inf,linetype="dotted", ...){
    list(geom_segment(aes(x = x_left, xend = x_right, y = y_top, yend = y_top),...),
         geom_segment(aes(x = x_right, xend = x_right, y = y_top, yend = y_bottom),...),
         geom_segment(aes(x = x_left, xend = x_right, y = y_bottom, yend = y_bottom),...),
         geom_segment(aes(x = x_left, xend = x_left, y = y_top, yend = y_bottom),...))
}


scDRcoexScatterPlot <- function(sub_ggData, inpConf,#inpdrX, inpdrY,
                                inp1, inp2, inpGrp,
                                inpsiz, inpcols, inpcol2 ,inpfsz, inpasp, inptxt,inpyrev, density_type){

    rat = (max(sub_ggData$val1) - min(sub_ggData$val1)) / (max(sub_ggData$val2) - min(sub_ggData$val2))

    # Do factoring
    ggCol = strsplit(inpConf[UI == inpGrp]$fCL, "\\|")[[1]]
    names(ggCol) = levels(sub_ggData$grpBy)
    gglvl = levels(sub_ggData$grpBy)[levels(sub_ggData$grpBy) %in% unique(sub_ggData$grpBy)]
    gglvl = switch(inpyrev,"As-it-is" = gglvl, "Reverse" = rev(gglvl))

    sub_ggData$grpBy = factor(sub_ggData$grpBy, levels = gglvl)

    colnames(sub_ggData) %<>% sub("grpBy",inpGrp,.)


    # Main plot
    pmain <- ggplot(sub_ggData, aes_string(x = "val1", y = "val2",color = inpGrp))+
        geom_point(size = inpsiz, shape = 16) +
        xlab(inp1)+ ylab(inp2)+
        sctheme(base_size = sList[inpfsz], XYval = inptxt) +
        scale_color_manual("", values = switch(inpcol2,
                                               "default" =ggCol,
                                               color_generator(inpcol2,length(gglvl))))
    # Marginal densities along x axis
    geom_x <- switch (density_type,
                      "density" =  list(geom_density(data = sub_ggData, aes_string(x = "val1", fill = inpGrp),
                                                     alpha = 0.5, size = 0.4)),
                      "density_ridges"= list(ggridges::geom_density_ridges(data = sub_ggData, aes_string(x = "val1", y = inpGrp,fill = inpGrp),
                                                                           alpha = 0.5, size = 0.4),
                                             scale_y_discrete(expand = expansion(mult = c(0.01, 1))))
    )
    xdens <- do.call(what = '+', args = list(axis_canvas(pmain, axis = "x"), geom_x))+
        scale_fill_manual("", values = switch(inpcol2,
                                              "default" =ggCol,
                                              color_generator(inpcol2,length(gglvl))))


    # Marginal densities along y axis
    # Need to set coord_flip = TRUE, if you plan to use coord_flip()
    geom_y <- switch (density_type,
                      "density" =  list(geom_density(data = sub_ggData, aes_string(x = "val2", fill = inpGrp),
                                                     alpha = 0.5, size = 0.4)),
                      "density_ridges"= list(ggridges::geom_density_ridges(data = sub_ggData, aes_string(x = "val2", y = inpGrp,fill = inpGrp),
                                                                           alpha = 0.5, size = 0.4),
                                             scale_y_discrete(expand = expansion(mult = c(0.01, 1))))
    )

    ydens <- do.call(what = '+', args = list(axis_canvas(pmain, axis = "y", coord_flip = TRUE), geom_y))+
        coord_flip() +
        scale_fill_manual("", values = switch(inpcol2,
                                              "default" =ggCol,
                                              color_generator(inpcol2,length(gglvl))))
    p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
    p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

    ggOut = ggdraw(p2)
    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }
    return(ggOut)
}


scDRcoexLeg <- function(inp1, inp2, inpcols, inpfsz){
    # Generate coex color palette
    cInp = strsplit(inpcols, "; ")[[1]]
    if(cInp[1] == "Red (Gene1)"){
        c10 = c(255,0,0)
    } else if(cInp[1] == "Orange (Gene1)"){
        c10 = c(255,140,0)
    } else {
        c10 = c(0,255,0)
    }
    if(cInp[2] == "Green (Gene2)"){
        c01 = c(0,255,0)
    } else {
        c01 = c(0,0,255)
    }
    c00 = c(217,217,217) ; c11 = c10 + c01
    nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2
    gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1)))
    gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid
    gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid
    gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
    gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
    gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
    gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
    gg = gg[, c("v1", "v2", "cMix")]

    # Actual ggplot
    ggOut = ggplot(gg, aes(v1, v2)) +
        geom_tile(fill = gg$cMix) +
        xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) +
        scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) +
        scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) +
        sctheme(base_size = sList[inpfsz], XYval = TRUE)
    return(ggOut)
}


colorRamp2D <- function(col1 = "blue", col2 = "green",col3 = "red", col0 = "snow1",
                        nTot= 100){
    rotate <- function(x) t(apply(x, 2, rev))
    mm <- tcrossprod(seq(1,0,length.out = nTot))
    tmp1 <- sapply(col2rgb(col1)/255, function(x) 1-mm*(1-x))
    tmp2 <- sapply(col2rgb(col2)/255, function(x) 1-rotate(rotate(mm))*(1-x))
    tmp3 <- sapply(col2rgb(col3)/255, function(x) 1-rotate(mm)*(1-x))
    tmp4 <- sapply(col2rgb(col0)/255, function(x) 1-rotate(rotate(rotate(mm)))*(1-x))

    tmp <- tmp1*tmp2*tmp3*tmp4
    gg <- data.table(v1 = rep((nTot-1):0, times= nTot),
                     v2 = rep(0:(nTot-1),  each= nTot),
                     cMix = rgb(tmp))
    return(gg)
}


# co-expression legend with 3 and more colors
# https://stackoverflow.com/questions/11070101/2d-color-gradient-plot-in-r/11103414#11103414
#' @example scDRcoexLeg3(inpcol= "Orange (Gene1); Blue (Gene2); Red (Both)")
scDRcoexLeg3 <- function(inp1 = "SLURP1", inp2 = "HBA1",inpcols, inpfsz,nTot= 100){
    cInp = strsplit(inpcols, "; ")[[1]] %>% gsub(" .*","",.) %>% tolower()

    gg = colorRamp2D(col1 = cInp[1], col2 = cInp[2],col3 = cInp[3], col0 = "snow2",nTot= nTot)
    ggOut = ggplot(gg, aes(v1, v2)) +
        geom_tile(fill = gg$cMix) +
        xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) +
        scale_x_continuous(breaks = c(0, nTot-1), label = c("low", "high")) +
        scale_y_continuous(breaks = c(0, nTot-1), label = c("low", "high")) +
        sctheme(base_size = sList[inpfsz], XYval = TRUE)
    return(ggOut)
}


scDRcoexNum <- function(ggData, inp1, inp2){
    if(nrow(ggData) == 0) return(data.table(`expression > 0`="none", nCells=0, percent=0))
    # Actual data.table
    ggData$express = "none"
    ggData[val1 > 0]$express = inp1
    ggData[val2 > 0]$express = inp2
    ggData[val1 > 0 & val2 > 0]$express = "both"
    ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
    ggData = ggData[, .(nCells = .N), by = "express"]
    ggData$percent = 100 * ggData$nCells / sum(ggData$nCells)
    ggData = ggData[order(express)]
    colnames(ggData)[1] = "expression > 0"
    return(ggData)
}


scDRcoexCor <- function(ggData, inp1, inp2){
    if(nrow(ggData) <= 4) return("must have >4 observations")
    # calculate correlation
    res = Hmisc::rcorr(x = ggData[,val1],y = ggData[,val2],type = "spearman")
    rcorr_res = c(round(res$r[2,1], digits = 5),
                  round(res$P[2,1], digits = 5))
    names(rcorr_res) = c("correlation","cor_p-value")
    # Actual data.table
    ggData$express = "none"
    ggData[val1 > 0]$express = inp1
    ggData[val2 > 0]$express = inp2
    ggData[val1 > 0 & val2 > 0]$express = "both"
    ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
    ggData = ggData[, .(nCells = .N), by = "express"]
    ggData$percent = 100 * ggData$nCells / sum(ggData$nCells)
    ggData = ggData[order(express)]
    colnames(ggData)[1] = "expression > 0"

    # calculate fisher.test
    if(nrow(ggData) <= 2) {
        FISH_res = c(NaN, NaN)
        names(FISH_res) = c("odds ratio","fisher_p-value")
    } else {
        if(nrow(ggData) == 3){
            if(!"none" %in% ggData$`expression > 0`) ggData = rbind(ggData,
                                                                    data.table(`expression > 0`="none", nCells=0, percent=0))
            if(!"both" %in% ggData$`expression > 0`) ggData = rbind(data.table(`expression > 0`="both", nCells=0, percent=0),
                                                                    ggData)
        }
        FISH = fisher.test(matrix(ggData[,nCells],nrow =2))
        FISH_res = c(round(FISH$estimate, digits = 5),
                     round(FISH$p.value, digits = 5))
        names(FISH_res)[2] ="fisher_p-value"
    }

    return(c(rcorr_res,FISH_res))
}


## Create a table for printing out cell name and gene expression values
#' @param plot_brush Rshiny input$plot_brush
tableBrush <- function(ggData, plot_brush){
    area = FALSE
    if(is.null(plot_brush) &
       file.exists("tempData/plot_brush_tmp.rds")) {
        plot_brush = readRDS(file = "tempData/plot_brush_tmp.rds")
    }
    if(!is.null(plot_brush)) {
        area = (ggData[,X] > plot_brush$xmin) &
            (ggData[,X] < plot_brush$xmax) &
            (ggData[,Y] > plot_brush$ymin) &
            (ggData[,Y] < plot_brush$ymax)
        # save corrdicates inputs
        if(!dir.exists("tempData")) dir.create("tempData")
        saveRDS(plot_brush,file = "tempData/plot_brush_tmp.rds")
        area = (ggData[,X] > plot_brush$xmin) &
            (ggData[,X] < plot_brush$xmax) &
            (ggData[,Y] > plot_brush$ymin) &
            (ggData[,Y] < plot_brush$ymax)
    }
    return(ggData[area])
}


# a custom table container
sketch = htmltools::withTags(table(
    class = 'display',
    thead(
        tr(
            th(colspan = 2, 'Spearman Correlation'),
            th(colspan = 2, "Fisher's Exact Test")
        ),
        tr(
            lapply(c("Correlation","p-val","odds ratio","p-val"), th)
        )
    )
))

# @author jan-glx on StackOverflow
# @references \url{https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2}
# @seealso \code{\link[ggplot2]{geom_violin}}
#
geom_split_violin <- function(
    mapping = NULL,
    data = NULL,
    stat = 'ydensity',
    position = 'identity',
    ...,
    draw_quantiles = NULL,
    trim = TRUE,
    scale = 'area',
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
    return(layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomSplitViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            trim = trim,
            scale = scale,
            draw_quantiles = draw_quantiles,
            na.rm = na.rm,
            ...
        )
    ))
}
# A split violin plot geom
#
#' @importFrom scales zero_range
#' @importFrom ggplot2 GeomPolygon
#' @importFrom grid grobTree grobName
#
# @author jan-glx on StackOverflow
# @references \url{https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2}
# @seealso \code{\link[ggplot2]{geom_violin}}
#
GeomSplitViolin <- ggproto(
    "GeomSplitViolin",
    GeomViolin,
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
        data$xminv <- data$x - data$violinwidth * (data$x - data$xmin)
        data$xmaxv <- data$x + data$violinwidth * (data$xmax - data$x)
        grp <- data[1, 'group']
        if (grp %% 2 == 1) {
            data$x <- data$xminv
            data.order <- data$y
        } else {
            data$x <- data$xmaxv
            data.order <- -data$y
        }
        newdata <- data[order(data.order), , drop = FALSE]
        newdata <- rbind(
            newdata[1, ],
            newdata,
            newdata[nrow(x = newdata), ],
            newdata[1, ]
        )
        newdata[c(1, nrow(x = newdata) - 1, nrow(x = newdata)), 'x'] <- round(x = newdata[1, 'x'])
        grob <- if (length(x = draw_quantiles) > 0 & !zero_range(x = range(data$y))) {
            stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
            quantiles <- QuantileSegments(data = data, draw.quantiles = draw_quantiles)
            aesthetics <- data[rep.int(x = 1, times = nrow(x = quantiles)), setdiff(x = names(x = data), y = c("x", "y")), drop = FALSE]
            aesthetics$alpha <- rep.int(x = 1, nrow(x = quantiles))
            both <- cbind(quantiles, aesthetics)
            quantile.grob <- GeomPath$draw_panel(both, ...)
            grobTree(GeomPolygon$draw_panel(newdata, ...), name = quantile.grob)
        }
        else {
            GeomPolygon$draw_panel(newdata, ...)
        }
        grob$name <- grobName(grob = grob, prefix = 'geom_split_violin')
        return(grob)
    }
)

# prepare data for violin / boxplot
scVioBoxData <- function(inpConf, inpMeta, inp1, inp2, inp3,
                         inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                         inpH5, inpGene){
    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID,
                         inpConf[grp == TRUE]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID),
                     with = FALSE]
    colnames(ggData)[1] = c("X")
    N = ncol(ggData)
    colnames(ggData)[(N-2):N] = c("sub1","sub2","sub3")
    # Load in either cell meta or gene expr
    if(inp3 %in% inpConf$UI) colnames(ggData) %<>% sub(inp3,"grpBy",.)
    if(inp2 %in% inpConf$UI){
        ggData$val = inpMeta[[inpConf[UI == inp2]$ID]]
    } else {
        h5file <- H5File$new(inpH5, mode = "r")
        h5data <- h5file[["grp"]][["data"]]
        ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=)))
        ggData[val < 0]$val = 0
        set.seed(101)

        h5file$close_all()
    }
    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
    ggData = ggData[, c("X","grpBy","val")]

    return(ggData)
}

scVioBoxDataNum <- function(ggData,inpsplit){
    ggData_num = switch(inpsplit,
                   "stack" =  group_by(ggData,X),
                   "split" =  group_by(ggData,X, grpBy)) %>%
        summarise(cell.number = n(), mean = mean(val), sd = sd(val))

    return(ggData_num)
}

# Plot violin / boxplot
scVioBox <- function(ggData, inpConf, inp1, inp2, inp3,
                     inptyp, inpsplit, inpsig,
                     inpcols, inppts,inpvioscale,
                     inperr,inpsiz, inpfsz,inpylim,
                     plabel, pcut, inppvalpos,
                     inpfrt, inpleg,inplegpos = "bottom"){

    # Do factoring
    temp1 <- doFactoring(ggData,inpConf, col = "X", inp1, inpcols)
    temp3 <- doFactoring(temp1$ggData,inpConf, col = "grpBy", inp3, inpcols)

    gglvl1 = temp1$gglvl; ggCol1 = temp1$ggCol
    gglvl3 = temp3$gglvl; ggCol3 = temp3$ggCol

    ggData = temp3$ggData

    #if data missing,  add rows with 0
    #Count = switch(inpsplit,
    #               "stack" =  group_by(ggData,X),
    #               "split" =  group_by(ggData,X, grpBy)) %>% summarise(n = n())

    #if(inpsig & any(Count$n < 3)){
    #    miss_data = switch(inpsplit,
    #                       "stack" =  as.data.frame(table(ggData$X)),
    #                       "split" =  as.data.frame(table(ggData$X,ggData$grpBy)))
    #    miss_data %<>% filter(Freq < 3)
    #    colnames(miss_data) = c("X","grpBy","val")
    #    miss_data$val = 3- miss_data$val
    #    miss_data = miss_data[rep(seq(nrow(miss_data)), miss_data$val),]
    #    miss_data$val = 0
    #    ggData %<>% rbind(miss_data)
        #ggData %<>% split(f = ggData$X) %>% lapply(function(obj) obj[order(obj$grpBy,obj$paired),])
        #ggData = ggData[gglvl1]
        #ggData %<>% bind_rows()
    #}

    # Actual ggplot
    ggOut = switch (inptyp,
                    "violin" = ggplot(ggData, aes(X, val)) +
                                        scale_fill_manual(values=switch(inpsplit,
                                                                        "stack"= ggCol1[gglvl1],
                                                                        "split" = ggCol3[gglvl3]))+
                                        xlab(paste0("\n",inp1))+
                                        ylab(paste0(inp2,"\n"))+
                                    switch(as.character(inpsplit == "split" & length(gglvl3) == 2),
                                           # https://github.com/kassambara/ggpubr/issues/266
                                           # only want to change the fill color for this aesthetic in the stats results
                                           "TRUE" = geom_split_violin(inherit.aes = FALSE,
                                                                      aes(X, val, fill = switch(inpsplit,
                                                                                        "stack"= X,
                                                                                        "split" = grpBy))),
                                           "FALSE" = geom_violin(inherit.aes = FALSE,scale = inpvioscale,
                                                                 aes(X, val, fill = switch(inpsplit,
                                                                                   "stack"= X,
                                                                                   "split" = grpBy)))),
                    "boxplot" =  ggboxplot(data = ggData,
                                           x = "X",
                                           y = "val",
                                           fill = switch(inpsplit,
                                                         "stack"= "X",
                                                         "split" = "grpBy"),
                                           #shape = ifelse(inppts,19,NA),
                                           palette = switch(inpsplit,
                                                            "stack"= ggCol1[gglvl1],
                                                            "split" = ggCol3[gglvl3]),
                                           xlab = paste0("\n",inp1),
                                           ylab = paste0(inp2,"\n"),
                                           outlier.shape = NA,
                                           line.color = "gray",
                                           bxp.errorbar = inperr),
                    "barplot" =  ggbarplot(data = ggData,
                                           x = "X",
                                           y = "val",
                                           fill =  switch(inpsplit,
                                                          "stack"= "X",
                                                          "split" = "grpBy"),
                                           position = switch(inpsplit,
                                                             "stack"= position_stack(),
                                                             "split" = position_dodge()),
                                           palette = switch(inpsplit,
                                                            "stack"= ggCol1[gglvl1],
                                                            "split" = ggCol3[gglvl3]),
                                           add = c(ifelse(inperr,"mean_sd","mean")),
                                           xlab = paste0("\n",inp1),
                                           ylab = paste0(inp2,"\n"))
                    )
    if(inpsig){
        stat.test <- switch(inpsplit,
                            "stack" = ggData,
                            "split" = group_by(ggData, X)) %>%
            wilcox_test(as.formula(paste("val","~",
                                         switch(inpsplit,
                                                "stack" = "X",
                                                "split" = "grpBy"))),
                        p.adjust.method = "bonferroni") %>%
            add_significance("p") %>%
            rstatix:::remove_ns(col = plabel,
                                signif.cutoff = switch(plabel,
                                                       "p" =  10^(-pcut),
                                                       "p.adj" =  10^(-pcut),
                                                       "p.signif" = "ns",
                                                       "p.adj.signif" = "ns")) %>%
            add_xy_position(x = "X", dodge = 0.8,step.increase = 0.02 * (1+inpylim),
                            scales = "free",
                            stack = FALSE)#ifelse(inpsplit == "stack",TRUE,FALSE) )
        stat.test$y.position = stat.test$y.position + (stat.test$y.position - min(stat.test$y.position)) * inpylim +
            mean(stat.test$y.position) * inppvalpos
        ggOut = ggOut + stat_pvalue_manual(
            stat.test,  label = plabel, tip.length = 0.01,size = inpsiz*5)
    }

    if(inppts) ggOut = ggOut + geom_jitter(size = ifelse(inppts,inpsiz,NA), shape = 19)
    if(inpylim != 0) {
        Ylim = max(ggOut$data$val * (1+inpylim), ifelse(inpsig,max(stat.test$y.position),0),na.rm = TRUE)
        ggOut = ggOut + ylim(NA,Ylim)
    }
    ggOut = ggOut +
        sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
        theme(legend.title = element_blank(),
              legend.position = switch(as.character(inpleg),
                                       "FALSE" = "none",
                                       "TRUE" = inplegpos),
              axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                            hjust = switch(as.character(inpfrt),
                                                           "0" = 0.5,
                                                           "30" = 1,
                                                           "45" = 1,
                                                           "90" = 0
                                            ),
                                            vjust = switch(as.character(inpfrt),
                                                           "0" = 0,
                                                           "30" = 1,
                                                           "45" = 1,
                                                           "90" = 0.5)
              )
        )
    return(ggOut)
}

# Plot proportion plot
scPropData <- function(inpConf, inpMeta, inp1, inp2,
                       inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                       inpnorm){

    # Prepare ggData
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID,
                         inpConf[UI == inp2]$ID,
                         inpConf[UI == "patient"]$ID,
                         inpConf[UI == "orig.ident"]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "grpBy","paired","orig.ident","sub1", "sub2","sub3")

    ggData = switch(inpnorm,
                    "Cell % by sample" = {
                        ggData = Subset(ggData, inpsub1_2)
                        ggData_total = ggData[, .(SampleCells = .N), by = "orig.ident"]
                        ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
                        ggData = ggData[, .(nCells = .N), by = c("X", "grpBy","paired","orig.ident")]
                        ggData[["SampleCells"]] = ggData_total[["SampleCells"]][match(ggData[["orig.ident"]],
                                                                                      ggData_total[["orig.ident"]])]
                        ggData[["pctCells"]] = ggData[["nCells"]] /ggData[["SampleCells"]] *100
                        ggData
                    },
                    "Cell % by Y-axis" = {
                        ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
                        ggData_total = ggData[, .(SampleCells = .N)]
                        ggData = ggData[, .(nCells = .N), by = c("X","grpBy")]
                        ggData[["pctCells"]] = ggData[["nCells"]] /ggData_total[["SampleCells"]] *100
                        ggData
                    },
                    {
                        #ggData_total = ggData[, .(SampleCells = .N), by = "X"]
                        #ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
                        #ggData = ggData[, .(nCells = .N), by = c("X", "grpBy","paired","orig.ident")]
                        #ggData[["SampleCells"]] = ggData_total[["SampleCells"]][match(ggData[["X"]],
                        #                                                              ggData_total[["X"]])]
                        #ggData[["pctCells"]] = ggData[["nCells"]] /ggData[["SampleCells"]] *100
                        #ggData

                        ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
                        ggData = ggData[, .(nCells = .N), by = c("X", "grpBy")]
                        ggData = ggData[, {tot = sum(nCells)
                        .SD[,.(pctCells = 100 * sum(nCells) / tot,
                               nCells = nCells), by = c("grpBy")]}, by = c("X")]
                        ggData
                    })
    return(ggData)
}


# Plot proportion plot
scProp <- function(ggData, inpConf, inp1, inp2,
                   inptyp, inpsplit, inpnorm, inpsig, inpAddLine, inpcols,inpyrange, inpflpxy,
                   inpflpx, inppts, inpsiz, inpfilter,
                   inpfsz, inplsz,inplab, inperr,
                   plabel, pcut, inppvalpos,
                   inpfrt, inpleg, inplegpos = "bottom"){

    inppair = grepl("^paired",inpsig)
    # Actual ggplot

    # Do factoring
    temp1 <- doFactoring(ggData,inpConf, col = "X", inp1, inpcols)
    temp2 <- doFactoring(temp1$ggData,inpConf, col = "grpBy", inp2, inpcols)

    gglvl1 = temp1$gglvl; ggCol1 = temp1$ggCol
    gglvl2 = temp2$gglvl; ggCol2 = temp2$ggCol

    ggData = temp2$ggData

    if(inptyp %in% c("barplot","boxplot","piechart","pairedplot")){

        ggOut = switch (inptyp,
                        "barplot" =  ggbarplot(data = ggData,
                                               x = "X",
                                               y = ifelse(inpnorm == "Cell Number","nCells","pctCells"),
                                               fill = "grpBy",
                                               position = switch(inpsplit,
                                                                  "stack"= position_stack(),
                                                                  "split" = position_dodge()),
                                               palette = ggCol2[gglvl2],
                                               add = c(ifelse(inperr & inpsplit == "split","mean_sd","mean")),
                                               ylab = ifelse(inpnorm == "Cell Number","Average Cell Number","Average Cell Proportion (%)"),
                                               xlab = inp1)+
                            geom_jitter(size = ifelse(inppts,inpsiz,NA), shape = 19),
                         "boxplot" =  ggboxplot(data = ggData,
                                                x = "X",
                                                y = ifelse(inpnorm == "Cell Number","nCells","pctCells"),
                                                fill = switch(inpsplit,
                                                              "stack"= "X",
                                                              "split" = "grpBy"),
                                                shape = ifelse(inppts,19,NA),
                                                palette = switch(inpsplit,
                                                                 "stack"= ggCol1[gglvl1],
                                                                 "split" = ggCol2[gglvl2]),
                                                ylab = ifelse(inpnorm == "Cell Number","Average Cell Number","Average Cell Proportion (%)"),
                                                xlab = inp1,
                                                outlier.shape = NA,
                                                line.color = "gray",
                                                bxp.errorbar = inperr)+
                            geom_jitter(size = ifelse(inppts,inpsiz,NA), shape = 19),
                        "piechart" =  {
                            ggData[["value"]] = switch(inpnorm,
                                                       "Cell Number" = ggData[["nCells"]],
                                                       round(ggData[["pctCells"]], digits = ifelse(inpfilter < 0.15,2,1)))
                            ggData %<>% group_by(X) %>%
                                summarise(value = sum(value))

                            #https://r-charts.com/part-whole/pie-chart-labels-outside-ggplot2/
                            df2 <- ggData %>%
                                mutate(csum = rev(cumsum(rev(value))),
                                       pos = value/2 + lead(csum, 1),
                                       pos = if_else(is.na(pos), value/2, pos)) %>%
                                filter(value > inpfilter)
                            df2$value %<>% paste0(ifelse(inpnorm == "Cell Number","","%"))
                            df2$label = paste(df2$X,df2$value)

                            ggplot(ggData, aes(x = "" , y = value, fill = X)) +
                                geom_col(width = 1, color = 1) +
                                coord_polar(theta = "y") +
                                xlab("") +ylab("")+
                                scale_fill_manual(name = "",
                                                  values = ggCol1[gglvl1]) +
                                switch(inplab,
                                       "No %" = NULL,
                                       "No labels"  = geom_text_repel(data = df2, aes(y = pos, label = value),
                                                                      color = "grey10",force = 3,
                                                                      size = lList[inpfsz], nudge_x = 1, nudge_y = 1, show.legend = FALSE),
                                       "black text" = geom_text_repel(data = df2, aes(y = pos, label = label),
                                                                      color = "grey10",force = 3,
                                                                      size = lList[inpfsz], nudge_x = 1, nudge_y = 1, show.legend = FALSE),
                                       "color labels" = geom_label_repel(data = df2, aes(y = pos, label = label),force = 3,
                                                                         size = lList[inpfsz], nudge_x = 1, nudge_y = 1, show.legend = FALSE)
                                )
                            },
                         "pairedplot" = {
                             y = ifelse(inpnorm == "Cell Number","nCells","pctCells")
                             not_y = ifelse(inpnorm == "Cell Number","pctCells","nCells")
                             ggData1 = ggData[, .(pctCells = 100 * sum(nCells) / SampleCells,
                                        nCells = sum(nCells)),
                                    by = c("X","paired","SampleCells")] %>%
                             pivot_wider(!c("SampleCells",not_y),
                                         names_from = "paired",values_from = "pctCells", values_fill = 0)  %>%
                            tibble::column_to_rownames("X") %>% t() %>% as.data.frame() %>% 
                             #pivot_longer(!X,names_to = "paired", values_to = "pctCells") %>% 
                             #pivot_wider(names_from = "X",values_from = "pctCells") %>%
                             dplyr::mutate("is_increasing" = ifelse(.[,(levels(ggData$X)[2])]> .[,(levels(ggData$X)[1])],
                                                                    "Increase","Decrease")) %>%
                                 tibble::rownames_to_column("patient") %>%
                                 pivot_longer(!c("is_increasing","patient"),names_to = inp1, values_to = y)
                             ggData1$is_increasing %<>% factor(levels = c("Increase","Decrease"))
                           #https://www.biostars.org/p/306348/
                             ggplot(ggData1, aes_string(x = inp1, y = y)) +
                                 geom_boxplot(aes_string(fill = inp1), alpha = 0.2, col = "grey")+
                                 geom_point(aes(color = patient),size = inpsiz*3)+
                                scale_colour_manual(values = c(color_generator(inpcols,length(unique(ggData$paired))),
                                                               color_generator("JCO")[c(1,4)]))+
                                 geom_line(aes(group = patient, color = is_increasing),
                                           size = ifelse(inpAddLine,sList[inplsz]/50,0))+
                                 
                                 ylab(ifelse(inpnorm == "Cell Number","Number of Cells","Proportion of Cancer-associated TCRs (%)"))+
                                 xlab(inp1)+
                                 stat_compare_means(method = "t.test", paired = TRUE,
                                                    size = sList[inpfsz]/4,step.increase = 0.5)
                             }
                        )
        #if data missing,  add rows with 0

        if(grepl("Wilcoxon",inpsig)){
            stat.test <- switch(inpsplit,
                                "stack" = ggData,
                                "split" = group_by(ggData, X)) %>%
                wilcox_test(as.formula(paste(y,"~",
                                             switch(inpsplit,
                                                    "stack" = "X",
                                                    "split" = "grpBy"))),
                            p.adjust.method = "bonferroni",
                            paired = inppair) %>%
                add_significance("p") %>%
                rstatix:::remove_ns(col = plabel,
                          signif.cutoff = switch(plabel,
                                                 "p" =  10^(-pcut),
                                                 "p.adj" =  10^(-pcut),
                                                 "p.signif" = "ns",
                                                 "p.adj.signif" = "ns")) %>%
                add_xy_position(x = "X", dodge = 0.8,step.increase = 0.12,
                                           scales = "fixed",
                                           stack = FALSE)#ifelse(inpsplit == "stack",TRUE,FALSE) )
            stat.test$y.position = stat.test$y.position + mean(stat.test$y.position) * inppvalpos
            ggOut = ggOut + stat_pvalue_manual(
                stat.test,  label = plabel, tip.length = 0.01,size = inpsiz*3)
        }

        if(grepl("chisq",inpsig)){
            if(inpnorm == "Cell % by sample"){
                df_prob <-  ggData %>%
                    group_by(X, grpBy) %>%
                    summarise(sum_SampleCells = sum(SampleCells), .groups = "drop") %>%
                    pivot_wider(names_from = "grpBy",values_from = "sum_SampleCells", values_fill = 0) %>%
                    tibble::column_to_rownames("X")
                prob <- rep(rowSums(df_prob)/2, ncol(df_prob))/(rowSums(df_prob)*ncol(df_prob))
            } else  {
                df = as.matrix(table(ggData$X, ggData$grpBy))
                prob = rep(1/2/ncol(df), 2*ncol(df))
            }
            stat.test <- ggData %>%
                group_by(X, grpBy) %>%
                summarise(sum_nCells = sum(nCells), .groups = "drop") %>%
                pivot_wider(names_from = "grpBy",values_from = "sum_nCells", values_fill = 0) %>% 
                tibble::column_to_rownames("X") %>% 
                chisq_test(p = prob, rescale.p = TRUE,simulate.p.value = TRUE) %>%
            add_significance("p") %>%
            rstatix:::remove_ns(col = plabel,
                                signif.cutoff = switch(plabel,
                                                       "p" =  10^(-pcut),
                                                       "p.adj" =  10^(-pcut),
                                                       "p.signif" = "ns",
                                                       "p.adj.signif" = "ns"))
            stat.test$y.position =  group_by(ggData,grpBy) %>%
                summarise(mean = mean(ifelse(inpnorm == "Cell Number",nCells,pctCells))) %>% .[,"mean"] %>% sum
            stat.test$group1 = levels(ggData$X)[1]
            stat.test$group2 = levels(ggData$X)[2]
            stat.test$y.position = stat.test$y.position + max(stat.test$y.position) * inppvalpos
            ggOut = ggOut + stat_pvalue_manual(
                stat.test,  label = plabel, tip.length = 0.01,size = inpsiz*3)
        }
    if(inpyrange != "1"){
        yrange  = layer_scales(ggOut)$y$range$range
        yrange[2] = yrange[2]* (1+as.numeric(inpyrange))
        ggOut = ggOut + ylim(yrange)
    }    
    if(inpflpxy){
        ggOut = ggOut + coord_flip()
    }
    if(inpflpx){
        ggOut = ggOut +  scale_x_discrete(limits = rev)
    }
    } else if(inptyp == "scatterplot"){
        ggData$X %<>% droplevels()
        ggData$index = plyr::mapvalues(ggData$X,from = levels(ggData$X),to = seq_len(length(levels(ggData$X))))
        ggData$index %<>% as.integer()

        colnames(ggData) %<>% sub("grpBy",inp2,.)
        ggOut = ggscatter(data = ggData,
                          x = "index",
                          y = ifelse(inpnorm == "Cell Number","nCells","pctCells"),
                          color = inp2,
                          #shape = ifelse(inppts,19,NA),# Extending the regression line
                          palette = ggCol2[gglvl2],
                          #cor.method = "spearman",
                          #cor.coef = TRUE,
                          #add = ifelse(inpAddLine,"reg.line","none"),
                          ylab = ifelse(inpnorm == "Cell Number","Number of Cells","Cell Proportion (%)"),
                          size = sList[inplsz]/5)+
            scale_x_continuous(name =inp1,
                               breaks = seq_len(length(levels(ggData$X))),
                               labels = levels(ggData$X))
        if(inpAddLine) ggOut = ggOut+ geom_line(aes_string(group = inp2, colour = inp2))

        if(inpflpxy){
            ggOut = ggOut + coord_flip()
        }
        if(inpflpx){
            ggOut = ggOut +  scale_x_continuous(limits = rev)
        }
    }

    ggOut = ggOut +
        cowplot::theme_cowplot(font_size = sList[inpfsz],font_family = "Avenir")+
        theme(legend.title = element_blank(),
              legend.position = switch(as.character(inpleg),
                                       "FALSE" = "none",
                                       "TRUE" = inplegpos)
              )
    if(inptyp == "piechart") {
        ggOut = ggOut +  theme(axis.line =  element_blank(),
                                                    axis.text.x = element_blank())

    } else {
        ggOut = ggOut +  theme(
                    axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                                  hjust = switch(as.character(inpfrt),
                                                                 "0" = 0.5,
                                                                 "30" = 1,
                                                                 "45" = 1,
                                                                 "90" = 0
                                                  ),
                                                  vjust = switch(as.character(inpfrt),
                                                                 "0" = 0,
                                                                 "30" = 1,
                                                                 "45" = 1,
                                                                 "90" = 0.5)))
    }
    return(ggOut)
}


# Get gene list
scGeneList <- function(inp, inpGene){
    geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])),
                          present = TRUE)
    geneList[!gene %in% names(inpGene)]$present = FALSE
    return(geneList)
}

# prepare data for scBubbHeat
scBubbHeatData <- function(inpConf, inpMeta, inp, inpGrp, inpPlt,
                           inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                           inpH5, inpGene, inpScl, inpRow, inpCol){
    # Identify genes that are in our dataset
    geneList = scGeneList(inp, inpGene)
    geneList = geneList[present == TRUE]
    shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
    shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

    # Prepare ggData
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData = inpMeta[, c("sampleID",
                         inpConf[grp == TRUE]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID), with = FALSE]
    N = ncol(ggData)
    colnames(ggData)[(N-2):N] = c("sub1","sub2","sub3")
    ggData$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    Colnames = colnames(ggData)
    for(iGene in geneList$gene){
        ggData[[iGene]] = h5data[inpGene[iGene],]
    }
    h5file$close_all()

    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
    ggData$grpBy %<>% droplevels()
    # Aggregate
    ggData = pivot_longer(ggData, !Colnames, names_to = "geneName", values_to = "val") %>%
        as.data.table()

    ggData$val = expm1(ggData$val)

    ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)),
                    by = c("geneName", "grpBy")]
    ggData$val = log1p(ggData$val)

    # Scale if required
    colRange = range(ggData$val)
    validate(
        need(!anyNA(colRange), paste(paste(ggData[["geneName"]][is.na(ggData$val)],collapse = ","),
                                     "has/have no expression!"))
    )
    if(inpScl){
        ggData[, val:= scale(val), keyby = "geneName"]
        colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
    }

    # hclust row/col if necessary
    ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
    tmp = ggMat$geneName
    ggMat = as.matrix(ggMat[, -1])
    rownames(ggMat) = tmp
    if(inpRow){
        hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat))))
        ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
    } else {
        ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
    }
    if(inpCol){
        hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
        ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
    }
    return(ggData)
}

# Plot gene expression bubbleplot / heatmap
scBubbHeat <- function(ggData,inp, inpGene,inpRow, inpCol,inpPlt,
                       inpcols,inpcolinv, inpfsz, inpfrt, inpsiz, inpasp,
                       inpleg, save = FALSE){
    # Identify genes that are in our dataset
    geneList = scGeneList(inp, inpGene)
    geneList = geneList[present == TRUE]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))

    # hclust row/col if necessary
    ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
    tmp = ggMat$geneName
    ggMat = as.matrix(ggMat[, -1])
    rownames(ggMat) = tmp
    if(inpRow){
        hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat))))
        ggRow = ggplot() + coord_flip() +
            geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
                               labels = unique(ggData$grpBy), expand = c(0, 0)) +
            scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                               labels = hcRow$labels$label, expand = c(0, 0.5)) +
            sctheme(base_size = sList[inpfsz]) +
            theme(axis.title = element_blank(), axis.line = element_blank(),
                  axis.ticks = element_blank(), axis.text.y = element_blank(),
                  axis.text.x = element_text(color="white", angle = 45, hjust = 1,vjust=0))
        ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
    } else {
        ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
    }
    if(inpCol){
        hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
        ggCol = ggplot() +
            geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
            scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                               labels = hcCol$labels$label, expand = c(0.05, 0)) +
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
                               labels = unique(ggData$geneName), expand=c(0,0)) +
            sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            theme(axis.title = element_blank(), axis.line = element_blank(),
                  axis.ticks = element_blank(), axis.text.x = element_blank(),
                  axis.text.y = element_text(color = "white"))
        ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
    }
    # Actual plot according to plottype
    cols =  color_generator(inpcols)
    if(inpcolinv) cols = rev(cols)
    if(inpPlt == "Bubbleplot"){
        # Bubbleplot
        ggOut = ggplot(ggData, aes(grpBy, geneName,fill = val, size = prop)) +
            geom_point(color = "black",pch=21) +
            scale_x_discrete(expand = c(0.05, 0)) +
            scale_y_discrete(expand = c(0, 0.5)) +
            scale_size_continuous("proportion", range = c(0, inpsiz),
                                  limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
            scale_fill_gradientn("expression", limits = colRange, colours = cols) +
            cowplot::theme_cowplot(font_size = sList[inpfsz])+
            #sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            guides(color = guide_colorbar(barwidth = 10)) +
            theme(axis.title = element_blank(), legend.box = "horizontal",
                  axis.text.x =    element_text(angle = as.numeric(inpfrt), hjust = 1,vjust = 0.5))
    } else {
        # Heatmap
        ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) +
            geom_tile() +
            scale_x_discrete(expand = c(0.05, 0)) +
            scale_y_discrete(expand = c(0, 0.5)) +
            scale_fill_gradientn("expression", limits = colRange, colours = cols) +
            cowplot::theme_cowplot(font_size = sList[inpfsz])+
            #sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            guides(fill = guide_colorbar(barwidth = 10)) +
            theme(axis.title = element_blank(), legend.box = "horizontal",
                  axis.text.x =    element_text(angle = as.numeric(inpfrt), hjust = 1,vjust = 0.5))
    }
    rat = (length(unique(ggData$grpBy)) + as.integer(inpRow)) / (length(unique(ggData$geneName)) + as.integer(inpCol))
    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }

    # Final tidy
    if(inpleg) {ggLeg = g_legend(ggOut)} else {
        ggLeg = ggplot() + geom_point() +theme_cowplot(line_size = 0)
    }
    ggOut = ggOut + theme(legend.position = "none",
                          axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                                        hjust = switch (as.character(inpfrt),
                                                                        "0" = 0.5,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0
                                                        ),
                                                        vjust = switch (as.character(inpfrt),
                                                                        "0" = 0,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0.5)
                          )
    )
    if(!save){
        if(inpRow & inpCol){ggOut =
            grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                         layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
        } else if(inpRow){ggOut =
            grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                         layout_matrix = rbind(c(1,3),c(2,NA)))
        } else if(inpCol){ggOut =
            grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                         layout_matrix = rbind(c(3),c(1),c(2)))
        } else {ggOut =
            grid.arrange(ggOut, ggLeg, heights = c(7,2),
                         layout_matrix = rbind(c(1),c(2)))
        }
    } else {
        if(inpRow & inpCol){ggOut =
            arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                        layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
        } else if(inpRow){ggOut =
            arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                        layout_matrix = rbind(c(1,3),c(2,NA)))
        } else if(inpCol){ggOut =
            arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                        layout_matrix = rbind(c(3),c(1),c(2)))
        } else {ggOut =
            arrangeGrob(ggOut, ggLeg, heights = c(7,2),
                        layout_matrix = rbind(c(1),c(2)))
        }
    }
    return(ggOut)
}


# Plot gene expression bubbleplot / heatmap
scCorBubbHeat <- function(inpConf, inpMeta, inpdf, inpPlt,
                       inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,
                       inpH5, inpGene, inpScl, inpRow, inpCol,
                       inpcols,inpcolinv, inpfsz, inpfrt, inpsiz, inpasp,
                       inpleg, save = FALSE){
    # Identify genes that are in our dataset
    geneList = unique(inpdf$gene[inpdf$gene %in% names(inpGene)])
    shiny::validate(need(length(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
    shiny::validate(need(length(geneList) > 1, "Please input at least 2 genes to plot!"))

    # Prepare ggData
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData = data.table()
    for(iGene in geneList){
        tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1_1]$ID, inpConf[UI == inpsub2_1]$ID), with = FALSE]
        N = ncol(tmp)
        colnames(tmp)[(N-1):N] = c("sub1","sub2")
        tmp$geneName = iGene
        tmp$val = h5data[inpGene[iGene],]
        ggData = rbindlist(list(ggData, tmp))
    }
    h5file$close_all()

    ggData = Subset(ggData, inpsub1_2, inpsub2_2)
    ggData$grpBy = paste(ggData$sub1,ggData$sub2, ggData$geneName)
    # Aggregate
    ggData = pivot_wider(ggData, names_from = grpBy, values_from = val)
    ggData[is.na(ggData)] = 0
    corr <- cor(ggData[5:ncol(ggData)],method = "spearman", use = "complete.obs")
    p.mat <- cor_pmat(ggData[5:ncol(ggData)],method = "spearman", use = "complete.obs")

    ggcorrplot(corr, hc.order = FALSE,show.diag = FALSE,
               method ="square", type = "lower", p.mat = p.mat)
    # Scale if required
    colRange = range(ggData$val)
    if(inpScl){
        ggData[, val:= scale(val), keyby = "geneName"]
        colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
    }
    ggcorrplot()
    # hclust row/col if necessary
    ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
    tmp = ggMat$geneName
    ggMat = as.matrix(ggMat[, -1])
    rownames(ggMat) = tmp
    if(inpRow){
        hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat))))
        ggRow = ggplot() + coord_flip() +
            geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
                               labels = unique(ggData$grpBy), expand = c(0, 0)) +
            scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                               labels = hcRow$labels$label, expand = c(0, 0.5)) +
            sctheme(base_size = sList[inpfsz]) +
            theme(axis.title = element_blank(), axis.line = element_blank(),
                  axis.ticks = element_blank(), axis.text.y = element_blank(),
                  axis.text.x = element_text(color="white", angle = 45, hjust = 1,vjust=0))
        ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
    } else {
        ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
    }
    if(inpCol){
        hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
        ggCol = ggplot() +
            geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
            scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                               labels = hcCol$labels$label, expand = c(0.05, 0)) +
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
                               labels = unique(ggData$geneName), expand=c(0,0)) +
            sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            theme(axis.title = element_blank(), axis.line = element_blank(),
                  axis.ticks = element_blank(), axis.text.x = element_blank(),
                  axis.text.y = element_text(color = "white"))
        ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
    }

    # Actual plot according to plottype
    cols =  color_generator(inpcols)
    if(inpcolinv) cols = rev(cols)
    if(inpPlt == "Bubbleplot"){
        # Bubbleplot
        ggOut = ggplot(ggData, aes(grpBy, geneName,fill = val, size = prop)) +
            geom_point(color = "black",pch=21) +
            scale_x_discrete(expand = c(0.05, 0)) +
            scale_y_discrete(expand = c(0, 0.5)) +
            scale_size_continuous("proportion", range = c(0, inpsiz),
                                  limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
            scale_fill_gradientn("expression", limits = colRange, colours = cols) +
            cowplot::theme_cowplot(font_size = sList[inpfsz])+
            #sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            guides(color = guide_colorbar(barwidth = 10)) +
            theme(axis.title = element_blank(), legend.box = "horizontal",
                  axis.text.x =    element_text(angle = as.numeric(inpfrt), hjust = 1,vjust = 0.5))
    } else {
        # Heatmap
        ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) +
            geom_tile() +
            scale_x_discrete(expand = c(0.05, 0)) +
            scale_y_discrete(expand = c(0, 0.5)) +
            scale_fill_gradientn("expression", limits = colRange, colours = cols) +
            cowplot::theme_cowplot(font_size = sList[inpfsz])+
            #sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
            guides(fill = guide_colorbar(barwidth = 10)) +
            theme(axis.title = element_blank(), legend.box = "horizontal",
                  axis.text.x =    element_text(angle = as.numeric(inpfrt), hjust = 1,vjust = 0.5))
    }
    rat = (length(unique(ggData$grpBy)) + as.integer(inpRow)) / (length(unique(ggData$geneName)) + as.integer(inpCol))
    if(inpasp == "Square") {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed") {
        ggOut = ggOut + coord_fixed()
    }

    # Final tidy
    if(inpleg) {ggLeg = g_legend(ggOut)} else {
        ggLeg = ggplot() + geom_point() +theme_cowplot(line_size = 0)
    }
    ggOut = ggOut + theme(legend.position = "none",
                          axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                                        hjust = switch (as.character(inpfrt),
                                                                        "0" = 0.5,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0
                                                        ),
                                                        vjust = switch (as.character(inpfrt),
                                                                        "0" = 0,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0.5)
                          )
    )

    return(ggOut)
}

# Plot Trajectory Curve lines
scCurvs <- function(inpConf, inpMeta, inp, inpGrp,
                    inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,
                    inpH5, inpGene,inplog, inpnorm,inpy_0,
                    inppts, inpcols, inplvl, inpfsz,inplsz, inpfrt,
                    inpleg,inplegpos = "bottom"){

    # Identify genes that are in our dataset
    geneList = scGeneList(inp, inpGene)
    geneList = geneList[present == TRUE]
    shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
    shiny::validate(need(nrow(geneList) >= 1, "Please input at least 1 genes to plot!"))

    # Prepare ggData
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData = data.table()
    for(iGene in geneList$gene){
        tmp = inpMeta[, c("sampleID", inpConf[grp == TRUE]$ID,
                          inpConf[UI == inpsub1_1]$ID, inpConf[UI == inpsub2_1]$ID), with = FALSE]
        N = ncol(tmp)
        colnames(tmp)[(N-1):N] = c("sub1","sub2")
        tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
        tmp$geneName = iGene
        tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=)))
        ggData = rbindlist(list(ggData, tmp))
    }
    h5file$close_all()

    ggData = Subset(ggData, inpsub1_2, inpsub2_2)

    ggData$grpBy %<>% droplevels()
    if (inpGrp %in% c(inpsub1_1, inpsub2_1)) {
        grp = which(c(inpsub1_1, inpsub2_1) %in% inpGrp)[1]
        if (is.null(list(inpsub1_2,inpsub2_2)[[grp]])) {
            GrpLvl = levels(pull(ggData, var = c(inpsub1_1, inpsub2_1)[grp]))
        } else GrpLvl = list(inpsub1_2,inpsub2_2)[[grp]]

        ggData$grpBy %<>% factor(levels = GrpLvl)
    }

    # Aggregate
    if(!inplog) ggData$val = expm1(ggData$val)

    #ggData = ggData[, .(mean_val = mean(val)),
    #                                by = c("geneName","grpBy"),allow.cartesian=TRUE]

    ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val",fun.aggregate = mean)

    # check if any empty group
    count_table = as.matrix(table(ggData$grpBy))
    empGroup = rownames(count_table)[count_table==0]
    if(length(empGroup)>0) for(g in empGroup) ggMat[,g] = 0
    setcolorder(ggMat, c("geneName", GrpLvl))

    ggData = melt(ggMat, id="geneName",variable.name = "grpBy", value.name = "mean_val")

    if(inpnorm){
        ggData = ggData[, .(grpBy = grpBy, max_val = max(mean_val), mean_val=mean_val),
                        by = c("geneName")]
        ggData = ggData[, .(val = mean_val /max_val),
                        by = c("geneName", "grpBy","max_val")]
        ggData[is.na(ggData)] = 0
        ggData$val = ggData$val * 100

    } else ggData$val = ggData$mean_val

    ggData$index = plyr::mapvalues(ggData$grpBy,
                                   from = levels(ggData$grpBy),
                                   to = seq_len(length(levels(ggData$grpBy))))
    ggData$index %<>% as.integer()

    # Do factoring
    ggData$geneName %<>% factor(levels = geneList$gene)
    gglvl = levels(ggData$geneName)
    ggCol = colorRampPalette(brewer.pal(12, "Paired"))(length(gglvl))

    # Actual ggplot
    ggOut = ggscatter(data = ggData,
                      x = "index", y ="val",
                      color = "geneName",
                      shape = ifelse(inppts,19,NA),# Extending the regression line
                      palette = switch(inpcols,
                                       "default" =ggCol,
                                       color_generator(inpcols,length(gglvl))),
                      ylab = ifelse(inpnorm,"Max expression level %",
                                    ifelse(inplog,"log mean UMI","mean UMI")),
                      size = 3,
                      #title = paste0(paste(geneList$gene, collapse = ", ")," expression"),
                      add = "loess",add.params =    list(size=sList[inplsz]/15),
                      conf.int.level = inplvl,
                      conf.int = TRUE)+
        theme(text = element_text(size=sList[inpfsz]),
              axis.text.x =    element_text(angle = as.numeric(inpfrt), hjust = 1))+
        scale_x_continuous(name =inpGrp, breaks = seq_len(length(levels(ggData$grpBy))),
                           labels = levels(ggData$grpBy))
    if(inpy_0) ggOut = ggOut + expand_limits(y = 0)

    ggOut = ggOut + theme(legend.position = switch (as.character(inpleg),
                                                    "FALSE" = "none",
                                                    "TRUE" = inplegpos),
                          axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                                        hjust = switch (as.character(inpfrt),
                                                                        "0" = 0.5,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0
                                                        ),
                                                        vjust = switch (as.character(inpfrt),
                                                                        "0" = 0,
                                                                        "30" = 1,
                                                                        "45" = 1,
                                                                        "90" = 0.5)
                          )
    )
    return(ggOut)
}

scFindMarkers <- function(inpConf, inpMeta,
                          inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                          inpident, inpident_1, inpident_2 = NULL){
    # Prepare ggData
    ggData = inpMeta[, c("sampleID",
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID,
                         inpConf[UI == inpident]$ID),
                     with = FALSE]
    colnames(ggData) = c("sampleID","sub1","sub2","sub3","Idents")
    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
    cells_1 = which(inpMeta$sampleID %in% ggData[Idents %in% inpident_1]$sampleID) - 1 %>%
        as.integer()
    if(is.null(inpident_2)) {
        gglvl = unique(as.character(ggData$Idents))
        inpident_2 = gglvl[-which(gglvl %in% inpident_1)]
    }
    cells_2 = which(inpMeta$sampleID %in% ggData[Idents %in% inpident_2]$sampleID) -1 %>%
        as.integer()
    # declared variable to the global environment
    assign("inpident_1", inpident_1, envir = globalenv())
    assign("inpident_2", inpident_2, envir = globalenv())
    assign("Cells1", cells_1, envir = globalenv())
    assign("Cells2", cells_2, envir = globalenv())
    if(!dir.exists("tempData")) dir.create("tempData")

    reticulate::py_run_file("scFindMarkers_ad.py")

    # prepare a info  to store different analysis conditions information

    info = paste0(paste(inpident_1,collapse = "_")," vs ",
                  paste(inpident_2,collapse = "_"), " in ",
                  paste(paste(inpsub1_2,collapse = "_"),
                        paste(inpsub2_2,collapse = "_"),
                        paste(inpsub3_2,collapse = "_"),collapse = "_"))
    write.table(info, file = "tempData/de_info.csv",quote = FALSE,row.names = F,col.names = F)
}

.checkIfIdentical <- function(inpsub1_2,inpsub2_2,inpsub3_2,inpident_1 = "",inpident_2 = "",
                              min_expr = 10,
                              file.name = "tempData/de_info.csv"){
    # prepare a info not to check if stored different analysis conditions is identical

    info = switch(file.name,
                  "tempData/de_info.csv" = paste0(
                                  paste(inpident_1,collapse = "_")," vs ",
                                  paste(inpident_2,collapse = "_"), " in ",
                                  paste(paste(inpsub1_2,collapse = "_"),
                                        paste(inpsub2_2,collapse = "_"),
                                        paste(inpsub3_2,collapse = "_"),collapse = "_")),
                  "tempData/cor_info.csv" = paste(c(
                                          paste(inpsub1_2,collapse = "_"),
                                          paste(inpsub2_2,collapse = "_"),
                                          paste(inpsub3_2,collapse = "_"),
                                          paste(min_expr,"%")),collapse = "_")
                  )
    if(!file.exists(file.name)) {
        return(FALSE)
    } else b = read.delim(file.name,header = FALSE)
    return(info == b$V1)
}



loadDeData <- function(h5File ="tempData/scFindMarkers.h5", key = "de",inprmgene = NULL){
    reticulate::py_run_string("import pandas as pd")
    pd <- reticulate::import("pandas")
    markers = pd$read_hdf(h5File, key = key)

    if(markers$gene[1] == toupper(markers$gene[1])){
        rm_gene = c("^MT-","^RPL","^RPS")
    } else rm_gene = c("^mt-","^Rpl","^Rps")
    # remove MT, RPS and PRL genes
    if(length(inprmgene) != 0){
        geneGrp = matrix(c("Mitochondrial (MT) genes",
                           "Ribosomal protein large subunit (RPL)",
                           "Ribosomal protein small subunit (RPS)",
                           rm_gene), nrow =3)
        rmgene = geneGrp[,1] %in% inprmgene %>% geneGrp[.,2] %>% paste(collapse = "|") %>%
            grep(markers$gene,value = T)
        markers = markers[!markers$gene %in% rmgene,]
    }

    return(markers)
}

#' laod correlation data by single gene
loadCorDataGene <- function(inpH5 = "tempData/scFindCor.h5",featureFile= "tempData/CorFeatures.rds", gene = ""){
    h5file <- H5File$new(inpH5, mode = "r")
    #print(h5file$ls(recursive=TRUE))

    features <- readRDS(featureFile) %>% .[["features"]]
    shiny::validate(need(gene %in% features, message = paste0("No gene is significantly correlated with ",gene,
                                                             ". Need lower the threshold or try other gene!")))

    geneIdx = which(features %in% gene[1])
    ggData <- data.table("correlation" = h5file[["cor"]]$read(args = list(geneIdx, quote(expr=))),
                         "p_val" = h5file[["pval"]]$read(args = list(geneIdx, quote(expr=))),
                         stringsAsFactors = F)
    h5file$close_all()

    ggData[["genes"]] = features
    ggData = ggData[-geneIdx,]
    ggData = ggData[!is.nan(ggData$p_val),]
    ggData = ggData[order(ggData$correlation,decreasing = T),]

    return(list("data" = ggData,
                "gene" = gene))
}

#' laod correlation matrix and pval matrix by multiple genes
loadCorPvalGenes <- function(inpH5 = "tempData/scFindCor.h5",
                             featureFile= "tempData/CorFeatures.rds",
                             gene = "", topN = 10, inpcutp = c("p_val_adj","p_val")[1],
                             corPosNeg = "both"){
    topN %<>% as.integer()
    h5file <- H5File$new(inpH5, mode = "r")
    #print(h5file$ls(recursive=TRUE))

    features <- readRDS(featureFile) %>% .[["features"]]
    shiny::validate(need(gene %in% features, message = paste0("No gene is significantly correlated with ",gene,
                                                              ". Need lower the threshold or try other gene!")))

    geneIdxs = which(features %in% gene)

    h5file <- H5File$new(inpH5, mode = "r")
    topNCorGenes = c()
    for(geneIdx in geneIdxs){
        topNCorGenes = switch(corPosNeg,
                              "positive only" = unique(c(topNCorGenes,tail(order(h5file[["cor"]][geneIdx,]),topN))),
                              "negative only" = unique(c(topNCorGenes,head(order(h5file[["cor"]][geneIdx,]),topN))),
                              "both" = unique(c(topNCorGenes,
                                                head(order(h5file[["cor"]][geneIdx,]),n = ceiling(topN/2)),
                                                tail(order(h5file[["cor"]][geneIdx,]),n = ceiling(topN/2))))
                              )
    }
    topNCorGenes = unique(c(geneIdxs,topNCorGenes))
    corMat = h5file[["cor"]][topNCorGenes,topNCorGenes]
    pvalMat = h5file[["pval"]][,topNCorGenes]
    if(inpcutp == "p_val_adj") pvalMat = apply(pvalMat, 2,function(x) p.adjust(x,method = "fdr"))
    pvalMat = pvalMat[topNCorGenes,]
    h5file$close_all()
    colnames(corMat) = rownames(corMat) =
        colnames(pvalMat) = rownames(pvalMat) = features[topNCorGenes]


    return(list("cor" = corMat, "pval" = pvalMat, "genes" = features[geneIdxs]))
}

# Adapted from http://pandas.pydata.org/pandas-docs/stable/io.html#io-external-compatibility

scFindCor <- function(inpConf, inpMeta,
                      inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                      min_expr = 10){
    # Prepare ggData
    ggData = inpMeta[, c("sampleID",
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID),
                     with = FALSE]
    colnames(ggData) = c("sampleID","sub1","sub2","sub3")
    ggData = Subset(ggData, inpsub1_2, inpsub2_2, inpsub3_2)
    cells = which(inpMeta$sampleID %in% ggData$sampleID) - 1 %>%
        as.integer()
    # declared variable to the global environment
    assign("Cells", cells, envir = globalenv())
    assign("Min_expr", min_expr/100, envir = globalenv())
    print("scFindCor start!!")
    if(!dir.exists("tempData")) dir.create("tempData")
    reticulate::py_run_file("scFindCor_ad.py")
    print("scFindCor end!!")

    # prepare a info  to store different analysis conditionsf information

    info = paste(c(
        paste(inpsub1_2,collapse = "_"),
        paste(inpsub2_2,collapse = "_"),
        paste(inpsub3_2,collapse = "_"),
        paste(min_expr,"%")),collapse = "_")
    write.table(info, file = "tempData/cor_info.csv",quote = FALSE,row.names = F,col.names = F)
}


Subset <- function(ggData, inpsub1_2, inpsub2_2 = NULL, inpsub3_2 = NULL){

    if(length(inpsub1_2) != 0 & length(inpsub1_2) != nlevels(ggData$sub1)){
        ggData = ggData[sub1 %in% inpsub1_2]
    }
    if(length(inpsub2_2) != 0 & length(inpsub2_2) != nlevels(ggData$sub2)){
        ggData = ggData[sub2 %in% inpsub2_2]
    }
    if(length(inpsub3_2) != 0 & length(inpsub3_2) != nlevels(ggData$sub3)){
        ggData = ggData[sub3 %in% inpsub3_2]
    }
    ggData = ggData[complete.cases(ggData)]
    return(ggData)
}

# VolcanoPlots to demonstrate Differential expressed genes
# https://zhuanlan.zhihu.com/p/82785739?utm_source=ZHShareTargetIDMore&utm_medium=social&utm_oi=642996063045423104
VolcanoPlots <- function(ggData, inp, inpGene,inpsub1_2,inpsub2_2,inpsub3_2,
                         inpcutp, inpcutpval, inpcutfc,inptop,
                         inpsort,inpcols, inputcolinv,inpalpha, inpsiz, inpfsz,
                         inpasp, inplab1,inplab2,inpleg, inplegpos = "bottom"){
    # Identify genes that are in our dataset
    geneList = scGeneList(inp, inpGene)
    geneList = geneList[present == TRUE]
    cols =  color_generator(inpcols,n = 10,alpha = inpalpha)[seq(1, 10,2)]
    if(inputcolinv) cols = rev(cols)
    cols.order =c("Most Down","Down","Stable","Up","Most Up")
    names(cols) = cols.order
    if(inplegpos == "right") cols.order = rev(cols.order)
    ggData[,paste0("log10_",inpcutp)] = -log10(ggData[,inpcutp])
    inf = is.infinite(ggData[,paste0("log10_",inpcutp)])
    ggData[inf,paste0("log10_",inpcutp)] = 337
    ggData[,"change"] = "Stable"
    inpcutpval %<>% as.numeric()
    inpcutfc %<>% as.numeric()

    ggData[ggData[,inpcutp] <= inpcutpval & ggData$avg_log2FC > 0,"change"] = "Up"
    ggData[ggData[,inpcutp] <= inpcutpval & ggData$avg_log2FC < 0,"change"] = "Down"

    ggData[ggData$avg_log2FC > inpcutfc & ggData[,inpcutp]  < inpcutpval,"change"] = "Most Up"
    ggData[ggData$avg_log2FC < -inpcutfc & ggData[,inpcutp]  < inpcutpval,"change"] = "Most Down"

    ggData$change %<>% factor(levels = cols.order)
    colnames(ggData)[grep("cluster",colnames(ggData))]="cluster"
    # 将需要标记的基因放置在单独的数组
    Up <- ggData[ggData$change %in% "Most Up",]
    Down <- ggData[ggData$change %in% "Most Down",]
    inptop %<>% as.integer()
    if(inpsort %in% c("p_val_adj","p_val")) {
        Up_gene_index <- rownames(Up)[Up[,inpsort] <= tail(head(sort(Up[,inpsort],decreasing = F),inptop),1)]
        Down_gene_index <- rownames(Down)[Down[,inpsort] <= tail(head(sort(Down[,inpsort],decreasing = F),inptop),1)]
    }
    if(inpsort == "avg_log2FC") {
        Up_gene_index <- rownames(Up)[Up[,inpsort] >= tail(head(sort(Up[,inpsort],decreasing = T),inptop),1)]
        Down_gene_index <- rownames(Down)[Down[,inpsort] <= tail(head(sort(Down[,inpsort],decreasing = F),inptop),1)]
    }
    # If too many gene with adj_p = 0
    if(length(Up_gene_index) >inptop) Up_gene_index = head(Up_gene_index,inptop)
    if(length(Down_gene_index) >inptop) Down_gene_index = tail(Down_gene_index,inptop)
    # prepare title
    cluster = stringr::str_split(ggData$cluster[1],patter = " vs.")[[1]]
    cluster = paste(rev(cluster),collapse = " <----    ----> ")
    subtitle = paste(c(paste(inpsub1_2,collapse = "."),
                       paste(inpsub2_2,collapse = "."),
                       paste(inpsub3_2,collapse = ".")),collapse = " ")
    if(subtitle == " ") subtitle = NULL
    ggOut<-ggplot(
        #设置数据
        ggData,
        mapping = aes_string(x = "avg_log2FC",
                             y = paste0("log10_",inpcutp),
                             fill = "change"))+
        geom_point(alpha= inpalpha,size=inpsiz,color = "black", pch=21)+
        ggtitle(label = paste(c(cluster,subtitle),collapse = " \n in "))

    # 辅助线
    ggOut = ggOut + geom_vline(xintercept=c(-inpcutfc,inpcutfc),lty=4,col="black",lwd=0.8)
    ggOut = ggOut + geom_hline(yintercept = -log10(inpcutpval),lty=4,col="black",lwd=0.8)

    # 坐标轴
    ggOut = ggOut + theme_bw()+
        sctheme(base_size = sList[inpfsz])+
        labs(x="log2(fold change)",
             y= paste("-log10 (",ifelse(inpcutp == "p_val_adj", "adjusted p-value","p-value"),")"))+

        # 图例
        theme(plot.title = element_text(hjust = 0.5),
              axis.title=element_text(size=sList[inpfsz]),
              legend.title = element_blank(),
              legend.text = element_text(size = sList[inpfsz]),
              legend.position = switch (as.character(inpleg),
                                      "FALSE" = "none",
                                      "TRUE" = inplegpos),
        )
    ggOut = ggOut + scale_fill_manual("",values=cols[unique(ggData$change)])

    if(inplab1 != "No labels" & length(c(Down_gene_index, Up_gene_index)) > 0){
        options(ggrepel.max.overlaps = Inf)
        ggOut = ggOut + switch(inplab1,
                               "black text" = geom_text_repel(data = ggData[c(Down_gene_index, Up_gene_index),],
                                                              aes_string(x = "avg_log2FC",
                                                                         y = paste0("log10_",inpcutp),
                                                                         label = "gene"),
                                                              colour = "grey10",
                                                              box.padding = unit(0.5, "lines"),
                                                              point.padding = unit(0.8, "lines"),
                                                              size = lList[inpfsz], seed = 42),
                               "black labels" = geom_label_repel(data = ggData[c(Down_gene_index, Up_gene_index),],
                                                                 aes_string(x = "avg_log2FC",
                                                                            y = paste0("log10_",inpcutp),
                                                                            label = "gene"),
                                                                 colour = "grey10",
                                                                 fill = "white",
                                                                 alpha = inpalpha,
                                                                 box.padding = unit(0.5, "lines"),
                                                                 point.padding = unit(0.8, "lines"),
                                                                 size = lList[inpfsz], seed = 42)
        )
    }
    if(inplab2 != "No labels" & any(geneList$gene %in% ggData$gene)){
        ggOut = ggOut + switch(inplab2,
                               "red text" = geom_text_repel(data = ggData[ggData$gene %in% geneList$gene,],
                                                              aes_string(x = "avg_log2FC",
                                                                         y = paste0("log10_",inpcutp),
                                                                         label = "gene"),
                                                              colour = "red",
                                                              force = 5,
                                                              box.padding = unit(0.5, "lines"),
                                                              point.padding = unit(0.8, "lines"),
                                                              size = lList[inpfsz], seed = 42),
                               "red labels" = geom_label_repel(data = ggData[ggData$gene %in% geneList$gene,],
                                                                 aes_string(x = "avg_log2FC",
                                                                            y = paste0("log10_",inpcutp),
                                                                            label = "gene"),
                                                                 colour = "red",
                                                                 fill = "white",
                                                                 force = 5,
                                                                 alpha = inpalpha,
                                                                 box.padding = unit(0.5, "lines"),
                                                                 point.padding = unit(0.8, "lines"),
                                                                 size = lList[inpfsz], seed = 42)
        )
    }
    rat = (max(ggData$avg_log2FC) - min(ggData$avg_log2FC)) /
        (max(ggData[,paste0("log10_",inpcutp)]) - min(ggData[,paste0("log10_",inpcutp)]))
    if(inpasp == "Square" & rat != Inf & rat != -Inf) {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed" & rat != Inf & rat != -Inf) {
        ggOut = ggOut + coord_fixed()
    }
    return(ggOut)
}

#' @param inp1 single gene input that Corplot will focus on.
#' @param inp single or multiple genes that manually input for visulization
#'
CorPlots <- function(ggData, inp1, inp = "", inpGene, cor_info,
                     inpcutp = c("p_val_adj","p_val")[1],
                     inpcutpval, inpcutfc,inptop,
                     inpcol1,inpcol2,
                     inpalpha, inpsiz, inpfsz,inpasp){
    # Identify genes that are in our dataset
    geneList = inp[inp %in% inpGene]

    inpcutpval %<>% as.numeric()
    inpcutfc %<>% as.numeric()

    if(inpcutp == "p_val_adj") ggData$p_val_adj = p.adjust(ggData$p_val,method = "fdr")
    ggData[[paste0("log10_",inpcutp)]] = -log10(ggData[[inpcutp]])
    ggData$idx = 1:nrow(ggData)
    # 将需要标记的基因放置在单独的数组
    inptop %<>% as.integer()
    Up_gene_index <- ggData %>% filter(correlation > inpcutfc) %>% slice_max(correlation, n = inptop) %>% .[["idx"]]
    Down_gene_index <- ggData %>% filter(correlation < -inpcutfc) %>% slice_min(correlation, n= inptop) %>% .[["idx"]]


    # If too many gene with adj_p = 0
    if(length(Up_gene_index) >inptop) Up_gene_index = head(Up_gene_index,inptop)
    if(length(Down_gene_index) >inptop) Down_gene_index = tail(Down_gene_index,inptop)
    # prepare title
    subtitle = cor_info
    if(subtitle == " ") subtitle = NULL

    # prepare selected label
    label.select = unique(c(ggData$genes[c(Up_gene_index,Down_gene_index)],geneList))
    label.select = label.select[label.select %in% ggData$genes]
    shiny::validate(need(length(label.select) > 0, "No gene pass filter. Lower the threshold!"))

    label_ggData = ggData[ggData$genes %in% label.select,]
    label_ggData[["mark"]] = inpcol1
    if(length(geneList) >0) label_ggData[label_ggData$genes %in% geneList, "mark"] = inpcol2
    options(ggrepel.max.overlaps = (inptop*2+length(geneList)))
    ggOut <- ggline(ggData,
                    x = "correlation",
                    y = paste0("log10_",inpcutp),
                    numeric.x.axis = TRUE,
                    ylab = paste("-log10 (",ifelse(inpcutp == "p_val_adj", "adjusted p-value","p-value"),")"),
                    xlab = "Spearman Correlation",
                    font.label = list(size = sList[inpfsz],
                                      face = "plain",
                                      color = label_ggData$mark),
                    label = "genes",             # Add point labels
                    label.select = label.select,
                    repel = TRUE,
                    max.overlaps = Inf,
                    title = paste(c(paste("gene correlated with",inp1),subtitle),collapse = " \n in "))+
        geom_point(data = label_ggData,
                   color = label_ggData$mark,
                   alpha= inpalpha,size=inpsiz,pch=21)+
        ggtitle(label = paste(c(paste("gene correlated with",inp1),subtitle),collapse = " \n in "))
    ggOut
    # 辅助线
    ggOut = ggOut + geom_vline(xintercept=c(-inpcutfc,inpcutfc),lty=4,col="black",lwd=0.8)
    ggOut = ggOut + geom_hline(yintercept = -log10(inpcutpval),lty=4,col="black",lwd=0.8)

    # 坐标轴
    ggOut = ggOut + theme_bw()+
        sctheme(base_size = sList[inpfsz])+

        # 图例
        theme(plot.title = element_text(hjust = 0.5),
              axis.title=element_text(size=sList[inpfsz])
        )
    rat = (max(ggData$correlation) - min(ggData$correlation)) /
        (max(ggData[[paste0("log10_",inpcutp)]]) - min(ggData[[paste0("log10_",inpcutp)]]))
    if(inpasp == "Square" & rat != Inf & rat != -Inf) {
        ggOut = ggOut + coord_fixed(ratio = rat)
    } else if(inpasp == "Fixed" & rat != Inf & rat != -Inf) {
        ggOut = ggOut + coord_fixed()
    }
    return(ggOut)
}


#'@param corMat correlation matrix of selected genes
#'@param pvalMat pvalue matrixo of selected genes
#'@param inp2 single or multiple genes input that CorNetwork will focus on.
#'
CorNetwork <- function(corMat, pvalMat, cor_info,
                     inpcutpval, inpcutfc,
                     corPosNeg = "positive only",
                     inpcols,inpcolinv,
                     linkDistance,linkWidth, inpfsz,showName,
                     inpleg= FALSE, inplegpos = "bottom"){

    inpcutpval %<>% as.numeric()
    inpcutfc %<>% as.numeric()

    #将矩阵中不符合条件的r值替换为0；

    if(corPosNeg == "positive only") corMat[which(corMat <= inpcutfc)]=0
    if(corPosNeg == "negative only") corMat[which(corMat >= -inpcutfc)]=0
    if(corPosNeg == "both") corMat[which(abs(corMat) <= inpcutfc)]=0

    corMat[which(pvalMat > inpcutpval)]=0
    # 删掉相关系数矩阵数据全都为0的行和列；
    # 对角线为1，rowSums 和 colSums 必然 >=1
    keep = which(rowSums(abs(corMat)) > 1)
    corMat <- corMat[names(keep),names(keep)]

    shiny::validate(need(!any(dim(corMat) == c(0,0)), "No gene pass filter. Either lower the correlation cut-off,
                                         or increase the p-value cut-off"))

    # Remove duplicate edges
    g1 <- graph.adjacency(corMat,weight=T,mode="undirected")

    g1 <- simplify(g1)

    # Find group membership
    #计算群体结构（short random walks）；
    wt <- cluster_walktrap(g1, steps = 6)
    members <- membership(wt)
    #使用默认颜色列表；
    #V(g1)$color  <- c$membership+1

    # Convert igraph to list for networkD3
    sj_list <- igraph_to_networkD3(g1, group = members)

    ############################################################################################
    #http://www.vesnam.com/Rblog/viznets6/

    # prepare link color


    # generate color scale based on corrlation value range.
    Range = range(sj_list$links$value, na.rm = T)
    n = ceiling(max(Range)*10)*4

    ggCol = switch(inpcols,
                   "default" = rep("#333",n*2),
                   color_generator(inpcols,n = n*2))
    if(inpcolinv) ggCol = rev(ggCol)

    linkcols = cut(sj_list$links$value, breaks = c(-n:n)/40)

    sj_list$links$color = plyr::mapvalues(linkcols,
                              from = levels(linkcols),
                              to = ggCol)
    sj_list$links$color %<>% droplevels()
    # Plot as a forceDirected Network
    sj_list$nodes$group = 1

    p <- forceNetwork(Links = sj_list$links, Nodes = sj_list$nodes,
                      #colourScale = JS('force.alpha(1); force.restart(); d3.scaleOrdinal(d3.schemeCategory20);'),
                      Source = 'source',Value = "value",
                      Target = 'target', NodeID = 'name', Group = 'group',
                      linkColour = sj_list$links$color,
                      fontFamily = "Avenir",
                      fontSize = 20,
                      #Nodesize =  "nodeBetweenness",
                      opacity = 1, #as.numeric(input$opacity)
                      opacityNoHover = ifelse(showName,5,0),#as.numeric(input$showName)
                      legend = inpleg, bounded = FALSE,
                      linkDistance = as.numeric(linkDistance),
                      linkWidth = as.numeric(linkWidth),#JS("function(d) { return d.value/5; }"),
                      #radiusCalculation = JS(" Math.sqrt(d.nodesize)+6"),
                      charge = -30,
                      zoom = FALSE
                      )

    return(p)
}

# Wrapper Functions for scRepertoire v1.3.4
scRepertoire <- function(inpConf, inpMeta, inp1, inp2,
                         inpsub1_1, inpsub1_2,inpsub2_1, inpsub2_2,inpsub3_1, inpsub3_2,
                         inptyp, inpscale, inpindex, inpcloneCall,
                         inplvls, inpcols,
                         inpflpxy, inppts, inpsiz, inpfsz, inplsz,
                         inpfrt, inpleg, inplegpos = "bottom"){
    # Prepare ggData
    cloneCall <- switch(inpcloneCall,
                        "TCR/Ig genes" = "CTgene",
                        "CDR3 nucleotide" = "CTnt",
                        "CDR3 amino acid" = "CTaa",
                        "TCR/Ig + CDR3 nucleotide" = "CTstrict")

    ggData = inpMeta[, c(inpConf[UI == inp1]$ID,
                         inpConf[UI == inp2]$ID,
                         inpConf[UI == "patient"]$ID,
                         inpConf[UI == cloneCall]$ID,
                         inpConf[UI == "tcr.frequency"]$ID,
                         inpConf[UI %in% inpsub1_1]$ID,
                         inpConf[UI %in% inpsub2_1]$ID,
                         inpConf[UI %in% inpsub3_1]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", inp2,"paired",cloneCall,"tcr.frequency","sub1", "sub2", "sub3")

    ggData =  Subset(ggData, inpsub1_2, inpsub2_2,inpsub3_2)


    # Do factoring

    temp1 <- doFactoring(ggData,inpConf, col = "X", inp1, inpcols)
    temp2 <- doFactoring(temp1$ggData,inpConf, col = inp2, inp2, inpcols)

    gglvl1 = temp1$gglvl; ggCol1 = temp1$ggCol
    gglvl2 = temp2$gglvl; ggCol2 = temp2$ggCol
    ggData = temp2$ggData
    rm(temp1);rm(temp2)

    # Actual ggplot
    ggOut = switch (inptyp,
                    "Unique Barplot" = scquantContig(ggData,cloneCall = cloneCall,
                                                     scale = inpscale,gglvl1,
                                                     cols= ggCol1[gglvl1]),

                    "Clone sizes distribution" = scabundanceContig(ggData,inptyp,
                                                              cloneCall = cloneCall, chain = "both",
                                                              group = "X",
                                                              cols= ggCol1[gglvl1],
                                                              inppts,inpsiz, inpfsz,inplsz),
                    "Cumulative clone sizes distribution" = scabundanceContig(ggData,inptyp,
                                                                         cloneCall = cloneCall, chain = "both",
                                                                         group = "X",
                                                                         cols= ggCol1[gglvl1],
                                                                         inppts,inpsiz, inpfsz,inplsz),
                    "Proportion Barplot" = scClonalProportion(ggData,Split = c(10, 100, 1000, 10000, 30000, 100000),
                                                              cloneCall = cloneCall, chain = "both",inpcols),
                    "Paired Diversity" = scClonalDiversity(ggData, inptyp,
                                                           cloneCall = cloneCall, chain = "both", scale=inpscale,
                                                          inpindex = inpindex,inp1, groupBy = inp2, x.axis = "X",
                                                          cols = ggCol1[gglvl1],
                                                          gglvl1, inppts,inpsiz, inpfsz, inplsz = inplsz,
                                                          n.boots = 100),
                    "unpaired Diversity" = scClonalDiversity(ggData, inptyp,
                                                             cloneCall = cloneCall, chain = "both", scale=inpscale,
                                                             inpindex = inpindex,inp1, groupBy = inp2, x.axis = "X",
                                                             cols = ggCol1[gglvl1],
                                                             gglvl = gglvl1, inppts,inpsiz,inpfsz, inplsz,
                                                             n.boots = 100),
                    "Paired scatter Clonotype" = ScatterClonotype(ggData, inptyp, inpConf,
                                                           cloneCall = cloneCall, chain = "both", scale=inpscale,
                                                           inp2, x.axis = gglvl1[1], y.axis = gglvl1[2],
                                                           inplvls,inpcols, inpsiz, inplsz),
                    "unpaired scatter Clonotype" = ScatterClonotype(ggData, inptyp, inpConf,
                                                                    cloneCall = cloneCall, chain = "both", scale=inpscale,
                                                                    inp2, x.axis = gglvl1[1], y.axis = gglvl1[2],
                                                                    inplvls,inpcols, inpsiz, inplsz),
                    "Novel & Expand & Contract Barplot_1" = ScClonotypeBar(ggData, inptyp, inpConf,
                                                                         cloneCall = cloneCall, chain = "both",
                                                                         x.axis = gglvl1[1], y.axis = gglvl1[2],
                                                                         inp2, gglvl2, ggCol2,inpscale,inpcols,
                                                                         inplvls,inppts,inpsiz, inplsz),
                    "Novel & Expand & Contract Barplot_2" = ScClonotypeBar(ggData, inptyp, inpConf,
                                                                         cloneCall = cloneCall, chain = "both",
                                                                         x.axis = gglvl1[1], y.axis = gglvl1[2],
                                                                         inp2, gglvl2, ggCol2,inpscale,inpcols,
                                                                         inplvls,inppts,inpsiz, inplsz)
    )

    if(inpflpxy){
        ggOut = ggOut + coord_flip()# + scale_x_discrete(limits = rev)
    }
    ggOut = ggOut +
        sctheme(base_size = sList[inpfsz])+
        theme(legend.position = switch (as.character(inpleg),
                                        "FALSE" = "none",
                                        "TRUE" = inplegpos),
              legend.title = element_blank(),
              axis.text.x =    element_text(angle = as.numeric(inpfrt),
                                            hjust = switch (as.character(inpfrt),
                                                            "0" = 0.5,
                                                            "30" = 1,
                                                            "45" = 1,
                                                            "90" = 0
                                            ),
                                            vjust = switch (as.character(inpfrt),
                                                            "0" = 0,
                                                            "30" = 1,
                                                            "45" = 1,
                                                            "90" = 0.5)
              )
        )

    return(ggOut)
}
#' @param ggData The product of combineTCR() into data.frame
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param group The column header used for grouping.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"

scquantContig <- function (ggData, cloneCall = "CTstrict", scale = FALSE,
                           gglvl, chain = "both",cols)
{
    df <- base::split(ggData, by = "X")
    df = df[gglvl]

    # Prepare ggData
    inpcloneCall <- switch(cloneCall,
                        "CTgene" = "TCR/Ig genes",
                        "CTnt" = "CDR3 nucleotide",
                        "CTaa" = "CDR3 amino acid",
                        "CTstrict" = "TCR/Ig + CDR3 nucleotide")
    x <- "values"
    labs <- "Samples"
    Con.df <- data.frame(matrix(NA, length(df), 3))
    colnames(Con.df) <- c("contigs", "values", "total")
    for (i in seq_along(df)) {
        Con.df[i, 1] <- length(unique(df[[i]][[cloneCall]]))
        Con.df[i, 2] <- names(df)[i]
        Con.df[i, 3] <- length(df[[i]][[cloneCall]])
    }
    if (scale == TRUE) {
        y <- "scaled"
        Con.df$scaled <- Con.df$contigs/Con.df$total * 100
        ylab <- paste("Percent of Unique",inpcloneCall,"Clonotype")
    }
    else {
        y <- "contigs"
        ylab <- paste("Unique",inpcloneCall,"Clonotype")
    }
    Con.df[, x] %<>% factor(levels = gglvl)
    plot <- ggplot(aes(x = Con.df[, x], y = Con.df[, y], fill = as.factor(Con.df[, x])),
                   data = Con.df) + stat_summary(geom = "errorbar",
                                                 fun.data = mean_se,
                                                 position = "dodge", width = 0.5) +
        labs(fill = labs) + stat_summary(fun = mean, geom = "bar",
                                         color = "black", lwd = 0.25) + theme_classic() + xlab("Samples") +
        ylab(ylab) + scale_fill_manual(values = cols)
    return(plot)
}

#' @param ggData The product of combineTCR() into data.frame
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param group The column header used for grouping.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
scabundanceContig <- function(ggData, inptyp, cloneCall = "CTstrict", chain = "both",
                              group = "X",cols,inppts,inpsiz, inpfsz = "Medium",inplsz = "Medium") {
    df <- split(ggData, f = ggData$X)

    Con.df <- NULL
    xlab <- "Clone size"
    names <- names(df)

    for (i in seq_along(df)) {
        if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
        }
        data1 <- df[[i]] %>% group_by(df[[i]][[cloneCall]]) %>%
            summarise(Abundance=n())
        colnames(data1)[1] <- cloneCall
        data1$values <- names[i]

        label <- df[[i]][1,..group]
        data1[,paste(group)] <- label
        Con.df<- rbind.data.frame(Con.df, data1) }
    Con.df <- data.frame(Con.df)
    fill <- group

    ylab <- switch (inptyp,
                    "Clone sizes distribution" = "Clone sizes distribution",
                    "Cumulative clone sizes distribution" = "Cumulative clone sizes distribution"
    )
    Con.df_group = Con.df %>% group_by(X,Abundance) %>% summarise(Count=n())
    Con.df_group %<>% group_by(X) %>%
        mutate(Pk = Count/sum(Count))
    if(inptyp == "Cumulative clone sizes distribution") Con.df_group %<>%
        arrange(desc(Abundance)) %>% group_by(X) %>%
        mutate(cumsum_Pk = cumsum(Pk))



    plot = ggplot(data = Con.df_group,
                  aes_string(x = "Abundance",
                             y = switch (inptyp,"Clone sizes distribution" = "Pk",
                                         "Cumulative clone sizes distribution" = "cumsum_Pk"),
                             color = "X"))+
        geom_point(shape = ifelse(inppts,19,NA),size = inpsiz)+
        scale_color_manual(values=cols)+
        geom_smooth(method = "lm", se=FALSE,size = lList[inplsz]/5) +
        xlab(xlab)+   ylab(ylab)+

        yscale("log10", .format = TRUE)+
        xscale("log10", .format = TRUE)+
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                     parse = TRUE, size = sList[inpfsz]/3)
    return(plot)
    }



#' @param ggData The product of combineTCR() into data.frame
#' @param split The cutpoints for the specific clonotypes.#'
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param group The column header used for grouping.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
scClonalProportion <- function(ggData,Split = c(10, 30, 100, 300,1000,3000, 10000,3000, 100000),
                               cloneCall = "CTstrict", chain = "both",inpcols) {
    df <- split(ggData, f = ggData$X)
    lvl = switch (class(ggData$X),
                  "factor" = levels(ggData$X),
                  "character" = unique(ggData$X)
    )
    Con.df <- NULL
    mat <- matrix(0, length(df), length(Split),
                  dimnames = list(names(df),paste0('[', c(1, Split[-length(Split)] + 1), ':', Split, ']')))
    if (chain != "both") {
        for (x in seq_along(df)) {
            df[[x]] <- off.the.chain(df[[x]], chain, cloneCall)
        }
    }
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    Cut <- c(1, Split[-length(Split)] + 1)
    for (i in seq_along(Split)) {
        mat[,i] <- vapply(df, function (x) sum(na.omit(x[Cut[i]:Split[i]])),
                          FUN.VALUE = numeric(1))
    }
    mat = mat[,colSums(mat) > 0]
    gglvl = colnames(mat)

    mat %<>% as.data.frame()
    mat %<>% tibble::rownames_to_column(var = "X")
    mat_melt <- pivot_longer(as.data.frame(mat),cols = starts_with("["))
    mat_melt$X %<>% factor(levels = lvl)
    mat_melt$name %<>% factor(levels = gglvl)
    if(inpcols == "default") inpcols = "colorBlind"
    plot <- ggplot(mat_melt, aes(x=as.factor(X), y=value, fill=name)) +
        geom_bar(stat = "identity", position="fill",
                 color = "black", lwd= 0.25) +
        scale_fill_manual(name = "",
                          values = color_generator(inpcols,length(gglvl))) +
        xlab("Samples") +
        ylab("TCR Repertoire %") +
        theme_classic()
    return(plot)
}

#Use to shuffle between chains
off.the.chain <- function(dat, chain, cloneCall) {
    chain1 <- toupper(chain) #to just make it easier
    if (chain1 %in% c("TRA", "TRD", "IGH")) {
        x <- 1
    } else if (chain1 %in% c("TRB", "TRG", "IGL")) {
        x <- 2
    } else {
        warning("It looks like ", chain, " does not match the available options for `chain = `")
    }
    dat[[cloneCall]] <- stringr::str_split(dat[[cloneCall]], "_", simplify = TRUE)[,x]
    return(dat)
}

#' @param ggData The product of combineTCR(), combineBCR(),  or expression2List().
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param groupBy Variable in which to group the diversity calculation
#' @param x.axis Additional variable in which to split the x.axis
#' @param exportTable Exports a table of the data into the global environment
#' in addition to the visualization
#' @param n.boots number of bootstraps to downsample in order to get mean diversity
#' @importFrom stringr str_sort str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
#' @author Andrew Malone, Nick Borcherding
scClonalDiversity <- function(ggData, inptyp,cloneCall = "CTstrict", chain = "both", scale=FALSE,
                              inpindex = "Shannon",inp1 = inp1,
                              groupBy = inp2,x.axis = "X", gglvl, cols,inppts,inpsiz,inpfsz,inplsz = "Medium",
                              n.boots = 100) {
    set.seed(101)
    ggData %<>% as.data.frame()
    df <- split(ggData, f = ggData$X)

    Min <- c()
    mat <- NULL
    mat_a <- NULL
    sample <- c()
    df <- bind_rows(df, .id = "element.names")
    df$group.element <- paste0(df[,groupBy], ".", df[,x.axis],".",df[,"paired"])
    group.element.uniq <- unique(df$group.element)
    df[,"group.element"] %<>% factor(levels = group.element.uniq)

    df <- split(df, f = df[,"group.element"])

    Min <- sapply(df,function(x) length(which(!is.na(unique(x[,cloneCall])))))
    keep = Min >= 10
    Min <- min(Min[keep],100) # need to specify a minimal sample_n to make results reproducible

    #Calculating diversity using Vegan R package
    #' @importFrom vegan diversity estimateR
    scdiversityCall <- function(data, index) {

        output <- switch (index,
                          "Shannon" = vegan::diversity(data[,"Freq"], index = "shannon"),
                          "Simpson" = vegan::diversity(data[,"Freq"], index = "simpson"),
                          "Inv.Simpson" = vegan::diversity(data[,"Freq"], index = "invsimpson"),
                          "Chao" = vegan::estimateR(data[,"Freq"])[2], #Chao,
                          "ACE" = vegan::estimateR(data[,"Freq"])[4] #ACE
        )
        return(output)
    }
    df = df[keep]
    for (i in seq_along(df)) {
        data <- as.data.frame(table(df[[i]][,cloneCall]))
        mat_a <- NULL
        sample <- c()

        for (j in seq(seq_len(n.boots))) {
            x <- sample_n(data, Min)
            sample <- scdiversityCall(x, index = inpindex)
            mat_a <- rbind(mat_a, sample)
        }
        mat_a[is.na(mat_a)] <- 0
        mat_a<- colMeans(mat_a)
        mat_a<-as.data.frame(t(mat_a))
        mat <- rbind(mat, mat_a)
    }
    colnames(mat) <- inpindex

    mat[,"grpBy"] <- stringr::str_split(group.element.uniq[keep], "[.]", simplify = TRUE)[,1]
    mat[,x.axis] <- stringr::str_split(group.element.uniq[keep], "[.]", simplify = TRUE)[,2]
    mat[,"paired"] <- stringr::str_split(group.element.uniq[keep], "[.]", simplify = TRUE)[,3]

    inppair = ifelse(inptyp == "Paired Diversity", TRUE,FALSE)

    #if no paired data avaible, remove
    if(inppair) {
        Count = mat %>% group_by(paired) %>% summarise(n = n())
        if(any(Count$n %% length(gglvl) != 0) ){
            keep = Count[Count$n %% length(gglvl) == 0,"paired"]
            mat %<>% filter(paired %in% keep$paired)
            mat %<>% split(f = mat$X) %>% lapply(function(obj) obj[order(obj$grpBy),])
            mat = mat[gglvl]
            mat %<>% bind_rows()
        }
    }

    #rownames(mat) <- names(df)

    mat_melt <- pivot_longer(as.data.frame(mat),cols = inpindex,
                             names_to = "variable") %>% as.data.frame()
    mat_melt[,"X"] %<>% factor(levels = gglvl)

    ggOut = switch (inptyp,
                         "Paired Diversity" = ggpaired(data = mat_melt,
                                                       x = "X",
                                                       y = "value",
                                                       color = "X",
                                                       shape = ifelse(inppts,19,NA),# Extending the regression line
                                                       point.size = inpsiz,
                                                       size = inpsiz/5,
                                                       palette = cols,
                                                       ylab = paste(inpindex, "Index Score"),
                                                       xlab = inp1,
                                                       add = ifelse(inppts,"point","none"),
                                                       line.color = "gray",
                                                       line.size = ifelse(inppair,sList[inplsz]/50,0)
                                                       ),
                         "unpaired Diversity" = ggboxplot(data = mat_melt,
                                                          x = "X",
                                                          y = "value",
                                                          fill = "X",
                                                          shape = ifelse(inppts,19,NA),# Extending the regression line
                                                          palette = cols,
                                                          ylab = paste(inpindex, "Index Score"),
                                                          xlab = inp1,
                                                          outlier.shape = NA,
                                                          line.color = "gray",
                                                          line.size = ifelse(inppair,sList[inplsz]/50,0))+
                             geom_jitter(size = ifelse(inppts,inpsiz,NA), shape = 19)
    )

    if(length(gglvl) ==3){
        comparisons <- list( gglvl[1:2], gglvl[2:3], gglvl[c(1,3)] )
        ggOut = ggOut +stat_compare_means(method = "wilcox.test",paired = inppair,
                                          comparisons = comparisons, size =sList[inpfsz]/4)
    }
    if(length(gglvl) ==2){
        comparisons <- list(gglvl[1:2])
        ggOut = ggOut +stat_compare_means(method = "wilcox.test",paired = inppair,
                                          comparisons = comparisons, size =sList[inpfsz]/4)
    }

    #if (exportTable == TRUE) { return(mat) }
    return(ggOut)
}


#' @param ggData The product of combineTCR() into data.frame
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param inp2 The column header used for grouping.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param x.axis name of the list element to appear on the x.axis
#' @param y.axis name of the list element to appear on the y.axis
#' use for size of dots
ScatterClonotype <- function (ggData, inptyp,inpConf,
                             cloneCall = "CTstrict", chain = "both", scale=inpscale,
                             inp2 = inp2, x.axis = NULL, y.axis = NULL,
                             inplvls,inpcols, inpsiz,inplsz){
    set.seed(101)
    colnames(ggData) %<>% sub(inp2, "grpBy",.)
    ggData = ggData[,c(cloneCall,"X", "grpBy", "tcr.frequency"),with = FALSE]
    ggData = pivot_wider(ggData, names_from = "X", values_from = "tcr.frequency",
                          values_fn = sum)
    ggData[is.na(ggData)] = 0

    ggData %<>% mutate(logFC = log1p(ggData[[y.axis]])-log1p(ggData[[x.axis]])) %>%
                 mutate(class = cut(logFC, breaks=c(-Inf, -inplvls, inplvls, Inf),
                                    labels=c("Contracted","Persistent","Expanded")))
    #ggData[ggData[[y.axis]] == 0, "class",with = FALSE]= "Persistent"
    ggData$class %<>% as.character()
    ggData[(ggData[[x.axis]] == 0 & ggData[["class"]] == "Expanded"), "class",with = FALSE]= "Novel"
    ggData$class %<>% factor(levels = c("Novel", "Persistent","Expanded","Contracted"))

    if(inptyp == "Paired scatter Clonotype") {
        gglvl = levels(ggData$class)
        ggCol = color_generator(inpcols,length(gglvl))
        names(ggCol) = gglvl
        ggCol["Persistent"] = "#CCC9C9"
    }

    if(inptyp == "unpaired scatter Clonotype") {
        gglvl = levels(ggData$grpBy)
        ggData$grpBy %<>% as.character()
        ggData[ggData[["class"]] %in% "Persistent", "grpBy"] = "Persistent"
        ggData$grpBy %<>% factor(c("Persistent", gglvl))

        temp2 <- doFactoring(ggData,inpConf, col = "grpBy", inp2, inpcols)
        gglvl = temp2$gglvl; ggCol = c("#CCC9C9",temp2$ggCol)
        names(ggCol) = gglvl
        ggData = temp2$ggData
    }

    plot <- ggscatter(ggData, x = x.axis, y = y.axis,
                      color = switch (inptyp,
                                      "Paired scatter Clonotype" = "class",
                                      "unpaired scatter Clonotype" = "grpBy"),
                      palette = ggCol[gglvl],
                      size = inpsiz,
                      xlab = paste("\nClone percentage in",x.axis,"(%)"),
                      ylab = paste("Clone percentage in",y.axis,"(%)\n")
                      )+
        theme(legend.title = element_blank()) +
            geom_abline(slope = 1, intercept = 0, alpha = 0.4, lty = 2, size =lList[inplsz]/10) +
        scale_x_continuous(trans = 'log10',#limits=c(0.0001, 10),
                           breaks=10^c(-4:1),
                           labels=c(0.0001,0.001,0.01,0.1,1,10),
                           oob = scales::squish_infinite) +
        scale_y_continuous(trans = 'log10',#limits=c(0.0001, 10),
                           breaks=10^c(-4:1),
                           labels=c(0.0001,0.001,0.01,0.1,1,10),
                           oob = scales::squish_infinite)
    return(plot)
}


ScClonotypeBar <- function(ggData, inptyp, inpConf,
                           cloneCall = cloneCall, chain = "both",
                           x.axis, y.axis,
                           inp2, gglvl2, ggCol2,inpscale,inpcols,
                           inplvls,inppts, inpsiz, inplsz){
    set.seed(101)
    colnames(ggData) %<>% sub(inp2, "grpBy",.)
    ggData = ggData[,c(cloneCall,"X", "grpBy", "tcr.frequency"),with = FALSE]
    ggData = pivot_wider(ggData, names_from = "X", values_from = "tcr.frequency",
                          values_fn = sum)
    ggData[is.na(ggData)] = 0

    ggData %<>% mutate(logFC = log1p(ggData[[y.axis]])-log1p(ggData[[x.axis]])) %>%
        mutate(class = cut(logFC, breaks=c(-Inf, -inplvls, inplvls, Inf),
                           labels=c("Contracted","Persistent","Expanded")))
    #ggData[ggData[[y.axis]] == 0, "class",with = FALSE]= "Persistent"
    ggData$class %<>% as.character()
    ggData[(ggData[[x.axis]] == 0 & ggData[["class"]] == "Expanded"), "class",with = FALSE]= "Novel"
    ggData$class %<>% factor(levels = c("Novel", "Persistent","Expanded","Contracted"))

    ggCol1 = color_generator(inpcols,length(levels(ggData$class)))
    names(ggCol1) = levels(ggData$class)
    ggCol1 = ggCol1[names(ggCol1) !="Persistent"]

    ggData %<>% data.table()
    ggData = ggData[, .(nCells = .N), by = c("class","grpBy")]
    ggData %<>% filter(class != "Persistent")
    ggData = switch(inptyp,
                     "Novel & Expand & Contract Barplot_1" = ggData[, {tot = sum(nCells)
                                                             .SD[,.(tot = tot,
                                                                    nCells = nCells,
                                                                    pctCells = 100 * sum(nCells) / tot),by = "grpBy"]},
                                                             by = "class"],
                     "Novel & Expand & Contract Barplot_2" =  ggData[, {tot = sum(nCells)
                                                             .SD[,.(tot = tot,
                                                                    nCells = nCells,
                                                                    pctCells = 100 * sum(nCells) / tot),by = "class"]},
                                                             by = "grpBy"]
                     )
    ggData$pctCells %<>% round(digits = 0)
    plot <- switch(inptyp,
                   "Novel & Expand & Contract Barplot_1" = ggbarplot(data = ggData,
                                                                   x = "class",
                                                                   y = ifelse(inpscale,"pctCells","nCells"),
                                                                   fill =  "grpBy",
                                                                   lab.size = inpsiz*2,
                                                                   lab.vjust = inpsiz/4,
                                                                   label = inppts, lab.pos = "in", lab.col = "white",
                                                                   width=0.95,
                                                                   position = position_stack(),
                                                                   palette = ggCol2[gglvl2],
                                                                   xlab = "",
                                                                   ylab = ifelse(inpscale,"Clone type percentage (%)\n",
                                                                                 "Number of Clone type\n")),
                   "Novel & Expand & Contract Barplot_2" = ggbarplot(data = ggData,
                                                   x = "grpBy",
                                                   y = ifelse(inpscale,"pctCells","nCells"),
                                                   fill = "class",
                                                   lab.size = inpsiz*2,
                                                   lab.vjust = inpsiz/4,
                                                   label = inppts, lab.pos = "in", lab.col = "white",
                                                   width=0.95,
                                                   position = position_stack(),
                                                   palette = ggCol1,
                                                   xlab = paste0("\n",inp2),
                                                   ylab = ifelse(inpscale,"Clone type percentage (%)\n",
                                                                 "Number of Clone type\n")))
    plot <- plot + theme(legend.title = element_blank())
    return(plot)
}

# strSplitClean(input$sc1n1inp1)
strSplitClean <- function(text){
    gsub(" ","",unlist(stringr::str_split(text, pattern = ",")))
}
