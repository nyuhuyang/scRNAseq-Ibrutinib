####################
## Density heatmap from Becht, E.Nat Biotechnol 2018 UMAP
## with Gaussian smoothing
####################
# Figure 6G
library(Seurat)
library(RColorBrewer)
library(magrittr)
library(smoothie)
library(fields)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

object = readRDS("data/OSU_SCT_20210821.rds")
umap = object@reductions$umap@cell.embeddings
 
Baseline_cellID <- colnames(object)[object$timepoint %in% "C1D-7"]
treated_cellID <- colnames(object)[object$timepoint %in% "C1D+1"]
N =1000
umap1_counts=log10(1+table(cut(umap[Baseline_cellID,1],breaks=N),
                           cut(umap[Baseline_cellID,2],breaks=N)))
umap2_counts=log10(1+table(cut(umap[treated_cellID,1],breaks=N),
                           cut(umap[treated_cellID,2],breaks=N)))
for(k in 8){
    print(k)
    umap1_counts = kernel2dsmooth(umap1_counts, kernel.type="gauss", nx=N, ny=N, sigma=k)
    umap2_counts = kernel2dsmooth(umap2_counts, kernel.type="gauss", nx=N, ny=N, sigma=k)
    
    jpeg(paste0(path,"Density heatmap_log10_sigma=",k,".jpeg"),height=1000,width=2000,res=300)
    par(mfrow=c(1,2))
    par(mar=c(1,1,3,3),bty="l")
    image(umap1_counts,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)),breaks=seq(0,max(umap1_counts),length.out=102),xaxt="n",yaxt="n")
    title(main="C1D-7",cex=0.4)
    x=seq((par("usr")[1]+4*par("usr")[2])/5,par("usr")[2],length.out=51)
    yb=rep((par("usr")[4]-par("usr")[3])*0.85+par("usr")[3],50)
    yt=rep((par("usr")[4]-par("usr")[3])*0.9,+par("usr")[3],50)
    rect(x[-length(x)],yb,x[-1],yt,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)),border=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)))
    text(x=mean(c(x[1],tail(x,1))),y=yt[1],pos=3,labels="Number of events\nper surface unit",xpd=NA,cex=0.8)
    text(x=x[1],y=yb[1],pos=1,labels="0",xpd=NA,cex=0.7)
    text(x=tail(x,1),y=yb[1],pos=1,labels="max(UMAP)",xpd=NA,cex=0.7)
    
    image(umap2_counts,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)),breaks=seq(0,max(umap2_counts),length.out=102),xaxt="n",yaxt="n")
    title(main="C1D+1",cex=0.4)
    rect(x[-length(x)],yb,x[-1],yt,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)),border=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)))
    text(x=mean(c(x[1],tail(x,1))),y=yt[1],pos=3,labels="Number of events\nper surface unit",xpd=NA,cex=0.8)
    text(x=x[1],y=yb[1],pos=1,labels="0",xpd=NA,cex=0.7)
    text(x=tail(x,1),y=yb[1],pos=1,labels="max(UMAP)",xpd=NA,cex=0.7)
    dev.off()
    
}
