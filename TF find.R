rm(list=ls())
options(stringsAsFactors = F)
source('functions_heatmap.R')
gs=read.csv('NRRA and NR induced gene.csv',row.names = 1)
markerGenes=read.csv('TF.csv',header = FALSE)[,1]
#library(stringr)
#markerGenes=str_to_title(markerGenes, locale = "en")
gene_subset <- as.matrix((gs[rownames(gs) %in% markerGenes,]))
gs=as.matrix(gene_subset)
write.csv(gs,file = "TF gene expression.csv")
annotation_colors = list(Cluster = c(Veh="#560047", RA="#a53bad", NR="#eb6bac", NRRA="#ffa8a0"))

cold <- colorRampPalette(c('#41b6c4','#253494','#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58'))
warm <- colorRampPalette(c('#fecc5c','#e31a1c','#800026','#800026','#800026','#800026','#800026','#800026'))
mypalette <- c(rev(cold(20)), warm(20))

#gs=t(scale(t(gs[cg,])))
gs=t(scale(t(gs)))
gs[gs>2]=2
gs[gs< -2]= -2
gs[1:4,1:4]

annotation_col = data.frame(
  Cluster = factor(rep(c("Veh","RA","NR","NRRA"),c(2,2,2,2))))
rownames(annotation_col) = colnames(gs)
library(pheatmap)
pheatmap(gs,fontsize = 4,
         #scale="row",
         cellheight = 4,
         show_colnames=FALSE, 
         show_rownames=TRUE, 
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         #cutree_rows = 6,
         method = "ward.D",
         color=viridis(10),
         annotation_col=annotation_col,
         annotation_colors=annotation_colors,
         #color=mypalette,
         # gaps_row=rowbreaks,
         border_color = FALSE,
         filename="NRRA and NR induced TF.pdf")