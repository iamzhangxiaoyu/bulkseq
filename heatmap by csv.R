#library(data.table)
gs=read.csv('DE_genes_groups copy.csv',row.names  =1)
gs=as.data.table(gs[!duplicated(gs$X), ])
#write.csv(gs,file = "DE_genes_groups.csv")
source("functions_heatmap.R")
gs=read.csv('NRRA induce(FC>1).csv',row.names = 1)
#cg=names(tail(sort(apply(gs,1,sd)),5000))

annotation_colors = list(Cluster = c(Veh="#560047", RA="#a53bad", NR="#eb6bac", NRRA="#ffa8a0"))

cold <- colorRampPalette(c('#41b6c4','#253494','#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58'))
warm <- colorRampPalette(c('#fecc5c','#e31a1c','#800026','#800026','#800026','#800026','#800026','#800026'))
mypalette <- c(rev(cold(20)), warm(20))

#gs=t(scale(t(gs[cg,])))
gs=t(scale(t(gs)))
gs[gs>2]=2
gs[gs< -1.5]= -1.5 
gs[1:4,1:4]

annotation_col = data.frame(
  Cluster = factor(rep(c("Veh","RA","NR","NRRA"),c(2,2,2,2))))
rownames(annotation_col) = colnames(gs)
library(pheatmap)
pheatmap(gs,fontsize = 8,
         #scale="row",
         show_colnames=FALSE, 
         show_rownames=FALSE, 
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
         filename="DE_NRRA induce(FC>1).pdf")
set.seed(123)
gene_clustering <- pheatmap::pheatmap(
  gs,fontsize = 8,
  #scale="row",
  show_colnames=FALSE, 
  show_rownames=FALSE, 
  cluster_cols=FALSE,
  cluster_rows=TRUE,
  method = "ward.D",
  color=viridis(10),
  annotation_col=annotation_col,
  annotation_colors=annotation_colors,
  #color=mypalette,
  # gaps_row=rowbreaks,
  border_color = FALSE,
  silent=TRUE
)

clusters <- cutree(gene_clustering$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
head(clustering)
table(clustering[,1])
write.csv(clustering, 
          file="hclust.csv")