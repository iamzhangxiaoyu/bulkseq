a=read.csv('46 Pluripotency genes.csv',header = FALSE)[,1]
library(stringr)
a=str_to_title(a, locale = "en")
write.csv(a,file='46 Pluripotency genes.csv')

gs=read.csv('DE_genes_groups.csv',row.names = 1)
a=read.table(paste0("combine.txt"),stringsAsFactors=FALSE)[,1]
length(a)
b <- a[which(a %in% rownames(gs))]
annotation_colors = list(Cluster = c(Veh="#560047", RA="#a53bad", NR="#eb6bac", NRRA="#ffa8a0"))

#c= gs[b,]
#write.csv(c, file = "165.csv")

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
source("functions_heatmap.R")
library(pheatmap)
pheatmap(gs[b,],fontsize = 6,
         #scale="row",
         show_colnames=FALSE, 
         show_rownames=FALSE, 
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         color=viridis(10),
         annotation_col=annotation_col,
         annotation_colors=annotation_colors,
         #color=mypalette,
         # gaps_row=rowbreaks,
         border_color = FALSE,
         filename="DE_undiff marker genes.pdf")
set.seed(123)
gene_clustering <- pheatmap::pheatmap(
  gs[b,],fontsize = 8,
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
          file="undiff_marker_4hclust.csv")