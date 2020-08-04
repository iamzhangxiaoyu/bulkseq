library(ggplot2)
library(ggrepel)
x <- read.csv("Veh vs NR.csv", row.names = 1)
pval_thr = 1e-3
fc_thr=log2(2)
# Identify the Up and Down regulated genes
x$isDE = "NO"

pos = which(x$logFC > fc_thr & x$FDR < pval_thr)
x$isDE[pos] <- "Up"

pos = which(x$logFC < -fc_thr & x$FDR < pval_thr)
x$isDE[pos] <- "Down"


UpGenes <-  rownames(subset(x, isDE=="Up"))
DownGenes <-  rownames(subset(x, isDE=="Down"))

#write.csv(x,file = "NRRA vs NR_DEG.csv")


m1= x$Sample_Veh
m2= x$Sample_NR


df <- data.frame(m1= m1,m2= m2)
df = log2(df+0.1)
rownames(df) <- rownames(x)


df$isDE= "NO"
df[UpGenes,"isDE"] <- "Up"
df[DownGenes,"isDE"] <- "Down"


markers_genes <- df[c("ENSMUSG00000007035","ENSMUSG00000020059","ENSMUSG00000022432","ENSMUSG00000036928","ENSMUSG00000005493",
                      "ENSMUSG00000051977","ENSMUSG00000068117","ENSMUSG00000022429","ENSMUSG00000021245","ENSMUSG00000085601",
                      "ENSMUSG00000086022","ENSMUSG00000028109",
                      "ENSMUSG00000051965","ENSMUSG00000018698","ENSMUSG00000034384","ENSMUSG00000041359","ENSMUSG00000056155",
                      "ENSMUSG00000066687","ENSMUSG00000024565","ENSMUSG00000021379","ENSMUSG00000024406","ENSMUSG00000000317",
                      "ENSMUSG00000013089","ENSMUSG00000025089"),]
markers_genes$gene <- rownames(markers_genes)

markers_genes$gene <- gsub("ENSMUSG00000007035","Msh5",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000020059","Sycp3",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000022432","Smc1b",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000036928","Stag3",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000005493","Msh4",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000051977","Prdm9",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000068117","Mei1",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000022429","Dmc1",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000021245","Mlh3",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000085601","Gm4969",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000086022","Rad51ap2",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000028109","Hormad1",markers_genes$gene)

markers_genes$gene <- gsub("ENSMUSG00000051965","Nanos2",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000018698","Lhx1",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000034384","Barhl2",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000041359","Tcl1",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000056155","Nanos3",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000066687","Zbtb16",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000024565","Sall3",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000021379","Id4",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000024406","Pou5f1",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000000317","Bcl6b",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000013089","Etv5",markers_genes$gene)
markers_genes$gene <- gsub("ENSMUSG00000025089","Gfra1",markers_genes$gene)

nbDown = sum(x$isDE=="Down")
nbUp = sum(x$isDE=="Up")
lbls = data.frame(m1=c(10,-2),
                  m2= c(-2,15),
                  isDE= c("Down","Up"),
                  count=c( paste0("Down:", nbDown), 
                           paste0("Up:", nbUp))
)




ggplot(df, aes(x=m1, y=m2,colour=isDE)) + geom_point(size=0.6) +
  theme_bw() + 
  scale_color_manual(values = c("royalblue","grey80","firebrick")) + 
  xlab("log2[Veh expression]") + ylab("log2[NR expression]") + 
  ggtitle("NR vs Veh") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_point(mapping = aes(x=m1, y=m2),markers_genes, color="black",size=2) + 
  geom_text_repel(data = markers_genes,mapping = aes(x=m1,y=m2,label=gene),color="black",segment.colour ="black", segment.size =0.5,force = 5) +
  geom_text_repel(data = lbls,mapping = aes(x=m1,y=m2,color=isDE,label=count)) +
  geom_abline(slope = 1,intercept = log2(1.2) * c(-1,1))
ggsave(filename = 'NR vs Veh.tiff',width = 12, height =8)
