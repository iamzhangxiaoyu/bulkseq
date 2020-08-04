##############################PCA#####################################
rm(list = ls())
options(stringsAsFactors = F)
gs=read.csv('DE_genes_groups.csv',row.names = 1)
dat=as.matrix(gs)
group_list = as.character(rep(c("Veh","RA","NR","NRRA"),c(2,2,2,2)))
table(group_list)
dat[1:4,1:4]
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list)
library("FactoMineR")
library("factoextra") 
library("ggplot2")
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')

##############################DEG#####################################
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图

bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
dev.off()
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
dim(dat)

