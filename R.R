rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor

setwd('/Users/zhangxiaoyu/Desktop/stra8 ko testis')
a=read.table('raw count.txt',row.names="gene_id",
             header = T ,sep = '\t')
dat=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),]
dat=log2(edgeR::cpm(dat)+1) 

hc=hclust(dist(t(dat)))
class(hc)
plot(hc,labels = FALSE)
clus = cutree(hc, 2) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
group_list= as.factor(clus) ##转换为因子属性
table(group_list) ##统计频数


colnames(dat) #取列名
library(stringr)
replicate=str_split(colnames(dat),'_',simplify = T)[,1] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(replicate)
n_g = apply(a,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）
# 这里我们定义， reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
# 而且大部分单细胞转录组技术很烂，通常超过75%的基因都没办法检测到。

df=data.frame(g=group_list,replicate=replicate,n_g=n_g) #新建数据框(细胞的属性信息)

##(样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)

df$all='all' #添加列，列名为"all"，没事意思，就是后面有需要
metadata=df
save(a,dat,df,file = 'input.Rdata')

#heatmap
load(file = 'input.Rdata')
a[1:4,1:4]
head(metadata) #head()函数显示操作前面的信息，默认前6行
metadata=df
group_list=metadata$g #'$'符，取列，取metadata矩阵的g列,取出层级聚类信息
table(group_list) ##这是全部基因集的聚类分组信息
cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
##画热图,针对top100的sd的基因集的表达矩阵,没有聚类分组
pheatmap(dat[cg,],show_colnames =F,show_rownames = F,
         filename = 'all_cells_top_1000_sd.png')
if(T){
  
  x=1:10;plot((x))
  scale(x);plot(scale(x))
  
  n=t(scale(t(dat[cg,]))) #scale()函数去中心化和标准化
  #对每个探针的表达量进行去中心化和标准化
  n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
  n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
  n[1:4,1:4]
  
  # pheatmap(n,show_colnames =F,show_rownames = F)
  
  ac=data.frame(g=group_list) #制作细胞（样本）分组矩阵
  rownames(ac)=colnames(n) ##ac的行名（样本名）等于n的列名（样本名）
  ##判断分组矩阵的行（样本数）和表达矩阵的列（样本数）是否相等
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,
           filename = 'all_cells_top_1000_sd_cutree1.png')
  
}

if(T){
  
  
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  
  ##这个聚类分组只是对top100的sd的基因集
  hc=hclust(dist(t(n))) 
  clus = cutree(hc, 2)
  group_list=as.factor(clus)
  table(group_list) ##这个聚类分组信息是针对top100的sd的基因集的，和针对全部基因集的分组结果不一样
  table(group_list,metadata$g) ## 其中 metadata$g 是前面步骤针对全部表达矩阵的层次聚类结果。
  
  ## 下面针对本次挑选100个基因的表达矩阵的层次聚类结果进行热图展示。
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,
           cellwidth = 10, cellheight = 0.5, fontsize = 8,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
           border_color = "NA",
           show_colnames =F,show_rownames = F,
           annotation_col=ac,
           filename = 'all_cells_top_1000_sd_cutree_2.png')
  dev.off() ##关闭画板
  
}


#PCA
dat_back=dat ##防止下面操作把数值搞坏的一个备份
dat=dat_back  ##表达矩阵数据
dat[1:4,1:4]
dat=t(dat)
dat=as.data.frame(dat) ##转换为数据框
dat=cbind(dat,group_list ) ##cbind()合并列（横向追加）;添加分组信息
dat[1:4,1:4]
## 表达矩阵可以随心所欲的取行列，基础知识需要打牢。
dat[1:4,12197:12199]
dat[,ncol(dat)] #ncol()列，返回列长值
table(dat$group_list)
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ) is removed
# before PCA analysis
## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #'-'表示“非”
fviz_pca_ind(dat.pca,repel =T,
             geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
             col.ind = dat$group_list, # color by groups 颜色组
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses 集中成椭圆
             legend.title = "Groups"
)
## 事实上还是有很多基因dropout非常严重。
ggsave('all_cells_PCA.png')

#gene number

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'input.Rdata')
a[1:4,1:4]
head(df) 

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
# 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
group_list=df$g
replicate=df$replicate
table(replicate)
n_g = apply(a,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）

### 绘图必须强推 ggpubr 
## 搜索到代码，然后修改即可出美图。
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'input_rpkm.Rdata')
#   n_g = apply(a,2,function(x) sum(x>0)) #统计每个样本有表达的有多少行（基因）

dat[1:4,1:4]
df=metadata
head(df) 
library(ggpubr)
ggviolin(df, x = "all", y = "n_g", fill = "all", 
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "replicate", y = "n_g", fill = "replicate",
         #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "g", y = "n_g", fill = "g", 
         add = "boxplot", add.params = list(fill = "white"))  + stat_compare_means()


#cell cycle
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'input.Rdata')
a[1:4,1:4]
head(df) 

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。

group_list=df$g
replicate=df$replicate
table(replicate)

a[1:6,1:6]
library(scran)
# https://mp.weixin.qq.com/s/nFSa5hXuKHrGu_othopbWQ
sce <- SingleCellExperiment(list(counts=dat)) 
#list() 创建列表

library(org.Mm.eg.db)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))

ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), 
                  keytype="SYMBOL", column="ENSEMBL")
#取探针名创建一个向量
#rownames(sce) 取行名（即实验检测到的基因）

if(F){
  assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
  save(assigned,file = 'cell_cycle_assigned.Rdata')
}
load(file = 'cell_cycle_assigned.Rdata')
head(assigned$scores)
table(assigned$phases)
draw=cbind(assigned$score,assigned$phases) #合并assigned$score列和assigned$phases列
colnames(draw)
attach(draw)
library(scatterplot3d)
scatterplot3d(G1, S, G2M, angle=20,
              color = rainbow(3)[as.numeric(as.factor(assigned$phases))],
              grid=TRUE, box=FALSE)
detach(draw)

library(pheatmap)
cg=names(tail(sort(apply(dat,1,sd)),100))
n=t(scale(t(dat[cg,])))
pheatmap(n,show_colnames =F,show_rownames = F)
library(pheatmap)
df$cellcycle=assigned$phases
ac=df
rownames(ac)=colnames(n)
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,
         filename = 'all_cells_top_100_sd_all_infor.png')
dev.off()
head(ac)
table(ac[,c(1,5)])


#DEG
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'input.Rdata')
a[1:4,1:4]
head(df)
table(df$g)
## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞数量）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。
exprSet=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] 
exprSet=exprSet[!grepl('ERCC',rownames(exprSet)),]
group_list=ifelse(df$g==1,'KO','WT')
table(group_list)
library(edgeR)
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('KO-WT'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='KO-WT', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
}

write.csv(DEG_limma_voom, file = "DEG_limma_voom.csv", quote = F, row.names = T)
boxplot(dat['Shank1',]~df$g)


group_list=ifelse(df$g==1,'KO','WT')
table(group_list)
library(edgeR)
do_limma_RNAseq <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('KO-WT'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='KO-WT', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
  return(DEG_limma_voom)
}
deg1=do_limma_RNAseq(exprSet,group_list)
boxplot(dat['Shank1',]~df$g)
boxplot(dat['Prelid1',]~df$g)
boxplot(dat[rownames(deg1)[1],]~df$g)

group_list=ifelse(df$g==2,'me','other');table(group_list)
deg2=do_limma_RNAseq(exprSet,group_list)
boxplot(dat[rownames(deg2)[1],]~df$g)
boxplot(dat[rownames(deg2)[2],]~df$g)

group_list=ifelse(df$g==3,'me','other');table(group_list)
deg3=do_limma_RNAseq(exprSet,group_list)

group_list=ifelse(df$g==4,'me','other');table(group_list)
deg4=do_limma_RNAseq(exprSet,group_list)

deg1=deg1[order(deg1$logFC,decreasing = T),]

deg2=deg2[order(deg2$logFC,decreasing = T),]
deg3=deg3[order(deg3$logFC,decreasing = T),]
deg4=deg4[order(deg4$logFC,decreasing = T),]

cg=c(head(rownames(deg1),18))
     head(rownames(deg2),18),
     head(rownames(deg3),18),
     head(rownames(deg4),18)
)
library(pheatmap)
g=df$g
mat=dat[cg,]
mat=mat[,order(g)]
ac=data.frame(group=g)
rownames(ac)=colnames(dat)
mat=mat[head(rownames(deg1),100),]

n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]

pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         annotation_col = ac)



plot(deg1$logFC,-log10(deg1$adj.P.Val))
with(deg1,plot( logFC,-log10( adj.P.Val)))
with(deg2,plot( logFC,-log10( adj.P.Val)))
with(deg3,plot( logFC,-log10( adj.P.Val)))
with(deg4,plot( logFC,-log10( adj.P.Val)))

diff1=rownames(deg1[abs(deg1$logFC)>3,])
diff2=rownames(deg2[abs(deg2$logFC)>3,])
diff3=rownames(deg3[abs(deg3$logFC)>3,])
diff4=rownames(deg4[abs(deg4$logFC)>3,])

library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
df <- bitr(diff1, fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
kk  <- enrichKEGG(gene         = df$ENTREZID,
                  organism     = 'mmu', 
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9)
head(kk)[,1:6]
dotplot(kk)




























