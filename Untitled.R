library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
a=read.csv('107 meiosis gene.csv',header = FALSE)[,1]
library(stringr)
a=str_to_title(a, locale = "en")
write.csv(a,file = "107 meiosis gene.csv")

dat=read.csv('DE_genes_groups.csv',row.names = 1)
dat=as.matrix(dat)
a=read.table(paste0("meiosis gene.txt"),stringsAsFactors=FALSE)[,1]
length(a)
b <- a[which(a %in% rownames(dat))]
length(b)
gene_subset <- colSums(log(dat[b, ]))/length(b)
gene_subset=t(gene_subset)
gene_subset<- as.data.frame(gene_subset)
write.csv(gene_subset,file = 'percent.meiosis.csv')


