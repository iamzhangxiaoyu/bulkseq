gs=read.csv('DE_genes_groups copy.csv',row.names = 1)
library(stringr)
batch=str_split(colnames(gs),'_',simplify = T)[,2] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(batch)
library(sva)
ex_b_sva = ComBat(dat=as.matrix(gs), 
                  batch=batch 
)
write.csv(ex_b_sva, file = "batch remove.csv")
