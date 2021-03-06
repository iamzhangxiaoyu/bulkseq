rm(list = ls()) 
options(stringsAsFactors = F)
library(DESeq2)
# 创建需要配对比较的列表
createList <- function(group=NULL) {
  
  tumorsam <- names(group)
  sampleList = list()
  treatsamList =list()
  treatnameList <- c()
  ctrlnameList <- c()
  
  #A-1: 类1 vs 其他
  sampleList[[1]] = tumorsam
  treatsamList[[1]] = intersect(tumorsam, names(group[group=="WT_Veh"])) # 亚型名称需要根据情况修改
  treatnameList[1] <- "WT_Veh" # 该亚型的命名
  ctrlnameList[1] <- "Others" # 其他亚型的命名
  
  #A-2: 类2 vs 其他
  sampleList[[2]] = tumorsam
  treatsamList[[2]] = intersect(tumorsam, names(group[group=="STRA8KO_Veh"]))
  treatnameList[2] <- "STRA8KO_Veh"
  ctrlnameList[2] <- "Others"
  
  
  #如果有更多类，按以上规律继续写
  
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}

# 配对DESeq2
twoclassDESeq2 <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) { # 循环读取每一次比较的内容
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]] 
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="") # 生成最终文件名
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { # 因为差异表达分析较慢，因此如果文件存在，在不覆盖的情况下（overwt=F）不再次计算差异表达
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    
    saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    
    # 差异表达过程，具体参数细节及输出结果解释，请参阅相关document
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula("~ Type")) # 设计矩阵仅包含亚型信息，若有批次效应请修改
    
    dds$Type <- relevel(dds$Type,ref = "control")
    
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    
    resData <- as.data.frame(res[order(res$padj),])
    resData$id <- rownames(resData)
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")
    #输出到文件
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}


# 读取read count表达矩阵
expr <- read.table("easy_input_counts copy.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
expr[1:3, 1:3]

# 读取亚型信息
subt <- read.table("easy_input_subtype copy.txt", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)

n.sub.label <- unique(subt$TCGA_Subtype) # 亚型名称
n.sub.label

n.sub <- length(table(subt$TCGA_Subtype)) # 亚型个数
n.sub


### 创建配对比较的列表信息 ###
group <- subt$TCGA_Subtype
names(group) <- rownames(subt) 
complist <- createList(group=group)

### 执行配对DESeq2差异表达 ###

# 差异表达分析过程比较慢请耐心等待，这里为了加速分析过程，只选取方差top25%的基因


twoclassDESeq2(res.path = ".", #所有配对差异表达结果都会输出在res.path路径下
               countsTable = expr[row.names(expr),intersect(colnames(expr),rownames(subt))],
               prefix = "SKCM", #文件名以SKCM开头
               complist = complist,
               overwt = F)


DEfiles <- c("SKCM_deseq2_test_result.WT_Veh_vs_Others.txt",
             "SKCM_deseq2_test_result.STRA8KO_Veh_vs_Others.txt")

degs.list <- list()
for (i in 1:n.sub) {
  degs <- read.table(DEfiles[i],sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
  head(degs)
  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(degs))
}
