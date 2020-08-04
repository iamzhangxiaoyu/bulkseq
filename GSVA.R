library(DESeq2)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
createList <- function(group=NULL) {
  
  tumorsam <- names(group)
  sampleList = list()
  treatsamList =list()
  treatnameList <- c()
  ctrlnameList <- c()
  
  #A-1: 类1 vs 其他
  sampleList[[1]] = tumorsam
  treatsamList[[1]] = intersect(tumorsam, names(group[group=="CD1WT_Veh"])) 
  treatnameList[1] <- "CD1WT_Veh" 
  ctrlnameList[1] <- "Others" 
  
  #A-2: 类2 vs 其他
  sampleList[[2]] = tumorsam
  treatsamList[[2]] = intersect(tumorsam, names(group[group=="CD1WT_T"]))
  treatnameList[2] <- "CD1WT_T"
  ctrlnameList[2] <- "Others"
  
  #A-3: 类3 vs 其他
  sampleList[[3]] = tumorsam
  treatsamList[[3]] = intersect(tumorsam, names(group[group=="WT_Veh"]))
  treatnameList[3] <- "WT_Veh"
  ctrlnameList[3] <- "Others"
  
  #A-4: 类4 vs 其他
  sampleList[[4]] = tumorsam
  treatsamList[[4]] = intersect(tumorsam, names(group[group=="WT_T"]))
  treatnameList[4] <- "WT_T"
  ctrlnameList[4] <- "Others"
  
  
  
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}


twoclassDESeq2 <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) { 
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]] 
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="") 
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { 
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    
    saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    
    
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula("~ Type")) 
    
    dds$Type <- relevel(dds$Type,ref = "control")
    
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    
    resData <- as.data.frame(res[order(res$padj),])
    resData$id <- rownames(resData)
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}
expr <- read.table("easy_input_counts.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
expr[1:3, 1:3]
subt <- read.table("easy_input_subtype.txt", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
n.sub.label <- unique(subt$TCGA_Subtype) 
n.sub.label
n.sub <- length(table(subt$TCGA_Subtype)) 
n.sub

group <- subt$TCGA_Subtype
names(group) <- rownames(subt) 
complist <- createList(group=group)




var <- apply(expr,1,sd)
index <- var > quantile(var)[4]

twoclassDESeq2(res.path = ".", 
               countsTable = expr[index,intersect(colnames(expr),rownames(subt))],
               prefix = "SKCM", 
               complist = complist,
               overwt = F)


DEfiles <- c("SKCM_deseq2_test_result.CD1WT_Veh_vs_Others.txt",
             "SKCM_deseq2_test_result.CD1WT_T_vs_Others.txt",
             "SKCM_deseq2_test_result.WT_Veh_vs_Others.txt",
             "SKCM_deseq2_test_result.WT_T_vs_Others.txt")

degs.list <- list()
for (i in 1:n.sub) {
  degs <- read.table(DEfiles[i],sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
  head(degs)
  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(degs))
}


library(clusterProfiler) 
library(GSVA)
library(pheatmap)
library(gplots)
subtype_specific_gsea <- function(msigdb=NULL,n.top=10,mode=c("up","down"),degs.list=NULL,subtype.label=NULL,nPerm.gsea=1000,minGSSize.gsea=10,maxGSSize.gsea=500,pvalueCutoff.gsea=1){
  
  MSigDB <- read.gmt(msigdb)
  GSEA.list <- top.gs <- list() 
  
  if(!is.element(mode, c("up", "dn"))) { stop("mode must be up or dn!\n") }
  
  for (i in 1:n.sub) {
    degs <- degs.list[[n.sub.label[i]]]
    geneList <- degs$log2FC; names(geneList) <- rownames(degs)
    geneList <- sort(geneList,decreasing = T) 
    
    cat(paste0("GSEA for ",subtype.label[i]," starts!\n"))
    GSEA.list[[subtype.label[i]]] <- GSEA(geneList = geneList,
                                          TERM2GENE=MSigDB,
                                          nPerm = nPerm.gsea,
                                          minGSSize = minGSSize.gsea,
                                          maxGSSize = maxGSSize.gsea,
                                          seed = T,
                                          verbose = F,
                                          pvalueCutoff = pvalueCutoff.gsea) 
    
    GSEA.dat <- as.data.frame(GSEA.list[[subtype.label[i]]])
    
    if(mode == "up") {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = T),] 
    } else {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = F),] 
    }
    
    
    write.table(GSEA.dat,paste0(subtype.label[[i]],"_degs_",mode,"_gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    
    
    top.gs[[subtype.label[i]]] <- rownames(GSEA.dat)[1:n.top] 
  }
  
  
  gs <- list()
  for (i in as.character(unlist(top.gs))) {
    gs[[i]] <- MSigDB[which(MSigDB$ont %in% i),"gene"]
  }
  
  return(list(mode=mode,top.gs=top.gs,gs=gs))
}
library("AnnotationDbi")
library("annotate")
library("stats4")
library("IRanges")
library("S4Vectors")




msigdfFile = "c5.all.v7.1.symbols.gmt"
n.top = 10
mode = "up" #"up"和"dn"二选一
gs.up <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)
gsva_gs.up <- gsva(as.matrix(expr), gs.up$gs, method="gsva") 
dim(gsva_gs.up)
gsva_gs.up_mean <- data.frame(row.names = rownames(gsva_gs.up)) 
for (i in n.sub.label) {
  gsva_gs.up_mean <- cbind.data.frame(gsva_gs.up_mean,
                                      data.frame(rowMeans(gsva_gs.up[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.up_mean) <- n.sub.label
jco <- c("#F2CCCC","#E6D8CF","#D5E3F0","#FDE7DA","#E2D6EC", "#CCEFDB")

annRows <- data.frame(subtype=rep(n.sub.label,each=n.top), names = unlist(gs.up$top.gs), stringsAsFactors = F)
annRows <- annRows[!duplicated(annRows$names),]; rownames(annRows) <- annRows$names # 倘若出现一条通路在>=2个亚型中上调，去掉重复值，这种情况在亚型较多的时候会发生


annColors <- list(subtype=c("CD1WT_Veh"=jco[1],"CD1WT_T"=jco[2],"WT_Veh"=jco[3],"WT_T"=jco[4]))

filename <- paste0("subtype_specific_top_",mode,"_gsea.pdf")
pheatmap(gsva_gs.up_mean[rownames(annRows),],
         cellwidth = 10, cellheight = 10,
         #color = bluered(64), #自定义颜色
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA, #如果想要边框，就去掉这行
         annotation_row = annRows[,"subtype",drop = F],
         annotation_colors = annColors,
         filename = filename)
mode = "dn"
gs.dn <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)
gsva_gs.dn <- gsva(as.matrix(expr), gs.dn$gs, method="gsva") 
gsva_gs.dn_mean <- data.frame(row.names = rownames(gsva_gs.dn)) 
for (i in n.sub.label) {
  gsva_gs.dn_mean <- cbind.data.frame(gsva_gs.dn_mean,
                                      data.frame(rowMeans(gsva_gs.dn[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.dn_mean) <- n.sub.label
annRows <- data.frame(subtype=rep(n.sub.label,each=n.top), names = unlist(gs.dn$top.gs), stringsAsFactors = F)
annRows <- annRows[!duplicated(annRows$names),]; rownames(annRows) <- annRows$names 


annColors <- list(subtype=c("CD1WT_Veh"=jco[1],"CD1WT_T"=jco[2],"WT_Veh"=jco[3],"WT_T"=jco[4]))

filename <- paste0("subtype_specific_top_",mode,"_gsea.pdf")
pheatmap(gsva_gs.dn_mean[rownames(annRows),],
         cellwidth = 10, cellheight = 10,
         #color = bluered(64), #自定义颜色
         border_color = NA, #如果想要边框，就去掉这行
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = annRows[,"subtype",drop = F],
         annotation_colors = annColors,
         filename = filename)

