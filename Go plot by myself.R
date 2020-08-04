library(ggplot2)
enrich <- read.csv("GO_term_BP_gene_ND vs NDRA_down.csv",
                     header=T,stringsAsFactors=F,comment.char="",quote="\"")  


enrich1 <- enrich[order(enrich$qvalue),]   #对富集结果按照qvalue进行从小到大排序，保证最显著的通路在前
enrich2 <- enrich1[1:10,]                            #这里画图只展示top10的通路
count <- as.numeric(unlist(strsplit(enrich2$GeneRatio,"/2185",fixed=T))) #提取每条通路里面差异表达的基因数
enrich3 <- data.frame(enrich2[,3],count,enrich2[,7])
colnames(enrich3) <- c("Description","count","qvalue")

p <- ggplot(data=enrich3,aes(x=Description,y=count,fill=qvalue))   #fill=qvalue  fill颜色填充，使用连续值qvalue
p1 <- p + geom_bar(stat="identity") + coord_flip()            #coord_flip()颠倒坐标轴
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=20))

p3 <- p2 + ylim(0,150) + scale_fill_gradient(low="red",high="blue")   #ylim(0,200) 更改横坐标的范围这里坐标轴颠                                                                                                                  倒了，虽然看起来是x轴，但其实是y轴

p4 <- p3 + scale_x_discrete(limits=rev(enrich3[,1])) +labs(x="",y="",title="Biological_Process")

#输出为png格式的图片
#png("S01_S03_S05_vs_S02_S04_S06_Biological_Process_enrich.png",width=680,height=480)
#print(p4)
#dev.off()

#输出为pdf的文件
pdf("ND vs NDRA_down.pdf",width=10, height = 10)
print(p4)
dev.off()
