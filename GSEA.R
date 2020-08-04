library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(cowplot)



# 要画的通路
geneSetID <- c("mmu04140", "mmu04142")

# 突出显示感兴趣的基因
selectedGeneID <- c("Atg5")

# 自定义足够多的颜色
mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

gsym.fc <- read.table("easy_input_rnk.txt", header = T)
dim(gsym.fc)
head(gsym.fc)
#把gene symbol转换为ENTREZ ID
#此处物种是人，其他物种的ID转换方法，请参考FigureYa52GOplot
gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

head(gsym.id)
dim(gsym.id)

#让基因名、ENTREZID、foldchange对应起来
gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
head(gsym.fc.id)

#按照foldchange排序
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
head(gsym.fc.id.sorted)

#获得ENTREZID、foldchange列表，做为GSEA的输入
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID
head(id.fc)

#查看clusterProfiler用法
browseVignettes("clusterProfiler")

#这一条语句就做完了KEGG的GSEA分析
kk <- gseKEGG(id.fc, organism = "mmu",nPerm  = 1000,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,)
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
write.table(sortkk,"gsea_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)

for (i in geneSetID) {
  gseaplot(kk, i)
  myGeneList <- enrichplot:::gsInfo(kk, i)
  row.names(myGeneList) <- gsym.fc$gsym
  myGeneList$id <- gsym.fc$ENTREZID 
  write.csv(myGeneList, paste0("gsea_genelist_", i, "_group1.csv"))
}

x <- kk
geneList <- position <- NULL ## to satisfy codetool

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,2)

# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  
  #scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) + #不画网格
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
#p.res


# 画rank
rel_heights <- c(1.5, .5, 1.5) # 上中下三个部分的比例

i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}
#head(gsdata)

p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + 
  scale_color_manual(values = mycol) + #用自定义的颜色
  
  theme_bw() +
  theme(panel.grid = element_blank()) + #不画网格
  
  theme(legend.position = "none",
        plot.margin = margin(t=-.1, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  #scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))
#p2


# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
#head(df2)

# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]
head(selectgenes)



p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  
  #scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw() +
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        panel.grid = element_blank()) +
  
  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
#p.pos


# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)
