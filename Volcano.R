#使用国内镜像安装包
#options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

setwd("C:/Users/xiaoyu/Desktop/RNAseq")

# 全部基因差异表达分析结果
x <- read.csv("easy_input_limma.csv", row.names = 1)
x$label<- rownames(x)
head(x)

# 突出展示感兴趣的基因
selectedGeneID <- read.csv("easy_input_selected.csv")
head(selectedGeneID)

# 提取感兴趣的基因的差异分析结果
x$gsym <- row.names(x)
selectgenes <- merge(selectedGeneID, x, by = "gsym")
head(selectgenes)

#plot_mode <- "classic" #经典版
plot_mode <- "advanced" #酷炫版

logFCcut <- 1.5 #log2-foldchange
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value

#for advanced mode
logFCcut2 <- 2.5
logFCcut3 <- 5
pvalCut2 <- 0.0001
pvalCut3 <- 0.00001

#置x，y軸的最大最小位置
xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))
ymin <- 0
ymax <- -log10(x$P.Value)[3] * 25.1

# 基因名的颜色，需大于等于pathway的数量，这里自定义了足够多的颜色
mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

if (plot_mode == "classic"){
  # 簡單的setting for color
  x$color_transparent <- ifelse((x$P.Value < pvalCut & x$logFC > logFCcut), "red", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "blue","grey"))
  # 簡單的setting for size
  size <- ifelse((x$P.Value < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # 複雜的的setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey", n1)
  names(cols)<- rownames(x)
  
  #不同阈值的点的颜色
  cols[x$P.Value < pvalCut & x$logFC >logFCcut]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # 複雜的的setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  #不同阈值的点的大小
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 6
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 6
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(P.Value), label = label, color = pathway)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10), #刻度线的位置
    labels = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10),
    limits = c(-11, 11) #x轴范围，两侧对称才好看
  ) +
  #或用下面这行：
  #xlim(c(xmin, xmax)) + 
  
  #画阈值分界线
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + #虚线的形状和粗细
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
  ) +
  theme(panel.grid=element_blank())

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}
p1

# 显示 logFC > n 的基因的基因名
n = 9
p1 + geom_text_repel(aes(x = logFC, y = -log10(P.Value), 
                         label = ifelse(logFC > n, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

# 突出显示候选基因
p2 <- p1 + 
  # 在感兴趣的基因外面画个黑色圈
  geom_point(data = selectgenes, alpha = 1, size = 4.6, shape = 1, 
             stroke = 1, #圈粗细
             color = "black") +
  
  # 显示感兴趣的基因的基因名
  scale_color_manual(values = mycol) + 
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, #不显示图例
                  size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  guides(color=guide_legend(title = NULL)) 

p2

# 显示pathway
np <- length(unique(selectgenes$pathway))
(labelsInfo <- data.frame(pathway = names(table(selectgenes$pathway)),
                          col = mycol[1:np]))

p2 + annotation_custom(tableGrob(labelsInfo$pathway, rows = c(rep("", np)), cols = "",
                                 theme = ttheme_minimal(base_colour = labelsInfo$col)),
                       ymin = ymax - 2, ymax = ymax, xmin = xmin - 1.5, xmax = xmin)

# 保存到PDF文件
if (plot_mode == "classic"){ggsave("volcano_classic.pdf", width=6,height=5)} else if (plot_mode == "advanced") {ggsave("Volcano_advanced.pdf",width=12,height=10)} else {stop("Unsupport mode")}