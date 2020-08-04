library(ggplot2)
library(dplyr)
library(plyr)
expr_df <- read.csv(file='easy_input_expr.csv',row.names = 1, 
                    header = TRUE, sep=",", stringsAsFactors = FALSE)
meta_df <- read.csv(file='easy_input_meta.csv', row.names = 1,
                    header = TRUE, sep=",",stringsAsFactors = FALSE)
#查看前3个基因在前4个sample中的表达矩阵
expr_df[1:3,1:4]
pca.results <- prcomp(expr_df, center = TRUE, scale. = FALSE)

#定义足够多的颜色，用于展示分组
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

pca.results <- read.csv(file='PCA.csv',row.names = 1, 
                        header = TRUE, sep=",", stringsAsFactors = FALSE)

pca.rotation <- pca.results$rotation
pca.rotation
pca.pv <- summary(pca.results)$importance[2,]
pca.pv

low_dim_df <- as.data.frame(pca.results$x[,c(1,2)])
low_dim_df$group <- meta_df$group
#查看前3行
low_dim_df
write.csv(pca.results, file = "PCA.csv")
add_ellipase <- function(p, x="PC1", y="PC2", group="group",
                         ellipase_pro = 0.95,
                         linetype="dashed",
                         colour = "black",
                         lwd = 2,...){
  obs <- p$data[,c(x, y, group)]
  colnames(obs) <- c("x", "y", "group")
  ellipse_pro <- ellipase_pro
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(obs, 'group', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
  })
  names(ell)[2:3] <- c('x', 'y')
  
  ell <- ddply(ell, .(group) , function(x) x[chull(x$x, x$y), ])
  p <- p + geom_polygon(data = ell, aes(x=x,y=y,group = group), 
                        colour = colour,
                        alpha = 1,fill = NA,
                        linetype=linetype,
                        lwd =lwd)
  return(p)
}


#计算坐标轴标签
pc1.pv <- paste0(round(pca.pv['PC1'],digits = 3) * 100, "%")
pc2.pv <- paste0(round(pca.pv['PC2'],digits = 3) * 100, "%")

#画出各个样本在二维空间的点
p <- ggplot(pca.results) + 
  geom_point(aes(x=PC1, y=PC2, shape=group), size=10, #点的大小
             #shape=c('1','2'),#点的形状
             alpha=0.5) +#设置点为半透明，出现叠加的效果
  #如果使用默认的颜色，就在下面这行前面加个#
  scale_color_manual(values = mycol[1:length(unique(meta_df$group))]) +
  #还能调整整体的颜色亮度
  #scale_colour_hue(l=45) + 
  theme_bw() + #去除背景色
  scale_shape(solid = FALSE)+
  
  
  #图例
  guides(color=guide_legend(title = NULL)) +
  theme(legend.background = element_blank(), #移除整体边框
        #图例的左上角置于绘图区域的左上角
        legend.position = c(0,1),legend.justification = c(0,1),
        legend.text = element_text(size=12)) + #字体大小
  
  #调整坐标轴标签
  xlab(paste0("PC1" )) + 
  ylab(paste0("PC2" )) 
p

ggsave('PCA.pdf',width = 5.5,height = 4)