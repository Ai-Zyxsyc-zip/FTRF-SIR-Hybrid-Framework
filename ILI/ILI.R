library(rstudioapi) #设置工作路径
library(readr)
library(dplyr)
library(randomForestSRC)
library(glmnet)
library(ggplot2)
library(ROCR)
library("KRLS")


# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)


#读取数据，去除第1行，和1, 5, 7, 8, 9, 10, 11, 12列的信息
Data <- read.csv("ILINet.csv", skip = 1)
Data <- Data[, -c(1, 5, 7, 8, 9, 10, 11, 12, 15)]
#Data
# 交换Data数据框中的第四列和第六列
Data <- Data[, c(1:3, 6, 5, 4)]

#统计每个地区出现的频率
all_region <- table(Data$REGION) %>% as.data.frame()
#all_region


#筛选频次大于700次的地区
tar_region <- all_region[which(all_region$Freq > 700),]
#tar_region
region_name <- tar_region$Var1 %>% as.character()
#region_name


#重新排列
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })
#datalist

#重新命名
covs     <- c("YEAR", "WEEK", "NUM..OF.PROVIDERS",
              "ILITOTAL", "X.UNWEIGHTED.ILI")


#受限于电脑性能的调整

#选择地区数量
num_region <- 2

tar_list <- lapply(1:num_region, function(k){
  tem <- datalist[[k]][, covs]
 
   # 随机抽取10%的数据
  index <- sample(1:nrow(tem), size = floor(0.2 * nrow(tem)))
  tem <- tem[index, ]
  
  colnames(tem) <-c(paste0('X', 1:4), 'Y')
  tem <- tem %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  
  return(tem)
})


for(kk in 1:num_region){
  simulate(nreps = 100, nsite = num_region, nvars = 4, site = kk)
}

kk <- 1
y <- read.csv(paste0('mse_of_ILI_', kk , ".csv")) 
y <- y[, -c(3, 8, 9)]
y <- as.matrix(y) %>% as.vector()
x <- rep(c('RF-tar', "FTRF", "SR", "GMA", "RF-agg", "TLRF", "TLKRR"), each = 50) 
x <- factor(x, levels = c('RF-tar', "FTRF", "SR", "GMA", "RF-agg", "TLRF", "TLKRR"))
Data <- data.frame(X = x, Y = y)
p <- ggplot(Data, aes(x = X, y = Y,fill=X)) + 
  labs(title = 'ILI', x = NULL, y = "MSE") +
  stat_boxplot(geom = "errorbar", width = 0.2, aes(x = X,y = Y)) +
  geom_boxplot(width=0.4, outlier.shape = NA) +
  scale_fill_manual(values=c("#7BDFF2","#B2F7EF", "#bdc0dc","#EFF7F6","#F7D6E0", "#F2B5D4", "#c6e114"))+
  scale_x_discrete(labels =  c('RF-tar', "FTRF", expression(SR^{'+'}), expression(GMA^{'+'}),  "RF-agg", "TL-RF", "TL-KRR"))+
  theme_bw() + 
  theme(axis.title.y=element_text(size=20),
        axis.text = element_text(size = 20)) +
  # scale_y_continuous(limits=c(0,2.2),breaks=seq(0,2,1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "mm")) +
  theme(axis.title=element_text(size= 20,  face = "bold"),
        axis.text.y = element_text(size = 19,  face = "bold"),
        axis.text.x = element_text(size = 19,  face = "bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA, color="black", size= 2,  linetype="solid")) + 
  theme(legend.position= "none")
p
ggsave(paste0("ILI.png"), p, width = 8, height = 6, dpi = 300)
ggsave(paste0("ILI_", kk, ".png"), p, width = 8, height = 6, dpi = 300)


