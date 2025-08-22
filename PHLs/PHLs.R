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


#读取数据，去除第1行，和1, 3列的信息
Data <- read.csv("WHO_NREVSS_Public_Health_Labs.csv", skip = 1)
Data <- Data[, -c(1, 3, 11 ,12 )]
#Data


# 删除第四列值为'X'的行
Data <- Data %>% filter(Data[, 4] != "X")
# 删除第4列到第11列的值都为0的行
Data <- Data %>% filter(!apply(Data[, 4:8], 1, function(x) all(x == 0)))
#数值型
Data <- Data %>%
  mutate_at(vars(2:ncol(Data)), ~ as.numeric(as.character(.)))


#统计每个地区出现的频率
all_region <- table(Data$REGION) %>% as.data.frame()
#all_region


#筛选频次大于8次的地区
tar_region <- all_region[which(all_region$Freq > 8),]
#tar_region
region_name <- tar_region$Var1 %>% as.character()
#region_name


#重新排列
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })
#datalist
#

# 合并每 5 个地区的数据为一组
datalist <- lapply(seq(1, length(datalist), by = 10), function(i) {
  # 获取每 5 个数据的子集
  region_group <- datalist[i:min(i + 4, length(datalist))]  # 防止超出范围
  # 将这些地区的数据合并成一个数据框
  do.call(rbind, region_group)
})

colnames <- lapply(datalist, function(df) colnames(df)[3:8])
#colnames
#colnames[[1]]
#重新命名
#covs     <-   c("TOTAL.SPECIMENS","A..2009.H1N1.","A..H3.","A..Subtyping.not.Performed.")
covs     <-   c("TOTAL.SPECIMENS", "BYam", "BVic", "B")
#covs

#受限于电脑性能的调整

#选择地区数量
num_region <- 3

tar_list <- lapply(1:num_region, function(k){
  tem <- datalist[[k]][, covs]
 
   # 随机抽取20%的数据
  index <- sample(1:nrow(tem), size = floor(1 * nrow(tem)))
  tem <- tem[index, ]
  
  colnames(tem) <-c(paste0('X', 1:3), 'Y')
  tem <- tem %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  
  return(tem)
})


for(kk in 1:num_region){
  simulate(nreps = 100, nsite = num_region, nvars = 3, site = kk)
}


kk <- 3
y <- read.csv(paste0('mse_of_PHLs_', kk , ".csv")) 
y <- y[, -c(3, 8, 9)]
y <- as.matrix(y) %>% as.vector()
x <- rep(c('RF-tar', "FTRF", "SR", "GMA", "RF-agg", "TLRF", "TLKRR"), each = 50) 
x <- factor(x, levels = c('RF-tar', "FTRF", "SR", "GMA", "RF-agg", "TLRF", "TLKRR"))
Data <- data.frame(X = x, Y = y)
p <- ggplot(Data, aes(x = X, y = Y,fill=X)) + 
  labs(title = 'PHLs_B', x = NULL, y = "MSE") +
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

ggsave(paste0("PHLs_", kk, ".png"), p, width = 8, height = 6, dpi = 300)


