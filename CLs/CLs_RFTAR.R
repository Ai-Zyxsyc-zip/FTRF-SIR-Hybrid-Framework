library(rstudioapi)
library(readr)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
set.seed(1234)
# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# 读取数据，去除第1行，和1, 6, 7, 8,9,10列的信息
Data <- read.csv("WHO_NREVSS_Clinical_Labs.csv", skip = 1)
Data <- Data[, -c(1, 6, 7, 8,9,10)]
# 删除第四列值为'X'的行
Data <- Data %>% filter(Data[, 4] != "X")

# 统计每个地区出现的频率
all_region <- table(Data$REGION) %>% as.data.frame()

# 筛选频次大于400次的地区
tar_region <- all_region[which(all_region$Freq > 400),]
region_name <- tar_region$Var1 %>% as.character()

# 重新排列数据列表
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })

# 选择前 num_region 个地区
num_region <- 2  # 假设选择2个地区

covs <- c("YEAR", "WEEK", "TOTAL.SPECIMENS")

tar_list <- lapply(1:num_region, function(k){
  tem <- datalist[[k]][, covs]
  
  # 随机抽取数据
  index <- sample(1:nrow(tem), size = floor(1 * nrow(tem)))
  tem <- tem[index, ]
  
  # 重命名列
  colnames(tem) <- c(paste0('X', 1:2), 'Y')
  
  # 转换为数字类型
  tem <- tem %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  
  return(tem)
})

set.seed(1234)
# 构建目标地区数据
data_fit <- tar_list[[1]]
data_fit <- data_fit[order(data_fit$X1, data_fit$X2), ]
data_fit$Y <- data_fit$Y

# 滑动窗口 + RF建模预测
window_size <- 4
step <- 4
n <- nrow(data_fit)

predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric())

for (start in seq(1, n - window_size * 2, by = step)) {
  train_indices <- start:(start + window_size - 1)
  test_indices <- (start + window_size):(start + 2 * window_size - 1)
  
  if (max(test_indices) > n) break
  
  train_data <- data_fit[train_indices, c("X1", "X2", "Y")]
  test_data <- data_fit[test_indices, c("X1", "X2")]
  actual_Y <- data_fit$Y[test_indices]
  
  rf_model <- tryCatch({
    rfsrc(Y ~ ., data = train_data, ntree = 500)
  }, error = function(e) {
    message("训练失败：", e$message)
    return(NULL)
  })
  
  if (is.null(rf_model)) next
  
  pred_rf <- predict(rf_model, newdata = test_data)
  pred_Y <- pred_rf$predicted
  
  predicted_all <- rbind(predicted_all,
                         data.frame(time = test_indices, I = pred_Y, obs = actual_Y))
}

# # 可视化
# cls <- ggplot(predicted_all, aes(x = time)) +
#   geom_point(aes(y = obs, color = "实际"), size = 2) +
#   geom_line(aes(y = I, color = "预测"), size = 1.2) +
#   scale_color_manual(values = c("预测" = "blue", "实际" = "red")) +
#   labs(
#     title = "Alabama RF-tar预测结果",
#     x = "时间/周",
#     y = "感染人数",
#     color = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = c(0.95, 0.95),
#     legend.justification = c("right", "top"),
#     legend.background = element_rect(fill = "white", color = "black")
#   )

# Visualization
cls <- ggplot(predicted_all, aes(x = time)) +
  geom_point(aes(y = obs, color = "Observed"), size = 2) +
  geom_line(aes(y = I, color = "Predicted"), size = 1.2) +
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  labs(
    title = "Alabama RF-tar Prediction Results",
    x = "Time/Week",
    y = "Number of Infected Individuals",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )



# 展示与保存图像
cls
ggsave("cls_rf_tar.png", plot = cls,
       width = 5, height = 5, dpi = 300, bg = "white")

# MSE和R²
mse <- mean((predicted_all$obs - predicted_all$I)^2)
r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)
cat("MSE:", mse, "\n")
cat("R^2:", r2, "\n")