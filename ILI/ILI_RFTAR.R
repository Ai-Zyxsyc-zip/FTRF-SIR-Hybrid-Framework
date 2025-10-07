library(randomForestSRC)
set.seed(1234)
# 滑动窗口预测

# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# 读取数据
Data <- read.csv("ILINet.csv", skip = 1)
Data <- Data[, -c(1, 5, 7, 8, 9, 10, 11, 12, 15)]
Data <- Data[, c(1:3, 6, 5, 4)]

# 统计频率大于700的地区
all_region <- table(Data$REGION) %>% as.data.frame()
tar_region <- all_region[which(all_region$Freq > 700),]
region_name <- tar_region$Var1 %>% as.character()
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })

#选择地区数量
num_region <- 2
covs <- c("YEAR", "WEEK", "NUM..OF.PROVIDERS", "ILITOTAL", "X.UNWEIGHTED.ILI")
tar_list <- lapply(1:num_region, function(k){
  tem <- datalist[[k]][, covs]
  colnames(tem) <- c(paste0('X', 1:4), 'Y')
  tem <- tem %>% mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  return(tem)
})

# 设置地区名称变量
region_label <- "Alabama"
window_size <- 4
step <- 4
n <- nrow(data_fit)

predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric(), R2 = numeric())

for (start in seq(1, n - window_size*2, by = step)) {
  train_indices <- start:(start + window_size - 1)
  test_indices <- (start + window_size):(start + 2*window_size - 1)
  
  if (max(test_indices) > n) break
  
  # 确保数据没有缺失值
  train_data <- data_fit[train_indices, c("X1", "X2", "X3", "X4")]
  train_Y <- data_fit$Y[train_indices]
  if (any(is.na(train_data)) || any(is.na(train_Y))) {
    next
  }
  
  test_data <- data_fit[test_indices, c("X1", "X2", "X3", "X4")]
  actual_Y <- data_fit$Y[test_indices]
  if (any(is.na(test_data)) || any(is.na(actual_Y))) {
    next
  }
  
  # 合并训练数据与目标变量
  train_data_with_y <- cbind(train_data, Y = train_Y)
  
  # 训练 RF 模型
  rf_model <- tryCatch({
    rfsrc(Y ~ X1 + X2 + X3 + X4, data = train_data_with_y, ntree = 500)
  }, error = function(e) {
    message("Error in training RF model: ", e$message)
    return(NULL)
  })
  
  # 检查是否模型训练成功
  if (is.null(rf_model)) {
    message("Skipping this iteration due to model training failure.")
    next
  }
  
  # 预测
  pred_rf <- predict(rf_model, newdata = test_data)
  pred_Y <- pred_rf$predicted
  
  # 计算 R²
  mean_actual_Y <- mean(actual_Y)
  ss_total <- sum((actual_Y - mean_actual_Y)^2)
  ss_residual <- sum((actual_Y - pred_Y)^2)
  R2 <- 1 - (ss_residual / ss_total)
  
  # 保存结果
  predicted_all <- rbind(predicted_all, data.frame(time = test_indices, I = pred_Y, obs = actual_Y, R2 = R2))
}

# # 绘图对象
# ili <- ggplot(predicted_all, aes(x = time)) +
#   geom_point(aes(y = obs, color = "实际观测"), size = 2) +
#   geom_line(aes(y = I, color = "预测"), size = 1.2) +
#   scale_color_manual(values = c("预测" = "blue", "实际观测" = "red")) +
#   labs(
#     title = paste0(region_label, " RF-tar预测结果"),
#     x = "时间/周",
#     y = "感染比例",
#     color = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = c(0.95, 0.95),
#     legend.justification = c("right", "top"),
#     legend.background = element_rect(fill = "white", color = "black")
#   )


# Plot object
ili <- ggplot(predicted_all, aes(x = time)) +
  geom_point(aes(y = obs, color = "Observed"), size = 2) +
  geom_line(aes(y = I, color = "Predicted"), size = 1.2) +
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  labs(
    title = paste0(region_label, " RF-tar Prediction Results"),
    x = "Time/Week",
    y = "Infection Proportion",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )



# 展示图片
ili
# 保存为白底 PNG 文件，文件名中也使用地区名
ggsave(
  filename = paste0("ili_rf_tar.png"),
  plot = ili,
  width = 5, height = 5, dpi = 300, bg = "white"
)

# 计算 MSE（均方误差）
mse <- mean((predicted_all$obs - predicted_all$I)^2)

# 输出 MSE
cat("MSE:", mse, "\n")

r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)
cat("R^2:", r2, "\n")
