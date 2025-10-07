set.seed(1234)  # 请根据需要修改随机种子

# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# 读取数据，去除第1行，和1,6,7,8,9,10,11,12,13,14列的信息
Data <- read.csv("WHO_NREVSS_Combined_prior_to_2015_16.csv", skip = 1)
Data <- Data[, -c(1,6,7,8,9,10,11,12,13,14)]

# 删除第4列值为'X'的行
Data <- Data %>% filter(Data[, 4] != "X")

# 统计每个地区出现的频率
all_region <- table(Data$REGION) %>% as.data.frame()

# 筛选频次大于200次的地区
tar_region <- all_region[which(all_region$Freq > 200), ]
region_name <- tar_region$Var1 %>% as.character()

# 重新排列
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })

# 选择地区数量
num_region <- 2

# 滑动窗口参数
window_size <- 4  # 窗口大小
step <- 4         # 步长

# 为每个选定地区进行数据预处理和模型训练
tar_list <- lapply(1:num_region, function(k) {
  tem <- datalist[[k]][, c("YEAR", "WEEK", "TOTAL.SPECIMENS")]
  
  # 随机抽取数据
  index <- sample(1:nrow(tem), size = floor(1 * nrow(tem)))
  tem <- tem[index, ]
  
  # 移除常数列
  constant_columns <- sapply(tem, function(col) length(unique(col)) == 1)
  tem <- tem[, !constant_columns]
  
  # 重命名列：最后一列命名为Y，其余列依次命名为X1, X2, X3, ...
  colnames(tem) <- c(paste0('X', 1:(ncol(tem) - 1)), 'Y')
  
  # 将所有列转换为数字类型
  tem <- tem %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  
  return(tem)
})

# 选择第一个地区的数据进行滑动窗口训练
data_fit <- tar_list[[1]]

# 按YEAR（X1）和WEEK（X2）排序
data_fit <- data_fit[order(data_fit$X1, data_fit$X2), ]

# 进行滑动窗口预测
predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric(), R2 = numeric())

for (start in seq(1, n - window_size*2, by = step)) {
  # 定义训练集和测试集索引
  train_indices <- start:(start + window_size - 1)
  test_indices <- (start + window_size):(start + 2*window_size - 1)
  
  # 如果测试集索引超出了数据集范围，则停止循环
  if (max(test_indices) > n) break
  
  # 确保训练数据和测试数据没有缺失值
  train_data <- data_fit[train_indices, c("X1", "X2")]  # 只选择X1和X2作为特征
  train_Y <- data_fit$Y[train_indices]
  if (any(is.na(train_data)) || any(is.na(train_Y))) {
    next
  }
  
  test_data <- data_fit[test_indices, c("X1", "X2")]  # 只选择X1和X2作为特征
  actual_Y <- data_fit$Y[test_indices]
  if (any(is.na(test_data)) || any(is.na(actual_Y))) {
    next
  }
  
  # 合并训练数据与目标变量
  train_data_with_y <- cbind(train_data, Y = train_Y)
  
  # 训练 RF 模型
  rf_model <- tryCatch({
    rfsrc(Y ~ X1 + X2, data = train_data_with_y, ntree = 10)
  }, error = function(e) {
    message("Error in training RF model: ", e$message)
    return(NULL)
  })
  
  # 如果模型训练失败，跳过该轮迭代
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
# cpt <- ggplot(predicted_all, aes(x = time)) +
#   geom_point(aes(y = obs, color = "实际观测"), size = 2) +
#   geom_line(aes(y = I, color = "预测"), size = 1.2) +
#   scale_color_manual(values = c("预测" = "blue", "实际观测" = "red")) +
#   labs(
#     title = "Alabama RF-tar预测结果",
#     x = "时间/周",
#     y = "感染数量",
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
cpt <- ggplot(predicted_all, aes(x = time)) +
  geom_point(aes(y = obs, color = "Observed"), size = 2) +
  geom_line(aes(y = I, color = "Predicted"), size = 1.2) +
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  labs(
    title = "Alabama RF-tar Prediction Results",
    x = "Time/Week",
    y = "Number of Infections",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )



# 展示图
cpt

# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
# 保存为白底 PNG 文件
ggsave("cpt_rf_tar.png", plot = cpt,
       width = 5, height = 5, dpi = 300, bg = "white")



# 计算 MSE（均方误差）
mse <- mean((predicted_all$obs - predicted_all$I)^2)

# 输出 MSE
cat("MSE:", mse, "\n")

# 计算 R²
r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)
cat("R^2:", r2, "\n")
