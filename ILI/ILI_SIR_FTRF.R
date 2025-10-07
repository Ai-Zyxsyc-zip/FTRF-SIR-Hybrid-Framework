
library(dplyr)
library(ggplot2)
library(deSolve)
library(randomForestSRC)
library(rstudioapi)

set.seed(1234)

# 1. 设置工作目录为当前脚本所在路径
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# 2. 读取和预处理数据
Data <- read.csv("ILINet.csv", skip = 1)
Data <- Data[, -c(1, 5, 7, 8, 9, 10, 11, 12, 15)]
Data <- Data[, c(1:3, 6, 5, 4)]

# 筛选出现频率大于700的地区
all_region <- table(Data$REGION) %>% as.data.frame()
tar_region <- all_region[which(all_region$Freq > 700),]
region_name <- tar_region$Var1 %>% as.character()
datalist <- lapply(region_name, function(v) { Data[Data$REGION == v, ] })

# 选取前两个地区数据
num_region <- 2
covs <- c("YEAR", "WEEK", "NUM..OF.PROVIDERS", "ILITOTAL", "X.UNWEIGHTED.ILI")
tar_list <- lapply(1:num_region, function(k){
  tem <- datalist[[k]][, covs]
  tem$X.UNWEIGHTED.ILI <- as.numeric(as.character(tem$X.UNWEIGHTED.ILI)) / 100
  colnames(tem) <- c(paste0('X', 1:4), 'Y')
  tem <- tem %>% mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}"))
  return(tem)
})

# 设置目标与辅助地区数据
target_data <- tar_list[[1]]  # 目标地区，比如 Alabama
aux_data <- tar_list[[2]]     # 辅助地区，比如 Alaska

target_data <- target_data[order(target_data$X1, target_data$X2), ]
aux_data <- aux_data[order(aux_data$X1, aux_data$X2), ]

# 3. 定义 SIR 微分方程模型
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# 4. 拟合辅助站点 SIR 模型并预测感染人数序列
fit_sir_predict <- function(obs_data, time_points) {
  I_obs <- obs_data$Y
  N <- 1
  I0 <- I_obs[1]
  S0 <- N - I0
  R0 <- 0
  init_state <- c(S = S0, I = I0, R = R0)
  
  loss_function <- function(params) {
    names(params) <- c("beta", "gamma")
    out <- ode(y = init_state, times = time_points, func = sir_model, parms = params)
    pred_I <- as.data.frame(out)$I
    sum((I_obs - pred_I)^2)
  }
  
  init_params <- c(beta = 0.2, gamma = 0.1)
  lower_bounds <- c(beta = 0, gamma = 0)
  upper_bounds <- c(beta = 2.0, gamma = 1.0)
  
  opt <- optim(par = init_params, fn = loss_function, method = "L-BFGS-B",
               lower = lower_bounds, upper = upper_bounds)
  
  best_params <- opt$par
  out <- ode(y = init_state, times = time_points, func = sir_model, parms = best_params)
  return(as.data.frame(out)$I)
}

# 5. 滑动窗口预测参数设置
window_size <- 4
step <- 4
n <- nrow(target_data)

predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric())

# 6. 滑动窗口循环：训练FTRF模型并预测
for (start in seq(1, n - window_size * 2, by = step)) {
  train_idx <- start:(start + window_size - 1)
  test_idx <- (start + window_size):(start + 2 * window_size - 1)
  if (max(test_idx) > n) break
  
  train_data <- target_data[train_idx, ]
  test_data <- target_data[test_idx, ]
  
  aux_window <- aux_data[train_idx, ]
  
  # 辅助站点SIR拟合与预测，作为辅助特征
  hatI_train <- fit_sir_predict(aux_window, seq_len(nrow(train_data)))
  hatI_test <- fit_sir_predict(aux_window, seq_len(nrow(test_data)))
  
  ftrf_train <- train_data
  ftrf_train$hatI <- hatI_train
  
  ftrf_test <- test_data
  ftrf_test$hatI <- hatI_test
  
  # 训练随机森林模型
  rf_model <- rfsrc(Y ~ ., data = ftrf_train, ntree = 100)
  
  # 预测感染人数
  pred_I <- predict(rf_model, newdata = ftrf_test)$predicted
  
  predicted_all <- rbind(predicted_all, data.frame(time = test_idx, I = pred_I, obs = test_data$Y))
}

# 7. 绘图展示预测结果
ili <- ggplot(predicted_all, aes(x = time)) +
  geom_point(aes(y = obs, color = "Observed"), size = 2) +
  geom_line(aes(y = I, color = "Predicted"), size = 1.2) +
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  labs(title = "Alabama SIR-FTRF Prediction Results", x = "Time/Week", y = "Infection Proportion", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white", color = "black"))

print(ili)

# 8. 保存图像
ggsave("ili_ftrf_sir.png", plot = ili, width = 5, height = 5, dpi = 300, bg = "white")

# 9. 计算并输出误差指标
mse <- mean((predicted_all$obs - predicted_all$I)^2)
r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)

cat("MSE:", mse, "\n")
cat("R^2:", r2, "\n")
