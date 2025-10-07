library(rstudioapi) # 设置工作路径
library(readr)
library(dplyr)
library(randomForestSRC)
library(glmnet)
library(ggplot2)
library(ROCR)
library(KRLS)
library(deSolve)
set.seed(1234)

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

#选择地区
data_fit <- tar_list[[1]]
data_fit <- data_fit[order(data_fit$X1, data_fit$X2), ]  # 先按YEAR（X1），再按WEEK（X2）排序
data_fit$Y <- data_fit$Y / 100

# SIR 模型函数
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# 参数拟合函数
fit_parameters <- function(par, data_time, data_I) {
  S0 <- par[1]
  beta <- par[2]
  gamma <- par[3]  # gamma也参与优化
  
  init <- c(S = S0, I = data_I[1], R = 1 - S0 - data_I[1])
  times <- data_time
  
  out <- tryCatch({
    ode(y = init, times = times, func = sir_model, parms = c(beta = beta, gamma = gamma))
  }, error = function(e) return(NULL))
  
  if (is.null(out)) return(Inf)
  
  pred_I <- out[, "I"]
  
  # 确保预测值没有 NA 或无穷大
  if (any(is.na(pred_I)) || any(is.infinite(pred_I))) {
    return(Inf)
  }
  
  # 计算误差
  return(sum((pred_I - data_I)^2))
}

# 更新参数范围
opt_result <- optim(par = c(0.99, 0.2, 0.1), 
                    fn = fit_parameters, 
                    data_time = train_time, 
                    data_I = train_I,
                    method = "L-BFGS-B", 
                    lower = c(0.01, 0.001, 0.01), 
                    upper = c(1, 1, 1))

# 滑动窗口预测
window_size <- 4
step <- 4
n <- nrow(data_fit)

predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric())

for (start in seq(1, n - window_size*2, by = step)) {
  train_indices <- start:(start + window_size - 1)
  test_indices <- (start + window_size):(start + 2*window_size - 1)
  
  if (max(test_indices) > n) break
  
  train_time <- 1:length(train_indices)
  train_I <- data_fit$Y[train_indices]
  
  opt_result <- optim(par = c(0.99, 0.2, 0.1), fn = fit_parameters, data_time = train_time, data_I = train_I,
                      method = "L-BFGS-B", lower = c(0.01, 0.001, 0.001), upper = c(1, 1, 1))
  
  S0_opt <- opt_result$par[1]
  beta_opt <- opt_result$par[2]
  gamma_opt <- opt_result$par[3]
  
  init <- c(S = S0_opt, I = train_I[length(train_I)], R = 1 - S0_opt - train_I[length(train_I)])
  test_time <- 1:length(test_indices)
  sir_output <- as.data.frame(ode(y = init, times = test_time, func = sir_model, parms = c(beta = beta_opt, gamma = gamma_opt)))
  
  pred_I <- sir_output$I
  actual_I <- data_fit$Y[test_indices]
  predicted_all <- rbind(predicted_all, data.frame(time = test_indices, I = pred_I, obs = actual_I))
}


# # 绘图对象
# ili <- ggplot(predicted_all, aes(x = time)) +
#   geom_point(aes(y = obs, color = "实际观测"), size = 2) +
#   geom_line(aes(y = I, color = "预测"), size = 1.2) +
#   scale_color_manual(values = c("预测" = "blue", "实际观测" = "red")) +
#   labs(
#     title = paste0(region_label, " SIR预测结果"),
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
    title = paste0(region_label, " SIR Prediction Results"),
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
  filename = paste0("ili_sir.png"),
  plot = ili,
  width = 5, height = 5, dpi = 300, bg = "white"
)


# 计算 MSE（均方误差）
mse <- mean((predicted_all$obs - predicted_all$I)^2)

# 输出 MSE
cat("MSE:", mse, "\n")
r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)
cat("R^2:", r2, "\n")