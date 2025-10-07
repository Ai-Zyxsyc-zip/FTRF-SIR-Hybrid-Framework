library(dplyr)
library(ggplot2)
library(deSolve)
library(randomForestSRC)

set.seed(1234)

# SIR 微分方程函数定义
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

fit_sir_predict <- function(obs_data, time_points) {
  I_obs <- obs_data$Y
  N <- 1
  I0 <- I_obs[1]
  S0 <- N - I0
  R0 <- 0
  init_state <- c(S = S0, I = I0, R = R0)
  
  # SSE损失函数
  loss_function <- function(params) {
    names(params) <- c("beta", "gamma")
    out <- ode(y = init_state, times = time_points, func = sir_model, parms = params)
    pred_I <- as.data.frame(out)$I
    sum((I_obs - pred_I)^2)
  }
  
  # 初始值 & 边界设置
  init_params <- c(beta = 0.2, gamma = 0.1)
  lower_bounds <- c(beta = 0, gamma = 0)
  upper_bounds <- c(beta = 2.0, gamma = 1.0)
  
  # 优化
  opt <- optim(par = init_params, fn = loss_function, method = "L-BFGS-B",
               lower = lower_bounds, upper = upper_bounds)
  
  # 预测
  best_params <- opt$par
  out <- ode(y = init_state, times = time_points, func = sir_model, parms = best_params)
  return(as.data.frame(out)$I)
}

# 基于某站点数据训练一个 SIR-RF 模型：输入(t, I)，输出 dI/dt / I
train_rf_sir_model <- function(data) {
  data <- data %>% arrange(X1, X2)
  I <- data$Y
  t <- 1:length(I)
  dI <- estimate_derivative(I)
  valid_idx <- which(!is.na(dI) & I > 0)
  df <- data.frame(t = t[valid_idx], I = I[valid_idx], dI_div_I = dI[valid_idx] / I[valid_idx])
  model <- rfsrc(dI_div_I ~ t + I, data = df, ntree = 100)
  return(model)
}

# 使用某模型预测在目标数据上的 dI/dt
predict_rf_sir <- function(model, I_vec, t_vec) {
  newdata <- data.frame(t = t_vec, I = I_vec)
  dI_div_I_pred <- predict(model, newdata)$predicted
  return(dI_div_I_pred * I_vec)
}

# 训练最终 FTRF 模型：输入(t, I, hatI'_2, ..., hatI'_k)，输出 I'
train_ftrf_model <- function(target_data, aux_derivs) {
  I <- target_data$Y
  t <- 1:length(I)
  dI <- estimate_derivative(I)
  valid_idx <- which(!is.na(dI))
  df <- data.frame(t = t[valid_idx], I = I[valid_idx], dI = dI[valid_idx])
  for (i in 1:length(aux_derivs)) {
    df[[paste0("hat_dI_", i)]] <- aux_derivs[[i]][valid_idx]
  }
  formula_str <- paste("dI ~ t + I +", paste0("hat_dI_", 1:length(aux_derivs), collapse = " + "))
  model <- rfsrc(as.formula(formula_str), data = df, ntree = 100)
  return(model)
}

# 使用 FTRF 预测未来的 I 序列（Euler）
predict_future_I <- function(ftrf_model, target_data, aux_derivs, steps = 10) {
  I_pred <- numeric(steps)
  I_pred[1] <- target_data$Y[nrow(target_data) - steps + 1]
  t0 <- nrow(target_data) - steps + 1
  for (i in 2:steps) {
    t_now <- t0 + i - 1
    row <- data.frame(t = t_now, I = I_pred[i - 1])
    for (j in 1:length(aux_derivs)) {
      row[[paste0("hat_dI_", j)]] <- aux_derivs[[j]][t_now]
    }
    dI_now <- predict(ftrf_model, newdata = row)$predicted
    I_pred[i] <- I_pred[i - 1] + dI_now
  }
  return(I_pred)
}

# 使用目标 & 辅助数据
target_data <- tar_list[[1]]  # Alabama
aux_data <- tar_list[[2]]     # Alaska

# 排序
target_data <- target_data[order(target_data$X1, target_data$X2), ]
aux_data <- aux_data[order(aux_data$X1, aux_data$X2), ]

# 滑动窗口参数
window_size <- 4
step <- 4
n <- nrow(target_data)

predicted_all <- data.frame(time = numeric(), I = numeric(), obs = numeric())

for (start in seq(1, n - window_size * 2, by = step)) {
  train_idx <- start:(start + window_size - 1)
  test_idx <- (start + window_size):(start + 2 * window_size - 1)
  if (max(test_idx) > n) break
  
  # 当前窗口的目标站点数据
  train_data <- target_data[train_idx, ]
  test_data <- target_data[test_idx, ]
  
  # 当前窗口的辅助站点数据（使用同一时间窗口）
  aux_window <- aux_data[train_idx, ]
  
  # 训练辅助站点 SIR 模型
  hatI_train <- fit_sir_predict(aux_window, seq_len(nrow(train_data)))
  hatI_test <- fit_sir_predict(aux_window, seq_len(nrow(test_data)))
  
  # 构造训练集与测试集（加辅助特征）
  ftrf_train <- train_data
  ftrf_train$hatI <- hatI_train
  
  ftrf_test <- test_data
  ftrf_test$hatI <- hatI_test
  
  # 模型训练 + 预测
  rf_model <- rfsrc(Y ~ ., data = ftrf_train, ntree = 100)
  pred_I <- predict(rf_model, newdata = ftrf_test)$predicted
  
  predicted_all <- rbind(predicted_all, data.frame(time = test_idx, I = pred_I, obs = test_data$Y))
}

# # 绘图
# cpt <- ggplot(predicted_all, aes(x = time)) +
#   geom_point(aes(y = obs, color = "实际"), size = 2) +
#   geom_line(aes(y = I, color = "预测"), size = 1.2) +
#   scale_color_manual(values = c("预测" = "blue", "实际" = "red")) +
#   labs(
#     title = "Alabama SIR-FTRF预测结果",
#     x = "时间/周", y = "感染数量", color = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = c(0.95, 0.95),
#     legend.justification = c("right", "top"),
#     legend.background = element_rect(fill = "white", color = "black")
#   )
# 

# Plotting
cpt <- ggplot(predicted_all, aes(x = time)) +
  geom_point(aes(y = obs, color = "Observed"), size = 2) +
  geom_line(aes(y = I, color = "Predicted"), size = 1.2) +
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +
  labs(
    title = "Alabama SIR-FTRF Prediction Results",
    x = "Time/Week", y = "Number of Infections", color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )


print(cpt)
# 获取当前脚本的目录，并设置为工作目录
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
# 保存为白底 PNG 文件
ggsave("cpt_ftrf_sir.png", plot = cpt,
       width = 5, height = 5, dpi = 300, bg = "white")
# 输出误差
mse <- mean((predicted_all$obs - predicted_all$I)^2)
r2 <- 1 - sum((predicted_all$I - predicted_all$obs)^2) / sum((predicted_all$obs - mean(predicted_all$obs))^2)
cat("MSE:", mse, "\n")
cat("R^2:", r2, "\n")
