
GMA  <- function(L, v, yhat, y, prior_w, lamda){
  # @para L       : iterations
  # @para v       : regularization parameter related to learning rate
  # @para yhat    : predicted value of response
  # @para y       : real response
  # @para prior_w : prior weight 
  # @para lamda   : regularization parameter of each model
  f_hat <- 0  
  Y_hat <- as.matrix(yhat) 
  Y     <- y
  err   <- c()
  final_w <- rep(0, ncol(Y_hat))
  for (l in 1:L) {
    initial_w <- rep(0, ncol(Y_hat))
    alpha     <- (l-1)/l;
    if(l == 1) mu <- 0 else if(l == 2) mu <-0.05 else mu <- v*(l-1)/l^2;
    if(l == 1) c  <- 1 else if(l == 2) c  <-0.25 else c <- 1/(20*v*(1-v)*(l-1));
    Q <- lapply(1:ncol(Y_hat), function(j)
    {
      sum((alpha * f_hat + (1 - alpha) * Y_hat[, j] - Y)^2) + 
        mu * sum(( f_hat - Y_hat[, j])^2) - lamda * c* log(prior_w[j])
    })
    id            <- which.min(unlist(Q))
    initial_w[id] <- 1 
    final_w       <- alpha * final_w + (1 - alpha) * initial_w
    f_hat         <- alpha * f_hat + (1 - alpha) * Y_hat[, id]
    err[l]        <- mean((c(f_hat) - Y)^2)
  }
  return(final_w)
}


gma_for_site <- function(new_trainset, new_testset, nsite)
{
  weight <- GMA(L = 100, v = 0.01, yhat = new_trainset[, paste0("Y", 1:nsite)], 
                y = new_trainset$Y, prior_w = 1, lamda = 0)
  pred  <- as.matrix(new_testset[, paste0("Y", 1:nsite)]) %*% weight
  error <- mean((pred - new_testset$Y)^2)
  return(error)
}

ftrf_for_site <- function(new_trainset, new_testset, site){
  y_position <- which(colnames(new_trainset) == 'Y')
  fml <- as.formula(paste("Y ~ ", paste0(colnames(new_trainset)[-c(y_position, y_position+site)], 
                                         collapse = "+")))
  model <- rfsrc(fml, new_trainset, ntree = 500, bootstrap = "by.root", samptype = "swr", block.size = 1)
  test_preds <- predict(model, new_testset)$predicted
  mse_for_site <- mean((new_testset$Y - test_preds)^2)
  return(mse_for_site)
}

sr_for_site <- function(new_trainset, new_testset, nsite)
{
  cv_sr   <- cv.glmnet(as.matrix(new_trainset[, paste0("Y", 1:nsite)]), new_trainset$Y, alpha = 0, family = "gaussian")
  fit_sr  <- glmnet(as.matrix(new_trainset[,paste0("Y", 1:nsite)]), new_trainset$Y, alpha = 0, family = "gaussian")
  pred_sr <- predict(fit_sr, newx = as.matrix(new_testset[, paste0("Y", 1:nsite)]), s = cv_sr$lambda.min)
  mse_for_site <- mean((new_testset$Y - pred_sr)^2)
  return(mse_for_site)
}


simulate <- function(nreps, nsite, nvars, site){
    MSE <- matrix(0, nrow = nreps, ncol = 10)
    for(nrep in 1:nreps){
      set.seed(nrep+100)
      index    <- sample(1:2, nrow(tar_list[[site]]), replace = TRUE)
      trainset <- tar_list[[site]][index==1,]
      testset  <- tar_list[[site]][index==2,]
      train_list <- c(list(trainset), tar_list[-site])
      fml      <- as.formula(paste('Y ~ ',  paste0("X", 1:nvars, collapse = '+')))
      ModelList0 <- lapply(train_list, function(x){
        rfsrc(fml, as.data.frame(x), block.size = 1, samptype = "swr", ntree = 1000)
      })
      
      ModelList1 <- lapply(ModelList0, function(model){
        tar_data <- train_list[[1]]
        y_pred   <- predict(model, tar_data)$predicted
        tar_data$Y <- tar_data$Y - c(y_pred)
        refit   <- rfsrc(fml, tar_data, samptype = "swr", ntree = 500)
        return(refit)
        })
      preds_for_site <- lapply(1:nsite, function(k){
          tar_data <- train_list[[1]]
          pred_0 <- predict(ModelList0[[k]], tar_data)$predicted
          pred_1 <- predict(ModelList1[[k]], tar_data)$predicted
          y_pred <- c(pred_0) + c(pred_1)
          return(y_pred)
        })
      preds_for_site <- do.call(cbind, preds_for_site)
      new_trainset   <- cbind(train_list[[1]], preds_for_site)
      colnames(new_trainset) <- c(paste0("X", 1:nvars), "Y", paste0("Y", 1:nsite))

      test_preds_for_site <- lapply(1:nsite, function(k){
          pred_0 <- predict(ModelList0[[k]], testset)$predicted
          pred_1 <- predict(ModelList1[[k]], testset)$predicted
          y_pred <- c(pred_0) + c(pred_1)
          return(y_pred)
        })
      test_preds_for_site <- do.call(cbind, test_preds_for_site)
      new_testset           <- cbind(testset, test_preds_for_site)
      colnames(new_testset) <- c(paste0("X", 1:nvars), "Y", paste0("Y", 1:nsite))

      # type 2
      preds_tr <- lapply(ModelList0, function(model){predict(model, trainset)$predicted})
      preds_tr <- do.call(cbind, preds_tr)
      aug_tr   <- cbind(trainset, preds_tr)
      colnames(aug_tr) <- c(paste0("X", 1:nvars), "Y", paste0("Y", 1:nsite))

      preds_te <- lapply(ModelList0, function(model){predict(model, testset)$predicted})
      preds_te <- do.call(cbind, preds_te)
      aug_te   <- cbind(testset, preds_te)
      colnames(aug_te) <- c(paste0("X", 1:nvars), "Y", paste0("Y", 1:nsite))
      
      # target learning
      y_pred    <- predict(ModelList0[[1]], testset)$predicted
      local_mse <- mean((y_pred - testset$Y)^2) 
      
      # Global learning
      all_data <- do.call(rbind, train_list)
      gl_fit   <- rfsrc(fml, all_data, samptype = "swr", ntree = 500)
      y_pred   <- predict(gl_fit, testset)$predicted
      gl_mse   <- mean((y_pred - testset$Y)^2) 
      
      #Stacked regression
      sr_mse_1  <- sr_for_site(new_trainset, new_testset, nsite = nsite)
      # sr_mse_2  <- sr_for_site(aug_tr, aug_te,  nsite = nsite)
      
      #GMA
      gma_mse_1 <- gma_for_site(new_trainset, new_testset, nsite = nsite)
      # gma_mse_2 <- gma_for_site(aug_tr, aug_te,  nsite = nsite)

      #FTRF1
      ftrf_mse1 <- ftrf_for_site(new_trainset, new_testset, site = 1)
      
      #FTRF2
      ftrf_mse2 <- ftrf_for_site(aug_tr, aug_te, site = 1)


      source_data <- do.call(rbind, train_list[-1])
      tar_data <- train_list[[1]]
      ## OTF-RF
      tlrf_fit   <- rfsrc(fml, source_data, samptype = "swr", ntree = 500)
      y_pred   <- predict(tlrf_fit, tar_data)$predicted
      tar_data$Y <- tar_data$Y - c(y_pred)
      tlrf_refit   <- rfsrc(fml, tar_data, samptype = "swr", ntree = 500)
      tlrf_pred <- c(predict(tlrf_fit, testset)$predicted) + c(predict(tlrf_refit, testset)$predicted)
      tlrf_mse   <- mean((tlrf_pred - testset$Y)^2) 

      # target learning
      bws <- c(0.8 * nvars, nvars, 1.25 * nvars, 1.5 * nvars, 1.8 * nvars)
      ta.cv.err <- c()
      for(jj in 1:5){
        ta.cv.err[jj] <- cv_bw(data = train_list[[1]], n_folds = 5, bw = bws[jj])
      }
      id.ta <- which.min(ta.cv.err)
      ta.tr.x <- train_list[[1]][, paste0("X", 1:nvars)]
      ta.tr.y <- train_list[[1]][, "Y"]
      ta.krr.fit  <- krls(X = ta.tr.x, y = ta.tr.y, sigma = bws[id.ta], print.level = 0)
      ta.te.x <- testset[, paste0("X", 1:nvars)]
      ta.te.y <- testset[, "Y"]
      ta.te.py <- predict(ta.krr.fit, ta.te.x)$fit
      ta.krr.mse <- mean((ta.te.y - c(ta.te.py))^2)

      # Global learning
      # all_data <- do.call(rbind, train_list)
      # gl.cv.err <- c()
      # for(jj in 1:5){
      #   gl.cv.err[jj] <- cv_bw(data = all_data, n_folds = 5, bw = bws[jj])
      # }
      # id.gl <- which.min(gl.cv.err)
      # gl.tr.x <- all_data[, paste0("X", 1:nvars)]
      # gl.tr.y <- all_data[, "Y"]
      # gl.krr.fit  <- krls(X = gl.tr.x, y = gl.tr.y, sigma = bws[id.gl], print.level = 0)
      # gl.te.x <- testset[, paste0("X", 1:nvars)]
      # gl.te.y <- testset[, "Y"]
      # gl.te.py <- predict(gl.krr.fit, gl.te.x)$fit
      # gl.krr.mse <- mean((gl.te.y - c(gl.te.py))^2)
      gl.krr.mse <- 0

      ## OTF-KRR
      #
      bws <- c(0.8 * nvars, nvars, 1.25 * nvars, 1.5 * nvars, 1.8 * nvars)
      so.cv.err <- c()
      for(jj in 1:5){
        so.cv.err[jj] <- cv_bw(data = source_data, n_folds = 5, bw = bws[jj])
      }
      id.so <- which.min(so.cv.err)
      so.tr.x <- source_data[, paste0("X", 1:nvars)]
      so.tr.y <- source_data[, "Y"]
      tlkrr_fit   <- krls(X = so.tr.x, y = so.tr.y, sigma = bws[id.so], print.level = 0)
      #
      tar_data <- train_list[[1]]
      ta.tr.x <- tar_data[, paste0("X", 1:nvars)]
      ta.tr.y <- tar_data$Y
      ta.te.x <- testset[, paste0("X", 1:nvars)]
      ta.te.y <- testset$Y
      ta.tr.py   <- predict(tlkrr_fit, ta.tr.x)$fit
      off_data <- tar_data
      off_data$Y <- ta.tr.y - c(ta.tr.py)
      ta.cv.err <- c()
      for(jj in 1:5){
        ta.cv.err[jj] <- cv_bw(data = off_data, n_folds = 5, bw = bws[jj])
      }
      id.ta <- which.min(ta.cv.err)
      tlkrr_refit   <- krls(X = ta.tr.x, y = ta.tr.y - c(ta.tr.py), sigma = bws[id.ta], print.level = 0)
      tlkrr_pred <- c(predict(tlkrr_fit, ta.te.x)$fit + c(predict(tlkrr_refit, ta.te.x)$fit))
      tlkrr_mse   <- mean((tlkrr_pred - ta.te.y)^2) 


      MSE[nrep, ] <- c(local_mse, ftrf_mse1, ftrf_mse2, sr_mse_1, gma_mse_1, gl_mse, tlrf_mse, ta.krr.mse,  gl.krr.mse, tlkrr_mse)
      write.csv(MSE, paste0('mse_of_PHLs_', site, ".csv"), row.names = FALSE)
    }
}

cv_bw <- function(data, n_folds, bw){
    cv_ids <- caret::createFolds(1 : nrow(data), k = n_folds)
    mse <- lapply(1:length(cv_ids), function(k){
        dat_tr <- data[-cv_ids[[k]], ]
        dat_te <- data[cv_ids[[k]], ]
        tr.x <- dat_tr[, grep('X', colnames(data), value = TRUE)]
        tr.y <- dat_tr[, "Y"]
        te.x <- dat_te[, grep('X', colnames(data), value = TRUE)]
        te.y <- dat_te[, "Y"]
        krr_fit <- krls(X = tr.x, y = tr.y, sigma = bw, print.level = 0)
        pred <- predict(krr_fit, te.x)$fit
        return(mean((pred - te.y)^2))
    })
    return(mean(unlist(mse)))
}
