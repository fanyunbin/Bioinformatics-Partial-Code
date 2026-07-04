# all ML algorithms
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  # for example: Enet [alpha=0.4]
  method = gsub(" ", "", method)
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  
  message("Run ", method_name, " algorithm for ", mode, "; ",
          ifelse(is.null(method_param), "", paste(method_param, collapse = ", ")), ";",
          " using ", ncol(Train_set), " Variables")
  
  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  
  if(!is.null(method_param)) {
    args = c(args, method_param)
  }
  
  obj <- tryCatch({
    do.call(what = paste0("Run", method_name),
            args = args)
  }, error = function(e) {
    message("Error in ", method_name, ": ", e$message)
    return(NULL)
  })
  
  if(mode == "Variable"){
    if(is.null(obj) || length(obj) == 0) {
      message("0 Variables retained; Using all variables\n")
      return(colnames(Train_set))
    } else {
      message(length(obj), " Variables retained;\n")
      return(obj)
    }
  } else {
    message("\n")
    return(obj)
  }
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  # 添加收敛控制参数以减少警告
  fit <- suppressWarnings(
    step(glm(formula = Train_label[[classVar]] ~ .,
             family = "binomial", 
             data = as.data.frame(Train_set),
             control = list(maxit = 100)),  # 增加最大迭代次数
         direction = direction, trace = 0)
  )
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = train(eval(parse(text = paste(classVar, "~."))), 
              data = data, 
              method="lda",
              trControl = trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  set.seed(123)
  
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial())
  
  n <- nrow(data)
  fold_vec <- sample(rep(1:10, length.out = n))
  folds <- matrix(0, nrow = n, ncol = 10)
  for (i in 1:10) {
    folds[fold_vec == i, i] <- 1
  }
  
  cvm <- cvrisk(fit, papply = lapply, folds = folds)
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(), 
                  control = boost_control(mstop = max(mstop(cvm), 40)))
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  if (!requireNamespace("plsRglm", quietly = TRUE)) {
    install.packages("plsRglm")
    library(plsRglm)
  }
  
  tryCatch({
    cv.plsRglm.res = plsRglm::cv.plsRglm(formula = Train_label[[classVar]] ~ ., 
                                 data = as.data.frame(Train_set),
                                 nt=10, verbose = FALSE)
    fit <- plsRglm::plsRglm(Train_label[[classVar]], 
                 as.data.frame(Train_set), 
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }, error = function(e) {
    message("Error in plsRglm: ", e$message)
    return(NULL)
  })
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  
  tryCatch({
    fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
                 data = cbind(Train_set, Train_label[classVar]),
                 ntree = 500, nodesize = rf_nodesize,  # 减少树的数量以加快速度
                 importance = T,
                 proximity = T,
                 forest = T)
    fit$subFeature = colnames(Train_set)
    
    if (mode == "Variable") {
      imp <- fit$importance
      if (!is.null(imp) && length(imp) > 0) {
        # 修复：如果所有重要性都是NA或0，返回所有特征
        if(all(is.na(imp)) || all(imp == 0, na.rm = TRUE)) {
          return(colnames(Train_set))
        }
        
        # 移除NA值
        imp_clean <- imp[!is.na(imp)]
        if(length(imp_clean) == 0) {
          return(colnames(Train_set))
        }
        
        # 选择重要性大于0的变量
        selected_vars <- names(imp_clean[imp_clean > 0])
        
        # 如果选出的变量太少，返回所有特征
        if (length(selected_vars) < 5) {
          return(colnames(Train_set))
        }
        
        return(selected_vars)
      } else {
        return(colnames(Train_set))
      }
    }
    
    if (mode == "Model") return(fit)
  }, error = function(e) {
    message("Error in RF: ", e$message)
    return(NULL)
  })
}

ExtractVar <- function(fit){
  if(is.null(fit)) {
    return(character(0))
  }
  
  Feature <- tryCatch({
    quiet(switch(
      EXPR = class(fit)[1],
      "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
      "glm" = names(coef(fit)),
      "svm.formula" = fit$subFeature,
      "train" = fit$coefnames,
      "glmboost" = names(coef(fit)[abs(coef(fit))>0]),
      "plsRglmmodel" = {
        if (!is.null(fit$Coeffs) && !is.null(rownames(fit$Coeffs))) {
          rownames(fit$Coeffs)[which(fit$Coeffs[,1] != 0)]
        } else {
          fit$subFeature
        }
      },
      "rfsrc" = {
        # 对于RF，返回所有特征
        fit$subFeature
      },
      "gbm" = {
        if (exists("summary.gbm")) {
          imp <- summary.gbm(fit, plotit = FALSE)
          if (!is.null(imp) && !is.null(imp$rel.inf)) {
            rownames(imp)[imp$rel.inf > 0]
          } else {
            fit$subFeature
          }
        } else {
          fit$subFeature
        }
      },
      "xgb.Booster" = fit$subFeature,
      "naiveBayes" = fit$subFeature,
      fit$subFeature
    ))
  }, error = function(e) {
    message("Error in ExtractVar: ", e$message)
    fit$subFeature
  })
  
  if(is.null(Feature) || length(Feature) == 0) {
    Feature <- fit$subFeature
  }
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  
  return(Feature)
}

# ... (其他函数保持不变，但您需要确保所有函数都存在)