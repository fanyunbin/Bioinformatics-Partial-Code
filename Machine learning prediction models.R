# 安装 BiocManager（如果未安装）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 从 Bioconductor 安装 ComplexHeatmap
BiocManager::install("ComplexHeatmap")

# 尝试加载
library(ComplexHeatmap)

# 引用包
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# 设置工作目录
setwd("C:\\Users\\fan'yun'bin\\Desktop\\生信\\机器学习\\OP\\ML模型")
source("refer.ML.R")

# 读取训练组的数据文件
Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr <- Train_data[,1:(ncol(Train_data)-1), drop=FALSE]
Train_class <- Train_data[,ncol(Train_data), drop=FALSE]

# 读取测试组的数据文件
Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr <- Test_data[,1:(ncol(Test_data)-1), drop=FALSE]
Test_class <- Test_data[,ncol(Test_data), drop=FALSE]
Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class <- Test_class[,c("Cohort", "Type")]

# 获取训练组和测试组的交集基因
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])

# 标准化数据
Train_set <- scaleData(data=Train_expr, centerFlags=T, scaleFlags=T) 
Test_set <- scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

# 读取机器学习方法的文件
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods <- methodRT$Model
methods <- gsub("-| ", "", methods)

# 准备机器学习模型的参数
classVar <- "Type"         # 设置分类的变量名
min.selected.var <- 3      # 降低基因数目的阈值（从5降到3）
Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

###################### 根据训练组数据构建机器学习模型 ######################
# 根据模型组合第一种机器学习方法筛选变量
# 根据模型组合第一种机器学习方法筛选变量
preTrain.var <- list()       # 用于保存各算法筛选的变量
set.seed(seed = 123)         # 设置种子

# 排除XGBoost作为预训练方法（如果存在）
preTrain.method <- preTrain.method[!grepl("XGBoost", preTrain.method)]

for (method in preTrain.method){
  tryCatch({
    preTrain.var[[method]] = RunML(method = method,          # 机器学习方法
                                   Train_set = Train_set,         # 训练组的基因表达数据
                                   Train_label = Train_class,    # 训练组的分类数据
                                   mode = "Variable",              # 选择运行模式(筛选变量)
                                   classVar = classVar)
  }, error = function(e) {
    message("Error in pre-training for method ", method, ": ", e$message)
    preTrain.var[[method]] <- colnames(Train_set)  # 出错时使用所有特征
  })
}

preTrain.var[["simple"]] <- colnames(Train_set)

# 根据模型组合第二种机器学习方法构建模型
model <- list()             # 初始化模型结果列表
set.seed(seed = 123)       # 设置种子
Train_set_bk = Train_set

# 排除XGBoost方法
methods_to_run <- methods[!grepl("XGBoost", methods)]

for (method in methods_to_run){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method_parts <- strsplit(method, "\\+")[[1]]
  
  if (length(method_parts) == 1) {
    method_parts <- c("simple", method_parts)
  }
  
  # 获取变量
  Variable = preTrain.var[[method_parts[1]]]
  if(is.null(Variable) || length(Variable) == 0) {
    Variable <- colnames(Train_set)
  }
  
  Train_set = Train_set_bk[, Variable, drop = FALSE]
  Train_label = Train_class
  
  # 添加错误处理
  tryCatch({
    model[[method_name]] <- RunML(method = method_parts[2],       # 机器学习方法
                                  Train_set = Train_set,         # 训练组的表达数据
                                  Train_label = Train_label,    # 训练组的分类数据
                                  mode = "Model",                 # 选择运行模式(构建模型)
                                  classVar = classVar)
  }, error = function(e) {
    message("Error in building model for ", method_name, ": ", e$message)
    model[[method_name]] <- NULL
  })
  
  # 如果某种机器学习方法筛选出的变量小于阈值，则该方法结果为空
  if(!is.null(model[[method_name]])) {
    vars <- tryCatch({
      ExtractVar(model[[method_name]])
    }, error = function(e) {
      colnames(Train_set)
    })
    
    if(length(vars) <= min.selected.var) {
      model[[method_name]] <- NULL
    }
  }
}

Train_set = Train_set_bk; rm(Train_set_bk)

# 移除NULL值
model <- model[!sapply(model, is.null)]

# 保存所有机器学习模型的结果
if(length(model) > 0) {
  saveRDS(model, "model.MLmodel.rds")
  cat("Successfully built", length(model), "models\n")
} else {
  stop("No valid models were built. Please check your data and method list.")
}

# 构建多变量逻辑回归模型（使用修复的版本）
FinalModel <- c("panML", "multiLogistic")[2]
if (FinalModel == "multiLogistic"){
  # 使用上面提供的修复版本
  logisticmodel <- list()
  
  for (method_name in names(model)) {
    fit <- model[[method_name]]
    if (!is.null(fit)) {
      vars <- tryCatch({
        ExtractVar(fit)
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(vars) && length(vars) > 0) {
        data_subset <- as.data.frame(Train_set[, vars, drop = FALSE])
        data_subset$response <- Train_class[[classVar]]
        
        if (length(vars) > 0) {
          formula_str <- paste("response ~", paste(vars, collapse = " + "))
          
          # 使用标准逻辑回归，但抑制警告
          suppressWarnings({
            tryCatch({
              tmp <- glm(as.formula(formula_str),
                         family = "binomial", 
                         data = data_subset)
              tmp$subFeature <- vars
              logisticmodel[[method_name]] <- tmp
            }, error = function(e) {
              # 如果失败，跳过这个模型
              NULL
            })
          })
        }
      }
    }
  }
  
  # 移除空值
  logisticmodel <- logisticmodel[!sapply(logisticmodel, is.null)]
  
  # 保存最终以多变量逻辑回归模型
  if (length(logisticmodel) > 0) {
    saveRDS(logisticmodel, "model.logisticmodel.rds")
    cat("Successfully built", length(logisticmodel), "logistic models\n")
  } else {
    message("Warning: No logistic models were built.")
  }
}

# 根据基因表达量计算每个样本的分类得分
model <- readRDS("model.MLmodel.rds")               # 使用各个机器学习模型的线性组合函数计算得分
methodsValid <- names(model)[!sapply(model, is.null)]  # 根据特征基因数目提取有效的模型

# 根据基因表达量预测样本的风险得分
# 根据基因表达量预测样本的风险得分
RS_list <- list()
for (method in methodsValid){
  cat("Processing method:", method, "\n")
  tryCatch({
    RS <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set, Test_set))
    cat("  Predicted", length(RS), "values for", method, "\n")
    RS_list[[method]] <- RS
  }, error = function(e) {
    cat("  ERROR in", method, ":", e$message, "\n")
    # 为出错的方法分配默认值
    RS_list[[method]] <- rep(0.5, nrow(rbind.data.frame(Train_set, Test_set)))
  })
}

# 检查结果
if(length(RS_list) > 0) {
  riskTab <- as.data.frame(t(do.call(rbind, RS_list)))
  riskTab <- cbind(id=row.names(riskTab), riskTab)
  write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)
} else {
  warning("No risk scores were calculated!")
}
riskTab <- as.data.frame(t(do.call(rbind, RS_list)))
riskTab <- cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)


methodsValid <- names(model)[!sapply(model, is.null)]

# 排除所有包含 plsRglm 的模型
methodsValid <- methodsValid[!grepl("plsRglm", methodsValid)]

cat("Using", length(methodsValid), "models (excluding plsRglm)\n")
cat("Valid methods:", paste(methodsValid, collapse = ", "), "\n")

# 找出具体出问题的模型
for (method in methodsValid){
  cat("Testing method:", method, "\n")
  
  tryCatch({
    fit <- model[[method]]
    data <- rbind.data.frame(Train_set, Test_set)
    new_data <- data[, fit$subFeature, drop = FALSE]
    
    # 测试 PredictClass
    label <- PredictClass(fit = fit, new_data = data)
    
    if (length(label) == nrow(data)) {
      cat("  OK: predicted", length(label), "values\n")
    } else {
      cat("  WARNING: expected", nrow(data), "values but got", length(label), "values\n")
      
      # 检查模型类型
      cat("  Model class:", class(fit)[1], "\n")
      
      # 直接测试原始预测
      cat("  Testing raw prediction...\n")
      pred_result <- tryCatch({
        quiet(switch(
          EXPR = class(fit)[1],
          "lognet"      = predict(fit, type = 'class', as.matrix(new_data)),
          "glm"         = predict(fit, type = 'response', as.data.frame(new_data)),
          "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = FALSE),
          "train"       = predict(fit, new_data, type = "raw"),
          "glmboost"    = predict(fit, type = "class", as.data.frame(new_data)),
          "plsRglmmodel" = plsRglm::predict.plsRglm(fit, newdata = as.data.frame(new_data), type = "response"),
          "rfsrc"        = predict(fit, as.data.frame(new_data)),
          "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)),
          "xgb.Booster" = predict(fit, as.matrix(new_data)),
          "naiveBayes" = predict(object = fit, type = "class", newdata = as.data.frame(new_data))
        ))
      }, error = function(e) {
        cat("  Raw prediction error:", e$message, "\n")
        NULL
      })
      
      if (!is.null(pred_result)) {
        cat("  Raw prediction structure:\n")
        print(str(pred_result))
      }
    }
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
  
  cat("\n")
}

















# 根据基因表达量预测样品的分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set, Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab <- cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

# 提取每种机器学习方法筛选到的变量(模型基因)
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- data.frame()
for (method in methodsValid) {
  vars <- ExtractVar(model[[method]])
  if (length(vars) > 0) {
    temp_df <- data.frame(features = vars, algorithm = method, stringsAsFactors = FALSE)
    fea_df <- rbind(fea_df, temp_df)
  }
}

write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# 计算每个模型的AUC值
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],      # 机器学习模型
                                Test_set = Test_set,        # 测试组的表达数据
                                Test_label = Test_class,    # 测试组的分类数据
                                Train_set = Train_set,      # 训练组的表达数据
                                Train_label = Train_class,  # 训练组的分类数据
                                Train_name = "Train",        # 训练组的标签
                                cohortVar = "Cohort",        # GEO的id
                                classVar = classVar)         # 分类变量
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab <- cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)

############################## 绘制AUC热图 ##############################
# 准备图形的数据
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

# 根据AUC的均值对机器学习模型进行排序
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),, drop = FALSE]

# 获取最优模型(训练组+测试组的AUC均值最大)
optimalModel <- rownames(AUC_mat)[1]
optimalModelGenes <- fea_df[fea_df$algorithm == optimalModel, "features"]
write.table(optimalModelGenes, file="optimalModelGenes.txt", sep="\t", quote=F, row.names=F, col.names=F)

# 设置热图注释的颜色
CohortCol <- c("red", "blue")
if(ncol(AUC_mat)>2){
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
}
names(CohortCol) <- colnames(AUC_mat)

# 绘制图形
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,    # AUC值的矩阵
                    avg_Cindex = avg_AUC,        # AUC均值
                    CohortCol = CohortCol,       # 数据集的颜色
                    barCol = "steelblue",        # 右侧柱状图的颜色
                    cellwidth = cellwidth, cellheight = cellheight,    # 热图每个格子的宽度和高度
                    cluster_columns = F, cluster_rows = F)      # 是否对数据进行聚类

# 输出热图
pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()


# 检查是否安装了 ComplexHeatmap
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  # 尝试从 Bioconductor 安装
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}

# 加载 ComplexHeatmap
library(ComplexHeatmap)

# 检查函数是否存在
if (exists("columnAnnotation")) {
  cat("columnAnnotation function exists\n")
} else {
  cat("columnAnnotation function does not exist\n")
}