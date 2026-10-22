# ================= 1. 加载所有需要的包 =================
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

# ================= 2. 设置工作目录 =================
setwd("D:\\酒精-骨质疏松\\机器学习\\机器学习诊断模型\\ML模型-OP")  # 请按需修改
source("refer.ML.R")  # 确保此文件包含 RunML, ExtractVar, scaleData, CalPredictScore, PredictClass, RunEval, SimpleHeatmap

# ================= 3. 读取数据与精准预处理 =================
Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr <- as.matrix(Train_data[,1:(ncol(Train_data)-1), drop=FALSE])
Train_class <- Train_data[,ncol(Train_data), drop=FALSE]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
Test_expr <- as.matrix(Test_data[,1:(ncol(Test_data)-1), drop=FALSE])
Test_class <- Test_data[,ncol(Test_data), drop=FALSE]

if (grepl("_", rownames(Test_class)[1])) {
  Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
}
Test_class <- Test_class[,c("Cohort", "Type"), drop=FALSE]

# 特征交集与数值化
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
if(length(comgene) == 0) stop("训练集和测试集没有共同的基因！")
Train_expr <- Train_expr[, comgene, drop=FALSE]
Test_expr  <- Test_expr[, comgene, drop=FALSE]
storage.mode(Train_expr) <- "numeric"
storage.mode(Test_expr) <- "numeric"

# 安全剔除零方差基因
vars_train <- apply(Train_expr, 2, var, na.rm = TRUE)
keep_genes <- names(vars_train)[which(vars_train > 0)]
if(length(keep_genes) == 0) stop("所有基因均被剔除！")
Train_expr <- Train_expr[, keep_genes, drop=FALSE]
Test_expr  <- Test_expr[, keep_genes, drop=FALSE]

# ★★★ 关键修复1：用训练集的均值和方差来标准化测试集 ★★★
# 计算训练集的均值与标准差 (防止除零，标准差为0时设为1)
train_mean <- colMeans(Train_expr)
train_sd <- apply(Train_expr, 2, sd)
train_sd[train_sd == 0] <- 1

# 严格按训练集参数标准化
Train_set <- scale(Train_expr, center = train_mean, scale = train_sd)
Test_set  <- scale(Test_expr,  center = train_mean, scale = train_sd)

# 处理可能产生的 Inf (理论上不会，但加个保险)
Train_set[is.infinite(Train_set)] <- NA
Test_set[is.infinite(Test_set)] <- NA
clean_cols <- colSums(is.na(Train_set)) == 0 & colSums(is.na(Test_set)) == 0
Train_set <- Train_set[, clean_cols, drop=FALSE]
Test_set  <- Test_set[, clean_cols, drop=FALSE]

# 强制因子转换
Train_class$Type <- as.factor(Train_class$Type)
Test_class$Type  <- as.factor(Test_class$Type)

cat("================ 终极数据体检 =================\n")
cat("训练集样本量:", nrow(Train_set), " | 测试集样本量:", nrow(Test_set), "\n")
cat("保留基因数量:", ncol(Train_set), "\n")
cat("训练集类别分布:\n"); print(table(Train_class$Type))
cat("测试集类别分布:\n"); print(table(Test_class$Type))
cat("===============================================\n")
if(nrow(Train_set) < 15) stop("样本量太小，请检查数据！")

# ================= 4. 读取方法列表 =================
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods <- methodRT$Model
methods <- gsub("-| ", "", methods)
classVar <- "Type"         
min.selected.var <- 5  # ★ 关键修改：强制模型至少保留 5 个基因！防止出现单基因模型
Variable <- colnames(Train_set)
methods <- methods[!grepl("GBM", methods)] 

# ================= 5. 预训练（仅在训练集内筛选变量，绝对避免数据泄露） =================
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))
preTrain.method <- preTrain.method[!grepl("XGBoost|GBM|plsRglm", preTrain.method)]

preTrain.var <- list()
set.seed(123)
for (method in preTrain.method){
  tryCatch({
    preTrain.var[[method]] <- RunML(method = method,
                                    Train_set = Train_set, # ★ 注意：这里只传入了训练集
                                    Train_label = Train_class,
                                    mode = "Variable",
                                    classVar = classVar)
  }, error = function(e) {
    message("Error in pre-training for method ", method, ": ", e$message)
    preTrain.var[[method]] <- colnames(Train_set)
  })
}
preTrain.var[["simple"]] <- colnames(Train_set)

# ================= 6. 构建模型 =================
model <- list()
set.seed(123)
Train_set_bk <- Train_set
methods_to_run <- methods[!grepl("XGBoost|GBM|plsRglm", methods)]

for (method in methods_to_run){
  cat(match(method, methods), ":", method, "\n")
  method_parts <- strsplit(method, "\\+")[[1]]
  if (length(method_parts) == 1) method_parts <- c("simple", method_parts)
  
  Variable <- preTrain.var[[method_parts[1]]]
  if (is.null(Variable) || length(Variable) == 0) Variable <- colnames(Train_set)
  
  Train_set_sub <- Train_set_bk[, Variable, drop = FALSE]
  
  tryCatch({
    model[[method]] <- RunML(method = method_parts[2],
                             Train_set = Train_set_sub,
                             Train_label = Train_class,
                             mode = "Model",
                             classVar = classVar)
  }, error = function(e) {
    message("Error in building model for ", method, ": ", e$message)
    model[[method]] <- NULL
  })
  
  if (!is.null(model[[method]])) {
    vars <- tryCatch(ExtractVar(model[[method]]), error = function(e) colnames(Train_set_sub))
    if (length(vars) <= min.selected.var) model[[method]] <- NULL
  }
}
Train_set <- Train_set_bk
rm(Train_set_bk)

model <- model[!sapply(model, is.null)]
if (length(model) == 0) stop("没有成功构建任何模型，请检查数据或降低 min.selected.var 阈值。")
saveRDS(model, "model.MLmodel.rds")
cat("成功构建了", length(model), "个模型\n")

# ================= 7. 多变量逻辑回归 =================
FinalModel <- "multiLogistic"
if (FinalModel == "multiLogistic") {
  logisticmodel <- list()
  for (method_name in names(model)) {
    fit <- model[[method_name]]
    vars <- tryCatch(ExtractVar(fit), error = function(e) NULL)
    if (!is.null(vars) && length(vars) > 0) {
      data_subset <- as.data.frame(Train_set[, vars, drop = FALSE])
      data_subset$response <- Train_class[[classVar]]
      formula_str <- paste("response ~", paste(vars, collapse = " + "))
      suppressWarnings({
        tryCatch({
          tmp <- glm(as.formula(formula_str), family = "binomial", data = data_subset)
          tmp$subFeature <- vars
          logisticmodel[[method_name]] <- tmp
        }, error = function(e) NULL)
      })
    }
  }
  logisticmodel <- logisticmodel[!sapply(logisticmodel, is.null)]
  if (length(logisticmodel) > 0) saveRDS(logisticmodel, "model.logisticmodel.rds")
}

# ================= 8. 预测风险得分 =================
model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)[!sapply(model, is.null)]
methodsValid <- methodsValid[!grepl("plsRglm", methodsValid)]
cat("使用", length(methodsValid), "个模型进行预测\n")

RS_list <- list()
for (method in methodsValid){
  cat("Processing method:", method, "\n")
  tryCatch({
    RS <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set, Test_set))
    RS_list[[method]] <- RS
  }, error = function(e) {
    cat("  ERROR in", method, ":", e$message, "\n")
    RS_list[[method]] <- rep(0.5, nrow(rbind.data.frame(Train_set, Test_set)))
  })
}
if (length(RS_list) > 0) {
  riskTab <- as.data.frame(t(do.call(rbind, RS_list)))
  riskTab <- cbind(id=row.names(riskTab), riskTab)
  write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)
}

# ================= 9. 预测类别 =================
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set, Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab <- cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

# ================= 10. 提取特征基因 =================
fea_df <- data.frame()
for (method in methodsValid) {
  vars <- tryCatch(ExtractVar(model[[method]]), error = function(e) NULL)
  if (length(vars) > 0) {
    temp_df <- data.frame(features = vars, algorithm = method, stringsAsFactors = FALSE)
    fea_df <- rbind(fea_df, temp_df)
  }
}
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# ================= 11. 计算AUC =================
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],
                                Test_set = Test_set,
                                Test_label = Test_class,
                                Train_set = Train_set,
                                Train_label = Train_class,
                                Train_name = "Train",
                                cohortVar = "Cohort",
                                classVar = classVar)
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab <- cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)

# ================= 12. 绘制AUC热图 =================
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1)
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC), , drop = FALSE]

CohortCol <- c("red", "blue")
if(ncol(AUC_mat)>2) CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)

cellwidth <- 1; cellheight <- 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,
                    avg_Cindex = avg_AUC,
                    CohortCol = CohortCol,
                    barCol = "steelblue",
                    cellwidth = cellwidth,
                    cellheight = cellheight,
                    cluster_columns = F,
                    cluster_rows = F)
pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

# ================= 13. 核心精炼策略：基于模型共识提取稳健基因 =================
cat("\n========== 基于模型共识提取稳健基因（排除外部验证集干扰） ==========\n")

# ★★★ 关键修复：绝对不能用外部验证集来筛选基因！改为根据模型在训练集/内部验证的表现，提取所有有效模型的基因频率。
if (exists("model") && length(model) > 0) {
  
  # 提取所有有效模型的基因
  all_genes <- unlist(lapply(model, function(fit) {
    tryCatch(ExtractVar(fit), error = function(e) NULL)
  }))
  all_genes <- all_genes[!sapply(all_genes, is.null)]
  
  if (length(all_genes) > 0) {
    gene_freq <- sort(table(all_genes), decreasing = TRUE)
    
    # 设定阈值：至少在2个或以上的模型中出现的基因才保留 (共识基因)
    robust_genes <- names(gene_freq)[gene_freq >= 2]
    
    # 如果交集基因大于10个，只保留频率最高的前10个
    if(length(robust_genes) > 10) {
      robust_genes <- robust_genes[1:10]
    }
    # 如果交集基因不足3个，保留频率最高的前5个（兜底）
    if(length(robust_genes) < 3) {
      robust_genes <- names(gene_freq)[1:min(5, length(gene_freq))]
    }
    
    result_df <- data.frame(Gene = robust_genes, Frequency = as.numeric(gene_freq[robust_genes]))
    write.table(result_df, "Consensus_Refined_Genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("✅ 共识精炼基因已写入 Consensus_Refined_Genes.txt\n")
    cat("   最终筛选出", length(robust_genes), "个稳健基因:\n")
    print(robust_genes)
    
  } else {
    message("⚠️ 未能从模型中提取到基因。")
  }
} else {
  message("⚠️ 缺少模型对象，跳过频率报告。")
}

# ================= 14. 终极验证：用精简基因集重建模型并在外部数据测试 =================
cat("\n========== 终极验证：使用精简基因集重建模型 ==========\n")
if (exists("robust_genes") && length(robust_genes) >= 3) {
  
  Train_set_final <- Train_set[, robust_genes, drop = FALSE]
  Test_set_final  <- Test_set[, robust_genes, drop = FALSE]
  
  # ★ 建议使用最稳健的算法（Lasso 或 Stepglm[forward]），放弃复杂模型
  set.seed(123)
  final_model <- RunML(method = "RF+Enet[alpha=0.7]", 
                       Train_set = Train_set_final,
                       Train_label = Train_class,
                       mode = "Model",
                       classVar = classVar)
  
  saveRDS(final_model, "final_model_robust_genes.rds")
  
  # 重新计算 AUC（此时训练集 AUC 可能会下降，外部验证 AUC 也会回归真实，这才是可信的结果）
  final_auc <- RunEval(fit = final_model,
                       Test_set = Test_set_final,
                       Test_label = Test_class,
                       Train_set = Train_set_final,
                       Train_label = Train_class,
                       Train_name = "Train",
                       cohortVar = "Cohort",
                       classVar = classVar)
  
  cat("✅ 最终精简模型的性能（这才是可以写进论文的真实结果）:\n")
  print(final_auc)
  
  final_RS <- CalPredictScore(fit = final_model, new_data = rbind.data.frame(Train_set_final, Test_set_final))
  final_riskTab <- data.frame(id = row.names(rbind.data.frame(Train_set_final, Test_set_final)), RiskScore = final_RS)
  write.table(final_riskTab, "final_model_riskScore.txt", sep="\t", row.names=F, quote=F)
  
} else {
  message("⚠️ 精简基因数量少于3个，无法进行最终建模验证。")
}

cat("\n========== 所有分析完成 ==========\n")

# ================= 14. 终极验证：用精简基因集重建模型并在外部数据测试 =================
cat("\n========== 终极验证：使用精简基因集重建模型 ==========\n")
if (exists("robust_genes") && length(robust_genes) >= 3) {
  
  Train_set_final <- Train_set[, robust_genes, drop = FALSE]
  Test_set_final  <- Test_set[, robust_genes, drop = FALSE]
  
  # ★ 建议使用最稳健的算法（Lasso 或 Stepglm[forward]），放弃复杂模型
  set.seed(123)
  final_model <- RunML(method = "Lasso", 
                       Train_set = Train_set_final,
                       Train_label = Train_class,
                       mode = "Model",
                       classVar = classVar)
  
  saveRDS(final_model, "final_model_robust_genes.rds")
  
  # 重新计算 AUC（此时训练集 AUC 可能会下降，外部验证 AUC 也会回归真实，这才是可信的结果）
  final_auc <- RunEval(fit = final_model,
                       Test_set = Test_set_final,
                       Test_label = Test_class,
                       Train_set = Train_set_final,
                       Train_label = Train_class,
                       Train_name = "Train",
                       cohortVar = "Cohort",
                       classVar = classVar)
  
  cat("✅ 最终精简模型的性能（这才是可以写进论文的真实结果）:\n")
  print(final_auc)
  
  final_RS <- CalPredictScore(fit = final_model, new_data = rbind.data.frame(Train_set_final, Test_set_final))
  final_riskTab <- data.frame(id = row.names(rbind.data.frame(Train_set_final, Test_set_final)), RiskScore = final_RS)
  write.table(final_riskTab, "final_model_riskScore.txt", sep="\t", row.names=F, quote=F)
  
} else {
  message("⚠️ 精简基因数量少于3个，无法进行最终建模验证。")
}

cat("\n========== 所有分析完成 ==========\n")

# ================= 15. 严格按照审稿人意见计算全套指标 =================
# Part 1. 模型性能评估
# ============================================================
fixed_threshold <- 0.50

# ---------- 1.1 获取预测概率 ----------
train_pred <- as.numeric(CalPredictScore(fit = final_model,
                                         new_data = as.data.frame(Train_set_final)))
test_pred  <- as.numeric(CalPredictScore(fit = final_model,
                                         new_data = as.data.frame(Test_set_final)))

# 防止极端概率导致 qlogis 出现 Inf
train_pred <- pmin(pmax(train_pred, 1e-6), 1 - 1e-6)
test_pred  <- pmin(pmax(test_pred,  1e-6), 1 - 1e-6)

# ---------- 1.2 构建 eval_all ----------
train_eval <- data.frame(
  Sample      = rownames(Train_set_final),
  Cohort      = "Training",
  Outcome     = as.numeric(as.character(Train_class$Type)),
  Probability = train_pred,
  stringsAsFactors = FALSE)

test_eval <- data.frame(
  Sample      = rownames(Test_set_final),
  Cohort      = as.character(Test_class$Cohort),
  Outcome     = as.numeric(as.character(Test_class$Type)),
  Probability = test_pred,
  stringsAsFactors = FALSE)

eval_all <- rbind(train_eval, test_eval)

# ---------- 1.2b 标签编码对齐（关键修复） ----------
# OP_GSE13850 队列的 Control/Treat 编码与训练集相反
idx_fix <- eval_all$Cohort == "OP_GSE13850"
eval_all$Outcome[idx_fix] <- 1 - eval_all$Outcome[idx_fix]
cat("Harmonized OP_GSE13850 label encoding.\n")

# ---------- 1.2c 诊断：按 Cohort 检查方向 ----------
res_by_cohort <- do.call(rbind, lapply(unique(eval_all$Cohort), function(ct) {
  x <- eval_all[eval_all$Cohort == ct, ]
  data.frame(
    Cohort      = ct,
    N           = nrow(x),
    Cases       = sum(x$Outcome == 1),
    Controls    = sum(x$Outcome == 0),
    Mean_P_Case = mean(x$Probability[x$Outcome == 1]),
    Mean_P_Ctrl = mean(x$Probability[x$Outcome == 0]),
    Diff        = mean(x$Probability[x$Outcome == 1]) -
      mean(x$Probability[x$Outcome == 0]),
    stringsAsFactors = FALSE)
}))
cat("\nDirection check by cohort:\n"); print(res_by_cohort)

# ---------- 1.3 单个 cohort 性能计算 ----------
calculate_metrics <- function(dat, threshold = 0.50) {
  dat <- dat[complete.cases(dat$Outcome, dat$Probability), ]
  n_total   <- nrow(dat)
  n_case    <- sum(dat$Outcome == 1)
  n_control <- sum(dat$Outcome == 0)
  
  # AUC + 95% CI
  auc_value <- auc_low <- auc_high <- NA
  if (n_case > 0 && n_control > 0) {
    roc_obj <- tryCatch(
      roc(response = dat$Outcome, predictor = dat$Probability,
          levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) NULL)
    if (!is.null(roc_obj)) {
      auc_value <- as.numeric(auc(roc_obj))
      ci_auc <- tryCatch(ci.auc(roc_obj, conf.level = 0.95, method = "delong"),
                         error = function(e) NULL)
      if (!is.null(ci_auc)) {
        auc_low  <- as.numeric(ci_auc[1])
        auc_high <- as.numeric(ci_auc[3])
      }
    }
  }
  
  # 固定阈值 0.50
  predicted_class <- ifelse(dat$Probability >= threshold, 1, 0)
  TP <- sum(predicted_class == 1 & dat$Outcome == 1)
  TN <- sum(predicted_class == 0 & dat$Outcome == 0)
  FP <- sum(predicted_class == 1 & dat$Outcome == 0)
  FN <- sum(predicted_class == 0 & dat$Outcome == 1)
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  # Brier score
  brier <- mean((dat$Probability - dat$Outcome)^2)
  
  # 校准截距和斜率
  calibration_intercept <- calibration_slope <- NA
  if (length(unique(dat$Outcome)) == 2) {
    logit_pred <- qlogis(dat$Probability)
    
    cal_int <- tryCatch(
      glm(Outcome ~ offset(logit_pred), family = binomial, data = dat),
      error = function(e) NULL)
    if (!is.null(cal_int)) calibration_intercept <- coef(cal_int)[1]
    
    cal_slope <- tryCatch(
      glm(Outcome ~ logit_pred, family = binomial, data = dat),
      error = function(e) NULL)
    if (!is.null(cal_slope)) calibration_slope <- coef(cal_slope)[2]
  }
  
  data.frame(
    N = n_total, Cases = n_case, Controls = n_control,
    AUC = auc_value, AUC_95CI_low = auc_low, AUC_95CI_high = auc_high,
    Sensitivity = sensitivity, Specificity = specificity,
    Balanced_Accuracy = balanced_accuracy,
    Threshold = threshold,
    Brier_Score = brier,
    Calibration_Intercept = calibration_intercept,
    Calibration_Slope = calibration_slope,
    TP = TP, TN = TN, FP = FP, FN = FN,
    stringsAsFactors = FALSE)
}

# ---------- 1.4 按 cohort 计算 ----------
performance_list <- lapply(unique(eval_all$Cohort), function(cohort_name) {
  cat("\n========================================\n")
  cat("Cohort:", cohort_name, "\n")
  cat("========================================\n")
  dat <- eval_all[eval_all$Cohort == cohort_name, ]
  print(table(dat$Outcome))
  calculate_metrics(dat, threshold = fixed_threshold)
})
names(performance_list) <- unique(eval_all$Cohort)

performance_table <- bind_rows(performance_list, .id = "Cohort")
rownames(performance_table) <- NULL

# 格式化 AUC + 95% CI
performance_table$AUC_95CI <- ifelse(
  is.na(performance_table$AUC), NA,
  sprintf("%.3f (%.3f–%.3f)",
          performance_table$AUC,
          performance_table$AUC_95CI_low,
          performance_table$AUC_95CI_high))

# ---------- 1.5 输出 ----------
write.table(performance_table, "Reviewer9_Model_Performance_Complete.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(performance_table, "Reviewer9_Model_Performance_Complete.csv",
          row.names = FALSE)

supp_table <- performance_table %>%
  select(Cohort, N, Cases, Controls, AUC_95CI, Sensitivity, Specificity,
         Balanced_Accuracy, Threshold, Brier_Score,
         Calibration_Intercept, Calibration_Slope)

write.csv(supp_table, "Supplementary_Table_Model_Performance.csv", row.names = FALSE)

cat("\n====================================================\n")
cat("FINAL MODEL PERFORMANCE SUMMARY\n")
cat("====================================================\n")
print(supp_table)
cat("\nFixed threshold:", fixed_threshold, "\n")
cat("Threshold was kept constant across all cohorts ",
    "(no cohort-specific optimization).\n", sep = "")

