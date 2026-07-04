# 加载必要的包
library(pROC)
library(dplyr)

# 设置工作目录（请根据实际情况调整）
setwd("E:\\酒精-骨质疏松\\机器学习\\机器学习诊断模型\\ML模型-RI")

# 1. 读取风险得分矩阵
riskTab <- read.table("model.riskMatrix.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 2. 读取测试集样本信息（Test_class）
# 注意：如果 Test_class 已在环境中，则跳过这步
Test_data <- read.table("data.test.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
Test_class <- Test_data[, ncol(Test_data), drop = FALSE]
Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(Test_class))
Test_class <- data.frame(
  Sample = rownames(Test_class),
  Cohort = Test_class$Cohort,
  Type = Test_class$Type,
  stringsAsFactors = FALSE
)

# 3. 取交集（确保样本一致）
common_samples <- intersect(rownames(riskTab), Test_class$Sample)
riskTab_sub <- riskTab[common_samples, , drop = FALSE]
Test_class_sub <- Test_class[match(common_samples, Test_class$Sample), ]

# 4. 获取模型列表和 cohort 列表
models <- colnames(riskTab_sub)
cohorts <- unique(Test_class_sub$Cohort)

# 5. 检查类别名称（例如 "Control" 和 "Treat" 可能不同）
unique_types <- unique(Test_class_sub$Type)
cat("检测到的类别名称：", paste(unique_types, collapse = ", "), "\n")
# 如果类别不是 "Control" 和 "Treat"，请相应地调整下面的 levels 参数
# 这里我们自动将第一个类别设为参考（较低风险），第二个设为较高风险
# 但为了稳妥，最好手动确认。假设您的类别是 "Control" 和 "Treat"（或 "0" 和 "1"）

# 我们采用自动检测方式：假设风险得分越高越可能属于第二个类别（按字母顺序或先出现的）
# 更好的做法是让 roc() 自动推断 direction，避免手动指定
# 我们首先尝试自动 direction，如果失败则尝试反向。

results <- list()

for (model in models) {
  for (cohort in cohorts) {
    idx <- which(Test_class_sub$Cohort == cohort)
    if (length(idx) < 2) next  # 样本数不足2，跳过
    true <- Test_class_sub$Type[idx]
    if (length(unique(true)) < 2) {
      message(paste("Skipping cohort", cohort, "because it has only one class."))
      next
    }
    pred <- riskTab_sub[idx, model]
    
    # 尝试计算 AUC，如果 direction 不对，让 roc 自动选择
    roc_obj <- tryCatch({
      roc(true, pred, quiet = TRUE, levels = unique(true))  # 自动检测方向
    }, error = function(e) {
      # 如果自动仍失败，尝试指定 direction = "<"
      roc(true, pred, quiet = TRUE, levels = unique(true), direction = "<")
    })
    
    # 如果仍然失败，跳过
    if (inherits(roc_obj, "try-error") || is.null(roc_obj)) {
      message(paste("Error in cohort", cohort, "for model", model))
      next
    }
    
    auc_val <- auc(roc_obj)
    ci_auc <- tryCatch({
      ci.auc(roc_obj, method = "delong")
    }, error = function(e) {
      # DeLong 可能因样本量小失败，尝试 bootstrap
      ci.auc(roc_obj, method = "bootstrap")
    })
    
    results[[paste(model, cohort, sep = "|")]] <- data.frame(
      Model = model,
      Cohort = cohort,
      AUC = as.numeric(auc_val),
      CI_lower = ci_auc[1],
      CI_upper = ci_auc[3],
      stringsAsFactors = FALSE
    )
  }
}

# 合并结果
if (length(results) > 0) {
  AUC_CI_df <- do.call(rbind, results)
  rownames(AUC_CI_df) <- NULL
  # 查看结果
  print(head(AUC_CI_df))
  # 保存
  write.csv(AUC_CI_df, "AUC_with_CI_all_models.csv", row.names = FALSE)
  cat("结果已保存至 AUC_with_CI_all_models.csv\n")
} else {
  cat("没有成功计算出任何 AUC，请检查数据。\n")
}