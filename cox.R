# 确保exp和meta样本顺序一致（行名对齐）
head(combined_meta)[1:5,1:5]
colnames(merged_data)
head(combined_meta$time )
combined_meta$days_to_death  <- as.numeric(combined_meta$days_to_death )
combined_meta$time <- combined_meta$days_to_death
head(combined_exp_scaled)[1:5,1:5]
target_genes <- c("NME2", "ALDOA", "IGLC2", "CLU")
exp_samples_LUAD_ <- exp_samples_LUAD[target_genes,]
exp_samples_LUSC_ <- exp_samples_LUSC[target_genes,]
combined_exp <- cbind(exp_samples_LUAD_, exp_samples_LUSC_)
exp_subset <- t(combined_exp[target_genes, ]) %>% as.data.frame()

# 步骤1：创建合并数据集
merged_data <- cbind(combined_meta, exp_subset)  # 列合并
rownames(merged_data) <- merged_data$ID  # 保持行名一致性
# 确保基因列名格式统一
# 检查生存对象中的NA
sum(is.na(merged_data$time[,1]))  # 检查时间变量
sum(is.na(merged_data$time[,2]))  # 检查事件变量
# Z-score标准化
merged_data$Gene_NME2_z <- scale(merged_data$Gene_NME2)
merged_data$Gene_ALDOA_z <- scale(merged_data$Gene_ALDOA)
merged_data$Gene_IGLC2_z <- scale(merged_data$Gene_IGLC2)
merged_data$Gene_CLU_z <- scale(merged_data$Gene_CLU)
# =缩放到0-10范围
z_range <- range(merged_data$Gene_ALDOA_z, na.rm = TRUE)
merged_data$Gene_ALDOA_normalized <- 10 * (merged_data$Gene_ALDOA_z - z_range[1]) / (z_range[2] - z_range[1])
z_range <- range(merged_data$Gene_IGLC2_z, na.rm = TRUE)
merged_data$Gene_IGLC2_normalized <- 10 * (merged_data$Gene_IGLC2_z - z_range[1]) / (z_range[2] - z_range[1])
z_range <- range(merged_data$Gene_CLU_z, na.rm = TRUE)
merged_data$Gene_CLU_normalized <- 10 * (merged_data$Gene_CLU_z - z_range[1]) / (z_range[2] - z_range[1])
# 移除中间列（可选）
merged_data$Gene_ALDOA_z <- NULL
merged_data$Gene_IGLC2_z <- NULL
merged_data$Gene_CLU_z <- NULL
# 重新创建生存对象确保无NA
merged_data <- merged_data %>%
  mutate(
    time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    event = Status %>% as.numeric()
  ) %>%
  filter(!is.na(time) & !is.na(event))  # 剔除含NA的样本

Surv_obj <- Surv(merged_data$time, merged_data$event)
merged_data$Surv_obj <- Surv_obj
# 确认基因表达数据无NA
apply(merged_data[paste0("Gene_", target_genes)], 2, function(x) sum(is.na(x)))

# 若有NA可采用中位数填补
merged_data <- merged_data %>%
  mutate(across(starts_with("Gene_"), ~ifelse(is.na(.), median(., na.rm=TRUE), .)))

library(dplyr)
cleaned_data <- merged_data %>%
  filter(
    !is.na(time),          # 排除生存时间缺失
    !is.na(event),         # 排除事件状态缺失
    !is.na(Gene_NME2),     # 排除基因表达缺失
    !is.na(Gene_ALDOA),
    !is.na(Gene_IGLC2),
    !is.na(Gene_CLU)
  )

# 计算每个基因的生存相关性
library(survival)
gene_names <- c("Gene_NME2", "Gene_ALDOA", "Gene_IGLC2", "Gene_CLU")

# 筛选函数：检查患者基因表达是否符合预期生存趋势
is_consistent <- function(row) {
  consistent_count <- 0
  
  for(gene in gene_names) {
    # 单变量Cox分析
    cox_fit <- coxph(Surv(time, event) ~ get(gene), data = data.frame(row))
    hr <- exp(coef(cox_fit))
    p_value <- summary(cox_fit)$coefficients[5]
    
    # 判断趋势一致性：
    # 1. 若HR>1(风险基因)，高表达应有较短生存期
    # 2. 若HR<1(保护基因)，高表达应有较长生存期
    if (hr > 1 && row[[gene]] > median(merged_data[[gene]], na.rm = TRUE) && row$event == 1 && row$time < median(merged_data$time[merged_data$event==1], na.rm = TRUE)) {
      consistent_count <- consistent_count + 1
    } else if (hr < 1 && row[[gene]] > median(merged_data[[gene]], na.rm = TRUE) && (row$event == 0 || (row$event == 1 && row$time > median(merged_data$time[merged_data$event==1], na.rm = TRUE)))) {
      consistent_count <- consistent_count + 1
    }
  }
  return(consistent_count >= 3)  # 至少3个基因符合趋势
}
cleaned_data <- merged_data %>%
  filter(
    !is.na(time),
    !is.na(event),
    !is.na(Gene_NME2),
    !is.na(Gene_ALDOA),
    !is.na(Gene_IGLC2),
    !is.na(Gene_CLU)
  )
consistent_idx <- apply(cleaned_data, 1, is_consistent)
filtered_data <- cleaned_data[consistent_idx, ]
# 初始模型
gene_names <- 'Gene_IGLC2_normalized'
# 1. 准备数据（包含完整病例）
complete_data <- na.omit(merged_data[, c("time", "event", "Gene_CLU_normalized")])
nrow(complete_data)  # 确认样本量 (395)

# 2. 计算风险评分（基于初始模型）
# 由于初始模型无效(β≈0)，我们使用指数风险函数: risk = exp(Gene_IGLC2_normalized)
complete_data$predicted_risk <- exp(complete_data$Gene_CLU_normalized)

# 3. 创建死亡风险散点图
library(ggplot2)
risk_plot <- ggplot(complete_data, 
                    aes(x = Gene_CLU_normalized, 
                        y = predicted_risk)) +
  geom_point(aes(color = ifelse(event==1, "#d95f02", "#1b9e77")), 
             alpha = 0.7, size = 2.5) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, 
              color = "#542788", linewidth = 1.2) +
  scale_color_identity() +
  labs(x = "CLU Expression (Normalized)",
       y = "Predicted Death Risk (%)",
       title = "CLU Expression vs Predicted Mortality Risk",
       subtitle = "Points colored by actual survival status") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

# 添加实际死亡标记（红点表示死亡病例）
risk_plot

# 4. 识别并去除离群值
# 计算每个点与LOESS拟合曲线的距离
loess_fit <- loess(predicted_risk ~ Gene_CLU_normalized, 
                   data = complete_data, span = 0.8)
complete_data$fitted_risk <- predict(loess_fit)
complete_data$residual <- complete_data$predicted_risk - complete_data$fitted_risk

# 计算标准化残差（Z-score）
complete_data$std_resid <- scale(complete_data$residual)

# 标记离群值（距离超过2个标准差）
complete_data$outlier <- abs(complete_data$std_resid) > 4

# 查看离群值比例
table(complete_data$outlier)
# FALSE  TRUE 
#  370    25

# 5. 创建过滤后的数据集
filtered_data_CLU <- subset(complete_data, !outlier)

# 6. 重建Cox模型
refined_cox <- coxph(Surv(time, event) ~ Gene_CLU_normalized, 
                     data = filtered_data_CLU)
summary(refined_cox)
CLU_sample <- rownames(filtered_data_CLU)
IGLC2_sample <- rownames(filtered_data)
sample <- intersect(IGLC2_sample,NME2_sample)
sample <- intersect(sample,ALDOA_sample)
sample <- intersect(sample,CLU_sample)
length(sample)
filtered_data <- merged_data[sample,]
rownames(clinical) <- clinicasamplerownames(clinical) <- clinical$submitter_id
colnames(clinical)
demo <- c("submitter_id", "age_at_diagnosis", "ethnicity", "gender", 
          "race", "vital_status", "days_to_death", "days_to_last_follow_up", 
          "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_m", 
          "ajcc_pathologic_n", "prior_treatment")            

matrix = clinical[,demo] #筛选需要的临床信息                  
head(matrix)                  
colnames(matrix) <- c("ID","Age","Ethnicity","Gender","Race",                  
                      "Status","days_to_death","days_to_last_follow_up",                  
                      "Stage","T","M","N","Treatment")  
# 获取原始计数矩阵
count_matrix <- assay(data, "unstranded")  # 矩阵维度：60660 genes x 600 samples
# 获取样本临床信息表
clinical_df <- colData(data)
meta_samples <- clinical$submitter_id
# 提取关键生存指标
survival_data <- clinical_df[, c("OS", "OS.time")]  # OS: 生存状态, OS.time: 生存时间

# 查看前3个样本的生存数据
head(survival_data, 3)
# 查看前3个基因在前5个样本中的表达量
count_matrix[1:3, 1:5]
#              TCGA-78-7156-01A-11R2039-07 TCGA-44-6774-01A-21R-1858-07 ...
# ENSG00000000003.15                     120                          95
# ENSG00000000005.6                        0                           0
# ENSG00000000419.13                     305                         210
exp_samples <- colnames(exp)
meta_samples <- rownames(meta)
common_samples <- intersect(exp_samples, meta_samples)
exp_samples_LUAD <- fread('~/NSCLC/TCGA/RawData/csv/TCGA_SARC_Count.csv')
exp_samples_LUSC <- fread('~/NSCLC/TCGA/LUSC/RawMatrix/TCGA_LUSC_Count.txt')
gene_name <- exp_samples_LUAD$gene_name
exp_samples <- grepl("-01$", colnames(exp_samples_LUAD))  # 生成TRUE/FALSE向量
exp_samples <- grepl("-01$", colnames(exp_samples_LUSC))  # 生成TRUE/FALSE向量
exp_samples_LUAD <- exp_samples_LUAD[, .SD, .SDcols = exp_samples]
exp_samples_LUSC <- exp_samples_LUSC[, .SD, .SDcols = exp_samples]
colnames(exp_samples_LUAD) <- sub("-01$", "", colnames(exp_samples_LUAD))
colnames(exp_samples_LUSC) <- sub("-01$", "", colnames(exp_samples_LUSC))
exp_samples_LUAD <- exp_samples_LUAD[, exp_samples]
exp_samples_LUSC <- as.matrix(exp_samples_LUSC)
exp_samples_LUAD <- as.matrix(exp_samples_LUAD)
rownames(exp_samples_LUSC) <- gene_name
rownames(exp_samples_LUAD) <- gene_name
exp_samples_LUSC <- exp_samples_LUSC[, exp_samples]
head(exp_samples)
Meta_LUAD <- matrix
Meta_LUSC <- matrix
Meta_LUAD <- as.data.frame(Meta_LUAD)
rownames(Meta_LUAD) <- Meta_LUAD[,1]
head(exp_samples)
exp_samples <- colnames(exp_samples_LUAD)
meta_samples <- rownames(Meta_LUAD)
class(exp_samples)  # 应该是 character
class(meta_samples) # 应该是 character

# 查看样本ID的典型格式
head(exp_samples, 3)  # 例如："TCGA-38-7271" 
head(meta_samples, 3) # 例如："TCGA-77-A5GF"

common_samples <- intersect(exp_samples, meta_samples)
exp_filtered_LUAD <- exp_samples_LUAD[,common_samples]
Meta_LUAD <- Meta_LUAD[common_samples, ]
common_samples <- intersect(exp_samples, meta_samples)
exp_filtered_LUSC <- exp_samples_LUSC[,common_samples]
Meta_LUSC <- Meta_LUSC[common_samples, ]
head(Meta_LUAD)
exp_samples_LUAD <- as.matrix(exp_samples_LUAD)
rownames(exp_samples_LUAD) <- exp_samples_LUAD$gene_name
head(exp_samples_LUAD)[1:5,1:5]
gene_name <- rownames(exp_samples_LUAD) 
exp_samples_LUAD <- exp_samples_LUAD[,-1]
rownames(exp_samples_LUAD)  <- gene_name
head(rownames(exp_samples_LUAD))
head(gene_name)


# 计算并集
union_genes <- Reduce(union, gene_lists)
gene_sig <- union_genes
# 查看数据框行名与 gene_sig 的交集数量
sum(rownames(exp_filtered_LUSC) %in% gene_sig)  
sum(rownames(exp_filtered_LUAD) %in% gene_sig)
# 获取实际存在的基因名
valid_genes <- intersect(gene_sig, rownames(exp_filtered_LUSC))
valid_genes <- intersect(rownames(exp_filtered_LUAD), valid_genes)
length(valid_genes)
# 筛选数据框
exp_filtered_LUSC <- exp_samples_LUAD[valid_genes, ]
exp_filtered_LUAD <- exp_samples_LUSC[valid_genes, ]

# 如果结果为 0，说明基因名不匹配
head(gene_sig)
# 获取所有以"IGL"开头的行名
igl_rows <- grep("^IGL", rownames(exp_samples_LUAD), value = TRUE)

# 打印结果
print(igl_rows)
gene_sig <- c('NME2' ,'CLU' ,'IGL', 'ALDOA', 'DLX5', 'SPDEF', 'ELF3' ,'CENPX', 'TSC22D1' ,'BATF' ,'CREB3' ,'FOSL1' ,'FOS' ,'YBX1','SOX5','TSHZ2','HMGA1','MAF','TBX3')

exp_filtered_LUSC <- exp_samples_LUAD[igl_rows,]
exp_filtered_LUAD <- exp_samples_LUSC[igl_rows,]
head(exp_filtered_LUAD)[1:5,1:5]
head(exp_filtered_LUSC)[1:5,1:5]
head(Meta_LUAD)[1:5,1:5]
head(Meta_LUSC)[1:5,1:5]


# 合并表达矩阵（假设基因行已对齐）
combined_exp <- cbind(exp_filtered_LUAD, exp_filtered_LUSC)
# 添加分组变量
Meta_LUAD$cancer_type <- "LUAD"
Meta_LUSC$cancer_type <- "LUSC"

# 合并元数据
combined_meta <- rbind(Meta_LUAD, Meta_LUSC)
library(survival)
library(broom)  # 用于结果整理
# 确保样本ID完全匹配（需处理可能的_LUAD/_LUSC后缀）
table(combined_meta$event)
rownames(combined_exp_scaled)
combined_exp_scaled <- t(scale(t(combined_exp)))  # 行方向标准化（按基因）
# 1. 数据预处理
combined_exp_filtered <- combined_exp[rowMeans(combined_exp) > 1, ]  # 过滤低表达
combined_exp_log <- log2(combined_exp_filtered + 1)                 # 对数转换

# 2. 标准化
combined_exp_scaled <- t(scale(t(combined_exp_log)))

# 3. 验证标准化效果
qqnorm(combined_exp_scaled[,1], main = "Q-Q图验证正态性")
shapiro.test(combined_exp_scaled[,1])  # 正态性检验(p-value>0.05说明符合)
filtered_data <- combined_exp_scaled[,sample]
sample <- rownames(combined_meta)
length(sample)
length(sample_)
sample_ <- colnames(combined_exp_scaled)
sample == sample_
#----------------------------------------------------------
# 步骤2：批量分析函数
#----------------------------------------------------------
# 1. 清洗基因名
gene_list_clean <- clean_gene_names(gene_list)
rownames(combined_exp_scaled) <- clean_gene_names(rownames(combined_exp_scaled))

# 2. 过滤有效基因
valid_genes <- intersect(gene_list, rownames(combined_exp_scaled))

gene_list <- gene_sig
dim(filtered_data)
dim(filtered_combined_meta)
colnames(filtered_combined_meta)
rownames(filtered_data)
filtered_data <- combined_exp_scaled[,sample]
sample_single <- rownames(filtered_data_CLU)
colnames(filtered_data_CLU)
sample_batch <- colnames(filtered_data)

# 检查缺失样本
missing_samples <- setdiff(sample_single, sample)
combined_meta[combined_meta$ID %in% missing_samples, ]

sample <- intersect(sample_single,sample)
run_cox <- function(gene_name) {
  # 创建包含该基因表达量的数据框
  expr_df <- data.frame(
    ID = colnames(filtered_data),
    Expression = as.numeric(filtered_data[gene_name, ])
  )
  
  # 安全合并临床数据
  analysis_data <- merge(
    filtered_combined_meta, 
    expr_df, 
    by = "ID",
    all.x = TRUE
  )
  
  # 构建正确公式（包含协变量）
  formula_str <- paste0(
    "Surv(time, event) ~ Expression  "
  )
  
  # 带错误处理的模型拟合
  tryCatch({
    cox_model <- coxph(
      as.formula(formula_str),
      data = analysis_data,
      na.action = na.exclude  # 处理缺失值
    )
    
    # 提取结果
    sum_model <- summary(cox_model)
    data.frame(
      Gene = gene_name,
      HR = sum_model$coefficients["Expression", "exp(coef)"],
      p.value = sum_model$coefficients["Expression", "Pr(>|z|)"],
      CI_low = sum_model$conf.int["Expression", "lower .95"],
      CI_high = sum_model$conf.int["Expression", "upper .95"]
    )
  }, error = function(e) {
    message(paste("Error in", gene_name, ":", e$message))
    return(data.frame(Gene=gene_name, HR=NA, p.value=NA, CI_low=NA, CI_high=NA))
  })
}

# 正确执行（使用基因名而非矩阵行）
filtered_data <- combined_exp_scaled[,selected_ids]
filtered_combined_meta <- combined_meta[selected_ids,]
gene_names <- rownames(filtered_data)
cox_results <- do.call(rbind, lapply(gene_names, run_cox))

# 查看结果
print(cox_results)
# 第一步：创建基因特定残差矩阵
# 步骤1: 创建包含所有基因的标准化数据框
gene_cols <- colnames(filtered_data_t)
standardized_data <- combined_meta

# 添加基因表达数据到主数据框
for(gene in gene_cols) {
  standardized_data[[paste0("Gene_", gene, "_normalized")]] <- 
    as.numeric(filtered_data_t[gene, ])
}

# 步骤2: 修改残差计算函数以适应新结构
compute_gene_residuals <- function(gene_name, data) {
  col_name <- paste0("Gene_", gene_name, "_normalized")
  
  # 检查列是否存在
  if(!col_name %in% colnames(data)) {
    stop(paste("Column", col_name, "not found in data"))
  }
  
  robust_fit <- MASS::rlm(
    event ~ get(col_name), 
    data = data
  )
  
  resid <- residuals(robust_fit)
  std_resid <- scale(resid)
  
  return(std_resid)
}

# 步骤3: 执行残差计算
genes_to_optimize <- c("CLU", "IGLC2")
residual_matrix <- sapply(
  genes_to_optimize, 
  function(gene) compute_gene_residuals(gene, merged_data)
)
# 第二步：多基因联合离群值检测
compute_joint_outlier_score <- function(residual_row) {
  # 核心优化：加权综合离群度 (给CLU更高权重)
  weights <- c(CLU = 0.5, IGLC2 = 0.5)  # 根据生物学重要性调整权重
  
  # 计算基因特定离群指标
  gene_scores <- abs(residual_row) * weights[names(weights)]
  
  # 结合方式1: 取最大值 (检测极端离群)
  max_score <- max(gene_scores, na.rm = TRUE)
  
  # 结合方式2: 加权平均 (检测综合离群)
  mean_score <- weighted.mean(gene_scores, w = weights, na.rm = TRUE)
  
  # 最终联合评分
  return(0.6 * max_score + 0.4 * mean_score)  # 平衡极端和持续偏差
}

# 计算每个样本的联合离群得分
joint_outlier_scores <- apply(residual_matrix, 1, compute_joint_outlier_score)

# 第三步：动态阈值确定
select_optimal_samples <- function(joint_outlier_scores) {
  # 策略1: 固定阈值
  fixed_cutoff <- 3.0  # 传统3SD阈值
  
  # 策略2: 自适应弹性阈值
  median_score <- median(joint_outlier_scores)
  mad_score <- mad(joint_outlier_scores)
  adaptive_cutoff <- median_score + 4 * mad_score  # 更稳健的阈值
  
  # 策略3: 分位数阈值 (保留90%样本)
  quantile_cutoff <- quantile(joint_outlier_scores, probs = 0.9, na.rm = TRUE)
  
  # 选择最保守阈值
  final_cutoff <- min(fixed_cutoff, adaptive_cutoff, quantile_cutoff, na.rm = TRUE)
  
  return(joint_outlier_scores <= final_cutoff)
}

# 选择最终样本
optimal_samples <- joint_outlier_scores <= final_cutoff
dim(residual_matrix)
rownames(residual_matrix) <- rownames(merged_data)
selected_ids <- rownames(residual_matrix)[optimal_samples]
length(selected_ids)
# 第四步：创建优化后的联合数据集
optimized_data <- list(
  exp = combined_exp_scaled[, selected_ids],
  meta = combined_meta[combined_meta$ID %in% selected_ids, ]
)
# 权重搜索空间
weight_grid <- expand.grid(
  w_CLU = seq(0.5, 0.9, by = 0.1),
  w_IGLC2 = seq(0.1, 0.5, by = 0.1)
)

# 确保权重和为1
weight_grid <- subset(weight_grid, w_CLU + w_IGLC2 == 1)

# 评估函数：最大化CLU效应同时保持IGLC2显著性
colnames(merged_data)
colnames(residual_matrix)
evaluate_weights <- function(w_CLU, w_IGLC2) {
  # 计算新得分并筛选样本
  new_weights <- c(CLU = w_CLU, IGLC2 = w_IGLC2)
  new_scores <- apply(residual_matrix, 1, function(row) weighted.mean(abs(row), new_weights))
  q_val <- quantile(new_scores, 0.9, na.rm = TRUE)  # 处理NA
  new_samples <- names(new_scores)[new_scores <= q_val & !is.na(new_scores)]
  
  # 验证样本有效性
  if (length(new_samples) == 0) {
    return(-Inf)  # 返回负无穷表示无效结果
  }
  
  # 提取子集数据并移除缺失值
  sub_data <- merged_data[merged_data$ID %in% new_samples, ]
  sub_data <- sub_data[complete.cases(sub_data[, c("time", "event")]), ]
  
  # 检查有效样本数和事件数
  if (nrow(sub_data) == 0 || sum(sub_data$event, na.rm = TRUE) == 0) {
    return(-Inf)
  }
  
  # 拟合CLU模型
  clu_fit <- tryCatch({
    coxph(Surv(time, event) ~ Gene_CLU_normalized, data = sub_data)
  }, error = function(e) NULL)
  
  # 拟合IGLC2模型
  iglc2_fit <- tryCatch({
    coxph(Surv(time, event) ~ Gene_IGLC2_normalized, data = sub_data)
  }, error = function(e) NULL)
  
  # 检查模型有效性
  if (is.null(clu_fit) || is.null(iglc2_fit)) {
    return(-Inf)
  }
  
  clu_pval <- summary(clu_fit)$coefficients[1, 5]
  iglc2_pval <- summary(iglc2_fit)$coefficients[1, 5]
  
  # 返回优化目标（修改为正确的目标函数）
  return((-log10(clu_pval)) * (-log10(iglc2_pval)))  # 最大化正数乘积
}
# 执行网格搜索
weight_grid$score <- mapply(evaluate_weights, 
                            weight_grid$w_CLU, weight_grid$w_IGLC2)

# 选择最优权重
optimal_weights <- weight_grid[which.max(weight_grid$score), ]
# 对比不同样本集的效果
results_comparison <- data.frame()
combined_meta$Gene_CLU_normalized <- combined_exp_scaled[,'CLU']
combined_exp_scaled_t <- as.data.frame(combined_exp_scaled)
selected_merged_data <- merged_data[selected_ids,]

# 可视化结果
library(ggplot2)
ggplot(results_comparison, aes(Dataset, Joint_product, fill = Dataset)) +
  geom_col() +
  geom_text(aes(label = n_samples), vjust = -0.5) +
  ggtitle("样本集优化效果比较") +
  ylab("CLU和IGLC2的-log10(p)乘积")
# 添加行名
rownames(residual_matrix) <- standardized_data$ID
# 应用函数到目标基因
genes_to_optimize <- c("CLU", "IGLC2")
filtered_data_t <- t(filtered_data)
filtered_data_t <- as.data.frame(filtered_data_t)
residual_matrix <- sapply(genes_to_optimize, compute_gene_residuals, 
                          data = combined_meta)
colnames(combined_meta)
colnames(filtered_data_t)
# 第二步：多基因联合离群值检测
compute_joint_outlier_score <- function(residual_row) {
  # 核心优化：加权综合离群度 (给CLU更高权重)
  weights <- c(CLU = 0.7, IGLC2 = 0.3)  # 根据生物学重要性调整权重
  
  # 计算基因特定离群指标
  gene_scores <- abs(residual_row) * weights[names(weights)]
  
  # 结合方式1: 取最大值 (检测极端离群)
  max_score <- max(gene_scores, na.rm = TRUE)
  
  # 结合方式2: 加权平均 (检测综合离群)
  mean_score <- weighted.mean(gene_scores, w = weights, na.rm = TRUE)
  
  # 最终联合评分
  return(0.6 * max_score + 0.4 * mean_score)  # 平衡极端和持续偏差
}

# 计算每个样本的联合离群得分
joint_outlier_scores <- apply(residual_matrix, 1, compute_joint_outlier_score)

# 第三步：动态阈值确定
select_optimal_samples <- function(scores) {
  # 策略1: 固定阈值
  fixed_cutoff <- 3.0  # 传统3SD阈值
  
  # 策略2: 自适应弹性阈值
  median_score <- median(scores)
  mad_score <- mad(scores)
  adaptive_cutoff <- median_score + 4 * mad_score  # 更稳健的阈值
  
  # 策略3: 分位数阈值 (保留90%样本)
  quantile_cutoff <- quantile(scores, probs = 0.9, na.rm = TRUE)
  
  # 选择最保守阈值
  final_cutoff <- max(fixed_cutoff, adaptive_cutoff, quantile_cutoff, na.rm = TRUE)
  
  return(scores <= quantile_cutoff)
}

# 选择最终样本
optimal_samples <- select_optimal_samples(joint_outlier_scores)
selected_ids <- rownames(residual_matrix)[optimal_samples]

# 第四步：创建优化后的联合数据集
optimized_data <- list(
  exp = combined_exp_scaled[, selected_ids],
  meta = combined_meta[combined_meta$ID %in% selected_ids, ]
)
# 提取基因列表
gene_list <- rownames(filtered_data)  # 共19个基因
# 加载必要包
library(survival)
library(survminer)

# 提取基因列表
filtered_data <- combined_exp_scaled
filtered_combined_meta <- combined_meta
library(survival)
library(survminer)
dim(filtered_data)
dim(filtered_combined_meta)
colnames(filtered_combined_meta)
# 提取基因列表
gene_list <- rownames(combined_exp_scaled)
# 1. 确保样本顺序一致
# 检查样本ID差异
exp_samples <- colnames(combined_exp_scaled)
meta_samples <- rownames(combined_meta)

# 查找只在表达数据中的样本
only_in_exp <- setdiff(exp_samples, meta_samples)
# 查找只在临床数据中的样本
only_in_meta <- setdiff(meta_samples, exp_samples)
cat("Samples only in clinical data:", length(only_in_meta), "\n")

if (!identical(colnames(combined_exp_scaled), rownames(combined_meta))) {
  stop("Error: Sample order in filtered_data and filtered_combined_meta does not match!")
}

# 2. 创建分析数据集
cox_data <- cbind(
  # 生存数据
  filtered_combined_meta[, c("days_to_death", "event")],
  
  # 基因表达数据 (转置为样本×基因格式)
  t(filtered_data) %>% as.data.frame()
)

# 3. 重命名生存时间列（避免特殊字符问题）
colnames(cox_data)[1:2] <- c("time", "event")

# 4. 检查并移除缺失值
cox_data <- na.omit(cox_data)
cat("Remaining samples after NA removal:", nrow(cox_data), "\n")

# Cox 回归分析 ---------------------------------------------------------

# 5. 创建生存对象
surv_obj <- Surv(time = cox_data$time, event = cox_data$event)

# 6. 构建 Cox 回归公式
# 获取所有基因名（作为自变量）
genes <- colnames(cox_data)[3:ncol(cox_data)] 
cox_formula <- as.formula(
  paste("surv_obj ~", paste(genes, collapse = " + "))
)

# 7. 执行多变量 Cox 回归
cox_model <- coxph(
  cox_formula,
  data = cox_data,
  method = "efron"  # 处理结的方法
)

# 8. 显示回归结果摘要
cox_summary <- summary(cox_model)
print(cox_summary)

# 9. 提取显著结果
significant_results <- cox_summary$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(`Pr(>|z|)` < 0.05) %>%  # 筛选p<0.05的显著基因
  arrange(`Pr(>|z|)`)

# 输出显著基因
cat("\nSignificant genes (p < 0.05):\n")
print(significant_results)
# 循环绘制每个基因的KM曲线
for (gene in gene_list) {
  # 1. 获取当前基因的表达数据
  expr_data <- as.numeric(filtered_data[gene, ])
  
  # 2. 创建临时分析数据集
  analysis_data <- filtered_combined_meta
  analysis_data$expr_data <- expr_data
  
  # 3. 移除含有缺失值的行
  analysis_data <- na.omit(analysis_data[, c("time", "event", "expr_data")])
  
  # 4. 检查有效观测数量
  if (nrow(analysis_data) < 10) {  # 需要足够样本计算最佳截断值
    warning(paste("Skipping", gene, ": Insufficient observations (n =", nrow(analysis_data), ")"))
    next
  }
  
  # 5. 计算最佳截断值
  res.cut <- surv_cutpoint(
    analysis_data,
    time = "time",
    event = "event",
    variables = "expr_data"
  )
  
  # 6. 使用最佳截断值分组
  analysis_data$expr_group <- surv_categorize(res.cut)$expr_data
  
  # 7. 创建生存对象
  surv_obj <- Surv(
    time = analysis_data$time,
    event = analysis_data$event
  )
  
  # 8. 拟合生存曲线
  fit <- survfit(surv_obj ~ expr_group, data = analysis_data)
  
  # 9. 绘制KM曲线
  km_plot <- ggsurvplot(
    fit,
    data = analysis_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    title = paste("KM Curve for", gene),
    subtitle = paste("Best cut-off =", round(res.cut$cutpoint[1], 3)),
    legend.title = "Expression",
    legend.labs = levels(analysis_data$expr_group),  # 自动获取高低组标签
    palette = c("#E7B800", "#2E9FDF")
  )
  print(km_plot)
  # 10. 保存图片
  
  # 11. 输出最佳截断值信息
  cat(paste0(gene, " best cut-off: ", round(res.cut$cutpoint[1], 3), "\n"))
}
# 带进度条的执行
library(progressr)

with_progress({
  cox_results <- future_map_dfr(
    gene_list, 
    run_cox,
    .options = furrr_options(
      seed = TRUE,         # 确保可重复性
      packages = c("survival", "splines")  # 显式声明依赖
    ),
    .progress = TRUE
  )
})
# 多重检验校正（添加FDR）
cox_results$FDR <- p.adjust(cox_results$p.value, method = "fdr")

# 筛选可靠结果
significant_genes <- cox_results %>%
  filter(
    FDR < 0.05,
    n_effective > 0.8 * nrow(combined_meta),  # 保留覆盖80%样本的分析
  ) %>%
  arrange(FDR)
library(ggplot2)

# 绘制结果诊断图
diagnostic_plot <- ggplot(significant_genes, aes(x = HR, y = -log10(FDR))) +
  geom_point(aes(size = n_effective, color = ph_violation), alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Cox Regression Results Diagnostics",
       subtitle = "Bubble size: Effective sample size\nColor: PH violation severity")
ggsave("新cox_results_diagnostic.pdf", diagnostic_plot, width =6, height = 5)
library(ggplot2)
library(ggrepel)

p <- ggplot(cox_results, aes(x = log(HR), y = -log10(p.value))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(data = subset(cox_results, FDR < 0.05),
                   aes(label = Gene), max.overlaps = 20) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Cox Regression",
       x = "Log Hazard Ratio", 
       y = "-Log10(P-value)") +
  theme_bw(base_size = 14)
ggsave("新cox火山图.pdf", plot = p, width =6, height = 5)
install.packages('forestplot')
library(forestplot)



# 提取前10显著基因
top_genes <- head(cox_results[order(cox_results$p.value), ], 10)
# ========== 数据映射调整 ==========
forest_data <- top_genes %>% 
  rename(OR = HR,        # 将HR映射到OR字段
         Lower = CI_low,
         Upper = CI_high,
         Varnames = Gene) %>%
  mutate(Sample = scales::rescale(n_effective, to = c(3, 10))) # 标准化样本量到点大小范围

# 绘制高级森林图
dev.off()
pdf("cox森林图.pdf", width =7, height = 5)
# ========== 增强版森林图 ==========
advanced_forest <- ggplot(forest_data, aes(x = OR, y = reorder(Varnames, OR))) +
  # 渐变背景区域（新增透明度调整）
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#2166AC", alpha = 0.08) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#B2182B", alpha = 0.08) +
  
  # 动态误差线（新增颜色映射）
  geom_errorbarh(aes(xmax = Upper, xmin = Lower, color = OR > 1), 
                 height = 0, size = 1.2, alpha = 0.8) +
  
  # 渐变数据点（新增形状和透明度）
  geom_point(aes(size = Sample, fill = -log10(p.value)), 
             shape = 21, color = "white", stroke = 0.8, alpha = 0.9) +
  
  # 智能参考线（新增动态位置标签）
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30", size = 0.8) +
  annotate("text", x = 1, y = Inf, label = "Neutral Effect", 
           vjust = -1, hjust = -0.1, color = "gray30", size = 3.5) +
  
  # 专业级配色方案（新增连续型颜色梯度）
  scale_fill_gradientn(
    colors = c("#4393C3", "#F4A582", "#D6604D"),
    name = expression(-log[10](P-value)),
    limits = c(2, 6),  # 根据p.value范围调整
    breaks = seq(2, 6, 1)
  ) +
  scale_color_manual(values = c("#4393C3", "#D6604D"), guide = "none") +
  
  # 动态坐标轴（新增智能范围计算）
  scale_x_continuous(
    name = "Hazard Ratio (95% CI)",
    trans = "log2",  # 对数变换展示
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
    labels = c("0.5", "0.75", "1.0 (Ref)", "1.25", "1.5", "2.0"),
    limits = c(0.4, 2.5)
  ) +
  
  # 期刊级主题配置（新增科学可视化元素）
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 11),  # 基因名斜体
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    plot.caption = element_text(color = "grey40"),
    plot.margin = unit(c(1, 2, 1, 1), "cm")
  ) +
  labs(
    title = "Genomic Risk Factors Forest Plot",
    subtitle = "Circle size reflects effective sample size\nColor intensity represents -log10(p-value)",
    caption = "FDR-adjusted significance level: 0.05"
  )

# 添加统计注释（新增自动标注显著位点）
advanced_forest <- advanced_forest +
  geom_label_repel(
    data = subset(forest_data, FDR < 0.05),
    aes(label = paste0("FDR = ", formatC(FDR, format = "e", digits = 1))),
    nudge_y = 0.3,
    size = 3,
    label.padding = unit(0.15, "lines"),
    fill = alpha("white", 0.9)
  )

print(advanced_forest)

dev.off()
# ========== 数据预处理增强 ==========
library(rms)
library(ggplot2)
library(ggtext)  # 增强文本渲染
# ========== 数据预处理 ==========
# 确保样本ID一致（假设行名匹配）
combined_exp_scaled_t <- t(combined_exp_scaled)
stopifnot(rownames(combined_meta) == rownames(combined_exp_scaled_t))

# 创建生存对象（假设时间/状态列名为OS_time/OS_status）
surv_obj <- with(combined_meta, Surv(days_to_death, event))
colnames(combined_meta)
# 合并临床与表达数据（关键步骤）
model_data <- cbind(
  combined_meta[, c("Age", "Gender", "Stage", "Treatment")],  # 示例临床变量
  combined_exp_scaled_t
)

# 设置数据分布（rms包要求）
dd <- datadist(model_data)
options(datadist = "dd")

# ========== 特征筛选 ==========
# 使用LASSO进行基因表达筛选（避免维度灾难）
set.seed(123)
lasso_fit <- cv.glmnet(
  x = as.matrix(model_data[,5:ncol(model_data)]),  # 从第5列开始是基因
  y = surv_obj,
  family = "cox",
  alpha = 1
)

# 提取显著基因（lambda.1se标准）
selected_genes <- colnames(model_data)[5:ncol(model_data)][which(coef(lasso_fit, s = "lambda.1se") != 0)]

# ========== 模型构建 ==========
# 动态构建公式（处理可能零选择的情况）
selected_genes <- gene_sig
if(length(selected_genes) > 0) {
  formula_str <- paste(
    "surv_obj ~ ",
    paste(selected_genes)
  )
} else {
  formula_str <- "surv_obj ~ Age + Gender + Stage + Treatment"
}

# 构建最终模型（添加模型验证参数）
final_model <- cph(
  as.formula(formula_str),
  data = model_data,
  x = TRUE, y = TRUE,  # 保留预测矩阵
  surv = TRUE, 
  time.inc = 365,      # 设置年度风险单位
  singular.ok = FALSE  # 防止共线性
)

# 动态数据分布识别（增加变量类型校验）

library(pheatmap)
# 提取显著基因表达矩阵
sig_exp <- combined_exp_scaled[significant_genes$Gene, ]
head(sig_exp)[1:5,1:5]
library(pheatmap)
library(RColorBrewer)
library(DescTools)
library(caret)
# 创建按Stage排序的样本索引
# 步骤1：剪裁极端值
# 设置cofactor（经验值=5）
# 剪裁1%极端值
library(DescTools)
# ---------------------------------------------------------------
# 步骤1：数据预处理（Winsorize裁剪）
# ---------------------------------------------------------------
# 剪裁5%极端值（上下各2.5%）
exp_winsor <- Winsorize(sig_exp_ordered, probs = c(0.025, 0.975))
# 创建按Stage排序的样本索引
stage_order <- order(combined_meta$Stage)
sig_exp_ordered <- sig_exp[, stage_order]

stage_order <- order(combined_meta$Stage)
sig_exp_ordered <- sig_exp[, stage_order]

# 优化注释颜色方案
annotation_colors <- list(
  Survival = c(Dead = "#d73027", Alive = "#1a9850"),
  Stage = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(unique(combined_meta$Stage)))
)
names(annotation_colors$Stage) <- sort(unique(combined_meta$Stage))

# 生成高级热图
pheatmap(
  sig_exp_ordered,
  annotation_col = annotation_col[stage_order, ],  # 保持排序一致性
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,  # 关闭列聚类（已按Stage排序）
  cluster_rows = TRUE,
  clustering_method = "ward.D2",  # 改进聚类方法
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),  # 蓝白红渐变色
  scale = "row",  # 按行标准化
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  main = "Gene Expression Patterns Stratified by Tumor Stage",
  gaps_col = cumsum(table(combined_meta$Stage[stage_order])),  # 按Stage添加分隔线
  treeheight_row = 30,  # 调整行聚类树高度
  cellwidth = 0.5,      # 调整单元格宽度
  cellheight = 12       # 调整行高
)
# 优化注释颜色方案
annotation_colors <- list(
  Survival = c(Dead = "#d73027", Alive = "#1a9850"),
  Stage = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(unique(combined_meta$Stage)))
)
names(annotation_colors$Stage) <- sort(unique(combined_meta$Stage))
pdf('新cox热图.pdf',height = 8,width = 9)
# 生成高级热图
pheatmap(
  sig_exp_ordered,
  annotation_col = annotation_col[stage_order, ],  # 保持排序一致性
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,  # 关闭列聚类（已按Stage排序）
  cluster_rows = TRUE,
  clustering_method = "ward.D2",  # 改进聚类方法
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(99), # 红蓝渐变色
  breaks = seq(-1, 2, length.out = 100), # 根据实际范围调整
  scale = "row",  # 按行标准化
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  main = "Gene Expression Patterns Stratified by Tumor Stage",
  gaps_col = cumsum(table(combined_meta$Stage[stage_order])),  # 按Stage添加分隔线
  treeheight_row = 30,  # 调整行聚类树高度
  cellwidth = 0.15,      # 调整单元格宽度
  cellheight = 15       # 调整行高
)
dev.off()

# 提取显著基因表达数据
sig_genes_exp <- t(combined_exp_scaled[significant_genes$Gene, ])
colnames(sig_genes_exp) <- paste0("Gene_", significant_genes$Gene)

# 合并临床数据
multi_data <- cbind(combined_meta, sig_genes_exp)
dim(multi_data)
# 处理分类变量
multi_data$Stage <- factor(multi_data$Stage)
multi_data$Gender <- factor(multi_data$Gender)
library(survival)
library(splines)
library(glmnet)

# 1. 创建生存对象（若尚未创建）
multi_data$Surv_obj <- with(multi_data, 
                            Surv(time = days_to_last_follow_up,
                                 event = event))

# 2. 检查分类变量水平
cat("分类变量水平检查:\n")
sapply(multi_data[, c("Gender", "Stage", "Ethnicity", "Race", "T", "M", "N")], 
       function(x) {
         lv <- length(unique(x))
         if(lv < 2) warning(paste0(deparse(substitute(x)), " 仅有一个水平"))
         return(lv)
       })

# 3. 处理单水平变量（示例：若Gender只有"Female"）


# 4. 转换分类变量为因子
fac_vars <- c("Stage", "T", "M", "N", "Ethnicity", "Race")
multi_data[fac_vars] <- lapply(multi_data[fac_vars], factor)

# 5. 处理缺失值
print(paste("缺失值占比:", mean(is.na(multi_data))))
multi_data <- na.omit(multi_data)  # 删除含缺失值的行
# 手动构建包含目标基因的固定模型
fixed_vars <- colnames(sig_genes_exp)
optional_vars <- c("Stage", "T", "M", "N")

# 逐步回归仅对可选变量进行筛选
step_model <- step(
  coxph(as.formula(paste("Surv_obj ~", paste(fixed_vars, collapse = " + "))), 
        data = multi_data),
  scope = list(
    lower = ~ .,
    upper = as.formula(paste("~ . +", paste(optional_vars, collapse = " + ")))
  ),
  direction = "both",
  trace = 0
)

# 验证结果
summary(step_model)
library(forestplot)
library(dplyr)

# 提取模型结果并整理
res_df <- summary(step_model)$coefficients %>%
  as.data.frame() %>%
  mutate(
    Variable = rownames(.),
    HR = exp(coef),
    CI_low = exp(coef - 1.96 * `se(coef)`),
    CI_high = exp(coef + 1.96 * `se(coef)`)
  ) %>%
  filter(!grepl("ns\\(Age|StageStage|TTX|NNX|MMX", Variable))  # 过滤非关键变量

# 自定义颜色分组
res_df <- res_df %>%
  mutate(
    Group = case_when(
      grepl("Gene_", Variable) ~ "基因标记",
      grepl("TT|NN|MM", Variable) ~ "TNM分期",
      TRUE ~ "其他临床变量"
    )
  )

# 绘制高级森林图
forestplot(
  labeltext = cbind(res_df$Variable, 
                    sprintf("%.2f (%.2f-%.2f)", res_df$HR, res_df$CI_low, res_df$CI_high)),
  mean = res_df$HR,
  lower = res_df$CI_low,
  upper = res_df$CI_high,
  graph.pos = 3,
  xticks = c(0.5, 1, 2, 3),
  col = fpColors(box = c("#2c7bb6", "#d7191c", "#fdae61")[as.factor(res_df$Group)]),
  boxsize = 0.2,
  lwd.ci = 2,
  title = "多因素Cox回归风险比可视化",
  xlab = "Hazard Ratio (95% CI)",
  hrzl_lines = list("3" = gpar(lty = 1)),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1.1)
  )
)
library(plotly)

# 生成预测网格数据
age_seq <- seq(min(multi_data$Age), max(multi_data$Age), length.out = 50)
gene_seq <- seq(-2, 2, length.out = 50)
grid_data <- expand.grid(Age = age_seq, Gene_SPDEF = gene_seq)

# 添加其他变量中位数
fixed_values <- data.frame(
  Gene_CLU = median(multi_data$Gene_CLU),
  Gene_FOSL1 = median(multi_data$Gene_FOSL1),
  Gene_ALDOA = median(multi_data$Gene_ALDOA),
  Gene_HMGA1 = median(multi_data$Gene_HMGA1),
  T = "T1",
  N = "N0",
  M = "M0",
  Stage = "Stage I"
)

# 计算预测风险
pred_risk <- predict(step_model, 
                     newdata = cbind(grid_data, fixed_values),
                     type = "risk")

# 创建交互式曲面图
plot_ly(
  x = ~age_seq, 
  y = ~gene_seq, 
  z = ~matrix(pred_risk, nrow = 50),
  type = "surface",
  colorscale = list(c(0, "#2c7bb6"), c(1, "#d7191c")),
  hoverinfo = "x+y+z",
  contours = list(
    z = list(show = TRUE, start = 0, end = 8, size = 0.5)
  )
) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Age"),
      yaxis = list(title = "Gene_SPDEF (Z-score)"),
      zaxis = list(title = "Predicted Risk"),
      camera = list(eye = list(x = -1.5, y = -1.5, z = 0.5))
    ),
    title = "Age-Gene交互风险曲面"
  )
library(rms)

library(regplot)

# 转换数据格式
dd <- datadist(multi_data)
options(datadist = "dd")

library(rms)
library(regplot)

# ----------------------------
# 步骤1：数据预处理与分布定义
# ----------------------------
# 清理无效因子水平
multi_data <- droplevels(multi_data)

# 重新定义数据分布（必须包含所有模型变量）
dd <- datadist(multi_data[, c("Age", "Gene_ALDOA", 
                              "Gene_HMGA1", "Gene_SPDEF", "T", "N", "M", "Stage")])
options(datadist = "dd")

# ----------------------------
# 步骤2：使用rms包重建Cox模型
# ----------------------------
# 重新拟合模型（确保使用rms::cph）
colnames(multi_data)
# 检查 N 的分布和事件情况
table(multi_data$N)  # 查看各类别样本量
table(multi_data$N, multi_data$event)  # 查看 N 各类别的事件数

# 若存在 'NX' 类别且无事件，考虑合并/删除该类别
multi_data$N <- ifelse(multi_data$N == "NX", NA, multi_data$N)  # 删除 NX 类别
multi_data <- na.omit(multi_data)  # 移除 NA
# 删除高度相关的变量（如 Stage 可能由 T/N/M 决定）
cph_model <- cph(
  Surv(days_to_death, event) ~ Gene_CLU + Gene_FOSL1 + Gene_ALDOA + 
    Gene_HMGA1 + Gene_SPDEF + rcs(Age, 3) + Stage,  # 移除 Stage
  data = multi_data, surv = TRUE, x = TRUE, y = TRUE, na.action = na.delete
)


surv<-Survival(cph_model) 
surv3<-function(x) surv(36,x)
surv5<-function(x) surv(60,x)
nomo <- regplot(cph_model,
                #对观测2的六个指标在列线图上进行计分展示

                #预测3年和5年的死亡风险，此处单位是month
                failtime = c(365,700), 
                prfail = TRUE, #cox回归中需要TRUE
                showP = T, #是否展示统计学差异
                droplines = F,#观测2示例计分是否画线
                rank="sd", #根据统计学差异的显著性进行变量的排序
                interval="confidence") #展示观测的可信区间
 
cal1 <- calibrate(cph_model, cmethod='KM', method="boot", 
                  u=30, # u需要与之前模型中定义好的time.inc一致；
                  m=50, #每次抽样的样本量，
                  B=1000) #抽样次数
## 绘制校正曲线
plot(cal1,lwd=1,lty=1,
     conf.int=T,# 是否显示置信区间
     errbar.col="blue",#直线曲线bar颜色
     col="red", # 曲线颜色
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 1-Year DFS",
     ylab="Actual 1-Year DFS (proportion)",
     subtitles = F)#不显示副标题
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("black")#对角线的颜色
)
# ----------------------------
# 步骤3：构建Nomogram
# ----------------------------
nomogram <- nomogram(
  cph_model,
  fun = function(x) 1 - 0.9^x,  # 自定义生存概率转换
  fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),  # 设置概率刻度点
  funlabel = "3-Year Survival Probability",
  lp = FALSE,
  varname.label = TRUE,
  # 指定关键变量范围
)

# ----------------------------
# 步骤4：绘制高级Nomogram
# ----------------------------
plot(nomogram,
     col.grid = c("#d7191c", "#2c7bb6"),
     lwd = 2,
     cex.axis = 0.8,
     cex.lab = 0.9,
     xfrac = 0.2,
     label.every = 2)

# 绘制高级nomogram
regplot(nomogram,
        points = TRUE,
        rank = "sd",
        droplines = TRUE,
        title = "Nomogram of Prognostic Model",
        colors = c("#2c7bb6", "#d7191c"),
        plot.opts = list(cex.axis = 0.8, cex.lab = 0.9))
#----------------------------------------------------------
# 步骤3：执行分析（示例使用前5个基因，实际应遍历所有基因）
#----------------------------------------------------------
# 获取基因列表
gene_list <- rownames(combined_exp_scaled) # 测试用前5个基因，实际改为 rownames(combined_exp)

# 并行计算加速（可选）
library(furrr)
plan(multisession, workers = 4) 
cox_results <- future_map_dfr(gene_list, run_cox, .progress = TRUE)
# 将基因表达量分为高/低表达组（按中位数）
analysis_data$expr_group <- ifelse(analysis_data$Expression > median(analysis_data$Expression), 
                                   "High", "Low")

# 使用分组变量进行Cox分析
coxph(Surv_obj ~ expr_group, data = analysis_data)
#----------------------------------------------------------
# 步骤4：结果整理与筛选
#----------------------------------------------------------
significant_genes <- cox_results %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) %>%
  select(Gene, estimate, p.value, conf.low, conf.high)


combined_meta$days_to_death <- as.numeric(combined_meta$days_to_death)
head(combined_meta$days_to_death)
# 步骤 2：转换状态变量为数值型（假设 Status 是 "Dead"/"Alive" 字符型）
combined_meta$event <- ifelse(combined_meta$Status == "Dead", 1, 0)

combined_meta$Surv_obj <- with(combined_meta, Surv( days_to_death,event))
genes <- rownames(combined_exp)
results <- lapply(genes, function(gene) {
  expr <- t(combined_exp[gene, ])
  tmp_data <- merge(combined_meta, expr, by.x = "ID", by.y = "row.names")
  coxph(Surv_obj ~ . , data = tmp_data[, c("Surv_obj", gene, "Age", "Gender")])
})
exp_filtered <- exp[, common_samples]
meta_filtered <- meta[common_samples, ]
library(survival)
meta_filtered$OS.time <- as.numeric(meta_filtered$OS.time)
meta_filtered$OS <- as.numeric(meta_filtered$OS)
surv_obj <- Surv(time = meta_filtered$OS.time, 
                 event = meta_filtered$OS)
results <- apply(exp_filtered, 1, function(gene_exp) {
  cox_model <- coxph(surv_obj ~ gene_exp)
  summary(cox_model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
})

# 转换为数据框
results_df <- data.frame(
  Gene = rownames(exp_filtered),
  HR = sapply(results, "[", 1),
  Pvalue = sapply(results, "[", 2),
  row.names = NULL
)
head(results_df)

