# 安装并加载必要的R包
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}
if (!require("limma")) {
  install.packages("limma")
  library(limma)
}

# ----------------------
# 1. 读取数据
# ----------------------
gene_path <- "C:/Users/吕杭/OneDrive/桌面/GENE.txt"
clinical_path <- "C:/Users/吕杭/OneDrive/桌面/clinical.txt"

# 读取基因数据
gene_data <- read.delim(gene_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# 读取临床数据
clinical_data <- read.delim(clinical_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(clinical_data) <- clinical_data$CGGA_ID  # 用CGGA_ID作为行名

# 打印所有临床列名（关键：让你确认实际列名）
cat("===== 临床数据所有列名如下 =====", "\n")
print(colnames(clinical_data))
cat("===============================", "\n")

# ----------------------
# 2. 基因表达数据预处理
# ----------------------
## 去除低表达基因
gene_mean <- rowMeans(gene_data)
gene_filtered <- gene_data[gene_mean >= 1, , drop=FALSE]
cat("去除低表达基因后：", nrow(gene_filtered), "个基因保留\n")

## log2标准化
if (max(gene_filtered) > 100) {
  gene_log2 <- log2(gene_filtered + 1)
  cat("已完成log2标准化\n")
} else {
  gene_log2 <- gene_filtered
  cat("数据已标准化，无需log2转换\n")
}

## 转置基因数据
gene_t <- as.data.frame(t(gene_log2))
cat("转置后基因数据：", nrow(gene_t), "个样本，", ncol(gene_t), "个基因\n")

# ----------------------
# 3. 临床数据预处理（修正：移除缺失的1p19q列）
# ----------------------
## 3.1 提取关键特征（移除缺失的1p19q_codeletion_status）
key_features <- c(
  "OS",  # 生存时间
  "Censor..alive.0..dead.1.",  # 生存状态（根据实际列名）
  "Grade",  # WHO分级
  "IDH_mutation_status"  # IDH突变状态（根据实际列名）
)

# 检查列名是否存在
missing_cols <- setdiff(key_features, colnames(clinical_data))
if (length(missing_cols) > 0) {
  stop("临床数据中缺少以下列：", paste(missing_cols, collapse=", "), "\n请根据上方打印的列名修正key_features！")
}

# 提取特征并修正列名
clinical_filtered <- clinical_data[, key_features, drop=FALSE]
colnames(clinical_filtered) <- c(
  "OS", 
  "Censor",  # 生存状态重命名
  "WHO_grade", 
  "IDH_status"  # IDH状态重命名
)

## 3.2 处理缺失值
missing_rate <- colMeans(is.na(clinical_filtered))
cat("临床特征缺失率：\n")
print(missing_rate)

for (col in colnames(clinical_filtered)) {
  if (is.numeric(clinical_filtered[[col]])) {
    clinical_filtered[[col]][is.na(clinical_filtered[[col]])] <- median(clinical_filtered[[col]], na.rm=TRUE)
  } else {
    freq <- table(clinical_filtered[[col]])
    clinical_filtered[[col]][is.na(clinical_filtered[[col]])] <- names(freq)[which.max(freq)]
  }
}

## 3.3 转换数据类型
clinical_filtered$Censor <- as.integer(clinical_filtered$Censor)  # 生存状态转为整数
clinical_filtered$WHO_grade <- as.factor(clinical_filtered$WHO_grade)
clinical_filtered$IDH_status <- as.factor(clinical_filtered$IDH_status)

# ----------------------
# 4. 样本匹配
# ----------------------
common_samples <- intersect(rownames(gene_t), rownames(clinical_filtered))
gene_matched <- gene_t[common_samples, , drop=FALSE]
clinical_matched <- clinical_filtered[common_samples, , drop=FALSE]

cat("样本匹配完成：共保留", length(common_samples), "个样本\n")

# ----------------------
# 5. 保存结果
# ----------------------
output_dir <- "C:/Users/吕杭/OneDrive/桌面/preprocessed_data"
if (!dir.exists(output_dir)) dir.create(output_dir)

save(gene_matched, clinical_matched, file=file.path(output_dir, "preprocessed_data.RData"))
write.csv(gene_matched, file=file.path(output_dir, "gene_matched.csv"))
write.csv(clinical_matched, file=file.path(output_dir, "clinical_matched.csv"))

cat("预处理完成！结果保存至：", output_dir, "\n")