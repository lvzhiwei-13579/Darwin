# 安装并加载必要的R包
if (!require("limma")) install.packages("limma")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(limma)      # 差异表达分析
library(pheatmap)   # 热图可视化
library(ggplot2)    # 火山图可视化
library(dplyr)      # 数据处理

# ----------------------
# 1. 加载数据并验证结构
# ----------------------
load("C:/Users/吕杭/OneDrive/桌面/preprocessed_data/preprocessed_data.RData")

# 恢复基因矩阵正确结构（行=基因，列=样本）
gene_original <- t(gene_matched)
cat("基因矩阵维度：", nrow(gene_original), "个基因 ×", ncol(gene_original), "个样本\n")
cat("临床数据样本数：", nrow(clinical_matched), "\n")

# 验证样本ID匹配
if (!all(colnames(gene_original) %in% rownames(clinical_matched))) {
  warning("部分样本ID在临床数据中不存在，将自动过滤")
  valid_samples <- intersect(colnames(gene_original), rownames(clinical_matched))
  gene_original <- gene_original[, valid_samples, drop = FALSE]
  clinical_matched <- clinical_matched[valid_samples, , drop = FALSE]
  cat("过滤后样本数：", length(valid_samples), "\n")
}

# ----------------------
# 2. 分组处理
# ----------------------
# 修正WHO分级名称（空格→下划线）
who_grade <- as.character(clinical_matched$WHO_grade)
who_grade_fixed <- gsub(" ", "_", who_grade)
clinical_matched$WHO_grade_fixed <- who_grade_fixed

# 筛选WHO III和IV级样本
keep_grades <- c("WHO_III", "WHO_IV")
keep_idx <- clinical_matched$WHO_grade_fixed %in% keep_grades
gene_subset <- gene_original[, keep_idx, drop = FALSE]
clinical_subset <- clinical_matched[keep_idx, , drop = FALSE]

# 最终样本数检查
cat("\n筛选后样本数：", ncol(gene_subset), "\n")
if (ncol(gene_subset) < 2) stop("样本数不足，无法进行分析")

# 定义分组因子
group_subset <- factor(
  clinical_subset$WHO_grade_fixed,
  levels = keep_grades
)
cat("分组分布：\n")
print(table(group_subset))

# ----------------------
# 3. 差异表达分析
# ----------------------
## 设计矩阵
design <- model.matrix(~ 0 + group_subset)
colnames(design) <- levels(group_subset)

## 差异计算
fit <- lmFit(gene_subset, design)
contrast_matrix <- makeContrasts(WHO_IV - WHO_III, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

## 提取差异基因（新增Direction列标记上下调）
degs <- topTable(fit2, coef = 1, number = Inf, adjust = "fdr") %>%
  mutate(Gene = rownames(.)) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%  # 标记上调/下调
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val)

cat("\n差异基因数量：", nrow(degs), "\n")
if (nrow(degs) == 0) {
  warning("未发现显著差异基因，使用宽松阈值")
  degs <- topTable(fit2, coef = 1, number = Inf, adjust = "fdr") %>%
    mutate(Gene = rownames(.)) %>%
    mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%  # 宽松阈值下也标记
    filter(abs(logFC) > 0.5, adj.P.Val < 0.1) %>%
    arrange(adj.P.Val)
  cat("宽松阈值下差异基因数量：", nrow(degs), "\n")
}

# ----------------------
# 4. 可视化（新增行注释+拆分热图）
# ----------------------
## 4.1 火山图（保持原代码不变）
volcano_data <- topTable(fit2, coef = 1, number = Inf, adjust = "fdr") %>%
  mutate(
    Significance = case_when(
      abs(logFC) > 1 & adj.P.Val < 0.05 ~ "显著差异",
      abs(logFC) > 0.5 & adj.P.Val < 0.1 ~ "潜在差异",
      TRUE ~ "无差异"
    )
  )

火山图 <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("red", "orange", "gray50")) +
  labs(x = "log2(表达倍数变化)", y = "-log10(校正P值)", title = "WHO IV vs III 差异基因火山图") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(火山图)

## 4.2 热图（去掉聚类树和标题）
# 选择基因
n_genes <- min(50, nrow(degs))
if (n_genes < 2) n_genes <- min(20, nrow(topTable(fit2, coef = 1, number = 100)))
top_degs <- rownames(topTable(fit2, coef = 1, number = n_genes))

# 提取并清洗热图数据
heatmap_raw <- gene_subset[top_degs, , drop = FALSE]
heatmap_clean <- heatmap_raw[rowSums(is.na(heatmap_raw)) == 0, ]  # 移除含NA的基因
if (nrow(heatmap_clean) < 2) stop("有效基因不足，无法绘制热图")

# 标准化（限制范围在-3到3之间）
heatmap_scaled <- t(scale(t(heatmap_clean)))
heatmap_scaled <- pmax(pmin(heatmap_scaled, 3), -3)  

# 准备行注释（标记基因上调/下调）
gene_direction <- degs %>% 
  select(Gene, Direction) %>% 
  filter(Gene %in% rownames(heatmap_clean))  # 仅保留热图中存在的基因

annotation_row <- data.frame(Sig = gene_direction$Direction)
rownames(annotation_row) <- gene_direction$Gene
annotation_row <- annotation_row[rownames(heatmap_clean), , drop = FALSE]  # 匹配热图基因顺序

# 分组注释配置（列注释）
annotation_col <- data.frame(Group = group_subset)
rownames(annotation_col) <- colnames(heatmap_scaled)
annotation_colors <- list(
  Group = c(WHO_III = "#2196F3", WHO_IV = "#FF5722"),  # 样本分组颜色
  Sig = c(Up = "#FF6666", Down = "#6666FF")           # 基因上下调颜色
)

# 绘制整体热图（含行注释，去掉聚类树和标题）
heatmap_plot <- pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,   # 列注释（样本分组）
  annotation_row = annotation_row,   # 行注释（基因上下调）
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,  # 关闭行聚类 → 左侧无黑色树状线
  cluster_cols = FALSE,  # 关闭列聚类 → 顶部无黑色树状线
  treeheight_row = 0,    # 行聚类树高度设为0（双重保险）
  treeheight_col = 0,    # 列聚类树高度设为0（双重保险）
  color = colorRampPalette(c("#2C7FB8", "white", "#E41A1C"))(100),
  border_color = NA      # 去掉单元格边框
  # 直接删除 main 参数，或设为 main = ""
)
print(heatmap_plot)

# 拆分上下调基因，分别绘制热图（避免拥挤）
up_genes <- gene_direction %>% filter(Direction == "Up") %>% pull(Gene)
down_genes <- gene_direction %>% filter(Direction == "Down") %>% pull(Gene)

# 绘制上调基因热图
if (length(up_genes) > 0) {
  pheatmap(
    heatmap_scaled[up_genes, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,  # 关闭行聚类
    cluster_cols = FALSE,  # 关闭列聚类
    treeheight_row = 0,
    treeheight_col = 0,
    border_color = NA,     # 去掉单元格边框
    color = colorRampPalette(c("#2C7FB8", "white", "#E41A1C"))(100)
  )
}

# 绘制下调基因热图
if (length(down_genes) > 0) {
  pheatmap(
    heatmap_scaled[down_genes, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,  # 关闭行聚类
    cluster_cols = FALSE,  # 关闭列聚类
    treeheight_row = 0,
    treeheight_col = 0,
    border_color = NA,   # 去掉单元格边框
    color = colorRampPalette(c("#2C7FB8", "white", "#E41A1C"))(100)
  )
}

# ----------------------
# 5. 保存结果（含拆分热图）
# ----------------------
output_dir <- "C:/Users/吕杭/OneDrive/桌面/DEGs_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(degs, file.path(output_dir, "筛选后差异基因列表.csv"), row.names = FALSE)

# 保存火山图
ggsave(
  file.path(output_dir, "火山图.png"),
  火山图,
  width = 10,
  height = 7,
  dpi = 300,
  type = "cairo-png"
)

# 保存整体热图
png(
  file.path(output_dir, "热图.png"),
  width = 1200,
  height = 1000,
  res = 200,
  type = "cairo-png"
)
print(heatmap_plot)
dev.off()

# 保存上调基因热图
if (length(up_genes) > 0) {
  png(
    file.path(output_dir, "上调基因热图.png"),
    width = 1200,
    height = 1000,
    res = 200,
    type = "cairo-png"
  )
  print(pheatmap(
    heatmap_scaled[up_genes, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    color = colorRampPalette(c("#2C7FB8", "white", "#E41A1C"))(100)
  ))
  dev.off()
}

# 保存下调基因热图
if (length(down_genes) > 0) {
  png(
    file.path(output_dir, "下调基因热图.png"),
    width = 1200,
    height = 1000,
    res = 200,
    type = "cairo-png"
  )
  print(pheatmap(
    heatmap_scaled[down_genes, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    color = colorRampPalette(c("#2C7FB8", "white", "#E41A1C"))(100)
  ))
  dev.off()
}

cat("\n所有结果已保存至：", output_dir, "\n")