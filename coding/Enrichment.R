# 第一步：设置网络超时和安装必要包
# 延长网络超时时间至120秒（默认60秒）
options(timeout = 120)

# 确保必要包已安装
if (!require("BiocManager")) install.packages("BiocManager")
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ----------------------
# 第二步：加载差异基因数据
# ----------------------
degs <- read.csv("C:/Users/吕杭/OneDrive/桌面/DEGs_results/筛选后差异基因列表.csv")
gene_list <- unique(as.character(degs$Gene))
cat("去重后基因数量：", length(gene_list), "\n")

# ----------------------
# 第三步：基因ID转换（忽略映射失败的警告）
# ----------------------
# 抑制部分基因映射失败的警告（不影响主要分析）
suppressWarnings({
  gene_df <- bitr(
    geneID = gene_list,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
})
cat("成功转换ID的基因数：", nrow(gene_df), "\n")
entrez_ids <- gene_df$ENTREZID

# ----------------------
# 第四步：GO富集分析（无需网络，可正常运行）
# ----------------------
go_enrich <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
go_results <- as.data.frame(go_enrich) %>% arrange(pvalue)
cat("GO显著富集条目数：", nrow(go_results), "\n")

# ----------------------
# 第五步：KEGG富集分析（解决连接问题）
# ----------------------
# 方案1：重试连接（最多3次）
kegg_enrich <- NULL
retry_count <- 0
max_retries <- 3

while (is.null(kegg_enrich) && retry_count < max_retries) {
  retry_count <- retry_count + 1
  cat("\n第", retry_count, "次尝试连接KEGG数据库...\n")
  
  tryCatch({
    kegg_enrich <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pAdjustMethod = "fdr",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE  # 优先使用在线数据
    )
  }, error = function(e) {
    cat("连接失败：", e$message, "\n")
  })
}

# 方案2：若仍失败，使用离线缓存数据（预先下载的KEGG通路信息）
if (is.null(kegg_enrich)) {
  warning("KEGG在线连接多次失败，尝试使用离线缓存数据（可能不是最新）")
  
  # 检查是否有本地缓存（若没有则自动下载一次）
  if (!file.exists("kegg_hsa_2023.RData")) {
    cat("正在下载KEGG离线缓存数据...\n")
    url <- "https://ndownloader.figshare.com/files/25950042"  # 第三方镜像的缓存数据
    download.file(url, "kegg_hsa_2023.RData", mode = "wb")
  }
  load("kegg_hsa_2023.RData")  # 加载离线通路数据（包含kegg_df数据框）
  
  # 使用离线数据进行富集分析
  kegg_enrich <- enricher(
    gene = entrez_ids,
    TERM2GENE = kegg_df[, c("pathway_id", "entrezid")],  # 通路-基因映射
    TERM2NAME = kegg_df[, c("pathway_id", "pathway_name")],  # 通路ID-名称映射
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05
  )
}

kegg_results <- as.data.frame(kegg_enrich) %>% arrange(pvalue)
cat("KEGG显著富集通路数：", nrow(kegg_results), "\n")

# ----------------------
# 第六步：可视化与保存结果
# ----------------------
output_dir <- "C:/Users/吕杭/OneDrive/桌面/Enrichment_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 保存GO结果
if (nrow(go_results) > 0) {
  write.csv(go_results, file.path(output_dir, "GO富集结果.csv"), row.names = FALSE)
  
  # 调整GO富集图参数（解决标签拥挤）
  go_plot <- dotplot(go_enrich, 
                     showCategory = 10,  # 减少显示条目数（原15→10）
                     split = "ONTOLOGY", 
                     font.size = 8,       # 整体字体缩小（原默认12→8）
                     label_format = 30) + # 标签最长30字符（超长截断）
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    ggtitle("GO富集分析") +
    # 旋转y轴标签（45度倾斜，右对齐）
    theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
          # 增加绘图区域宽度（按需调整）
          plot.width = unit(15, "cm"))
  
  print(go_plot)
  # 保存时增大宽度（如15cm→20cm，按需调整）
  ggsave(file.path(output_dir, "GO富集图.png"), 
         go_plot, 
         width = 20,  # 宽度增加（原12→20）
         height = 10, 
         dpi = 300)
}

# 保存KEGG结果
if (nrow(kegg_results) > 0) {
  write.csv(kegg_results, file.path(output_dir, "KEGG富集结果.csv"), row.names = FALSE)
  kegg_plot <- barplot(kegg_enrich, showCategory = 15) + ggtitle("KEGG富集分析")
  print(kegg_plot)
  ggsave(file.path(output_dir, "KEGG富集图.png"), kegg_plot, width = 10, height = 8, dpi = 300)
}

cat("\n分析完成！结果保存至：", output_dir, "\n")
cat("若KEGG结果仍为空，可能是网络问题，建议稍后重试或检查网络设置\n")
