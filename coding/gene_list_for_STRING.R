# 导出基因列表到文本文件（确保在之前的代码之后运行）
output_dir <- "C:/Users/吕杭/OneDrive/桌面/Enrichment_results"
genes_file <- file.path(output_dir, "gene_list_for_STRING.txt")

# 保存基因列表（每个基因占一行）
writeLines(all_genes, con = genes_file)

cat("基因列表已保存至：", genes_file, "\n")
cat("共", length(all_genes), "个基因\n")