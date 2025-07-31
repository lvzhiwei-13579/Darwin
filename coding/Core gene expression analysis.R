# 第一步：重置Bioconductor仓库设置（解决仓库替换问题）
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}
BiocManager::repositories()  # 显示当前仓库设置
BiocManager::install(ask = FALSE, update = FALSE)  # 重置仓库

# 第二步：分步安装必要包（优先解决依赖问题）
# 1. 安装基础包
required_packages <- c("igraph", "dplyr", "ggplot2", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}

# 2. 安装Bioconductor包（使用国内镜像加速）
BiocManager::install(
  c("org.Hs.eg.db", "biomaRt"),
  site_repository = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/",  # 清华镜像
  update = FALSE,
  ask = FALSE
)

# 第三步：验证安装
library(org.Hs.eg.db)  # 关键：人类基因ID数据库
library(biomaRt)
library(igraph)
library(dplyr)

# 第四步：测试biomaRt连接（使用国内可访问的服务器）
tryCatch({
  ensembl <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "asia.ensembl.org"  # 亚洲服务器，速度更快
  )
  cat("biomaRt连接成功！\n")
}, error = function(e) {
  warning("biomaRt连接失败，将使用离线ID转换方案\n")
})
# 第一步：安装并加载必要的程序包
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("biomaRt", "org.Hs.eg.db"), update = FALSE)  # 新增org.Hs.eg.db用于ID转换

required_packages <- c("igraph", "dplyr", "ggplot2", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(igraph)
library(dplyr)
library(ggplot2)
library(readr)
library(biomaRt)
library(org.Hs.eg.db)  # 人类基因ID转换数据库

# ----------------------
# 第二步：ID转换核心函数（Entrez ID -> 基因符号）
# ----------------------
# 将你的Entrez ID列表转换为基因符号（匹配STRING转换结果）
convert_entrez_to_symbol <- function(entrez_ids) {
  # 使用org.Hs.eg.db数据库进行转换（比biomaRt更稳定）
  symbols <- mapIds(
    x = org.Hs.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",  # 输入类型：Entrez ID
    column = "SYMBOL",     # 输出类型：基因符号
    multiVals = "first"    # 多个匹配时取第一个
  )
  # 移除转换失败的NA值
  symbols <- symbols[!is.na(symbols)]
  return(as.character(symbols))
}

# ----------------------
# 第三步：加载数据
# ----------------------
# 加载KEGG富集结果和差异基因数据（你的基因列表是Entrez ID）
kegg_results <- read.csv("C:/Users/吕杭/OneDrive/桌面/Enrichment_results/KEGG富集结果.csv")
degs <- read.csv("C:/Users/吕杭/OneDrive/桌面/DEGs_results/筛选后差异基因列表.csv")

# 定义STRING数据路径
string_file <- "C:\\Users\\吕杭\\OneDrive\\桌面\\Enrichment_results\\string_interactions.txt"
if (!file.exists(string_file)) {
  stop(paste("未找到文件，请确认路径：", string_file))
}

# 读取STRING数据（按空格分割，跳过表头行）
interactions <- read.table(
  string_file,
  sep = "",
  header = FALSE,
  skip = 1,  # 跳过第一行重复表头
  stringsAsFactors = FALSE
)
colnames(interactions) <- c("protein1", "protein2", "combined_score")

# 筛选置信度阈值（可适当降低，如300）
confidence_threshold <- 400
interactions <- interactions %>% 
  filter(combined_score > confidence_threshold)
cat("加载的互作关系数（过滤后）：", nrow(interactions), "\n")

# ----------------------
# 第四步：处理基因列表（将Entrez ID转换为基因符号）
# ----------------------
# 提取KEGG通路中的Entrez ID并转换为基因符号
significant_pathways <- kegg_results %>%
  filter(p.adjust < 0.05, !is.na(geneID))

if (nrow(significant_pathways) == 0) {
  significant_pathways <- kegg_results %>%
    filter(p.adjust < 0.1, !is.na(geneID))
  warning("显著通路较少，已放宽阈值至p.adjust < 0.1")
}

# 提取Entrez ID列表（假设geneID格式为"7143/3678/..."）
genes_list <- lapply(
  significant_pathways$geneID, 
  function(x) unlist(strsplit(as.character(x), "/"))
)
all_entrez <- unique(unlist(genes_list))
cat("显著通路中的Entrez ID总数：", length(all_entrez), "\n")
cat("前5个Entrez ID示例：", paste(head(all_entrez, 5), collapse = ", "), "\n")

# 转换为基因符号
all_symbols <- convert_entrez_to_symbol(all_entrez)
cat("转换后的基因符号总数：", length(all_symbols), "\n")
cat("前5个基因符号示例：", paste(head(all_symbols, 5), collapse = ", "), "\n")

# 若转换后数量不足，使用差异基因的Entrez ID转换
if (length(all_symbols) < 10) {
  warning("通路基因数量不足，使用差异基因进行分析")
  all_entrez <- unique(as.character(degs$Gene))  # 假设degs$Gene是Entrez ID
  all_symbols <- convert_entrez_to_symbol(all_entrez)
}

# ----------------------
# ----------------------
# 第五步：STRING的ENSP ID转换为基因符号（修复版）
# ----------------------
# 转换STRING中的protein1和protein2为基因符号（优化连接+分批处理）
convert_ensp_to_symbol <- function(ensp_ids, batch_size = 500) {
  # 1. 清理ENSP ID（去除物种前缀9606.）
  ensp_clean <- gsub("^9606\\.", "", ensp_ids)
  unique_ensp <- unique(ensp_clean)  # 去重减少转换量
  cat("需要转换的unique ENSP ID数量：", length(unique_ensp), "\n")
  
  # 2. 连接Ensembl服务器（使用https+亚洲服务器，稳定性更高）
  tryCatch({
    ensembl <- useMart(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      host = "https://asia.ensembl.org"  # 关键：添加https://，使用亚洲服务器
    )
  }, error = function(e) {
    # 备用服务器：若亚洲服务器失败，尝试欧洲服务器
    ensembl <- useMart(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      host = "https://www.ensembl.org"
    )
  })
  
  # 3. 分批转换（避免服务器负载过高）
  n_batches <- ceiling(length(unique_ensp) / batch_size)
  conv_list <- list()
  
  for (i in 1:n_batches) {
    cat("转换批次", i, "/", n_batches, "...\n")
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(unique_ensp))
    batch_ensp <- unique_ensp[start_idx:end_idx]
    
    # 重试机制：若单次失败，重试2次
    for (retry in 1:3) {
      tryCatch({
        batch_conv <- getBM(
          attributes = c("ensembl_peptide_id", "hgnc_symbol"),
          filters = "ensembl_peptide_id",
          values = batch_ensp,
          mart = ensembl
        )
        conv_list[[i]] <- batch_conv
        break  # 成功则跳出重试循环
      }, error = function(e) {
        if (retry < 3) {
          cat("第", retry, "次重试...\n")
          Sys.sleep(2)  # 等待2秒后重试
        } else {
          warning("批次", i, "转换失败，部分ID可能无法匹配\n")
          conv_list[[i]] <- data.frame(
            ensembl_peptide_id = batch_ensp,
            hgnc_symbol = NA
          )
        }
      })
    }
  }
  
  # 合并所有批次结果并构建转换字典
  conversion <- do.call(rbind, conv_list)
  conv_dict <- setNames(conversion$hgnc_symbol, conversion$ensembl_peptide_id)
  
  # 匹配原始输入的所有ID（包括重复）
  symbols <- sapply(ensp_clean, function(x) {
    ifelse(x %in% names(conv_dict), conv_dict[[x]], NA)
  })
  
  return(as.character(symbols))
}

# 转换并清洗STRING数据（保持原逻辑，使用修复后的转换函数）
interactions <- interactions %>%
  mutate(
    gene1 = convert_ensp_to_symbol(protein1),
    gene2 = convert_ensp_to_symbol(protein2)
  ) %>%
  filter(!is.na(gene1), !is.na(gene2)) %>%  # 移除转换失败的行
  mutate(
    gene1 = toupper(gene1),
    gene2 = toupper(gene2)
  )

# ----------------------
# 第六步：匹配基因列表（基因符号级别）
# ----------------------
filtered_interactions <- interactions %>%
  filter(gene1 %in% all_symbols | gene2 %in% all_symbols)

cat("与基因符号匹配的互作关系数：", nrow(filtered_interactions), "\n")

if (nrow(filtered_interactions) == 0) {
  # 最后的备选方案：降低置信度阈值
  warning("仍无匹配，尝试降低置信度阈值至300")
  interactions <- read.table(
    string_file,
    sep = "",
    header = FALSE,
    skip = 1,
    stringsAsFactors = FALSE
  )
  colnames(interactions) <- c("protein1", "protein2", "combined_score")
  interactions <- interactions %>% filter(combined_score > 300) %>%  # 降低阈值
    mutate(
      gene1 = convert_ensp_to_symbol(protein1),
      gene2 = convert_ensp_to_symbol(protein2)
    ) %>%
    filter(!is.na(gene1), !is.na(gene2)) %>%
    mutate(gene1 = toupper(gene1), gene2 = toupper(gene2))
  
  filtered_interactions <- interactions %>%
    filter(gene1 %in% all_symbols | gene2 %in% all_symbols)
  
  cat("降低阈值后匹配的互作关系数：", nrow(filtered_interactions), "\n")
  if (nrow(filtered_interactions) == 0) {
    stop("仍无匹配，请检查：1.基因列表是否为人类基因；2.STRING数据是否正确下载")
  }
}

# ----------------------
# 第七步：构建PPI网络（解决函数冲突问题）
# 强制转换为data.frame并明确使用dplyr函数
filtered_interactions <- as.data.frame(filtered_interactions)
network_data <- filtered_interactions %>%
  dplyr::select(gene1, gene2, combined_score) %>%  # 明确指定dplyr::select
  dplyr::rename(from = gene1, to = gene2, weight = combined_score)  # 明确指定dplyr::rename

# 构建网络
g <- graph_from_data_frame(network_data, directed = FALSE)
V(g)$degree <- degree(g)

# 识别核心基因
core_threshold <- quantile(V(g)$degree, 0.8, na.rm = TRUE)
V(g)$is_core <- V(g)$degree > core_threshold
core_genes <- names(V(g)[V(g)$is_core])
cat("核心基因数量：", length(core_genes), "\n")


# 第八步：可视化和保存结果
output_dir <- "C:/Users/吕杭/OneDrive/桌面/Enrichment_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 网络可视化参数
V(g)$color <- ifelse(V(g)$is_core, "#E41A1C", "#377EB8")
V(g)$size <- ifelse(V(g)$is_core, 12, 6)
V(g)$label <- ifelse(V(g)$is_core, names(V(g)), NA)

png(
  file.path(output_dir, "PPI网络_核心基因标注.png"),
  width = 1200,
  height = 1000,
  res = 300
)
par(mar = c(1, 1, 3, 1))
plot(
  g,
  layout = layout_with_fr,
  vertex.label.cex = 0.9,
  vertex.label.color = "black",
  edge.width = E(g)$weight / 200,
  edge.color = "gray50",
  main = "蛋白质互作网络（红色为核心基因）"
)
dev.off()

# 保存结果文件
write.csv(
  data.frame(gene_symbol = core_genes),
  file.path(output_dir, "核心基因列表.csv"),
  row.names = FALSE
)

write.csv(
  network_data,
  file.path(output_dir, "蛋白质互作关系_处理后.csv"),
  row.names = FALSE
)

cat("分析完成！结果保存至：", output_dir, "\n")