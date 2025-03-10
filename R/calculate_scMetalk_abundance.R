#' 根据代谢酶基因表达计算代谢物丰度
#' 
#' @param seurat_obj Seurat对象（需包含cell_type列）
#' @param enzyme_file 代谢酶文件路径（包含Metabolite和Gene列）
#' @param method 基因整合方法（"mean"或"sum"）
#' @return 代谢物丰度矩阵（行：代谢物，列：细胞类型）
#' @export
calculate_scMetalk_abundance <- function(seurat_obj, group.by = "cell_type",
                                           species = "mouse",method = "mean") {

  # 读取代谢酶关系表
  enzyme_df <- load_metabolic_data("enzyme", species)
  if (!all(c("Metabolite", "Gene") %in% colnames(enzyme_df))) {
    stop("代谢酶文件必须包含'Metabolite'和'Gene'列。")
  }
  
  # 获取细胞类型平均表达矩阵
  avg_expr <- data.frame(Seurat::AverageExpression(seurat_obj, 
                                        group.by = group.by)$RNA)
  
  # 初始化代谢物丰度矩阵
  cell_types <- colnames(avg_expr)
  metabolites <- unique(enzyme_df$Metabolite)
  metabolite_abundance <- matrix(
    0, 
    nrow = length(metabolites), 
    ncol = length(cell_types),
    dimnames = list(metabolites, cell_types)
  )
  
  # 遍历每个代谢物
  for (metabolite in metabolites) {
    # 提取相关代谢酶基因
    genes <- unique(enzyme_df$Gene[enzyme_df$Metabolite == metabolite])
    valid_genes <- intersect(genes, rownames(avg_expr))
    
    if (length(valid_genes) == 0) {
      warning(paste("代谢物", metabolite, "的代谢酶基因在数据中均未找到。"))
      next
    }
    
    # 计算基因表达得分
    scores <- avg_expr[valid_genes, ]
    if (method == "mean") {
      score <- colMeans(scores)
    } else if (method == "sum") {
      score <- colSums(scores)
    } else {
      stop("method参数必须是'mean'或'sum'。")
    }
    
    # 更新丰度矩阵
    metabolite_abundance[metabolite, ] <- score
  }
  
  # 将丰度矩阵添加到Seurat对象的元数据中
  seurat_obj@misc$metabolite_abundance <- metabolite_abundance
  return(seurat_obj)
}
