#' 计算细胞间代谢通讯强度矩阵（按Transporter和Receptor对分别计算）
#'
#' @param seurat_obj Seurat对象（需包含cell_type列）
#' @param species 种属类别（"human"或"mouse"）
#' @param method 基因整合方法（"mean"或"sum"）
#' @return 包含输出细胞、输入细胞、代谢物、转运体、受体和通讯强度的数据框
#' @export
calculate_scMetalk_communication <- function(seurat_obj,method = "mean", species = "human") {
  
  transporter_df  <- load_metabolic_data("transporter", species)
  receptor_df <- load_metabolic_data("receptor", species)
  metabolite_abundance <- seurat_obj@misc$metabolite_abundance
  
  if (!all(c("Metabolite", "Gene") %in% colnames(transporter_df)) ||
      !all(c("Metabolite", "Gene") %in% colnames(receptor_df))) {
    stop("转运体和受体文件必须包含'Metabolite'和'Gene'列。")
  }
  
  # 获取细胞类型平均表达矩阵
  avg_expr <- Seurat::AverageExpression(seurat_obj, 
                                        group.by = "cell_type", 
                                        assays = "RNA")$RNA
  
  # 初始化结果数据框
  result <- data.frame(
    sender = character(),
    receiver = character(),
    metabolite = character(),
    transporter = character(),
    receptor = character(),
    strength = numeric(),
    stringsAsFactors = FALSE
  )

  # 遍历每个代谢物
  metabolites <- rownames(metabolite_abundance)
  for (metabolite in metabolites) {
    # 提取相关转运体和受体基因
    transporter_genes <- unique(transporter_df$Gene[transporter_df$Metabolite == metabolite])
    receptor_genes <- unique(receptor_df$Gene[receptor_df$Metabolite == metabolite])
    
    # 如果没有找到相关基因，跳过该代谢物
    if (length(transporter_genes) == 0 || length(receptor_genes) == 0) {
      warning(paste("代谢物", metabolite, "的转运体或受体基因未找到。"))
      next
    }
    
    # 遍历每个转运体和受体对
    for (transporter in transporter_genes) {
      for (receptor in receptor_genes) {
        # 计算转运体和受体表达得分
        transporter_score <- .calculate_gene_score(avg_expr, transporter, method)
        receptor_score <- .calculate_gene_score(avg_expr, receptor, method)
        
        # 遍历所有细胞类型对
        cell_types <- colnames(metabolite_abundance)
        for (sender in cell_types) {
          for (receiver in cell_types) {
            # 计算通讯强度
            strength <- metabolite_abundance[metabolite, sender] * 
              transporter_score[sender] * 
              receptor_score[receiver]
            
            # 添加到结果数据框
            result <- rbind(result, data.frame(
              sender = sender,
              receiver = receiver,
              metabolite = metabolite,
              transporter = transporter,
              receptor = receptor,
              strength = strength,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  return(result)
}

# 辅助函数：计算单个基因的表达得分
.calculate_gene_score <- function(expr_matrix, gene, method) {
  if (!gene %in% rownames(expr_matrix)) {
    return(setNames(rep(0, ncol(expr_matrix)), colnames(expr_matrix)))
  }
  
  if (method == "mean") {
    return(expr_matrix[gene, ])
  } else if (method == "sum") {
    return(expr_matrix[gene, ])
  } else {
    stop("method参数必须是'mean'或'sum'。")
  }
}