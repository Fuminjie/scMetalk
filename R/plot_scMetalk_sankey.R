#' 使用桑基图进行细胞间沟通的可视化
#'
#'#' @param comm_result 通讯强度矩阵（来自calculate_metabolic_communication函数）
#' @param seurat_obj Seurat对象（用于获取基因表达量）
#' @param top_n 展示的Top N通讯流
#' @export
plot_scMetalk_alluvial <- function(comm_result, seurat_obj, top_n = 10) {
  # 检查输入
  if (!"strength" %in% colnames(comm_result)) {
    stop("comm_result必须包含'strength'列。")
  }
  if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
    stop("Seurat对象必须包含'cell_type'列。")
  }

  # 获取Top N通讯流
  top_comm <- comm_result[order(-comm_result$strength), ][1:top_n, ]

  # 获取细胞类型平均表达矩阵
  avg_expr <- Seurat::AverageExpression(seurat_obj,
                                        group.by = "cell_type",
                                        assays = "RNA")$RNA

  # 构建冲积图数据
  plot_data <- top_comm %>%
    dplyr::select(metabolite, sender, transporter, receptor, receiver, strength) %>%
    dplyr::mutate(
      metabolite = paste0("Metabolite: ", metabolite),
      sender = paste0("Sender: ", sender),
      transporter = paste0("Transporter: ", transporter),
      receptor = paste0("Receptor: ", receptor),
      receiver = paste0("Receiver: ", receiver)
    )

  # 绘制冲积图
  library(ggalluvial)
  ggplot(plot_data,
         aes(axis1 = metabolite,
             axis2 = sender,
             axis3 = transporter,
             axis4 = receptor,
             axis5 = receiver,
             y = strength)) +
    geom_alluvium(aes(fill = metabolite), width = 1/12) +  # 根据代谢物映射颜色
    geom_stratum(width = 1/12, fill = "grey80", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Metabolite", "Sender", "Transporter", "Receptor", "Receiver")) +
    theme_void() +
    labs(title = "Top N Metabolic Communication Flows",
         x = "Metabolite, Cell Types and Genes",
         y = "Communication Strength",
         fill = "Metabolite") +  # 图例标题
    theme(legend.position = "bottom") +
    scale_fill_brewer(palette="Paired")
}
