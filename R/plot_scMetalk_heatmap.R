#' 绘制代谢物丰度热图
#'
#' @param seurat_obj 包含代谢物丰度的Seurat对象
#' @export
plot_scMetalk_heatmap <- function(seurat_obj) {
  if (!"metabolite_abundance" %in% names(seurat_obj@misc)) {
    stop("请先运行calculate_metabolite_abundance()。")
  }
  mat <- seurat_obj@misc$metabolite_abundance
  mat = mat[-which(is.nan(mat)), ]
  pheatmap::pheatmap(
    mat, scale = "row",
    color = viridis::viridis(50),
    main = "Metabolite Abundance by Cell Type",
    clustering_method = "average"
  )
}
