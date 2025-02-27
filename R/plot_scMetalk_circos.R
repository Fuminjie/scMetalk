#' 使用cirocs图可视化细胞间的代谢通讯
#'
#' @param comm_data 数据框
#' @return 可视化circos图
#' @export
plot_scMetalk_circos <- function(
    comm_data,          # 必须参数：通讯数据框
    top_n = 20,         # 显示前N条连接
    cell_colors = c("#A6CEE3", "#CAB2D6", "#B2DF8A",  "#FB9A99", "#FDBF6F","#1F78B4", "#6A3D9A", "#33A02C", "#E31A1C", "#FF7F00"),  # 细胞类型颜色梯度
    metabolite_palette = "Set1",          # 代谢物颜色调色板
    type_colors = c(transporter="#E31A1C", receptor="#3498db"), # 基因类型颜色
    track_heights = c(3, 10),  # 轨道高度(mm)
    arrow_params = list(  # 箭头参数
      type = "curved",
      length = 0.15,
      width_factor = 1,
      col = "black"
    ),
    legend_position = c(-0.2, 0.95),  # 图例起始位置
    legend_cex = 0.7                 # 图例字体大小
) {
  # 加载必要包
  require(circlize)
  require(RColorBrewer)

  # 数据预处理
  comm_data <- comm_data[1:top_n, ]
  cell_group <- unique(c(comm_data$sender, comm_data$receiver))

  # 创建颜色映射
  cell_col <- setNames(
    colorRampPalette(cell_colors)(length(cell_group)),
    cell_group
  )

  # 代谢物颜色
  metabolite_colors <- if(metabolite_palette %in% rownames(brewer.pal.info)) {
    brewer.pal(name = metabolite_palette, n = length(unique(comm_data$metabolite)))
  } else {
    colorRampPalette(c("#e74c3c", "#f1c40f"))(length(unique(comm_data$metabolite)))
  }
  names(metabolite_colors) <- unique(comm_data$metabolite)

  # 基因数据准备
  genes <- c(
    structure(comm_data$transporter, names = comm_data$sender),
    structure(comm_data$receptor, names = comm_data$receiver)
  )
  genes <- genes[!duplicated(paste(names(genes), genes))]
  genes <- genes[order(names(genes))]

  # 构建弦图数据
  df_circos <- data.frame(
    from = paste(comm_data$sender, comm_data$transporter),
    to = paste(comm_data$receiver, comm_data$receptor)
  )

  # 标准化强度参数
  strength_norm <- comm_data$strength / max(comm_data$strength)

  # 绘制弦图
  chordDiagram(
    df_circos,
    order = paste(names(genes), genes),
    transparency = 0.2,
    preAllocateTracks = list(
      list(track.height = uh(track_heights[1], "mm")),
      list(track.height = uh(track_heights[2], "mm"))
    ),
    annotationTrack = NULL,
    directional = 1,
    link.arr.type = arrow_params$type,
    link.arr.length = arrow_params$length,
    link.arr.width = arrow_params$width_factor * strength_norm,
    link.arr.col = arrow_params$col,
    link.lwd = arrow_params$width_factor * strength_norm,
    link.largest.ontop = TRUE,
    col = metabolite_colors[comm_data$metabolite],
    direction.type = 'arrows'
  )

  # 高亮扇区
  highlight_sectors <- function() {
    # 基因类型高亮
    for(g in unique(genes)) {
      cells <- names(genes)[genes == g]
      gene_type <- ifelse(g %in% comm_data$transporter, "transporter", "receptor")
      for(c in cells) {
        highlight.sector(
          sector.index = paste(c, g),
          track.index = 2,
          col = adjustcolor(type_colors[gene_type], alpha.f = 0.5),
          text = g, cex = 0.5,
          text.vjust = 0.5,
          niceFacing = TRUE,
          lwd = 1
        )
      }
    }

    # 细胞类型高亮
    for(c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(
        sector.index = paste(c, gene),
        track.index = 1,
        col = cell_col[c],
        text = c,
        text.vjust = -1,
        niceFacing = TRUE,
        lwd = 1
      )
    }
  }
  highlight_sectors()

  # 添加图例
  add_legends <- function() {
    circos.clear()
    par(mar = c(1, 6, 1, 1), new = TRUE, xpd = TRUE)
    plot.new()

    line_height <- 0.1
    start_y <- legend_position[2]

    # 细胞类型图例
    legend(
      x = legend_position[1], y = start_y,
      title = "Cell Types",
      legend = names(cell_col),
      fill = cell_col,
      border = NA, bty = "n",
      cex = legend_cex,
      title.adj = 0
    )

    # 基因类型图例
    legend(
      x = legend_position[1],
      y = start_y - length(cell_col)*line_height,
      title = "Gene Types",
      legend = names(type_colors),
      fill = type_colors,
      border = NA, bty = "n",
      cex = legend_cex,
      title.adj = 0
    )

    # 代谢物图例
    legend(
      x = legend_position[1],
      y = start_y - length(cell_col)*line_height - 3*line_height,
      title = "Metabolites",
      legend = names(metabolite_colors),
      fill = metabolite_colors,
      border = NA, bty = "n",
      cex = legend_cex,
      title.adj = 0
    )

    par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
  }
  add_legends()
}
