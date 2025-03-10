#' 加载物种特异性代谢数据
#' @param data_type 数据类型：transporter/enzyme/receptor
#' @param species 物种：human/mouse
#' @return 请求的数据框
#' @export
load_metabolic_data <- function(data_type = c("transporter", "enzyme", "receptor"),
                                species = c("human", "mouse")) {
  
  # 参数验证
  data_type <- match.arg(data_type)
  species <- match.arg(species)
  
  # 构造数据对象名
  obj_name <- sprintf("%s_%s", data_type, species)
  
  # 获取数据路径
  data_path <- system.file(package = "scMetalk")
  target_file <- file.path(data_path, "data", paste0(obj_name, ".rda"))
  
  # 检查文件存在性
  if (!file.exists(target_file)) {
    stop(sprintf("数据文件 %s 不存在于路径：%s", obj_name, data_path))
  }
  
  # 加载数据到环境
  env <- new.env()
  load(target_file, envir = env)
  
  # 提取数据对象
  if (!exists(obj_name, envir = env)) {
    stop(sprintf("文件 %s 中未找到对象 %s", target_file, obj_name))
  }
  
  get(obj_name, envir = env)
}
