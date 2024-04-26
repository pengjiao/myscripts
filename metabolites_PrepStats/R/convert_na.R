
convert_na <- function(this.dir, sheet.name, chose.type, chose.comment, chose.color) {
  if(is.null(this.dir)){
    stop("Error: please provide the directory which excel files existed.")
  }

  file.paths <- list.files(this.dir, all.files = T, full.names = T, recursive = T)
  file.paths <- file.paths[grepl("xlsx$", file.paths)]

  file.names <- list.files(this.dir, all.files = T, full.names = F, recursive = T)
  file.names <- file.names[grepl("xlsx$", file.names)]

  char.split <- str_split(this.dir, "/", simplify = T)
  edited.path <- paste0(paste0(char.split[-length(char.split)], collapse = "/"), "/edited_metabolite")
  dir.create(edited.path)

  for (i in 1:length(file.paths)) {
    file.path <- file.paths[i]
    file.name <- file.names[i]

    sheet.names <- excel_sheets(file.path)

    if(is.null(sheet.name)){
      print(paste0("Please provide the sheet in an excel spreadsheet: the number or the name of ", file.name))
      break
    }
    if(!is.na(as.numeric(sheet.name))){
      cells <- cells[cells$sheet == sheet.names[sheet.name],]
    }
    if(sheet.name %in% sheet.names){
      data.input <- read_excel(file.path, sheet = sheet.name)
      formats <- xlsx_formats(file.path)
      cells <- xlsx_cells(file.path, sheets = sheet.name)
      cells <- cells[cells$sheet == sheet.name,]
    }

    cell.comments <- unique(cells$comment)
    cell.colors <- unique(formats$local$fill$patternFill$fgColor$rgb)

    if(is.null(chose.type) | chose.type == "comments"){
      if(is.na(cell.comments)){
        data <- data.input
      }
      else{
      if(is.null(chose.comment)){
        print(paste0("Plsease provide the right comment which will be converted to NAs: ", file.name))
        break
        }
      if(!chose.comment %in% cell.comments){
        print(paste0("Plsease provide the right comment which will be converted to NAs: ", file.name))
      break
        }
      else{
        cc <- cells$comment %in% chose.comment
        sM <- Matrix::sparseMatrix(cells$row, cells$col, x=cc)
        M <- as.matrix(sM[1:nrow(data.input),1:ncol(data.input)])
        data <- data.input
        data[M] <- NA
      }
      }
    }
    if(chose.type == "colors"){
      if(is.na(cell.colors)){
        data <- data.input
      }
      else{
      if(is.null(chose.color)){
        print(paste0("Plsease provide the right comment which will be converted to NAs: ", file.name))
        break
      }
      if(!chose.color %in% cell.colors){
        print(paste0("Plsease provide the right comment which will be converted to NAs: ", file.name))
        break
        }
      else{
        cc <- cells$local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == chose.color)
        sM <- Matrix::sparseMatrix(cells$row, cells$col, x=cc)
        M <- as.matrix(sM)
        data <- data.input
        data[M] <- NA
      }
      }
    }
    openxlsx::write.xlsx(data, file = paste0(edited.path, "/", file.name))
  }
}
