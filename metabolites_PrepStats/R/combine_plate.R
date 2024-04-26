
combine_plate <- function(this.dir){

  if(is.null(this.dir)){
    stop("Error: please provide the directory which excel files existed.")
  }

  file.names <- list.files(this.dir, full.names = T)
  lc.files <- file.names[grepl("LC", file.names)]
  fia.files <- file.names[grepl("FIA", file.names)]

  lc.input <- lapply(lc.files, function(x) data.frame(read_excel(x, col_types = "text", col_names = F)))

  df.list <- list()
  for (i in 1:length(lc.input)) {
    df <- data.frame(lc.input[[i]])
    sample.ident <- which(apply(df, 2, function(x) any(grepl("Quant 500_Cal1|Sample Identification", x))))
    bio.id <- which(apply(df, 1, function(x) any(grepl("Bio ID", x))))
    name.id <- which(apply(df, 1, function(x) any(grepl("Trigonelline", x))))

    if(length(bio.id) != 0 & length(bio.id) != 0){
      info.df <- df[c(bio.id, name.id),]
    }
    if(length(name.id) != 0){
      col.select <- which(grepl("Trigonelline", df[min(name.id),]))
      df.select <- data.frame(sample=as.character(df[,sample.ident]), df[,col.select:ncol(df)])
      df.reshape <- df.select[!is.na(df.select$sample),]

      df.preped <- df.reshape
      df.list[[i]] <- data.frame(plate_ID=paste0("plate", i), df.preped)
      colnames(df.list[[i]])[-c(1:2)] <- as.character(df.select[min(name.id),])[-1]
    }
    else{
      if(length(bio.id)!=0){
        col.select <- grepl("HMDB", df[bio.id,])
        df.select <- data.frame(sample=as.character(df[,sample.ident]), df[,col.select])
        df.reshape <- df.select[!is.na(df.select$sample),]

        colname.id <- min(which(apply(df.select, 1, function(x) any(grepl("HMDB", x)))))
        df.preped <- df.reshape
        df.list[[i]] <- data.frame(plate_ID=paste0("plate", i), df.preped)
        colnames(df.list[[i]])[-c(1:2)] <- as.character(df.select[colname.id,])[-1]
      }
    }
  }
  df.combined <- do.call(rbind, df.list)
  df.filtered <- df.combined[!grepl("Quant 500_|^PBS$|Sample Identification|\\.\\.[0-9]", df.combined$sample),]
  df.filtered[,-c(1:2)] <- apply(df.filtered[,-c(1:2)], 2, as.numeric)
  combined.lc <- df.filtered

  fia.input <- lapply(fia.files, function(x) read_excel(x, col_types = "text", col_names = F))

  df.list <- list()
  for (i in 1:length(fia.input)) {
    df <- data.frame(fia.input[[i]])
    sample.ident <- which(apply(df, 2, function(x) any(grepl("Quant 500_Cal1|Sample Identification", x))))
    bio.id <- which(apply(df, 1, function(x) any(grepl("Bio ID", x))))
    name.id <- which(apply(df, 1, function(x) any(grepl("^C0$", x))))

    if(length(bio.id) != 0 & length(bio.id) != 0){
      info.df <- df[c(bio.id, name.id),]
    }
    if(length(name.id) != 0){
      col.select <- which(grepl("^C0$", df[min(name.id),]))
      df.select <- data.frame(sample=as.character(df[,sample.ident]), df[,col.select:ncol(df)])
      df.reshape <- df.select[!is.na(df.select$sample),]

      df.preped <- df.reshape
      df.list[[i]] <- data.frame(plate_ID=paste0("plate", i), df.preped)
      colnames(df.list[[i]])[-c(1:2)] <- as.character(df.select[min(name.id),])[-1]
    }
    else{
      if(length(bio.id)!=0){
        col.select <- grepl("HMDB", df[bio.id,])
        df.select <- data.frame(sample=df[,sample.ident], df[,col.select])
        df.reshape <- df.select[!is.na(df.select$sample),]

        colname.id <- min(which(apply(df.select, 1, function(x) any(grepl("HMDB", x)))))
        df.preped <- df.reshape
        df.list[[i]] <- data.frame(plate_ID=paste0("plate", i), df.preped)
        colnames(df.list[[i]])[-c(1:2)] <- as.character(df.select[colname.id,])[-1]
      }
    }
  }
  df.combined <- do.call(rbind, df.list)
  df.filtered <- df.combined[!grepl("Quant 500_|^PBS$|Sample Identification|\\.\\.[0-9]", df.combined$sample),]
  df.filtered[,-c(1:2)] <- apply(df.filtered[,-c(1:2)], 2, as.numeric)
  combined.fia <- df.filtered
  return(list(combined.lc=combined.lc, combined.fia=combined.fia))
  }
