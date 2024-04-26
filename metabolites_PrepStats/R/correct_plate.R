
correct_plate <- function(input.list){
  data.type <- "combined.lc"
  original.df <- input.list[[data.type]]
  control.df <- original.df[grepl("control", original.df$sample, ignore.case = T),]
  sample.df <- original.df[!grepl("control", original.df$sample, ignore.case = T),]

  na.percent <- control.df %>%
    group_by(plate_ID) %>%
    summarise_if(is.numeric, funs(sum(is.na(.))/n() != 1))
  name.select <- names(na.percent[,-1])[apply(na.percent[,-1], 2, all)]
  control.filtered <- control.df[,c("plate_ID", "sample", name.select)]
  sample.filtered <- sample.df[,c("plate_ID", "sample", name.select)]

  plate.mean <- control.filtered %>%
    group_by(plate_ID) %>%
    summarise_if(is.numeric, mean, na.rm = T)
  qc.mean <- apply(plate.mean[, -1], 2, mean)
  correct.coeff <- data.frame(qc.mean/t(plate.mean[,-1]))
  colnames(correct.coeff) <- plate.mean$plate_ID

  df.final <- NULL
  for (i in colnames(correct.coeff)) {
    df.not.correct <- sample.filtered[sample.filtered$plate_ID == i,]
    df.correct <- data.frame(t(t(df.not.correct[,-c(1:2)]) * correct.coeff[,i]), check.names = F)
    df.final <- rbind(df.final, data.frame(plate_ID = df.not.correct$plate_ID,
                           sample = df.not.correct$sample,
                           df.correct,
                           check.names = F))
  }
  corrected.lc <- df.final

  data.type <- "combined.fia"
  original.df <- input.list[[data.type]]
  control.df <- original.df[grepl("control", original.df$sample, ignore.case = T),]
  sample.df <- original.df[!grepl("control", original.df$sample, ignore.case = T),]

  na.percent <- control.df %>%
    group_by(plate_ID) %>%
    summarise_if(is.numeric, funs(sum(is.na(.))/n() != 1))
  name.select <- names(na.percent[,-1])[apply(na.percent[,-1], 2, all)]
  control.filtered <- control.df[,c("plate_ID", "sample", name.select)]
  sample.filtered <- sample.df[,c("plate_ID", "sample", name.select)]

  plate.mean <- control.filtered %>%
    group_by(plate_ID) %>%
    summarise_if(is.numeric, mean, na.rm = T)
  qc.mean <- apply(plate.mean[, -1], 2, mean)
  correct.coeff <- data.frame(qc.mean/t(plate.mean[,-1]))
  colnames(correct.coeff) <- plate.mean$plate_ID

  df.final <- NULL
  for (i in colnames(correct.coeff)) {
    df.not.correct <- sample.filtered[sample.filtered$plate_ID == i,]
    df.correct <- data.frame(t(t(df.not.correct[,-c(1:2)]) * correct.coeff[,i]), check.names = F)
    df.final <- rbind(df.final, data.frame(plate_ID = df.not.correct$plate_ID,
                                           sample = df.not.correct$sample,
                                           df.correct,
                                           check.names = F))
  }
  corrected.fia <- df.final
  return(list(corrected.lc=corrected.lc, corrected.fia=corrected.fia))
}
