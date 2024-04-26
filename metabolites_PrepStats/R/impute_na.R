

impute_na <- function(input.list, na.cutoff, impute.method, k.num){

  if(is.null(na.cutoff)){
    na.cutoff <- 0.3
  }
  else{
    if(na.cutoff > 1 | na.cutoff < 0 | !is.numeric(na.cutoff)){
      stop("Please provide the correct NAs percentage cutoff: between 0 and 1.")
    }
  }

  if(is.null(impute.method)){
    impute.method <- "knn"
  }
  else{
    if(!impute.method %in% c("knn", "mean", "median")){
      stop("Please provide the three imputation methods: knn, mean and median.")
    }
  }

  data.type <- "corrected.lc"
  correct.df <- input.list[[data.type]]

  na.percent <- colSums(is.na(correct.df))/nrow(correct.df)
  filter.df <- correct.df[,na.percent <= na.cutoff]

  if(impute.method == "knn"){
    impute.log <- VIM::kNN(log2(filter.df[,-c(1:2)]), k = k.num, imp_var = F)
    impute.df <- 2^impute.log
  }

  if(impute.method == "mean"){
    mean.col <- apply(filter.df[,-c(1:2)], 2, mean, na.rm = T)
    impute.na <- filter.df[,-c(1:2)]
    for (i in 1:ncol(impute.na)) {
      impute.na[,i][is.na(impute.na[,i])] <- mean.col[i]
    }
    impute.df <- impute.na
  }

  if(impute.method == "median"){
    median.col <- apply(filter.df[,-c(1:2)], 2, median, na.rm = T)
    impute.na <- filter.df[,-c(1:2)]
    for (i in 1:ncol(impute.na)) {
      impute.na[,i][is.na(impute.na[,i])] <- median.col[i]
    }
    impute.df <- impute.na
  }

  impute.lc <- cbind(filter.df[,1:2], impute.df)

  data.type <- "corrected.fia"
  correct.df <- input.list[[data.type]]

  na.percent <- colSums(is.na(correct.df))/nrow(correct.df)
  filter.df <- correct.df[,na.percent <= na.cutoff]

  if(impute.method == "knn"){
    impute.log <- VIM::kNN(log2(filter.df[,-c(1:2)]), k = k.num, imp_var = F)
    impute.df <- 2^impute.log
  }

  if(impute.method == "mean"){
    mean.col <- apply(filter.df[,-c(1:2)], 2, mean, na.rm = T)
    impute.na <- filter.df[,-c(1:2)]
    for (i in 1:ncol(impute.na)) {
      impute.na[,i][is.na(impute.na[,i])] <- mean.col[i]
    }
    impute.df <- impute.na
  }

  if(impute.method == "median"){
    median.col <- apply(filter.df[,-c(1:2)], 2, median, na.rm = T)
    impute.na <- filter.df[,-c(1:2)]
    for (i in 1:ncol(impute.na)) {
      impute.na[,i][is.na(impute.na[,i])] <- median.col[i]
    }
    impute.df <- impute.na
  }

  impute.fia <- cbind(filter.df[,1:2], impute.df)

  return(list(imputed.lc = impute.lc, imputed.fia = impute.fia))

}
