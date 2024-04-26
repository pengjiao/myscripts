## calculation (semi) partial correlation in this case
##  strictly speaking, one needs to understand the difference between these two apparoches.
## partial correlation
# resY <- residuals(lm(Y ~ C1 + C2))
# resX <- residuals(lm(X ~ C1 + C2))
# round(cor(resX, resY), 03)

## semi-partial correlation version 1
# resX <- residuals(lm(X ~ C1+C2))
# cor(Y, resX)

# semi-partial correlation version 2
# sum1 <- summary(lm(Y ~ X+C1+C2))
# sum2 <- summary(lm(Y ~ C1+C2))
# round(sqrt(sum1$r.squared - sum2$r.squared), 03)


## ppcor::spcor.test(meta.noAD2$MMSE, immune.reordered.noAD2$T, meta.noAD2$Gender.y, method = "spearman")
library(tidyverse)
library(ppcor)
## semi-partial correlation between MMSE and LCMS controlling on sex
## semi-partial correlation between MMSE and LCMS controlling on sex
## semi-partial correlation between MMSE and LCMS controlling on sex
## semi-partial correlation between MMSE and LCMS controlling on sex
## semi-partial correlation between MMSE and LCMS controlling on sex

dta <- c()
for (i in colnames(lcms.noAD2[, -1])) {
  i_resutl <- ppcor::pcor.test(meta.noAD2[, "MMSE"], lcms.noAD2[, i], meta.noAD2[, "Gender.y"], method = "spearman")
  dta <- rbind(dta, c(i, (unlist(i_resutl))))
}
dta <- as.data.frame(dta)
colnames(dta)[1] <- "feature"
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_lcms_cor_MMSE_adjusted_on_sex.xlsx")

##
## semi-partial correlation between MMSE and FIA controlling on sex
dta <- c()
for (i in colnames(fia.noAD2[, -1])) {
  i_resutl <- ppcor::pcor.test(meta.noAD2[, "MMSE"], fia.noAD2[, i], meta.noAD2[, "Gender.y"], method = "spearman")
  dta <- rbind(dta, c(i, (unlist(i_resutl))))
}
dta <- as.data.frame(dta)
colnames(dta)[1] <- "feature"
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_fia_cor_MMSE_adjusted_on_sex.xlsx")


##
## semi-partial correlation between MMSE and Immune controlling on sex
dta <- c()
for (i in colnames(immune.reordered.noAD2[, -1])) {
  i_resutl <- ppcor::pcor.test(meta.noAD2[, "MMSE"], immune.reordered.noAD2[, i], meta.noAD2[, "Gender.y"], method = "spearman")
  dta <- rbind(dta, c(i, (unlist(i_resutl))))
}
dta <- as.data.frame(dta)
colnames(dta)[1] <- "feature"
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_immune_cor_MMSE_adjusted_on_sex.xlsx")



## semi-partial correlation between MRI and Immune controlling on sex
## semi-partial correlation between MRI and Immune controlling on sex
## semi-partial correlation between MRI and Immune controlling on sex
## recheck the order of samples
stopifnot(mri$sample == immune.mri$sample)
stopifnot(mri$sample == metabolites.mri$sample)

### sth wrong with one of the value in eTIV, missing, need to be filled out by mean value
mri[172, "eTIV"] <- mean(mri[, "eTIV"], na.rm = TRUE)
##
dta <- c()
for (i in colnames(immune.mri[, -1])) {
  for (j in colnames(mri[, -1])) {
    i_resutl <- ppcor::pcor.test(mri[, j], immune.mri[, i], meta.mri[, "Gender.y"], method = "spearman")
    dta <- rbind(dta, c(i, j, (unlist(i_resutl))))
  }
}
dta <- as.data.frame(dta)
colnames(dta)[1:2] <- c("feature1", "feature2")
### I think we should padjust for each cycle above not here with overall p values
### but I leave it as it is
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_immune_cor_MRI_adjusted_on_sex.xlsx")



## semi-partial correlation between MRI and Metabolites controlling on sex
## semi-partial correlation between MRI and Metabolites controlling on sex
## semi-partial correlation between MRI and Metabolites controlling on sex
### sth wrong with one of the value in eTIV, missing, need to be filled out by mean value
##
##
lcms.mri <- metabolites.mri[, c(1:69)]
fia.mri <- metabolites.mri[, -c(2:69)]


dta <- c()
for (i in colnames(lcms.mri[, -1])) {
  for (j in colnames(mri[, -1])) {
    i_resutl <- ppcor::pcor.test(mri[, j], lcms.mri[, i], meta.mri[, "Gender.y"], method = "spearman")
    dta <- rbind(dta, c(i, j, (unlist(i_resutl))))
  }
}
dta <- as.data.frame(dta)
colnames(dta)[1:2] <- c("feature1", "feature2")
### I think we should padjust for each cycle above not here with overall p values
### but I leave it as it is
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_lcms_cor_MRI_adjusted_on_sex.xlsx")


dta <- c()
for (i in colnames(fia.mri[, -1])) {
  for (j in colnames(mri[, -1])) {
    i_resutl <- ppcor::pcor.test(mri[, j], fia.mri[, i], meta.mri[, "Gender.y"], method = "spearman")
    dta <- rbind(dta, c(i, j, (unlist(i_resutl))))
  }
}
dta <- as.data.frame(dta)
colnames(dta)[1:2] <- c("feature1", "feature2")
### I think we should padjust for each cycle above not here with overall p values
### but I leave it as it is
dta$padj <- stats::p.adjust(dta$p.value, method = "BH")
openxlsx::write.xlsx(dta, file = "nc_scd_amci_ad1_fia_cor_MRI_adjusted_on_sex.xlsx")
