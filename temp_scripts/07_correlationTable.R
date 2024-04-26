##
library(tidyverse)
library(correlation)

### correlation table by variables from two data frames but the samples are in the same order!!!!
### you must make sure the order of samples from two data frames are the same before running the code below

results <- correlation::correlation(data = meta.noAD2[,c("MMSE","Age.y","Gender.y")], data2 = lcms.noAD2[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_LCMS_noControlling.xlsx")

results <- correlation::correlation(data = meta.noAD2[,c("MMSE","Age.y","Gender.y")], data2 = fia.noAD2[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_FIA_noControlling.xlsx")

results <- correlation::correlation(data = meta.noAD2[,c("MMSE","Age.y","Gender.y")], data2 = immune.reordered.noAD2[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_Immune_noControlling.xlsx")



#### below are the codes for making calcualtion on either male or female subjects
########################################################
########################################################
########################################################  dataset preparation
########################################################
########################################################
solna_m <- dplyr::filter(meta.noAD2[,c("sample","MMSE","Age.y","Gender.y")], Gender.y == "1")
solna_f <- dplyr::filter(meta.noAD2[,c("sample","MMSE","Age.y","Gender.y")], Gender.y == "2")

##
##################################### matched meta
reorder_idx8 <- match(solna_m$sample, lcms.noAD2$sample)
# Reordering and saving the output
lcms.noAD2.male <- lcms.noAD2[reorder_idx8, ]
rownames(lcms.noAD2.male) <- NULL
## recheck the order
stopifnot(lcms.noAD2.male$sample == solna_m$sample)
############################################################


##################################### matched meta
reorder_idx9 <- match(solna_m$sample, fia.noAD2$sample)
# Reordering and saving the output
fia.noAD2.male <- fia.noAD2[reorder_idx9, ]
rownames(fia.noAD2.male) <- NULL
## recheck the order
stopifnot(fia.noAD2.male$sample == solna_m$sample)
############################################################


##################################### matched meta
reorder_idx10 <- match(solna_m$sample, immune.reordered.noAD2$sample)
# Reordering and saving the output
immune.reordered.noAD2.male <- immune.reordered.noAD2[reorder_idx10, ]
rownames(immune.reordered.noAD2.male) <- NULL
## recheck the order
stopifnot(immune.reordered.noAD2.male$sample == solna_m$sample)
############################################################



##
##################################### matched meta
reorder_idx11 <- match(solna_f$sample, lcms.noAD2$sample)
# Reordering and saving the output
lcms.noAD2.female <- lcms.noAD2[reorder_idx11, ]
rownames(lcms.noAD2.female) <- NULL
## recheck the order
stopifnot(lcms.noAD2.female$sample == solna_f$sample)
############################################################


##################################### matched meta
reorder_idx12 <- match(solna_f$sample, fia.noAD2$sample)
# Reordering and saving the output
fia.noAD2.female <- fia.noAD2[reorder_idx12, ]
rownames(fia.noAD2.female) <- NULL
## recheck the order
stopifnot(fia.noAD2.female$sample == solna_f$sample)
############################################################


##################################### matched meta
reorder_idx13 <- match(solna_f$sample, immune.reordered.noAD2$sample)
# Reordering and saving the output
immune.reordered.noAD2.female <- immune.reordered.noAD2[reorder_idx13, ]
rownames(immune.reordered.noAD2.female) <- NULL
## recheck the order
stopifnot(immune.reordered.noAD2.female$sample == solna_f$sample)
############################################################
############################################################
############################################################  below are the calculations
############################################################
############################################################
############################################################

## only male
results <- correlation::correlation(data = solna_m[,c("MMSE","Age.y")], data2 = lcms.noAD2.male[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_LCMS_noControlling_onlyMale.xlsx")

results <- correlation::correlation(data = solna_m[,c("MMSE","Age.y")], data2 = fia.noAD2.male[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_FIA_noControlling_onlyMale.xlsx")

results <- correlation::correlation(data = solna_m[,c("MMSE","Age.y")], data2 = immune.reordered.noAD2.male[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_Immune_noControlling_onlyMale.xlsx")


## only female
results <- correlation::correlation(data = solna_f[,c("MMSE","Age.y")], data2 = lcms.noAD2.female[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_LCMS_noControlling_onlyFemale.xlsx")

results <- correlation::correlation(data = solna_f[,c("MMSE","Age.y")], data2 = fia.noAD2.female[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_FIA_noControlling_onlyFemale.xlsx")

results <- correlation::correlation(data = solna_f[,c("MMSE","Age.y")], data2 = immune.reordered.noAD2.female[,-1], method = "spearman",  p_adjust = "BH")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_Immune_noControlling_onlyFemale.xlsx")

