##########################################################################################
##################################  (NC + SCD1 + aMCI + AD1 + AD2) ####################
##########################################################################################

## get meta
meta <- SixHosp_AD_319subjects_Paper_Meta_Omics_AbPET_MRI_ATN$Meta_MMSE_PET |> data.frame()

##########################################################################################
################################## Metabolites Preprocess ################################
##########################################################################################
## get the metabolites
metabolites <- SixHosp_AD_319subjects_Paper_Meta_Omics_AbPET_MRI_ATN[["Metabolite"]]

## sometimes, classical data frame might work better compared to the newer tibble format
metabolites <- cbind(sample = rownames(metabolites), metabolites)

## lcms until Choline before C0, while the rest is all fia
# lcms <- metabolites[,1:69]
# fia <- metabolites[,-c(1:69)]

## it seems the "sample" column in meta contains the appropriate IDs for all 319 subjects
stopifnot(sum(!(is.element(metabolites$sample, meta$sample))) == 0)

## since we do many calculation below, so it is better to reorder the order of samples based on the order in meta
## Saving indices for how to reorder `omics` to match `meta`
reorder_idx <- match(meta$sample, metabolites$sample)
# Reordering and saving the output
metabolites.reordered <- metabolites[reorder_idx, ]
rownames(metabolites.reordered) <- NULL
## recheck the order
stopifnot(metabolites.reordered$sample == meta$sample)
###################
## lcms until Choline before C0, while the rest is all fia
lcms <- metabolites.reordered[, 1:69]
fia <- metabolites.reordered[, -c(2:69)]
rownames(lcms) <- NULL
rownames(fia) <- NULL
##########################################################################################
################################## Immune cell Preprocess ################################
##########################################################################################
## get the immune cells
immune <- SixHosp_AD_319subjects_Paper_Meta_Omics_AbPET_MRI_ATN[["Immune"]]

## sometimes, classical data frame might work better compared to the newer tibble format
immune <- cbind(sample = rownames(immune), immune)

## it seems the "sample" column in meta contains the appropriate IDs for all 319 subjects
stopifnot(sum((is.element(immune$sample, meta$sample))) == nrow(meta))

## since we do many calculation below, so it is better to reorder the order of samples based on the order in meta
## Saving indices for how to reorder `omics` to match `meta`
reorder_idx <- match(meta$sample, immune$sample)
# Reordering and saving the output
immune.reordered <- immune[reorder_idx, ]
rownames(immune.reordered) <- NULL
## recheck the order
stopifnot(immune.reordered$sample == meta$sample)

##########################################################################################
################################## MRI (AD2 has no MRI data) #############################
##########################################################################################

## get MRI
mri <- SixHosp_AD_319subjects_Paper_Meta_Omics_AbPET_MRI_ATN[["MRI"]]
rownames(mri) <- NULL
##
colnames(mri)[1] <- "sample"
##
## it seems the "sample" column in meta contains the appropriate IDs for all MRI subjects
stopifnot(sum((is.element(mri$sample, meta$sample))) == nrow(mri))
stopifnot(sum((is.element(mri$sample, metabolites$sample))) == nrow(mri))
stopifnot(sum((is.element(mri$sample, immune$sample))) == nrow(mri))

## we need to prepare the MRI-matched immune and metabolites data frame

##################################### MRI-matched immune cells
reorder_idx2 <- match(mri$sample, immune$sample)
# Reordering and saving the output
immune.mri <- immune[reorder_idx2, ]
rownames(immune.mri) <- NULL
## recheck the order
stopifnot(immune.mri$sample == mri$sample)
##################################### MRI-matched metabolites
reorder_idx3 <- match(mri$sample, metabolites$sample)
# Reordering and saving the output
metabolites.mri <- metabolites[reorder_idx3, ]
rownames(metabolites.mri) <- NULL
## recheck the order
stopifnot(metabolites.mri$sample == mri$sample)
##
##################################### MRI-matched meta
reorder_idx4 <- match(mri$sample, meta$sample)
# Reordering and saving the output
meta.mri <- meta[reorder_idx4, ]
rownames(meta.mri) <- NULL
## recheck the order
stopifnot(meta.mri$sample == mri$sample)


##########################################################################################
################################## (NC + SCD1 + aMCI + AD1) ####################
##########################################################################################


## we need to exclude AD2
meta.noAD2 <- subset(meta, Group != "AD_2")
lcms.noAD2 <- lcms[lcms$sample %in% meta.noAD2$sample, ]
fia.noAD2 <- fia[fia$sample %in% meta.noAD2$sample, ]
immune.reordered.noAD2 <- immune.reordered[immune.reordered$sample %in% meta.noAD2$sample, ]
##
## recheck the order of samples
stopifnot(meta.noAD2$sample == immune.reordered.noAD2$sample)
stopifnot(meta.noAD2$sample == lcms.noAD2$sample)
stopifnot(meta.noAD2$sample == fia.noAD2$sample)


#########################################################################################
##################################  (SCD1 + aMCI + AD1) ####################
##########################################################################################


## we need to exclude AD2 and NC
meta.noNCAD2 <- subset(meta.noAD2, Group != "NC")
lcms.noNCAD2 <- lcms.noAD2[lcms.noAD2$sample %in% meta.noNCAD2$sample, ]
fia.noNCAD2 <- fia.noAD2[fia.noAD2$sample %in% meta.noNCAD2$sample, ]
immune.reordered.noNCAD2 <- immune.reordered.noAD2[immune.reordered.noAD2$sample %in% meta.noNCAD2$sample, ]
##
## recheck the order of samples
stopifnot(meta.noNCAD2$sample == immune.reordered.noNCAD2$sample)
stopifnot(meta.noNCAD2$sample == lcms.noNCAD2$sample)
stopifnot(meta.noNCAD2$sample == fia.noNCAD2$sample)

## we need to prepare the MRI-matched immune and metabolites data frame

## remove NC group
which(is.element(mri$sample, meta.noNCAD2$sample))

mri.noNCAD2 <- mri[which(is.element(mri$sample, meta.noNCAD2$sample)), ]
stopifnot(mri.noNCAD2$sample %in% meta.noNCAD2$sample)


##################################### MRI-matched immune cells
reorder_idx5 <- match(mri.noNCAD2$sample, immune$sample)
# Reordering and saving the output
immune.mri.noNCAD <- immune[reorder_idx5, ]
rownames(immune.mri.noNCAD) <- NULL
## recheck the order
stopifnot(immune.mri.noNCAD$sample == mri.noNCAD2$sample)
##################################### MRI-matched metabolites
reorder_idx6 <- match(mri.noNCAD2$sample, metabolites$sample)
# Reordering and saving the output
metabolites.mri.noNCAD2 <- metabolites[reorder_idx6, ]
rownames(metabolites.mri.noNCAD2) <- NULL
## recheck the order
stopifnot(metabolites.mri.noNCAD2$sample == mri.noNCAD2$sample)
##
########################################
lcms.mri.noNCAD2 <- metabolites.mri.noNCAD2[, c(1:69)]
fia.mri.noNCAD2 <- metabolites.mri.noNCAD2[, -c(2:69)]
##################################### MRI-matched meta
reorder_idx7 <- match(mri.noNCAD2$sample, meta$sample)
# Reordering and saving the output
meta.mri.noNCAD2 <- meta[reorder_idx7, ]
rownames(meta.mri.noNCAD2) <- NULL
## recheck the order
stopifnot(meta.mri.noNCAD2$sample == mri.noNCAD2$sample)
############################################################

