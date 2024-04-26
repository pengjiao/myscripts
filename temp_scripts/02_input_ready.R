## subjects excluding NC and AD2
temp <- dplyr::left_join(lcms.noNCAD2, fia.noNCAD2)
temp <- dplyr::left_join(temp, immune.reordered.noNCAD2)
temp <- dplyr::left_join(temp, meta.noNCAD2)

## MRI subjects excluding NC and AD2 (AD2 has no MRI data)
temp2 <- dplyr::left_join(lcms.mri.noNCAD2, fia.mri.noNCAD2)
temp2 <- dplyr::left_join(temp2, immune.mri.noNCAD)
temp2 <- dplyr::left_join(temp2, meta.mri.noNCAD2)
temp2 <- dplyr::left_join(temp2, mri.noNCAD2)

## MRI subjects
temp3 <- dplyr::left_join(lcms.mri, fia.mri)
temp3 <- dplyr::left_join(temp3, immune.mri)
temp3 <- dplyr::left_join(temp3, meta.mri)
temp3 <- dplyr::left_join(temp3, mri)

## subjects excluding AD2
temp4 <- dplyr::left_join(lcms.noAD2, fia.noAD2)
temp4 <- dplyr::left_join(temp4, immune.reordered.noAD2)
temp4 <- dplyr::left_join(temp4, meta.noAD2)

save(lcms, lcms.mri, lcms.mri.noNCAD2, lcms.noAD2, lcms.noAD2.female, lcms.noAD2.male, lcms.noNCAD2,
     fia, fia.mri, fia.mri.noNCAD2, fia.noAD2,fia.noAD2.female, fia.noAD2.male, fia.noNCAD2,
     immune, immune.mri,immune.mri.noNCAD,immune.reordered, immune.reordered.noAD2, immune.reordered.noAD2.female,
     immune.reordered.noAD2.male, immune.reordered.noNCAD2, meta, meta.mri, meta.mri.noNCAD2, meta.noAD2, meta.noNCAD2,
     metabolites, metabolites.mri, metabolites.mri.noNCAD2, metabolites.reordered, file = "AD_paper_319subjects_20221019.RData")


