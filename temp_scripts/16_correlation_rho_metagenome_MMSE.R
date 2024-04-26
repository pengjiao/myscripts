library(phyloseq)
library(microbiome)
library(correlation)
###
metagenome.meta <- microbiome::meta(SixHosp_AD_Paper_Meta_Omics_AbPET_MRI_ATN_16S_metagenome_metadata_unified[["metagenome"]][["structure_phyloseq_relative_sum_to_100"]])
metagenome.otu <- phyloseq::otu_table(SixHosp_AD_Paper_Meta_Omics_AbPET_MRI_ATN_16S_metagenome_metadata_unified[["metagenome"]][["structure_phyloseq_relative_sum_to_100"]])
metagenome.otu <- metagenome.otu|>t()|>as.data.frame()
stopifnot(rownames(metagenome.otu) == rownames(metagenome.meta))

### correlation


results <- correlation::correlation(data = metagenome.meta[metagenome.meta$Group!="AD_2",c("MMSE","Age")], data2 = metagenome.otu[which(metagenome.meta$Group!="AD_2"),], method = "spearman",  p_adjust = "none")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_metagenome_noControlling.xlsx")

results <- correlation::correlation(data = metagenome.meta[metagenome.meta$Group!="AD_2",c("MMSE","Age")], data2 = metagenome.otu[which(metagenome.meta$Group!="AD_2"),], method = "spearman",  p_adjust = "fdr")
openxlsx::write.xlsx(as.data.frame(results), file = "nc_scd_mci_ad1_MMSE_metagenome_noControlling_fdr.xlsx")


### the speed of correlation is so slow we change to cor_test
DATA <- merge(metagenome.meta[metagenome.meta$Group!="AD_2",c("MMSE","Age","Gender")],metagenome.otu[which(metagenome.meta$Group!="AD_2"),], by = 0)
DATA.male <- subset(DATA, Gender == 1)
DATA.female <- subset(DATA, Gender == 2)

DATA.male.2 <- DATA.male[,-c(1,3,4)]
DATA.female.2 <- DATA.female[,-c(1,3,4)]

results.male <- rstatix::cor_test(DATA.male.2, vars = "MMSE", method = "spearman")
results.female <- rstatix::cor_test(DATA.female.2, vars = "MMSE", method = "spearman")



############### plot
input.mmse <- results.male

#input.mmse.sub <- dplyr::filter(input.mmse, abs(rho) > 0.2 & p < 0.05)
input.mmse.sub <- dplyr::filter(input.mmse, abs(cor) > 0.1 & p < 0.05)


## add one more column for the plus or minus sign of rho

input.mmse.sub <- input.mmse.sub%>%dplyr::mutate(direction = ifelse(cor > 0, "+","-"))


###################################  plot below ####################################

## one can plot the results based on the coefficients/difference
p <- ggpubr::ggbarplot(input.mmse.sub,
                       x = "var2", y = "cor",
                       fill = "direction",
                       color = "white",
                       palette = c("#1b9e77", "#d95f02"),
                       sort.val = "asc",
                       sort.by.groups = FALSE,
                       x.text.angle = 90,
                       ylab = "rho",
                       xlab = "",
                       rotate = TRUE,
                       title = "",
                       ggtheme = theme_pubclean()
)
# +
#   geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
#     width = 0.2, size = 0.6,
#     position = position_dodge(0.05), color = "black"
#   )
p <- ggpar(p, legend = "none")
# p <- p+ylim(-1,1.5)
p <- p + geom_hline(yintercept = 0, linetype = 1, color = "black", size = 0.2)
p <- p + geom_hline(yintercept = 0.2, linetype = 2, color = "red", size = 0.2)
p <- p + geom_hline(yintercept = -0.2, linetype = 2, color = "red", size = 0.2)

p <- p + theme(axis.text.x = element_text(angle = 0, hjust= 0.5))
p
## assign each output figure
# assign(paste0(i, ".NC.barplot"), p)
ggsave(filename = "metagenome_male.pdf")





