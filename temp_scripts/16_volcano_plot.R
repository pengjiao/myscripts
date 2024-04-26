################################ load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(correlation)
  library(EnhancedVolcano)
})

################################ calculate the correlation

### correlation table by variables from two data frames but the samples are in the same order!!!!
### you must make sure the order of samples from two data frames are the same before running the code below

### padj is usually used but we do not do this in this case for plotting the violin plot


results.cor.lcms <- correlation::correlation(
  data = meta.noAD2[, c("MMSE", "Gender.y")],
  data2 = lcms.noAD2[, -1],
  method = "spearman",
  p_adjust = "none"
)
openxlsx::write.xlsx(as.data.frame(results.cor.lcms), file = "nc_scd_mci_ad1_MMSE_LCMS_spearman.xlsx")

results.cor.fia <- correlation::correlation(
  data = meta.noAD2[, c("MMSE", "Gender.y")],
  data2 = fia.noAD2[, -1],
  method = "spearman",
  p_adjust = "none"
)
openxlsx::write.xlsx(as.data.frame(results.cor.fia), file = "nc_scd_mci_ad1_MMSE_FIA_spearman.xlsx")

results.cor.immune <- correlation::correlation(
  data = meta.noAD2[, c("MMSE", "Gender.y")],
  data2 = immune.reordered.noAD2[, -1],
  method = "spearman",
  p_adjust = "none"
)
openxlsx::write.xlsx(as.data.frame(results.cor.immune), file = "nc_scd_mci_ad1_MMSE_Immune_spearman.xlsx")


################################ Volcano (we only take MMSE correlation)
input <- read_excel("nc_scd_mci_ad1_MMSE_LCMS_spearman.xlsx") |> subset(Parameter1 == "MMSE")
res <- input
EnhancedVolcano(res,
  lab = res$Parameter2,
  x = "rho",
  y = "p",
  xlim = c(min(res$rho, na.rm = TRUE) - 0.1, max(res$rho, na.rm = TRUE) + 0.1),
  ylim = c(0, max(-log10(res$p), na.rm = TRUE) + 2),
  xlab = "Spearman's rho",
  title = "",
  caption = "",
  max.overlaps = Inf,
  legendPosition = "none",
  subtitle = "",
  col = c("grey", "forestgreen", "royalblue", "red2"),
  colAlpha = 1.0,
  shape = 16,
  pCutoff = 0.05,
  FCcutoff = 0.20,
  pointSize = 2.0,
  labSize = 5.0,
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = "black",
  gridlines.minor = TRUE,
  gridlines.major = FALSE,
  border = "partial"
)

ggsave("nc_scd_mci_ad1_MMSE_LCMS_spearman_volcanoplot.pdf", width = 5, height = 5)
