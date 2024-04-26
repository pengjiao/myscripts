library(tidyverse)
library(modelbased)
library(see)
library(patchwork)
library(ggpubr)
#### to plot the difference between male and female upon the association between immune cells and age
### we plot for immune cells, but one can plot the other omics as you want
### the only things you need to change are the input file name. e.g "lcms","fia", and the line 12
### mutate_each, variable selection.
## you just need to copy the whole script to a new R script file and change to "lcms"
stopifnot(immune.reordered$sample == meta$sample)
##
temp5 <- dplyr::left_join(immune.reordered, meta[, c("sample", "Group", "Gender.y", "Age.y")])
## we need to log transform each feature to make it more normal distribution.
## whenever you do lm, you need to make sure your "y" is normal like, otherwise, do some transformation on it
## whenever you do lm, you need to make sure your "y" is normal like, otherwise, do some transformation on it
temp5.log <- mutate_each(temp5, list(~ log10(. + 0.0001)), T:ncMo.cMO)
## gender as factor level
temp5.log <- temp5.log %>% mutate(across(Gender.y, as.factor))
##
DF <- NULL
for (i in names(temp5.log)[2:42]) {
  model <- lm(as.formula(paste0(i, "~ Group + Gender.y * Age.y")), data = temp5.log)
  # compute contrasts between groups
  contrasts <- estimate_contrasts(model, contrast = "Group", p_adjust = "none", standardize = TRUE)
  contrasts <- contrasts %>% mutate(feature = i)
  DF <- rbind(DF, contrasts)
}

## filter out those rows without NC
DF.nc.for <- DF %>% dplyr::filter(Level2 == "NC")
DF.nc.rev <- DF %>% dplyr::filter(Level1 == "NC")

## we need to reorder the difference
DF.nc.rev <- DF.nc.rev[, c(2, 1, 3:ncol(DF.nc.rev))]
colnames(DF.nc.rev)[1:2] <- colnames(DF.nc.for[1:2])
##
# DF.nc.rev$Difference <- DF.nc.rev$Difference * (-1)
# DF.nc.rev$CI_low <- DF.nc.rev$CI_low * (-1)
# DF.nc.rev$CI_high <- DF.nc.rev$CI_high * (-1)

DF.nc.rev[, c("Difference", "CI_low", "CI_high")] <- DF.nc.rev[, c("Difference", "CI_low", "CI_high")] * (-1)

## combine

DF.nc2 <- dplyr::bind_rows(DF.nc.for, DF.nc.rev)
## in case we need padjust value for other purposes
DF.nc3 <- DF.nc2 %>%
  group_by(Level1) %>%
  mutate(padj = p.adjust(p, method = "BH"))

## we can plot based on raw p values

feat.sub <- subset(DF.nc3, p < 0.05)
feat.sig <- unique(feat.sub$feature)

## take all the sig features

feat.sig.df <- DF.nc3[DF.nc3$feature %in% feat.sig, ] %>% arrange(feature)

##
## we also categorize p values
feat.sig.df$pval.grp <-
  ifelse(feat.sig.df$p < 0.05,
    "<0.05",
    ">0.05"
  )

##
feat.sig.df$pval.grp <- factor(feat.sig.df$pval.grp, levels = c("<0.05", ">0.05"))
## we make the plot
for (i in unique(feat.sig.df$Level1)) {
  p <- ggdotchart(subset(feat.sig.df, Level1 == i),
    x = "feature", y = "Difference",
    color = "pval.grp",
    sorting = "none",
    rotate = TRUE,
    dot.size = 3.0,
    #  add = "segment",
    title = paste0(i, "-NC"),
    y.text.col = FALSE,
    xlab = "",
    ylab = "Difference (relative to NC)",
    ggtheme = theme_pubr()
  )
  p <- ggpar(p, legend = "none")
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "blue")
  p <- p + geom_errorbar(aes(ymin = CI_low, ymax = CI_high), size = 0.2, width = 0.2)
  p <- p + scale_color_manual(values = c("red", "grey", "grey"))
  ## assign each output figure
  assign(paste0(i, "-NC_dotplot"), p)
}

## make the combined plot
p <- SCD1_graph + aMCI_graph + AD_1_graph + AD_2_graph + plot_layout(nrow = 1)
## save the plot
ggsave(filename = "immune_group_difference_adjusted_by_sex_age.pdf", width = 14, height = 8)
