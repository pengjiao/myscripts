library(tidyverse)
library(modelbased)
library(see)
library(patchwork)
#### to plot the difference between male and female upon the association between lcms cells and age
### we plot for lcms cells, but one can plot the other omics as you want
### the only thing you need to change is the input file name. e.g "lcms","fia"
## you just need to copy the whole script to a new R script file and change to "lcms"
stopifnot(lcms$sample == meta$sample)
##
temp5 <- dplyr::left_join(lcms, meta[, c("sample", "Group", "Gender.y", "Age.y")])
## we need to log transform each feature to make it more normal distribution.
## whenever you do lm, you need to make sure your "y" is normal like, otherwise, do some transformation on it
## whenever you do lm, you need to make sure your "y" is normal like, otherwise, do some transformation on it
temp5.log <- mutate_each(temp5, funs(log10(. + 0.0001)), TMAO:Xanthine)
## gender as factor level
temp5.log <- temp5.log %>% mutate(across(Gender.y, as.factor))
##
pstlist <- list()
for (i in names(temp5.log)[2:42]) {
  model <- lm(as.formula(paste0(i, "~ Gender.y * Age.y")), data = subset(temp5.log, Group != "AD_2"))
  # model <- lm(as.formula(paste0( i ,"~ Gender.y * Age.y")), data = subset(temp5.log, Group=="AD_2"))
  # model <- lm(as.formula(paste0( i ,"~ Gender.y * Age.y")), data = subset(temp5.log, Group=="AD_1"))
  # model <- lm(as.formula(paste0( i ,"~ Gender.y * Age.y")), data = subset(temp5.log, Group=="aMCI"))
  # model <- lm(as.formula(paste0( i ,"~ Gender.y * Age.y")), data = subset(temp5.log, Group=="SCD1"))
  # model <- lm(as.formula(paste0(i, "~ Gender.y * Age.y")), data = subset(temp5.log, Group == "NC"))
  # Recompute contrasts with a higher precision (for a smoother plot)
  contrasts <- estimate_contrasts(model, contrast = "Gender.y", at = "Age.y", length = 20, standardize = TRUE)
  # Add Contrast column by concatenating, here we need to check first and know the contrast.
  contrasts$Contrast <- paste(contrasts$Level1, "-", contrasts$Level2)
  contrasts$Contrast <- "Male - Female"
  # Plot
  p <- ggplot(contrasts, aes(x = Age.y, y = Difference, )) +
    # Add line and CI band
    geom_line(aes(color = Contrast), size = 2) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Contrast), alpha = 0.2) +
    # Add line at 0, indicating no difference
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("Age") +
    ylab("Difference% (men relative to women)") +
    ggplot2::theme_bw() +
    ggtitle(label = i) +
    theme(legend.position = "none")
  pstlist[[i]] <- p
}
##
wrap_plots(pstlist)
ggsave(filename = "nc_scd_amci_ad1_lcms.pdf", width = 16, height = 18)
# ggsave(filename = "nc_lcms.pdf", width = 16, height = 18)
# ggsave(filename = "scd_lcms.pdf", width = 16, height = 18)
# ggsave(filename = "amci_lcms.pdf", width = 16, height = 18)
# ggsave(filename = "ad1_lcms.pdf", width = 16, height = 18)
# ggsave(filename = "ad2_lcms.pdf", width = 16, height = 18)
