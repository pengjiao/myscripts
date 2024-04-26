library(tidyverse)
library(ggstatsplot)
library(ggsignif)
library(see)


temp0 <- dplyr::left_join(lcms, fia)
temp0 <- dplyr::left_join(temp0, immune.reordered)
temp0 <- dplyr::left_join(temp0, meta)

## metabolites or immune cells
## lcms, until Choline
## df <- temp0[, c(2:70, 508)]
## immune cells
## df <- temp0[, c(449:489, 508)]
## fia
# df <- temp0[,c(71:448,508)]

# converting to factor
temp0$Group <- factor(temp0$Group, levels = c("NC", "SCD1", "aMCI", "AD_1","AD_2"))

# set up the color
mycol <- see::metro_colors()[1:nlevels(temp0$Group)] |> unname()

## I just show one variable plot, but one can run a loop for whatever you like to make the plot
# creating the base plot
## save the plots into a list
pst <- list()

for (i in names(temp0)[2:489]) {
  
  p <- ggbetweenstats(temp0,
                      x = "Group",
                      y = !!i,
                      title = i,
                      pairwise.comparisons = FALSE,
                      results.subtitle = FALSE,
                      centrality.plotting = FALSE,
                      xlab = "",
                      ylab = "",
                      point.args = list(
                        position = ggplot2::position_jitterdodge(dodge.width = 0.6),
                        alpha = 1.0, size = 1.2, stroke = 0
                      ),
                      violin.args = list(width = 0.5, alpha = 0.2),
                      ggplot.component = list(scale_color_manual(values = mycol))
  )
  
  # using `pairwise_comparisons()` function to create a data frame with results
  set.seed(123)
  df <- pairwise_comparisons(temp0, Group, !!i, p.adjust.method = "BH") %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = base::c)) %>%
    dplyr::arrange(group1) %>%
    dplyr::mutate(asterisk_label = case_when(
      p.value > 0.05 ~ "ns",
      p.value < 0.05 & p.value > 0.01 ~ "*",
      p.value < 0.01 & p.value > 0.001 ~ "**",
      p.value < 0.001 & p.value > 0.0001 ~ "***",
      p.value < 0.0001 ~ "****"
    )) %>%
    dplyr::filter(asterisk_label != "ns")

  ### it is very likely that there are no siginicant comparison after fdr, then we just there is no need to 
  ## show the asterisk, therefore we do two kinds of plots - one with asterisk if any, one without it.
  
  if (nrow(df)>0) {
    
  # adding pairwise comparisons using `{ggsignif}` package
  p <- p +
    ggsignif::geom_signif(
      comparisons = df$groups,
      annotations = df$asterisk_label,
      test = NULL,
      step_increase = 0.1,
      na.rm = TRUE
    )
  p <- p+ ggplot2::theme_minimal() + theme(legend.position="none")
  p
  pst[[i]] <- p
  ## assign each output figure
  ## assign(paste0(i, ".mriplot"), p)
  ggsave(filename = paste0(i,"_boxplot.pdf"), width = 4, height = 5)
  
  }
 else {
   p 
   pst[[i]] <- p
   ## assign each output figure
   ## assign(paste0(i, ".mriplot"), p)
   ggsave(filename = paste0(i,"_boxplot.pdf"), width = 4, height = 5)
   
 }
 }
# wrap_plots(pst, nrow = 2)
# ggsave(filename = "mri_boxplot.pdf", width = 10, height = 8)