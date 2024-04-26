library(tidyverse)
library(ggpubr)
library(patchwork)

## the original dataset
temp4 <- dplyr::left_join(lcms.noAD2, fia.noAD2)
temp4 <- dplyr::left_join(temp4, immune.reordered.noAD2)
temp4 <- dplyr::left_join(temp4, meta.noAD2)

## metabolites or immune cells
## lcms
## df <- temp4[,c(2:70,508:510)]

## immune cells
df <- temp4[, c(449:489, 508:510)]

## fia
## df <- temp4[,c(71:448,508:510)]


ps.lst <- df |>
  group_split(Group) |>
  set_names(df$Group |> unique() |> sort())

### we must reorder the list to make sure the first element in the list is the reference 
ps.lst <- ps.lst[c("NC", "SCD1", "aMCI", "AD_1")]

## we remove the Group column
ps.lst2 <- map(ps.lst, ~ (.x |> dplyr::select(-Group)))



##############################################################


### save the above regression results into a data frame

for (nm in names(ps.lst2)) {
  ## 
  DF <- NULL
  ## we select the columns we need for the fit (controlling on sex and age)
  col.nm <- ps.lst2[[nm]] |>
    dplyr::select(!c("Gender.y", "Age.y")) |>
    colnames()
  ##
  for (i in col.nm) {
    ## we need to scale the data before lm, then we got standarized coefficients
    reg <- lm(as.formula(paste(i, " ~ Gender.y * Age.y")), data = as.data.frame(scale(ps.lst2[[nm]])))
    m <- reg |>
      broom::tidy() |>
      mutate(labtest = i)
    DF <- bind_rows(DF, m) |>
      dplyr::filter(term != "(Intercept)")
  }
  ## save the output
  openxlsx::write.xlsx(DF, paste0(nm, ".xlsx"))
  
  ## pretty the output for the plot, here you can subset either sex or age depending on your interest
 # stk <- subset(DF, term == "Age.y")
  stk <- subset(DF, term == "Gender.y")
  stk <- as.data.frame(stk)
  stk <- dplyr::arrange(stk, term)
  stk$padj <- p.adjust(stk$p.value, method = "BH")
  ## we categoruze p adjust values
  stk$padj.grp <-
    ifelse(stk$padj <= 0.05,
      "< 0.05",
      ifelse(stk$padj <= 0.1 & stk$padj > 0.05,
        "0.05 - 0.1",
        "> 0.1 "
      )
    )
  ## we also categorize p values
  stk$pval.grp <-
    ifelse(stk$p.value <= 0.05,
      "< 0.05",
      ifelse(stk$p.value <= 0.1 & stk$p.value > 0.05,
        "0.05 - 0.1",
        "> 0.1 "
      )
    )
## we make the plot
  p <- ggdotchart(stk,
    x = "labtest", y = "estimate",
    color = "pval.grp",
    sorting = "none", 
    rotate = TRUE, 
    dot.size = 3,
    add = "segment",
    y.text.col = FALSE,
    ## ylab = "Effect (old relative to young)",
    ylab = "Effect (women relative to men)",
    ggtheme = theme_pubr()
  )

  p <- p + geom_hline(yintercept = 0, linetype = 1, color = "red")
## assign each output figure 
  assign(paste0(nm, "_graph"), p)
}
## make the combined plot
NC_graph + SCD1_graph + aMCI_graph + AD_1_graph + plot_layout(nrow = 1)
## save the plot
## ggsave(filename = "immune_age_effect.pdf", width = 18, height = 14)
ggsave(filename = "immune_sex_effect.pdf", width = 18, height = 10)
