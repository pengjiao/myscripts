library(tidyverse)
library(ggstatsplot)

temp4 <- dplyr::left_join(lcms.noAD2, fia.noAD2)
temp4 <- dplyr::left_join(temp4, immune.reordered.noAD2)
temp4 <- dplyr::left_join(temp4, meta.noAD2)


for (i in names(temp4)[2:489]) {
  grouped_ggscatterstats(
    ## we do not want dummy number to represent sex in the figure, so we generate a new variable Sex
    data = temp4 %>% dplyr::mutate(Sex = ifelse(Gender.y == "1", "Male", "Female")),
    x = MMSE,
    y = !!i,
    type = "nonparametric",
    grouping.var = Sex,
    marginal = FALSE,
    xlab = "MMSE",
    ggtheme = ggplot2::theme_minimal()
  )
  ggsave(filename = paste0(i, "_cor_w_MMSE2_by_Sex.pdf"), width = 11, height = 4.5)
}
