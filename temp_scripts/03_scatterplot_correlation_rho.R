##
library(tidyverse)
library(ggpubr)

temp4 <- dplyr::left_join(lcms.noAD2, fia.noAD2)
temp4 <- dplyr::left_join(temp4, immune.reordered.noAD2)
temp4 <- dplyr::left_join(temp4, meta.noAD2)


for (i in names(temp4)[2:489]) {
  b <- ggplot(
    data = temp4 %>% dplyr::mutate(Sex = ifelse(Gender.y == "1", "Male", "Female")),
    aes_string(y = i, x = "MMSE")
  )
  c <- b + geom_point(aes(color = as.factor(Sex)), size = 1, shape = 1, alpha = 0.3) +
    geom_rug(aes(color = as.factor(Sex))) +
    geom_smooth(aes(color = as.factor(Sex)),
      method = lm,
      se = TRUE, fullrange = TRUE
    ) +
    scale_color_manual(values = c("red", "blue")) +
    ggpubr::stat_cor(aes(color = as.factor(Sex)), label.x = 3, method = "spearman", cor.coef.name = "rho")

  c + guides(color = guide_legend(title = "Sex")) +
    theme_bw()

  ggsave(filename = paste0(i, "_cor_w_MMSE_by_sex.pdf"), width = 4, height = 4)
}
