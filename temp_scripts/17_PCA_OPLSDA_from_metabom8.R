library(metabom8)
library(tidyverse)
## merge 
df <- dplyr::left_join(meta[,c("sample","Group")], lcms)
## we need to get sampleID to the row names of the data frame, it is needed for merging below
df2 <- tibble::column_to_rownames(df, var = "sample")
## 
input <- subset(df2, Group == "NC"| Group == "AD_1")
### you need to transform your matrix before any PCA or (O)PLSDA analysis
### you need to transform your matrix before any PCA or (O)PLSDA analysis
### you need to transform your matrix before any PCA or (O)PLSDA analysis

input[,-1] <- log2(input[,-1] + 0.0001)


### define X and Y

X <- input[, -1]
Y <- as.factor(input[, "Group"])

# Perform PCA
pca_model <- metabom8::pca(X = X, pc = 2, scale = "UV", center = TRUE)

# Plot PCA results: scores of the first two components
#


#######################################  PCA score plot ########################
p1 <- plotscores(
  obj = pca_model, pc = c(1, 2),
  an = list(Group = Y),
  title = "PCA - Scores plot"
)
p1 <- p1 + see::scale_color_metro_d() +
  ggplot2::theme(legend.position = "right")

p1
ggsave(filename = "pca_score_plot_NC_AD1.pdf")

#######################################  OPLS-DA ########################

# Train O-PLS model
opls.model=opls(X,Y, 
                center = TRUE,
                scale = "UV", 
                plotting = TRUE)

opls.model@summary
ggsave(filename = "oplsda_model_summary_plots_NC_AD1.pdf")

####################################### outlier detection ################
distX=dmodx(mod =opls.model, plot=TRUE)

ggsave(filename = "oplsda_distance_plot_checking_outliers_NC_AD1.pdf")

########################## OPLS-DA score plot ##########################

# Plot OPLS scores
p2 <- plotscores(obj=opls.model, an=list(
  Group=Y),         
  title='OPLS - Scores plot',  
  cv.scores = T)                 
p2 <- p2 + see::scale_color_metro_d() +
  ggplot2::stat_ellipse(geom = "polygon",level = 0.95, alpha = .1)+
  ggplot2::theme(legend.position = "right")

p2
ggsave(filename = "oplsda_score_plot_NC_AD1.pdf")


p3 <- opls_perm(opls.model, n = 200, plot = TRUE, mc = TRUE)

ggsave(filename = "oplsda_permutation_NC_AD1.pdf")



