library(ropls)
library(dplyr)
library(ggpubr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(ggmap)
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

##
## OPLS-DA
opls.res <- ropls::opls(
  x = input[,-1],
  y = as.factor(input[,"Group"]),
  predI = 1,
  orthoI = 1,
  algoC = c("default", "nipals", "svd")[2],
  crossvalI = 7,
  log10L = FALSE,
  permI = 100,
  scaleC = c("none", "center", "pareto", "standard")[4],
  subset = NULL,
  plotSubC = NA,
  fig.pdfC = c("none", "interactive", "myfile.pdf")[3],
  info.txtC = c("none", "interactive", "myfile.txt")[3]
)
## get predictive VIP
oplsda.vip <- ropls::getVipVn(opls.res)
## data frame
vip.df <- data.frame( feature = names(oplsda.vip), score = unname(oplsda.vip))
vip.df.top <- dplyr::arrange(vip.df, desc(score)) %>% top_n(20)

ggpubr::ggdotchart(vip.df.top, x = "feature", y = "score",
           sorting = "descending",                        
           add = "segments", 
           color = "red",
           xlab = "",
           ylab = "VIP score",
           rotate = TRUE,
           ggtheme = theme_pubr()                        
)
ggsave(filename = "oplsda_model_VIP_plots__inhouse_NC_AD1.pdf")

## get the score p1 and o1 for plot
opls.score <- data.frame(p1 = getScoreMN(opls.res), o1= getScoreMN(opls.res, orthoL = TRUE))
## merge by sampleID to get group info
opls.score2 <- merge(opls.score, df2, by = 0)

################# plot

p1 <- ggplot2::ggplot(opls.score2, aes(x = p1, y = o1, color = Group)) +
  geom_point() +
  stat_ellipse() + 
  xlab("predictive component (t1)") +
  ylab("orthogonal component (to1)") +
  ggplot2::geom_vline(xintercept = 0, size = .5) +
  ggplot2::geom_hline(yintercept = 0, size = .5) +
  ggplot2::theme_light()

d <- opls.res

p1 + labs(caption = paste0(
  "R2X    R2Y    Q2    pR2Y    pQ2\n",
  d@summaryDF$`R2X(cum)`, "    ",
  d@summaryDF$`R2Y(cum)`, "    ",
  d@summaryDF$`Q2(cum)`, "    ",
  d@summaryDF$pR2Y, "    ",
  d@summaryDF$pQ2
))

ggsave(filename = "oplsda_model_score_plots_inhouse_NC_AD1.pdf")


## get the loading score p1 and o1 for plot
opls.loading <- data.frame(p1 = getLoadingMN(opls.res), o1= getLoadingMN(opls.res, orthoL = TRUE))

################# plot

p2 <- ggplot2::ggplot(opls.loading, aes(x = p1, y =o1,label = rownames(opls.loading))) +
  geom_point() +
  ggrepel::geom_text_repel() +
  xlab("predictive component loading (t1)") +
  ylab("orthogonal component loading (to1)") +
  ggplot2::geom_vline(xintercept = 0, size = .5) +
  ggplot2::geom_hline(yintercept = 0, size = .5) +
  ggplot2::theme_light()
 
d <- opls.res

p2 + labs(caption = paste0(
  "R2X    R2Y    Q2    pR2Y    pQ2\n",
  d@summaryDF$`R2X(cum)`, "    ",
  d@summaryDF$`R2Y(cum)`, "    ",
  d@summaryDF$`Q2(cum)`, "    ",
  d@summaryDF$pR2Y, "    ",
  d@summaryDF$pQ2
))
ggsave(filename = "oplsda_model_loading_plots_inhouse_NC_AD1.pdf")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
############################## the plot below adapted from others ########################
#########################################################################################
############################## run the function code first then run on our model #################################
#########################################################################################
#########################################################################################
#########################################################################################

### https://github.com/Aariq/chemhelper
### https://raw.githubusercontent.com/zzzsssyyy1995/plot_oplsda/master/plot_oplsda.R
### https://github.com/LoosC/systemsseRology
### https://github.com/JeffreyMolendijk/ml-viewR
library(tidyverse)
library(ggmap)
plot.opls <- function(x,y,label=F){
  solna <- as.data.frame(cbind(getScoreMN(x),getScoreMN(x,orthoL=T)[,1]))
  colnames(solna) <- c("t1","ot1")
  solna %<>% rownames_to_column("Label") %>%
    dplyr::mutate(Groups=factor(y))
  p1 <- ggplot(solna,aes(x=t1,y=ot1,label=Label,color=Groups,fill=Groups))+
    stat_ellipse(show.legend = F,geom="polygon",segments = 51,
                 alpha=0.2,size=0.8,level=0.95,type="norm")+
    geom_point(size=3,position = position_jitter(height = 0.5,width = 1,seed = 2))+
    geom_hline(yintercept = 0,linetype=2)+
    geom_vline(xintercept = 0,linetype=2)+
    # scale_color_manual(values = c("#2b3d4f","#18bc9c"))+
    # scale_fill_manual(values = c("#2b3d4f","#18bc9c"))+
    labs(caption=paste0("R2X    R2Y    Q2    pR2Y    pQ2\n",
                        x@summaryDF$`R2X(cum)`,"    ",
                        x@summaryDF$`R2Y(cum)`,"    ",
                        x@summaryDF$`Q2(cum)`,"    ",
                        x@summaryDF$pR2Y,"    ",
                        x@summaryDF$pQ2))+
    xlab(paste0("T score[1] (",x@modelDF$R2X[1]*100,"%)"))+
    ylab("Orthogonal T score[1]")+
    guides(color=guide_legend(override.aes=list(size=3)))+
    theme_minimal()+
    theme(axis.text=element_text(size=12,face = "bold",colour = "black"),
          axis.title=element_text(size=12,face = "bold",colour = "black"),
          axis.ticks.length=unit(0, "cm"),
          plot.caption=element_text(size=12,face = "bold",colour = "black"),
          legend.title = element_blank(),
          legend.text=element_text(size=12,face = "bold",colour = "black"),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold",colour = "black"),
          plot.margin = unit(rep(1,4),"lines"),
    )
  if(label==T){
    p1<- p1+geom_text(show.legend = F)
  }
  library(aplot)
  library(ggpubr)
  xdensity <- ggplot(solna, aes(t1, fill=Groups)) +
    geom_density(alpha=.5,adjust=2) +
    # scale_fill_manual(values = c("#2b3d4f","#18bc9c"))+
    theme(legend.position = "none")+
    theme_nothing()
  xdensity
  # Marginal density plot of y (right panel)
  ydensity <- ggplot(solna, aes(ot1, fill=Groups)) +
    geom_density(alpha=.5,adjust=2) +
    # scale_fill_manual(values = c("#2b3d4f","#18bc9c"))+
    theme(legend.position = "none")+
    coord_flip()+
    ggmap::theme_nothing()
  ydensity
  
  p1 %<>%
    insert_top(xdensity,height = 0.2) %>%
    insert_right(ydensity,width = 0.2)
  p1
}
#################################################### run the code above first then run the code below
plot.opls(opls.res, y = as.factor(input[,"Group"]))
ggsave(filename = "oplsda_model_score_plots_adapated_1_NC_AD1.pdf")

#################  
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
############################## the second code adapted from others ########################
#########################################################################################
############################## run the function code first then run on our model #################################
#########################################################################################
#########################################################################################
#########################################################################################
### https://github.com/Aariq/chemhelper
### https://raw.githubusercontent.com/zzzsssyyy1995/plot_oplsda/master/plot_oplsda.R
### https://github.com/LoosC/systemsseRology
### https://github.com/JeffreyMolendijk/ml-viewR

###
## source this file if you would like to plot score and loading plots from ropls
## adapted from https://github.com/JeffreyMolendijk/ml-viewR/blob/master/R/view.R

gg_circle <- function(rx, ry, xc, yc, color = "black", fill = NA, ...) {
  x <- xc + rx * cos(seq(0, pi, length.out = 100))
  ymax <- yc + ry * sin(seq(0, pi, length.out = 100))
  ymin <- yc + ry * sin(seq(0, -pi, length.out = 100))
  annotate("ribbon", x = x, ymin = ymin, ymax = ymax, color = color, fill = fill, size = 0.5, linetype = 2, ...)
}

ropls.plot <- function(d, plottype = "score", xvar, yvar, hotelling = FALSE, ellipse = FALSE, col.var = NULL, col.pca = NULL) {
  N <- nrow(d@scoreMN)
  sm <- d@modelDF
  y <- d@suppLs$yMCN
  hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)
  
  # Define scores
  if (length(d@orthoScoreMN) == length(d@scoreMN)) {
    ## in case of OPLS-DA
    score <- data.frame(d@scoreMN, d@orthoScoreMN)
  } else {
    ## PCA, PLS-DA
    score <- data.frame(d@scoreMN)
  }
  
  # Define loadings
  if (length(d@orthoLoadingMN) == length(d@loadingMN)) {
    loading <- data.frame(d@loadingMN, d@orthoLoadingMN)
  } else {
    loading <- data.frame(d@loadingMN)
  }
  
  # plotting scores
  if (plottype == "score") {
    if (d@typeC == "PCA") {
      p <- ggplot(score, aes_string(x = xvar, y = yvar, col = col.pca)) +
        geom_point(size = 2)
    } else {
      p <- ggplot(score, aes_string(x = xvar, y = yvar, col = "y")) +
        geom_point(size = 2)
    }
    
    p <- p + xlab(paste(xvar, " (", sm[xvar, "R2X"] * 100, "%)")) + ylab(paste(yvar, " (", sm[yvar, "R2X"] * 100, "%)"))
    
    p <- p + labs(caption = paste0(
      "R2X    R2Y    Q2    pR2Y    pQ2\n",
      d@summaryDF$`R2X(cum)`, "    ",
      d@summaryDF$`R2Y(cum)`, "    ",
      d@summaryDF$`Q2(cum)`, "    ",
      d@summaryDF$pR2Y, "    ",
      d@summaryDF$pQ2
    ))
    p <- p + geom_hline(yintercept = 0, linetype = 1, color = "black", size = 0.1)
    p <- p + geom_vline(xintercept = 0, linetype = 1, color = "black", size = 0.1)
    p <- p + ggplot2::theme_minimal() 
    
    if (hotelling) {
      p <- p + gg_circle(
        rx = sqrt(as.numeric(stats::var(score %>% dplyr::select(xvar))) * hotFisN),
        ry = sqrt(as.numeric(stats::var(score %>% dplyr::select(yvar))) * hotFisN),
        xc = 0, yc = 0
      )
    }
    
    if (ellipse) {
      p <- p + stat_ellipse(
        geom = "polygon", alpha = 0.3, linetype = "blank",
        aes_string(fill = "y"), type = "norm"
      )
    }
  }
  
  # plotting loadings
  if (plottype == "loading") {
    p <- ggplot(loading, aes_string(x = xvar, y = yvar, col = col.var)) +
      geom_point(size = 2)
    p <- p + xlab(paste(xvar, " (", sm[xvar, "R2X"] * 100, "%)")) + ylab(paste(yvar, " (", sm[yvar, "R2X"] * 100, "%)"))
    p <- p + geom_hline(yintercept = 0, color = "gray") + geom_vline(xintercept = 0, color = "gray")
    p <- p + ggplot2::theme_minimal() 
  }
  
  return(p)
}
################################### run the code above and then run on the code below
#########################################################################################
#########################################################################################
########################
#Load dataset
data.x <- input[,-1]
data.y <- as.factor(input[,"Group"])
#Make models
ropls.pca = ropls::opls(x = data.x,info.txtC = "none", fig.pdfC = FALSE)
ropls.plsda = ropls::opls(x = data.x, y = data.y, info.txtC = "none", fig.pdfC = FALSE)
ropls.oplsda = ropls::opls(x = data.x, y = data.y, predI = 1, orthoI = 1, info.txtC = "none", fig.pdfC = FALSE)



### Visualizing results

ropls.plot(ropls.pca, xvar = "p1", yvar = "p2", hotelling = TRUE)
ropls.plot(ropls.pca, xvar = "p1", yvar = "p2", hotelling = TRUE, ellipse = FALSE, col.pca = data.y)
ropls.plot(ropls.pca, xvar = "p1", yvar = "p2", hotelling = TRUE, ellipse = TRUE, col.pca = data.y)
ropls.plot(ropls.pca, plottype = "loading", xvar = "p1", yvar = "p2", hotelling = TRUE, ellipse = FALSE)
ropls.plot(ropls.opls, xvar = "p1", yvar = "o1", hotelling = TRUE)
ropls.plot(ropls.opls, xvar = "p1", yvar = "o1", hotelling = TRUE, plottype = "loading")
ropls.plot(ropls.plsda, xvar = "p1", yvar = "p2", hotelling = TRUE)
ropls.plot(ropls.plsda, xvar = "p1", yvar = "p2", hotelling = TRUE, ellipse = TRUE) 
ropls.plot(ropls.oplsda, xvar = "p1", yvar = "o1", hotelling = TRUE)
ropls.plot(ropls.oplsda, xvar = "p1", yvar = "o1", hotelling = TRUE, ellipse = TRUE) 


