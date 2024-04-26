library(tidyverse)
library(NetCoMi)

## original input
temp4 <- dplyr::left_join(lcms.noAD2, fia.noAD2)
temp4 <- dplyr::left_join(temp4, immune.reordered.noAD2)
temp4 <- dplyr::left_join(temp4, meta.noAD2)


## metabolites or immune cells
## lcms
df <- temp4[, c(2:70, 508)]
## immune cells
## df <- temp4[,c(449:489,508)]

## fia
## df <- temp4[,c(71:448,508)]


ps.lst <- df %>%
  group_split(Group) %>%
  set_names(df$Group %>% unique() %>% sort())

### we must reorder the list to make sure the first element in the list is the reference
ps.lst <- ps.lst[c("NC", "SCD1", "aMCI", "AD_1")]

## we remove the Group column
ps.lst2 <- map(ps.lst, ~ (.x %>% dplyr::select(-Group)))

## we do correlation
ps.lst.cor <- map(ps.lst2, ~ (.x %>% cor()))


################################################################
################################################################
for (m in names(ps.lst.cor)) {
  if (m == names(ps.lst.cor)[1]) {
    net_single <- netConstruct(
      data = ps.lst.cor[[m]],
      seed = 1234,
      verbose = 3,
      sparsMethod = "t-test",
      sampleSize = nrow(ps.lst2[[m]]),
      adjust = "BH",
      dataType = "correlation"
    )

    ## analyze the network and detect the communities (not necessarily) within the network
    props_single <- netAnalyze(net_single,
      centrLCC = TRUE,
      clustMethod = "cluster_louvain"
    )

    ## summary
    ## summary(props_single, numbNodes = 5L)

    ########## plot the network
    ## separate layout
    pdf(file = paste0("network_single_layout_grp_", m, ".pdf"), height = 5, width = 5)
    set.seed(123456)
    p1 <- plot(props_single,
      sameLayout = FALSE,
      rmSingles = "TRUE",
      labelScale = TRUE,
      nodeColor = "cluster",
      nodeSize = "fix",
      highlightHubs = TRUE
    )
    dev.off()
    #### we need all the other two groups below show the same layout as this one
    my_layout <- p1$layout$layout1
  } else {
    #############################
    ###############################################################################
    ## Network construction and analysis of the result groups with layout of group1
    #############################################################################
    ################################################################

    ##
    ## Network construction and analysis
    net_single2 <- netConstruct(
      data = ps.lst.cor[[m]],
      seed = 1234,
      verbose = 3,
      sparsMethod = "t-test",
      sampleSize = nrow(ps.lst2[[m]]),
      adjust = "BH",
      dataType = "correlation"
    )

    ## analyze the network and detect the communities (not necessarily) within the network
    props_single2 <- netAnalyze(net_single2,
      centrLCC = TRUE,
      clustMethod = "cluster_louvain"
    )

    ## summary
    summary(props_single2, numbNodes = 5L)

    ########## plot the network
    ## separate layout
    pdf(file = paste0("network_single_layout_grp_", m, ".pdf"), height = 5, width = 5)
    set.seed(123456)
    p2 <- plot(props_single2,
      layout = my_layout,
      labelScale = TRUE,
      nodeColor = "cluster",
      nodeSize = "fix",
      highlightHubs = TRUE
    )
    dev.off()
  }
}
