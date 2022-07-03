library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(Matrix)
library(matrixStats)
library(dplyr)
library(stringr)



MAINDIR = "E:/lxy_pro/003"   #"E:/others/003_scr"
setwd(MAINDIR)
source("RFunxTreePlot.R")


# ===================================================================
#        setting direcories 
# ===================================================================

# DATADIR = sprintf("test_data/devTree/%s", "20201014") # "20200707"
# # DATADIR = "E:/Users/xyliu/data003/amph/tree"
# date_tag = c("20210115", "20201014", "20200707")[1]
# DATADIR = sprintf("test_data/devTree/%s", date_tag) # "20200707"
DATADIR = 'test_data/devTree/20210115'
figdir = file.path(DATADIR, 'figs')
if(!dir.exists(figdir)){dir.create(figdir, recursive = T)}


# ===================================================================
#         loading 1) tree-structure and 2) group-expressions
# ===================================================================

message("loading edges of developmental tree...")
# tailtag = c("", "-formal", "-labeled", "-rmvE15")[2]     #"-rmvE15"
# fn_tree_struct = sprintf("%s/tree_struct%s.csv", DATADIR, tailtag)
fn_tree_struct = "test_data/devTree/20210115/tree_struct-formal.csv"
df_struct0 = read.csv(fn_tree_struct, 
                      na.strings = '', 
                      as.is = T) 
rownames(df_struct0) = df_struct0$node
str(df_struct0)
# > str(df_struct0)
# 'data.frame':	93 obs. of  8 variables:
#   $ node       : chr  "B_0" "B_1" "B_2" "G3_1" ...
# $ parent     : chr  "32cell" "32cell" "32cell" "B_1" ...
# $ label      : chr  "B_0" "B_1" "B_2" "G3_1" ...
# $ stage      : chr  "B" "B" "B" "G3" ...
# $ stage_int  : int  6 6 6 7 7 7 7 7 7 7 ...
# $ stage_name : chr  "B" "B" "B" "G3" ...
# $ node_name  : chr  "B_0" "B_1" "B_2" "G3_1" ...
# $ parent_name: chr  "32cell" "32cell" "32cell" "B_1" ...


message("loading average expressions and expresson proportions...")
df_expr = read.csv(sprintf("%s/avg_expr_all.csv", DATADIR), row.names = 1)
df_exprprop = read.csv(sprintf("%s/expr_prop_all.csv", DATADIR), row.names = 1)
Genes = rownames(df_expr)

# ============== prepare for tree plot ===============
source("RFunxTreePlot.R")

#####[1] re-order the nodes for tree-branch order
nodes_early = c(paste0("B_", c(1, 0, 2)), 
                paste0("G3_", c(1, 5, 2, 3, 4, 0, 6)),
                paste0("G4_", c(2, 1, 3, 5, 4, 0, 6)),
                paste0("G5_", c(4, 1, 8, 5, 7, 3, 6, 2, 0, 9)))
nodes_late = subset(df_struct0, ! stage %in% c("B", "G3", "G4", "G5"))$node
nodes_ordered = c(nodes_early, nodes_late)
df_struct = df_struct0[nodes_ordered, ]
head(df_struct)


# =====================[ stages ]=========================
source("RFunxTreePlot.R")

treedt = makeTreeDataStage(df_struct, stg_levels = StageNames)
treedt@data
ggtree(treedt) + geom_nodelab()  # test

##### formal code ======================
root_lb = "32cell"

p_struct = plotStageTreeLR1(treedt, anno_list = LineageAnno,
                        root_lb = root_lb,
                        size_tree = 0.6, color_tree = 'grey30',
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        filename = sprintf("%s/devTree-basic.pdf", figdir),
                        offset.text = 0.16
)


################################################################################
# labeled style
################################################################################

source("RFunxTreePlot.R")
xmax = 9
stages = c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
# xticks = c(stages, rep("", xmax - length(stages)))

p_fig1 = plotStageTreeLR1(treedt, anno_list = LineageAnnoDetailed,
                          root_lb = root_lb,
                          size_node = 5.5,
                          size_tree = 0.6, color_tree = 'black',#'grey30',
                          alpha_point = 1,
                          size_tiplab=2, width = 6.5, height = 4.5,
                          xmax = 14,
                          filename = sprintf("%s/devTree-light.pdf", figdir),
                          offset.text = 0.16,
                          default_colr = 'grey20',
                          fontsize=3.2,
)
p_fig1
p_fig11 = p_fig1 + geom_nodelab(
  aes(x = x, 
      label=str_split_fixed(label, '_', 2)[, 2]), 
  color='black', #'grey5',
  size=3) +
  geom_tiplab(
    aes(x = x - 0.12, 
        label=str_split_fixed(label, '_', 2)[, 2]), 
    color='black', #'grey5',
    size=3) +
  scale_x_discrete(limits=c(1:xmax), labels = stages) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 150, 6, 6)) +
  theme(legend.position="left")
p_fig11
ggsave(filename = sprintf("%s/devTree-light-labeled0.pdf", figdir),
       width = 6.5, height = 4.5,
       plot = p_fig11)

hlightClade(p_fig11, anno_list=LineageAnno, 
            anno_color_list=LineageColor, #list(), 
            # barsize=1, 
            extend=0.4,
            alpha = 0.1,
            # align = T, offset = 0.4,
            default_colr = 'grey90'
) -> p_fig111 # with colored shadow
ggsave(filename = sprintf("%s/devTree-light-labeled.pdf", figdir),
       width = 6.5, height = 4.5,
       plot = p_fig111)






