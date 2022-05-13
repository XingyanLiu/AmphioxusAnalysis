# plot trees for the webserver
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

MAINDIR = "E:/lxy_pro/003"
setwd(MAINDIR)
source("RFunxTreePlot.R")

# ===================================================================
#        setting directories
# ===================================================================

DATADIR = sprintf("test_data/devTree/%s", "20201014") # "20200707"
# DATADIR = "E:/Users/xyliu/data003/amph/tree"
date_tag = c("20210115", "20201014", "20200707")[1]
DATADIR = sprintf("test_data/devTree/%s", date_tag) # "20200707"

figdir = sprintf("%s/figs%s", DATADIR, '')
if(!dir.exists(figdir)){dir.create(figdir, recursive = T)}

### annotations of genes
Annos = readRDS("resources/Annos.rds")

# ===================================================================
#         loading 1) tree-structure and 2) group-expressions
# ===================================================================
message("loading edges of developmental tree...")
df_struct0 = read.csv(sprintf("%s/tree_struct-formal.csv", DATADIR),
                      na.strings = '', 
                      as.is = T) # c("parent", "node", "label"))
rownames(df_struct0) = df_struct0$node
str(df_struct0)

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

##### basic structure
tree0 = makePhyloFromDf(df_struct)

########################[ jpeg ]#######################
source("RFunxTreePlot.R")

figdir_ALL = "F:/amph_exprTree_jpeg"
all_gids = rownames(df_expr)
length(all_gids)
.gid = all_gids[1]

for (.gid in all_gids){
  .gnm = gsub('/', '_', x=Annos[.gid, 'gene_short_name'])
  .fn_save = sprintf("%s/exprTree-%s(%s).jpeg", figdir_ALL, .gid, .gnm)
  
  if (file.exists(.fn_save)){
    sprintf("skipped for %s(%s) has already been plotted.", .gid, .gnm) %>% message()
    next
  }
  plotExprTreeLR_full(tree0, .gid, df_expr, df_exprprop,
                      anno_list = LineageAnnoDetailed,
                      annos=Annos, cln = "tt",
                      color_nodelb = 'grey5',
                      size_tree = 1.,
                      size_nodelb = 5,
                      size_node_range=c(3, 12.5),
                      size_xticks = 15,
                      nodelb_offset = - 0.,
                      tiplb_offset = - 0.1,
                      xticks=c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0"),
                      height=600,
                      filename = .fn_save,
                      fontsize=5 # lineage annotations
                      ) -> p
  
}


