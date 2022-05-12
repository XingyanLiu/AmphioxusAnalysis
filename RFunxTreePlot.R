# BiocManager::install("ggtree")
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(RColorBrewer)


### =================================================
###         Global settings
### =================================================

TREE_W = 7.5 # width of tree
TREE_H = 6 # height of tree
ladderiz = F

COLORS_tab20 = c('#1f77b4',  '#aec7e8',
                 '#ff7f0e',  '#ffbb78',
                 '#2ca02c',  '#98df8a',
                 '#d62728',  '#ff9896',
                 '#9467bd',  '#c5b0d5',
                 '#8c564b',  '#c49c94',
                 '#e377c2',  '#f7b6d2',
                 '#7f7f7f',  '#c7c7c7',
                 '#bcbd22',  '#dbdb8d',
                 '#17becf',  '#9edae5')
COLORS_plasma = c(
  '#f0f921', '#fdc328', '#f89441','#e56b5d','#cb4679',
  '#a82296', '#7d03a8', '#4b03a1', '#0d0887'
)

LineageAnno = list(
  G3_0 = "Epithelial ectoderm",
  G3_1 = "Endoderm",
  G3_2 = "Mesoderm",
  G3_3 = "Notochord",
  G3_4 = "Neural ectoderm",
  G3_5 = "Uncertain",
  G3_6 = "Primordial germ cells"
)

LineageAnnoDetailed = list(
  G3_0 = "Epithelial ectoderm",      #
  G3_3 = "Notochord",                #
  # Neural ectoderm
  L0_6 = "Posterior neural precursor",
  L0_9 = "Anterior neural precursor",
  L0_10 = "Differentiated neurons",
  # Mesoderm
  G6_1 = "Ventral mesoderm\nposteriorto-pharynx",
  N3_6 = "Pharyngeal mesoderm",
  L0_17 = "Tail bud mesoderm",
  # 7-5
  L0_7 = "Uncertain",                #
  # Endoderm
  L0_8 = "Mid-gut",
  L0_3 = "Fore-gut",
  L0_5 = "Fore-gut",
  L0_0 = "Hind-gut",
  L0_15 = "Pharynx",
  # PGC
  L0_18 = "Primordial germ cells"    #
)

StageNameMap = list(
  E6 = "B",
  E7 = "G3",
  E8 = "G4",
  E9 = "G5",
  E10 = "G6",
  E11 = "N0",
  E12 = "N1",
  E13 = "N3",
  E14 = "L0"
)

StageNames = c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")

# CMAP_lineage = COLORS_tab20[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)]
CMAP_lineage = COLORS_LINE = c(
  '#596e79', '#b3b3b3', '#c7b198', # B1-B3
  '#40bad5', '#984ea3', '#36622b', 
  '#035aa6', '#fcbf1e', '#af0404', 
  '#dd7631'
)
COLORS_stg = list(
  colorRampPalette(brewer.pal(11, 'Spectral'))(9),
  COLORS_plasma
  )[[2]]

### =================================================
###         helper functions
### =================================================

makeListColor = function(lst_names, colrs=COLORS_tab20){
  clst = as.list(colrs[1:length(lst_names)])
  names(clst) = lst_names
  clst
}


LineageColor = makeListColor(names(LineageAnno), CMAP_lineage[-c(1:3)])



nodeFromPlt = function(plt, node_lb){
  subset(plt$data, label == node_lb)$node
}

diyColors = function(cmap_name="RdPu", n=7, bg="#F0F0F0"){
  # `cmap_name` must be a validate name in brewer.pal()
  colors = brewer.pal(n, cmap_name)
  if(!is.null(bg)){colors[1] = bg}
  colors
}




makeTreeDataStage = function(df, parent="parent", node="node", label="label",
                             stg_levels = levels(df$stage),
                             cols_ext = NULL #c("label", "stage", "stage_int")
                             ){
  cols_tree = c(parent, node, label)
  tb_tree = as_tibble(df[, cols_tree]) %>% filter(!is.na(parent))
  colnames(tb_tree) = c("parent", "node", "label")
  if (is.null(cols_ext)){
    cols_ext = colnames(df)[! colnames(df) %in% cols_tree[-3]]
  }
  tb_extdata = as_tibble(df[, cols_ext]) #
  tb_extdata$stage = factor(tb_extdata$stage, levels = stg_levels, ordered = T)
  tree0 = as.phylo(tb_tree)
  tree0$edge.length = 1
  # prepare for stage-tree
  treedt = full_join(x = as_tibble(tree0),
                     y = tb_extdata, # y = tb_data_expr,
                     by='label') %>% as.treedata()
  treedt
}


makePhyloFromDf = function(df, parent="parent", node="node", label="label", 
                           na.rmv=T){
  # basic tree structure
  tb_tree = as_tibble(df[, c(parent, node, label)]) 
  if (na.rmv){
    tb_tree = tb_tree %>% filter(!is.na(parent))
  }
  colnames(tb_tree) = c("parent", "node", "label")
  tree0 = as.phylo(tb_tree)
  tree0$edge.length = 1
  tree0
}



makeTreeDataExpr = function(tree0, df_expr, df_exprprop=NULL, gene_name){
  # combine basic structure with expressions
  if (! gene_name %in% rownames(df_expr)){
    df_expr[gene_name, ] = 0
  }
  tb_expr = as_tibble(t(df_expr[gene_name, ]))
  if(!is.null(df_exprprop)){
    if (! gene_name %in% rownames(df_expr)){
      df_expr[gene_name, ] = 0
    }
    tb_expr$proportion = t(df_exprprop[gene_name, ])[, 1]
  }
  tb_expr$label = colnames(df_expr)
  # treedt_expr = full_join(x = as_tibble(tree0),
  #                         y = tb_expr,
  #                         by='label') %>% as.treedata() 
  treedt_expr = left_join(x = as_tibble(tree0),
                          y = tb_expr,
                          by='label') %>% as.treedata() 
  treedt_expr
}


### used for annotation of lineages
labelClade = function(plt, anno_list=LineageAnno, 
                      anno_color_list=list(), 
                      barsize=1, extend=0.4,
                      align = T, offset = 0.4,
                      default_colr = 'grey30',
                      fontsize = 3.88,
                      ...){
  for (node_lb in names(anno_list)){
    cld_lb = anno_list[[node_lb]]
    colr = ifelse(is.null(anno_color_list[[cld_lb]]), 
                  default_colr, 
                  anno_color_list[[node_lb]])
    plt = plt + geom_cladelabel(node = nodeFromPlt(plt, node_lb), #soffset.text = 
                                label = cld_lb, 
                                barsize = barsize,
                                fontsize = fontsize,
                                color = colr,
                                extend = extend,
                                align = align, 
                                offset = offset, 
                                ...)
    
  }
  plt
  
}
hlightClade = function(plt, anno_list=LineageAnno, 
                      anno_color_list=list(), 
                      # barsize=1, 
                      extend=0.4,
                      alpha = 0.5, 
                      # align = T, offset = 0.4,
                      default_colr = 'grey30',
                      # fontsize = 3.88,
                      ...){
  for (node_lb in names(anno_list)){
    cld_lb = anno_list[[node_lb]]
    colr = ifelse(is.null(anno_color_list[[node_lb]]), 
                  default_colr, 
                  anno_color_list[[node_lb]])
    plt = plt + geom_hilight(node = nodeFromPlt(plt, node_lb), #soffset.text = 
                             fill = colr, #"steelblue", 
                             alpha = alpha, 
                             extend = extend,
                             ...)
    
  }
  plt
  
}


### =================================================
###         plotting functions
### =================================================





# displaying expression levels and proportions
plotExprTree = function(treedt_expr, gene_name, title=gene_name,
                        colors_expr=diyColors("RdPu"),
                        size_tree = 1.,
                        size_tiplab = 3,
                        filename=NULL,
                        width=TREE_W, height=TREE_H, ...){
  mapping_expr = parse(text = sprintf("aes(color=%s, size=proportion)", gene_name))
  p_expr = ggtree(treedt_expr, size=size_tree, color='lightgrey',
                  #mapping = eval(mapping_expr), 
                  # continuous = TRUE, # for coloring branches, might cause error
                  ladderize = ladderiz, ...) +
    layout_dendrogram() +
    geom_point(eval(mapping_expr), alpha=1) +
    scale_color_gradientn(colors = colors_expr) +
    geom_tiplab(aes(label=label), 
                color='grey', 
                angle=90,
                hjust=1, offset=-0.2,
                align=TRUE, size=size_tiplab) + 
    # scale_x_reverse() + coord_flip() +
    ggtitle(title) +
    # geom_nodelab(aes(label=label, x = branch + 0.1, y=y+0.2), color='grey', size=3) +
    # geom_nodelab(eval(mapping_expr_nd), size=3) +
    theme(legend.position="right") 
  if (!is.null(filename)){
    ggsave(filename = filename, 
           plot=p_expr, 
           width=width, height=height)
  }
  p_expr
}

# displaying expression levels and proportions
plotExprTreeLR = function(treedt_expr, gene_name, title=gene_name,
                        colors_expr=diyColors("RdPu"),
                        color_tree='lightgrey',
                        xmax=NA,
                        anno_list = NULL,
                        anno_offset = 0.8,
                        size_tree = 1.,
                        size_tiplab = 3,
                        size_node_range = c(1, 6),
                        tiplab = TRUE,
                        filename=NULL,
                        width=TREE_W * 1.2, 
                        height=TREE_H,
                        ...){
  mapping_expr = parse(text = sprintf("aes(color=%s, size=proportion)", gene_name))
  p_expr = ggtree(treedt_expr, size=size_tree, color=color_tree,
                  #mapping = eval(mapping_expr), 
                  # continuous = TRUE, # for coloring branches, might cause error
                  ladderize = ladderiz) +
    geom_point(eval(mapping_expr), alpha=1) +
    scale_color_gradientn(colors = colors_expr) +
    scale_size(range=size_node_range) #c(5, 30)
  if (tiplab){
    p_expr <- p_expr + 
      geom_tiplab(aes(label=label), 
                  color='grey', 
                  align=TRUE, 
                  size=size_tiplab) 
  }
    # scale_x_reverse() + coord_flip() +
  if (! is.null(title)){
    p_expr <- p_expr + ggtitle(title)
  }
    
  p_expr <- p_expr + theme(legend.position="right") 
  
  if (! is.null(anno_list)){
    # leave more space for tip labels
    xmax = ifelse(is.na(xmax), max(p_expr$data$x) + 3 + anno_offset, xmax)
    p_expr = p_expr + xlim(0, xmax)
    
    p_expr = labelClade(p_expr, anno_list = anno_list, 
                        offset = anno_offset,
                        ...)
  }
  if (!is.null(filename)){
    ggsave(filename = filename, 
           plot=p_expr, 
           width=width, height=height)
  }
  p_expr
}

plotExprTreeLR_full = function(
  tree0, gid, 
  df_expr, df_exprprop, 
  tt=NULL,
  annos=NULL, cln='tt',
  colors_expr=diyColors("RdPu"),
  color_tree='grey70',
  color_nodelb = 'grey20',
  color_anno_default = 'grey10',
  xmax=9,
  xticks=NULL,#c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
  anno_list = NULL, # for `labelClade``
  anno_offset = 0.2, #old=0.8
  nodelb_offset = - 0.35,
  tiplb_offset = - 0.45,
  size_tree = 1.,
  size_node_range = c(1, 6),
  size_nodelb = 3,
  size_xticks = 3,
  extend = 0.25,
  plot.margin = margin(6, 150, 6, 6),
  filename=NULL,
  wscale = 1.2, hscale = 1,
  # width=TREE_W * 1.2, 
  # height=TREE_H, # ignored if JPEG format
  ... # for `labelClade`: fontsize = 3,
){
  # making tree data
  treedt_expr = makeTreeDataExpr(tree0, df_expr, df_exprprop, gid)
  if(is.null(tt)){
    tt = ifelse(is.null(annos), gid, subset(annos, ID == gid)[[cln]])
    print(sprintf("title of figure: %s", tt))
    }
  p_expr = plotExprTreeLR(treedt_expr, gid, 
                      title = tt, 
                      color_tree=color_tree,
                      size_tree = size_tree,
                      size_node_range=size_node_range,
                      tiplab = FALSE,
                      xmax=NA,
                      anno_list = anno_list, 
                      anno_offset = anno_offset,
                      extend=extend,
                      default_colr = color_anno_default,
                      ...)
  # xmax = 9
  # stages = c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
  # xticks = c(stages, rep("", xmax - length(stages)))
  p_expr = p_expr + geom_nodelab(
    aes(x = x + nodelb_offset, 
        label=str_split_fixed(label, '_', 2)[, 2]), 
    color=color_nodelb, size=size_nodelb) +
    geom_tiplab(
      aes(x = x + tiplb_offset, 
          label=str_split_fixed(label, '_', 2)[, 2]), 
      color=color_nodelb, size=size_nodelb) +
    scale_x_discrete(limits=c(1: xmax), labels = xticks) +
    coord_cartesian(clip = 'off') + # avoiding texts from being cut
    theme_tree2(plot.margin=plot.margin) +
    theme(legend.position="left", 
          axis.text.x= element_text(size=size_xticks)#family, face, colour, 
          )
  
  if (!is.null(filename)){
    if (str_ends(filename, '.jpeg')){
      w = 900# * wscale
      h = 600 * hscale
      print(filename)
      jpeg(filename, width = w, height = h)
      print(p_expr)
      dev.off()
      
    }else{
      w = 7.5 * wscale
      h = 6 * hscale
      ggsave(filename = filename, plot=p_expr, 
             width=w, height=h)
    }
  }
  p_expr
  
}


plotLabeledExprTree = function(treedt_expr, gene_name, 
                               title=gene_name,
                               colors_expr=diyColors("RdPu"),
                               color_tree = 'lightgrey',
                               size_tree = 1.,
                               size_label = 2.4, 
                               xmax=NA,
                               anno_list = NULL,
                               anno_offset = 0.4,
                               filename=NULL,
                               width=TREE_W * 1.2, 
                               height=TREE_H, 
                               ...){
  
  mapping_expr_lb = parse(text = sprintf("aes(label=label,fill=%s)", gene_name))
  p_expr_lb = ggtree(treedt_expr, size=size_tree, color=color_tree,
                  #mapping = eval(mapping_expr), 
                  # continuous = TRUE, # for coloring branches, might cause error
                  ladderize = ladderiz) +
    scale_color_gradientn(colors = colors_expr) +
    geom_label(mapping = eval(mapping_expr_lb), color='black', size=size_label) +
    scale_fill_gradientn(colors = colors_expr) +
    ggtitle(title) +
    # scale_x_reverse() + coord_flip() +
    # theme_tree2(fgcolor='lightgrey', legend.position="right")
    theme(legend.position="right")
  
  if (! is.null(anno_list)){
    # leave more space for tip labels
    xmax = ifelse(is.na(xmax), max(p_expr_lb$data$x) + 3.8, xmax)
    p_expr_lb = p_expr_lb + xlim(0, xmax)
    
    p_expr_lb = labelClade(p_expr_lb, anno_list = anno_list, 
                        offset = anno_offset,
                        ...)
  }
  if (!is.null(filename)){
    ggsave(plot=p_expr_lb, filename =filename,
           width=width, height=height)
  }
  p_expr_lb
}


plotLabeledStageTree = function(treedt, 
                                size_tree=1.5,
                                size_label=2.5,
                                colors_stg = COLORS_stg,
                                # colors_stg = brewer.pal(11, 'Spectral'),
                                anno_list = NULL,
                                xmax=NA,
                                filename=NULL,
                                width=TREE_W, height=TREE_H, ...){
  # === filled node labels ===
  p = ggtree(treedt, size=size_tree, color='lightgrey', ladderize = ladderiz) + #, layout = 'slanted')
    geom_label2(aes(label=label, fill = stage_int, subset= label != "root"), 
               alpha = 1, color='black', 
               size=size_label, 
               na.rm = T) +
    # scale_color_manual(values = colors_stg) +
    scale_fill_gradientn("stage", colors = colors_stg) +
    # geom_tiplab(aes(label=label), color='darkgrey', align=TRUE, hjust = -0.5, size=3) + 
    theme_tree(legend.position="right") 
  
  if (! is.null(anno_list)){
    # leave more space for tip labels
    xmax = ifelse(is.na(xmax), max(p$data$x) + 2.5, xmax)
    p = p + xlim(0, xmax)
    
    p = labelClade(p, anno_list = anno_list, ...)
  }
  if (!is.null(filename)){
    ggsave(plot=p, filename = filename,
           width=width, height=height)
  }
  p
  
}





plotStageTree = function(treedt, 
                         size_tree=1.5,
                         size_node=4.5,
                         size_tiplab=2.,
                         colors_stg = COLORS_stg,
                         # n_level = 15,
                         # colors_stg=colorRampPalette(brewer.pal(11, 'Spectral'))(n_level),
                         filename=NULL,
                         width=TREE_W, height=TREE_H, ...){
  
  p = ggtree(treedt, size=size_tree, color='lightgrey', 
             ladderize = ladderiz, ...) + 
    layout_dendrogram() +
    geom_point(aes(color=stage), size = size_node, na.rm = T, alpha=1) +
    scale_color_manual(values = colors_stg) +
    geom_tiplab(aes(label=label), 
                color='darkgrey', 
                angle=90,
                hjust=1, offset=-0.2,
                align=TRUE, size=size_tiplab) + 
    theme_tree(legend.position="right", plot.margin=margin(1, 20, 1, 1))# +
    # xlim(NA, 1) + # useless
    # scale_x_reverse() + coord_flip() + ### used when `layout_denrogram` does nnt work (R3.5)
    # theme(legend.position="right") #,plot.margin = unit(c(0.01, 0.01, 0.05, 0.01), units = "npc"),
  
  if (!is.null(filename)){
    ggsave(filename = filename, 
           plot=p, 
           width=width, height=height)
  }
  p
}


plotStageTreeLR = function(treedt, 
                         size_tree=1.5,
                         size_node=4.5,
                         size_tiplab=2.,
                         xmax=NA,
                         anno_list = NULL,
                         root_lb = "root",
                         colors_stg = COLORS_stg,
                         # n_level = 15,
                         # colors_stg=colorRampPalette(brewer.pal(11, 'Spectral'))(n_level),
                         filename=NULL,
                         width=TREE_W, height=TREE_H, ...){
  
  p = ggtree(treedt, size=size_tree, color='lightgrey', 
             ladderize = ladderiz) + 
    geom_point2(aes(color=stage, subset= label != root_lb), 
                size = size_node, na.rm = T, alpha=1) +
    scale_color_manual(values = colors_stg) +
    theme_tree(legend.position="left", plot.margin=margin(1, 20, 1, 1)) 
  # theme(legend.position="right") #,plot.margin = unit(c(0.01, 0.01, 0.05, 0.01), units = "npc"),
  
  if (! is.null(anno_list)){
    # leave more space for tip labels
    xmax = ifelse(is.na(xmax), max(p$data$x) + 4, xmax)
    p = p + xlim(0, xmax)
    
    p = labelClade(p, anno_list = anno_list, default_colr = 'gray50',...)
  }
  
  if (!is.null(filename)){
    ggsave(filename = filename, 
           plot=p, 
           width=width, height=height)
  }
  p
}
plotStageTreeLR1 = function(treedt, 
                           size_tree=1.5,
                           size_node=4.5,
                           size_tiplab=2.,
                           alpha_point=1,
                           xmax=NA,
                           anno_list = NULL,
                           root_lb = "root",
                           color_tree = 'lightgrey',
                           colors_stg = COLORS_stg,
                           # n_level = 15,
                           # colors_stg=colorRampPalette(brewer.pal(11, 'Spectral'))(n_level),
                           filename=NULL,
                           width=TREE_W, height=TREE_H, ...){
  
  p = ggtree(treedt, size=size_tree, color=color_tree, 
             ladderize = ladderiz) + 
    geom_point2(aes(fill=stage, subset= label != root_lb), 
               color=color_tree, 
               shape=21,
               size = size_node, 
               na.rm = T, alpha=alpha_point) +
    scale_fill_manual(values = colors_stg) +
    theme_tree(legend.position="left", plot.margin=margin(1, 20, 1, 1)) 
  # theme(legend.position="right") #,plot.margin = unit(c(0.01, 0.01, 0.05, 0.01), units = "npc"),
  
  if (! is.null(anno_list)){
    # leave more space for tip labels
    xmax = ifelse(is.na(xmax), max(p$data$x) + 4, xmax)
    p = p + xlim(0, xmax)
    
    p = labelClade(p, anno_list = anno_list, ...)
  }
  
  if (!is.null(filename)){
    ggsave(filename = filename, 
           plot=p, 
           width=width, height=height)
  }
  p
}




### =================================================
###         Wrapper of plotting-functions 
### =================================================

plotTreesBoth = function(genes, tree0, df_expr, df_exprprop,
                         annos=NULL, cln="tt",
                         fname=NULL, width=TREE_W, height=TREE_H){
  plts = sapply(genes, function(gid){
    tt = ifelse(is.null(annos), gid, subset(annos, ID == gid)[[cln]])
    treedt_expr = makeTreeDataExpr(tree0, df_expr, df_exprprop, gid)
    p1 = plotExprTree(treedt_expr, gid, title = tt)
    p2 = plotLabeledExprTree(treedt_expr, gid, title = tt)
    list(p1, p2)
  })
  if (! is.null(fname)){
    print(fname)
    pdf(fname, width = width, height = height)
    invisible(lapply(plts, print))
    dev.off()
  }
  
  plts
}


plotTreesBothLR = function(genes, tree0, df_expr, df_exprprop,
                         annos=NULL, cln="tt",
                         anno_list_lineage = NULL,
                         fname=NULL, 
                         width=TREE_W *1.2, height=TREE_H,
                         ...){
  
  plts = sapply(genes, function(gid){
    tt = ifelse(is.null(annos), gid, subset(annos, ID == gid)[[cln]])
    print(sprintf("title of figure: %s", tt))
    treedt_expr = makeTreeDataExpr(tree0, df_expr, df_exprprop, gid)
    p1 = plotExprTreeLR(treedt_expr, gid, 
                        title = tt, 
                        anno_list = anno_list_lineage, 
                        ...)
    p2 = plotLabeledExprTree(treedt_expr, gid, 
                             title = tt, 
                             anno_list = anno_list_lineage, 
                             ...)
    list(p1, p2)
  }) # simplified
  if (! is.null(fname)){
    print(fname)
    pdf(fname, width = width, height = height)
    invisible(lapply(plts, print))
    dev.off()
  }
  
  plts # a list of N*2 length-lists (dot and labeled expression tree)
}


plotTreesBothLR_jpeg = function(genes, tree0, df_expr, df_exprprop,
                           annos=NULL, cln="tt",
                           anno_list_lineage = NULL,
                           fname=NULL, 
                           width=900, 
                           height=900,
...){
  
  plts = sapply(genes, function(gid){
    tt = ifelse(is.null(annos), gid, subset(annos, ID == gid)[[cln]])
    print(sprintf("title of figure: %s", tt))
    treedt_expr = makeTreeDataExpr(tree0, df_expr, df_exprprop, gid)
    p1 = plotExprTreeLR(treedt_expr, gid, 
                        title = tt, 
                        anno_list = anno_list_lineage,
                        ...)
    p2 = plotLabeledExprTree(treedt_expr, gid, 
                             title = NULL, #tt, 
                             anno_list = anno_list_lineage,
                             ...)
    list(p1, p2)
  }) # simplified
  # if (! is.null(fname)){
  #   lapply(plts, function(p){
  #     message(fname)
  #     jpeg(fname, width = width, height = height)
  #     invisible(print(p))
  #     dev.off()
  #   })
  if (! is.null(fname)){
    print(fname)
    jpeg(fname, width = width, height = height)
    invisible(print(cowplot::plot_grid(plotlist = plts, ncol = 1)))
    dev.off()
  }
  
  plts # a list of N*2 length-lists (dot and labeled expression tree)
}














