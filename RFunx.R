library("Seurat")
library("stringr")
source("COLOR_SETS.R")
COLORS = COLORS_z26

# DATADIR_main = "E:/Users/xyliu/data003/amph"
dir_anno = "E:/Users/xyliu/data003/amph"
if (dir.exists("E:/Users/xyliu/data003/amph")){
  dir_anno = "E:/Users/xyliu/data003/amph"
}else{
  dir_anno = "E:/lxy_pro/003/AmphioxusDevTree/resources"
}
  
Annos = readRDS(sprintf("%s/Annos.rds", dir_anno))
str(Annos)

StageOrd = c(
  E1 = "2cell",
  E2 = "4cell",
  E3 = "8cell",
  E4 = "32cell",
  E5 = "256cell",
  E6 = "B",
  E7 = "G3",
  E8 = "G4",
  E9 = "G5",
  E10 = "G6",
  E11 = "N0",
  E12 = "N1",
  E13 = "N3",
  E14 = "L0",
)

getGeneNames = function(gids,
                        anno_df=Annos, 
                        key="gene_short_name", 
                        from_srt=T){
  if (from_srt){
    gids = gsub("-", "_", gids)
  }
  anno_df[gids, key]
  
}

shortMultiNames = function(nms, sep=",", etc = "..", colla_sep=sep){
  sapply(nms, function(nm){
    nm_split = str_split(nm, sep, simplify = T)
    if(length(nm_split) >= 3){
      nm_split = nm_split[1, c(1, 2)]
      return(paste(c(nm_split, etc), collapse = colla_sep))
    }else{
      return(nm)
    }
  }, USE.NAMES = F)
  
}

make_foo_from_mapping <- function(mapping, # list, or a named vactor
                                  self_default = TRUE, 
                                  default_val = NULL,
                                  print_info = FALSE) {
  if (self_default){
    foo = function(x){
      if (! x %in% names(mapping)){
        return(x) # for invalid key, return itself as the value 
      }else{
        if(print_info){message(sprintf("TEST: changing %s to %s", x, mapping[x]))}
        return(mapping[x])
      }
    }
  }else{
    foo = function(x){
      if (! x %in% names(mapping)){
        return(default_val) # for invalid key, return the given default value 
      }else{
        return(mapping[x])
      }
    }
  }
  
  return(foo)
}

change_names = function(nms, mapping, ...){
  ## mns should be a vector
  # mapping should be a list or a named vector
  # if (is.list(mapping)){
  #   mapping = as.vector(mapping)
  # }
  # nms_new = mapping[nms]
  nms_new = sapply(
    as.vector(nms), 
    make_foo_from_mapping(mapping, ...),
    simplify=T, USE.NAMES = F
  ) %>% as.vector()
  
  if (! is.null(names(nms))){
    names(nms_new) <- names(nms)
  }
  return(nms_new)
}



WrapperDimPlot = function(obj, 
                          reduction='umap',
                          group.by = 'tree.ident',
                          dir_fig = NULL,
                          sname="temp",
                          label=T,
                          ...){
  plt = DimPlot(obj, reduction = 'umap', group.by = group.by, label = label,...)
  if(!is.null(dir_fig)){
    filename = sprintf("%s/%s_%s_%s.pdf", dir_fig, reduction, sname, group.by)
    ggsave(filename = filename, 
           plot = plt, 
           width = 5, height = 4)
    message(sprintf("figure saved into: %s", filename))
  }
  return(plt)
}


WrapperDotPlot = function(obj, 
                          genes_dot,
                          groupby=NULL,
                          gene_labs = NULL,
                          tt = NULL,
                          color_range = c("lightgrey", "mediumblue"),
                          x_rotate=90,
                          transpose = F,
                          dir_fig = NULL,
                          sname="temp",
                          w=NULL, h=NULL,
                          wscale=1., hscale=1.,
                          ...){
  message("Dot Feather plot...")
  # genes_dot = unique(markers_selected$gene)
  
  n_grps = nlevels(Idents(obj))
  
  plt_dot = DotPlot(obj, features = genes_dot, group.by = groupby,
                    cols = color_range, ...) + 
    theme(axis.text.x=element_text(angle=x_rotate, hjust=1)) #+ coord_flip()
  if (transpose){
    plt_dot = plt_dot + coord_flip()
    if(is.null(w)){w = min(5 + 0.16 * n_grps, 49.5)}
    if(is.null(h)){h = min(4 + 0.12 *length(genes_dot), 49.5)}
  }else{
    if(is.null(w)){w = min(6 + 0.15 *length(genes_dot), 49.5)}
    if(is.null(h)){h = min(4.8 + 0.16 * n_grps, 49.5)}
  }
  w = w * wscale
  h = h * hscale
  if (! is.null(gene_labs)){
    plt_dot = plt_dot + scale_x_discrete(labels = rev(gene_labs))
  }
  if (! is.null(tt)){
    plt_dot = plt_dot + ggtitle(tt)
  }
  if(!is.null(dir_fig)){
    filename = sprintf("%s/dot_%s.pdf", dir_fig, sname)
    ggsave(filename = filename, 
           plot = plt_dot, width = w,
           height = h)
    message(sprintf("figure saved into: %s size=(%.2f, %.2f)", filename, w, h))
  }
  return(plt_dot)
}




WrapperFeaturePlot = function(obj,
                              genes_to_plot,
                              n_col = NULL,
                              color_range = c("lightgrey", "mediumblue"),
                              dir_fig = NULL,
                              sname="temp",
                              ...) {
  message("Feather plot...")
  # genes_to_plot = unique(top_n(markers_selected, 1, pct.1)$gene)
  n = length(genes_to_plot)
  if(is.null(n_col)){
    n_col = ifelse(n %% 4 >= 3 || n >= 10, yes = 4, no = 3)
  }
  n_row = ceiling(n / n_col)

    plt_marker = FeaturePlot(obj, features = genes_to_plot, ncol = n_col, 
                             cols = color_range, ...)#, min.cutoff = "q8")
  if(!is.null(dir_fig)){
    filename = sprintf("%s/scatter_markers_%s.pdf", dir_fig, sname)
    ggsave(filename = filename, 
           plot = plt_marker, 
           width = 4 * n_col, height = 4*n_row)
    message(sprintf("figure saved into: %s", filename))
  }
  return(plt_marker)
}

WrapperEbowPlot = function(obj, 
                           n_pcs_use=NULL,
                           dir_fig = NULL,
                           sname="temp",
                           ...){
  # WrapperEbowPlot(obj, n_pcs_use, dir_fig, sname, ndims = n_pcs_calc)
  plt <- ElbowPlot(obj, ...) 
  if(!is.null(n_pcs_use)){
    plt = plt + geom_vline(xintercept = n_pcs_use, color='blue')
  }
  if(!is.null(dir_fig)){
    filename = sprintf("%s/PCs_%s.pdf", dir_fig, sname)
    ggsave(filename = filename, 
           plot = plt, width = 4, height = 3)
    message(sprintf("figure saved into: %s", filename))
  }
  return(plt)
  
}



LoadCounts = function(dir_data, sname,
                      tag="_afterQC",
                      fn_raw_rds = sprintf("%s/%s%s.rds", dir_data, sname, tag),
                      dn_raw_mtx = sprintf("%s/%s%s_mtx", dir_data, sname, tag)
                      ){
  
  if(file.exists(fn_raw_rds)){
    message(sprintf("Loading data from %s", fn_raw_rds))
    cnt = readRDS(fn_raw_rds)
  }else{
    
    message(sprintf("Loading data from %s", dn_raw_mtx))
    cnt = Seurat::Read10X(dn_raw_mtx, gene.column = 1)
    
    message("backup .rds file...")
    saveRDS(cnt, fn_raw_rds)
  }
  return(cnt)
}


LoadGeneAnnotations = function(maindir_data=NULL, as.gr=TRUE){
  if(is.null(maindir_data)){
    maindir_data = "E:/Users/xyliu/data003/amph"
  }
  
  # load gene annotations (NOTE: better use `read.delim()` )
  fn_gene_annos = file.path(maindir_data, "gene_annos_merged_srt_bed.tsv")
  genes_df <<- read.delim(fn_gene_annos, header=TRUE, sep='\t', as.is=T)
  message("`gene_df` added to the namespace")
  print(levels(as.factor(genes_df$chr)))
  str(genes_df)
  if(as.gr){
    suppressMessages(library("GenomicRanges", quietly =TRUE))
    genes_gr <<- GRanges(
      seqnames = as.factor(genes_df$chr),
      ranges = IRanges(genes_df$start, genes_df$end),
      name = genes_df$ID_srt,
      # NOTE: seurat replaced '_' in the IDs by '-'
      strand = genes_df$strand)
    message("`gene_gr` added to the namespace")
  }
  
}


GenesPerCell = function(mat){
  genes_per_cell = Matrix::colSums(mat > 0)
  print("Summary of the number of genes per cell:")
  print(summary(genes_per_cell))
  return(genes_per_cell)
}

CountsPerCell = function(mat){
  counts_per_cell = Matrix::colSums(mat)
  print("Summary of the number of counts per cell:")
  print(summary(counts_per_cell))
  return(counts_per_cell)
}


SaveMM = function(mat, dirname){
  if(! dir.exists(dirname)){dir.create(dirname)}
  message("saving .rds file...")
  saveRDS(mat, file = file.path(dirname, "rawconts.rds"))
  message("saving .mtx file...")
  Matrix::writeMM(obj = mat, file = file.path(dirname, "matrix.mtx"))
  message("saving names...")
  write.table(colnames(mat),
              file = file.path(dirname, "barcodes.tsv"),
              col.names = F, row.names = F, quote = F)
  write.table(rownames(mat),
              file = file.path(dirname, "genes.tsv"),
              col.names = F, row.names = F, quote = F)
  print("The matrix has been saved into directory:")
  print(dirname)
  
}





RoughQC = function(dn, name="genes", max_cut=100,
                   plot_ranks = TRUE, 
                   dn_new=NULL){
  print(dn)
  
  # BUSpaRse::read_count_output
  # Matrix::readMM(file.path(fname, "gene.mtx"))
  res_mat <- BUSpaRse::read_count_output(dn, 
                               name = name, tcc = FALSE)
  
  print(sprintf("shape of the RAW matrix: %d (genes), %d (cells)", 
                dim(res_mat)[1], dim(res_mat)[2]))
  
  cell_tot_counts <- Matrix::colSums(res_mat)
  print(summary(cell_tot_counts))
  
  # Compute barcode rank
  # barcodeRanks {DropletUtils}
  bc_rank <- barcodeRanks(res_mat)
  
  ### --------Plot----------
  if(plot_ranks){
    plt <- qplot(bc_rank$total, bc_rank$rank, geom = "line") +
    geom_vline(xintercept = bc_rank$knee, color = "blue", linetype = 2) +
    geom_vline(xintercept = bc_rank$inflection, color = "green", linetype = 2) +
    # annotate("text", y = 1000, x = 1.5 * c(bc_rank$knee, bc_rank$inflection),
    #          label = c("knee", "inflection"), color = c("blue", "green")) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Barcode rank", x = "Total UMI count")
  plot(plt)
  
  }
  
  ## --- Filter the matrix
  cutoff <- as.integer(min(max_cut, bc_rank$inflection))
  print(sprintf("Min-gene cutoff: %d (inflection: %.1f)", cutoff, bc_rank$inflection))
  
  ## filter cells
  res_mat <- res_mat[, cell_tot_counts > cutoff]
  
  ## filter genes
  gene_tot_counts <- Matrix::rowSums(res_mat)
  res_mat <- res_mat[gene_tot_counts > 0, ]
  
  ## print("shape of the filered matrix:")
  print(sprintf("shape of the filered matrix: %d (genes), %d (cells)", 
                dim(res_mat)[1], dim(res_mat)[2]))
  print("Counts summary:")
  print(summary(Matrix::colSums(res_mat)))
  
  ## --- Prepare and save the filtered matrix
  if (is.null(dn_new)){
    dn_new <- file.path(dn, sprintf("nCcut_%d", cutoff))
  }
  SaveMM(res_mat, dn_new)

  
  return(res_mat)
}














