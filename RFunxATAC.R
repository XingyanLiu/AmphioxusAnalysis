library(stringr)
source("COLOR_SETS.R")
COLORS = COLORS_z26

DATADIR_main = "E:/Users/xyliu/data003/amph"
RNAStages = paste0("E", c(1:15))
ATACStages = c("blastula","G3", "G6",  "N1", "N3", "L0")

StageMap = list(
  E6 = "blastula",
  E7 = "G3",
  E10 = "G6",
  E12 = "N1",
  E13 = "N3",
  E14 = "L0"
)


Stages_polyT = paste0("E", c(1:8, 10))


LineageOrdAll = c(
  E6_0 = 'B_0',
  E6_1 = 'B_1',
  E6_2 = 'B_2',
  E6_3 = "Primordial germ cells",
  
  E7_6 = "Primordial germ cells",
  E7_0 = "Epithelial ectoderm",
  E7_4 = "Neural ectoderm",
  E7_3 = "Notochord",
  E7_2 = "Mesoderm",
  E7_5 = "Tail bud stem cells",
  E7_1 = "Endoderm"
)
LineageOrd = LineageOrdAll[-c(1:3)]

LineageColors = c(
  "B_0" = '#596e79',
  "B_1" = '#b3b3b3',
  "B_2" = '#c7b198',
  "Primordial germ cells" = '#40bad5',    # 7-6
  "Epithelial ectoderm" = '#984ea3',      # 7-0
  "Neural ectoderm" = '#36622b',          # 7-4
  "Notochord" = '#035aa6',                # 7-3
  "Mesoderm" = '#fcbf1e',                 # 7-2
  "Unassigned" = '#af0404',      # 7-5
  "Endoderm" = '#dd7631'                  # 7-1
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

ColorSets = list(
  lineage = LineageColors,
  RNAcluster = COLORS_tab20,
  refined_group = COLORS_tab20
)

# =============== functions for names =================

StageSname = function(stg, polyT = Stages_polyT){
  sn = ifelse(stg %in% polyT, paste0(stg, '_polyT'), stg)
  sn
}

ShowColors = function(colors=COLORS){
  orders = seq(colors)
  y = rep(1, length(colors))
  plot(orders, y, col=colors, pch=19)
}

# ========================= functions for cicero GA ===========================

prepareAnno4GA = function(gene_anno){
  # ( preparing data )
  #[ Add a column for the pData table indicating the gene if a peak is a promoter ]#
  # Create a gene annotation set that only marks the transcription start sites of 
  # the genes. We use this as a proxy for promoters.
  # To do this we need the first exon of each transcript
  pos <- subset(gene_anno, strand == "+")
  pos <- pos[order(pos$start),] 
  # remove all but the first exons per transcript
  pos <- pos[!duplicated(pos$transcript),] 
  # make a 1 base pair marker of the TSS
  pos$end <- pos$start + 1 
  
  neg <- subset(gene_anno, strand == "-")
  neg <- neg[order(neg$start, decreasing = TRUE),] 
  # remove all but the first exons per transcript
  neg <- neg[!duplicated(neg$transcript),] 
  neg$start <- neg$end - 1
  
  gene_annotation_sub <- rbind(pos, neg)
  # Make a subset of the TSS annotation columns containing just the coordinates 
  # and the gene name
  gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]
  # Rename the gene symbol column to "gene"
  names(gene_annotation_sub)[4] <- "gene"
  
  gene_annotation_sub
}

peakMetaFromNames = function(pknames, sep1=":", sep2="-", sep_new="_"){
  tmp = str_split(pknames, sep1, simplify = TRUE)
  chrs = tmp[, 1]
  sites = str_split(tmp[, 2], sep2, simplify = TRUE)
  peaknames1 <- paste(chrs, sites[,1], sites[,2], sep=sep_new)
  head(peaknames1)
  # site_name chromosome bp1 bp2
  peak_meta = data.frame(site_name = peaknames1, 
                         chromosome = chrs, 
                         bp1 = as.numeric(sites[,1]), 
                         bp2 = as.numeric(sites[,2]), 
                         row.names = peaknames1) 
  peak_meta
}

makeMatSummary = function(mat, row_names = c("n_features", "n_counts")){
  summ = list()
  Matrix::colSums(mat > 0) -> summ$gs_per_cell
  Matrix::colSums(mat) -> summ$cnts_per_cell
  summ$summdf = data.frame(rbind(summary(summ$gs_per_cell), summary(summ$cnts_per_cell)), 
                           row.names = row_names)
  summ
}

histMatSummary = function(summ_list,
                          fn_pdf = "temp_mat_summary_hist.pdf",
                          xg = "features per cell",
                          xc = "counts per cell",
                          colr = "lightblue"){
  pdf(fn_pdf, width = 6, height = 6)
  par(mfcol = c(2, 2))
  hist(summ_list$gs_per_cell, xlab = xg, main = xg, col = colr)
  hist(log10(summ_list$gs_per_cell), xlab = paste(xg, "(log10)"), 
       main = xg, col = colr)
  hist(summ_list$cnts_per_cell, xlab = xc, main = xc, col = colr)
  hist(log10(summ_list$cnts_per_cell), xlab =  paste(xc, "(log10)"), 
       main = xc, col = colr)
  dev.off()
  print(fn_pdf)
}


processCDS4ATAC = function(cds){
  message("detecting genes...")
  cds <- detect_genes(cds)
  message("estimating size factors...")
  cds <- estimate_size_factors(cds)
  message("preprocessing...")
  cds <- preprocess_cds(cds, method = "LSI")
  message("computing UMAP...")
  cds <- reduce_dimension(cds, reduction_method = 'UMAP', 
                          preprocess_method = "LSI")
  message("Done!")
  cds
}

runCustomedCicero = function(cds, 
                             genomic_coords,
                             window_size = 5e+05,
                             distance_constraint = 250000,
                             s = 0.75,
                             sample_num = 100
                             ){
  print("Starting Cicero")
  print("Calculating distance_parameter value")
  distance_parameters <- estimate_distance_parameter(
    cds, sample_num=sample_num, genomic_coords = genomic_coords,
    window = window_size,
    s = s,
    distance_constraint = distance_constraint)
  mean_distance_parameter <- mean(unlist(distance_parameters))
  print("Running models")
  cicero_out = generate_cicero_models(
    cds,
    distance_parameter = mean_distance_parameter,
    s = s, 
    window = window_size, 
    genomic_coords = genomic_coords
  )
  print("Assembling connections")
  conns <- assemble_connections(cicero_out)
  conns
}



# ================================================================================

clusterATAC = function(obj,
                       use_pcs=2:30, 
                       nneigh=20, reso=0.6){
  obj <- FindNeighbors(obj, reduction="lsi",
                       k.param = nneigh, dims=use_pcs)
  obj <- FindClusters(obj, resolution = reso)
  obj
}


processATAC = function(obj, n_pcs=30, 
                       use_pcs=2:n_pcs, min.cutoff="q0",
                       cluster=FALSE, 
                       nneigh=20, reso=0.6){
  library("Signac")
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = min.cutoff)
  obj <- RunSVD(
    object = obj,
    assay = 'peaks',
    n = n_pcs,
    reduction.key = 'LSI_',
    reduction.name = 'lsi'
  )
  obj <- RunUMAP(object = obj, reduction = 'lsi', dims = use_pcs)
  if (cluster){

    obj = clusterATAC(obj, use_pcs = use_pcs, nneigh = nneigh, reso = reso)
  }
  obj
}


makeSimpleGeneActivity = function(peak_cnt, fname_gtf, upstream=1500, ...){

  chroms = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10",
             "Chr11","Chr12","Chr13","Chr14","Chr15","Chr1620","Chr17","Chr18","Chr19")
  ga_raw <- CreateGeneActivityMatrix(peak.matrix = peak_cnt, 
                                     annotation.file = fname_gtf,
                                     seq.levels = chroms,
                                     upstream = upstr, 
                                     verbose = TRUE,...)
  ga_raw
}


orderStr = function(fct, as.char=T){
  # if is.character()
  tmp = as.numeric(as.character(fct))
  lvs = as.character(c(min(tmp): max(tmp)))
  if (as.char){
    fct = factor(as.character(fct), levels = lvs)
  }else {
    fct = factor(tmp, levels = lvs)
  }

  fct
}


peakMetaFromNames = function(pknames, sep1=":", sep2="-", sep_new="_"){
  tmp = str_split(pknames, sep1, simplify = TRUE)
  chrs = tmp[, 1]
  sites = str_split(tmp[, 2], sep2, simplify = TRUE)
  peaknames1 <- paste(chrs, sites[,1], sites[,2], sep=sep_new)
  head(peaknames1)
  # site_name chromosome bp1 bp2
  peak_meta = data.frame(site_name = peaknames1, 
                         chromosome = chrs, 
                         bp1 = as.numeric(sites[,1]), 
                         bp2 = as.numeric(sites[,2]), 
                         row.names = peaknames1) 
  peak_meta
}



normalize_log1first = function(mat, s = 300){
  mat = log1p(mat)
  coln = colnames(mat)
  csums = Matrix::colSums(mat)
  if (is.null(s)){
    s = median(csums)
    message(sprintf("Using median of the column-sums as the scale: %.2f", s))
  }
  sf = s / csums
  mat = mat %*% Matrix::sparseMatrix(seq(sf), seq(sf), x=sf)
  colnames(mat) = coln
  mat
}

NormalizeLog = function(obj, assay_name="RNA", 
                       new_assay_name = assay_name, s=300){
  if (new_assay_name == assay_name){
    obj = SetAssayData(obj, slot = "data", 
                       new.data = normalize_log1first(
                         GetAssayData(obj, assay = assay_name, slot = "counts"),
                         s=s), 
                       assay = assay_name)
  }else{
    new_assay = GetAssay(obj, assay = assay_name)
    new_assay = SetAssayData(new_assay, slot = "data", 
                             new.data = normalize_log1first(
                               GetAssayData(new_assay, slot = "counts"),
                               s=s))
    obj[[new_assay_name]] = new_assay
    message(sprintf("Setting a new assay: %s", new_assay_name))
    
  }
  
  obj
}


# ==============================================================================

# KNN classifier
labelByKNN = function(obj, key_knn = "harmony", col_trans="refined_group",
                      col_split = "tech", ref_groups = c("scRNA-seq"),
                      k=21){
  # obj: merged object
  library("class")
  
  ### L2-normalization of feature (e.g. PC) space.
  X = Reductions(obj, key_knn)@cell.embeddings %>% 
    apply(1, function(x){x / norm(x, "2")}) %>% 
    t()
  
  y = obj@meta.data[[col_trans]]
  is.ref = obj@meta.data[[col_split]] %in% ref_groups
  X_train = X[is.ref, ]
  X_test = X[! is.ref, ]
  y_train = factor(y[is.ref])
  
  print(table(y_train))
  
  ### train a knn-voting classifier
  message("begin KNN classification")
  transfered_lbs = knn(train = X_train, test = X_test, cl = y_train, 
                       k=k, prob=T)
  
  y[! is.ref] = as.character(transfered_lbs) # necessary to be transformed 
  print(table(transfered_lbs))
  
  res_knn = list(all_labels = y, 
                 transfered_lbs= transfered_lbs, probs = attr(transfered_lbs, "prob"))
  return(res_knn)
}



# =======================[ plotting functioons ]========================

histPredictScores = function(scores, score_cut = 0.4,
                             fname="hist_predict_scores_tmp.pdf",
                             width = 4, height = 3,
                             main = "prediction scores"){
  pdf(fname, width = width, height = height)
  hist(
    scores,
    xlab="prediction score",
    col="lightblue",
    xlim=c(0, 1),
    main=main
  )
  abline(v=score_cut, col="red", lwd=2, lty=2)
  dev.off()
}

plotContinMat = function(mat, 
                         filename="continency_mat_tmp.pdf",
                         width=5, height=5,
                         norm_row=T, 
                         cluster_cols = F, 
                         cluster_rows = F,...){
  if (norm_row){
    colnm = colnames(mat)
    mat = mat %>% apply(1, function(x){x / sum(x)}) %>% t()
  }
  col_ord = order(apply(mat, 1, which.max))
  colnames(mat) = colnm

  pheatmap::pheatmap(mat[col_ord, ], 
                     cluster_cols = cluster_cols, 
                     cluster_rows = cluster_rows,
                     filename = filename,
                     width=width, height = height, ...)

}





















