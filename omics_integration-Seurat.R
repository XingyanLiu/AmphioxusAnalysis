library("stringr")
library("dplyr")
library("Matrix")
library("ggplot2")
library("GenomicRanges")
library("Seurat")
library("cowplot")
# library("cicero")
library("Signac")
library("harmony")


setwd("D:/Users/xyliu/003")
source("RFunxATAC.R")


# ========================================================

DATADIR = "E:/Users/xyliu/data003/amph/ATAC"
# DATADIR = "E:/others/003_atac"
datadir_ga = sprintf("%s/GA/%s", DATADIR, "0625-default")

# ========================================================

### parameters for gene activity matrix
ga_type = c("simple", "cicero")[1]
norm_mtd = c("norm", "log1f")[1]
upstr = 1500          # only used when `ga_tupe == "simple"`
recalculate_ga = F

mtd_trans = sprintf("anchor_%s_%s%s", ga_type, norm_mtd, "")

### result directory
resdir = sprintf("%s/results/%s-%s", DATADIR, "0627", mtd_trans)
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)

print(resdir)
# ========================================================
# source("RFunxATAC.R")
NPCs_trans = list(
  E6 = 25, # 10 (old)
  E7 = 30, # 15 (old)
  E8 = 30, # 20
  E9 = 30,
  E10 = 30,
  E11 = 30,
  E12 = 30,
  E13 = 40,
  E14 = 40
)
stages = names(StageMap)#[-1]
stages = c( "E13", "E14")

for (stg in stages){

  sn_rna = StageSname(stg)
  sn_atac = StageMap[[stg]]
  n_pcs = NPCs_trans[[stg]]

  if (stg == "E6"){sn_rna = paste0(sn_rna, "_sub1")}
  message(sprintf("processing stage %s & %s", stg, sn_atac))

  # ==== RNA-seq data ====
  ### way 1: something wrong...
  obj_rna = readRDS(sprintf("%s/RNA/obj_rna_%s.rds", DATADIR, sn_rna))
  hvgs = VariableFeatures(obj_rna, assay = "RNA")

  # str(obj_rna@meta.data)

  # ========================================================
  #     ATAC-seq: 0. preprocess
  # ========================================================


  #### 0.2: load processed Seurat object =====
  obj_atac = readRDS(sprintf("%s/QC/processed_%s.rds", DATADIR, sn_atac))

  # ========================================================
  #     ATAC-seq: 1. Gene activity
  # ========================================================

  ### Gene Activities ======
  if (ga_type == "simple"){

    fname_ga_simple = sprintf("%s/GA/simple_ga_raw_%s_up%d.rds",
                              DATADIR, sn_atac, upstr)
    if (recalculate_ga | !file.exists(fname_ga_simple)){
      fname_gtf = sprintf("%s/genes_cellranger_new.gtf", DATADIR)
      ga_raw = makeSimpleGeneActivity(peak_cnt = obj_atac[['peaks']]@counts,
                                      fname_gtf, upstr)

      message(sprintf("saving (simple) gene activaties (upstream = %d)...\n\t%s",
                      upstr, fname_ga_simple))
      saveRDS(ga_raw, fname_ga_simple)
    }else{
      ga_raw = readRDS(fname_ga_simple)
    }


  }else if (ga_type == "cicero"){

    fname_ga = sprintf("%s/cicero_ga_raw_%s.rds", datadir_ga, sn_atac)
    message(sprintf("using prepared Cicero gene activity from:\n\t%s", fname_ga))
    ga_raw = readRDS(fname_ga)
  }
  assay_name_ga = "GA"
  obj_atac[[assay_name_ga]] = CreateAssayObject(counts = ga_raw)

  if (norm_mtd == "norm"){
    # scale_fct = 1000
    scale_fct = median(obj_atac@meta.data[[paste0("nCount_", assay_name_ga)]])
    message(sprintf("NOTE: using standard normalization, scale.factor = %.1f", scale_fct))
    obj_atac = NormalizeData(obj_atac, assay=assay_name_ga,
                             scale.factor = scale_fct)

  }else{
    obj_atac = NormalizeLog(obj_atac, assay_name_ga, s = NULL)
  }


  # ========================================================
  #     Integration : anchors, labels
  # ========================================================
  # n_pcs = 30 # has been set before
  name_trans = "refined_group" #"stg_leiden"
  rna_assay = "RNA"
  only_polyT = TRUE
  use_pcs = 2: min(n_pcs, 30) # use the same number of PCs as scRNA-seq data

  rna_assay_impute = rna_assay
  k.filter = ifelse(ga_type == "simple", 300, 200)
  k.weight = 30 # 50 by default

  if (only_polyT){
    message("NOTE: using only `polyT` to transfer labels")
    obj_rna_ref = subset(obj_rna, subset= primer == 'polyT')
  }else{
    obj_rna_ref = obj_rna
  }

  if (stg == "E6"){
    message("(E6) using only cells with `pseudo_batch == 1`")
    obj_rna_ref = subset(obj_rna, subset = pseudo_batch == "1")
  }


  transfer.anchors <- FindTransferAnchors(
    reference = obj_rna_ref,
    query = obj_atac,
    features = hvgs,
    reference.assay = rna_assay,
    query.assay = "GA",
    reduction = "cca",
    npcs = n_pcs, # 30 by default
    k.filter = k.filter, # 200 by default
  )

  #### scRNA-seq based annotation ====
  celltype.predictions <- TransferData(
    anchorset = transfer.anchors,
    refdata = as.factor(obj_rna_ref@meta.data[[name_trans]]),
    weight.reduction = obj_atac[["lsi"]],
    dims = use_pcs,
    k.weight = k.weight
    )

  obj_atac <- AddMetaData(obj_atac, metadata = celltype.predictions)
  # to make the colors match
  obj_atac$predicted.id <- factor(obj_atac$predicted.id,
                                  levels = levels(obj_rna@meta.data[[name_trans]]))



  ##### group counts =======
  group_cnts = rbind(table(obj_rna@meta.data[[name_trans]]),
                     table(obj_atac$predicted.id))
  rownames(group_cnts) = c("scRNA-seq", "scATAC-seq")
  write.csv(group_cnts, sprintf("%s/group_counts_%s_%s.csv", resdir, sn_rna, sn_atac))
  print(group_cnts)

  # ============= save labels ===============
  # mtd_trans = sprintf("knn_%s_%s", ga_type, norm_mtd)
  dflb = data.frame(self_cluster = obj_atac$seurat_clusters)
  dflb[[mtd_trans]] = obj_atac$predicted.id
  write.csv(dflb,
            sprintf("%s/transferred_labels_%s_%s.csv", resdir, stg, sn_atac))

  contin_mat = table(dflb)
  colnames(contin_mat) = paste(stg, colnames(contin_mat), sep="_")

  # source("RFunxATAC.R")
  plotContinMat(contin_mat,
                filename=sprintf("%s/contin_mat_%s_%s.pdf", figdir, stg, sn_atac))

  # ==========================================
  #  prediction scores
  score_cut = 0.4

  histPredictScores(obj_atac$prediction.score.max,
                    score_cut = score_cut,
                    fname = sprintf("%s/predict_scores_%s_%s.pdf", figdir, stg, sn_atac),
                    main = sprintf("prediction scores (%s -> %s)", stg, sn_atac))


  table(obj_atac$prediction.score.max > score_cut)


  # ====================Vis======================

  ####### Vis (separately) =====
  p1 <- DimPlot(obj_atac, reduction = "umap", cols = COLORS,
                group.by = "predicted.id", label = TRUE, repel = TRUE) +
    ggtitle("scATAC-seq cells") +
    NoLegend()
  p2 <- DimPlot(obj_rna, group.by = name_trans,  cols = COLORS,
                label = TRUE, repel = TRUE) +
    ggtitle("scRNA-seq cells") +
    NoLegend()
  p12 = plot_grid(p1, p2)
  ggsave(filename = sprintf("%s/vis_sep_%s_%s.pdf", figdir, stg, sn_atac),
         plot = p12,
         width=10, height = 5)


  # ========================================================
  #     Integration : transfer RNA counts (inputation)
  # ========================================================

  refdata <- GetAssayData(
    object = obj_rna_ref,
    assay = rna_assay_impute,
    slot = "data"
  );

  imputation <- TransferData(
    anchorset = transfer.anchors,
    refdata = refdata,
    weight.reduction = obj_atac[["lsi"]],
    dims = use_pcs,
    k.weight = k.weight,
  )

  obj_atac[[rna_assay_impute]] = imputation

  rm(imputation) # free memory
  rm(refdata) # free memory
  # obj_atac.filtered <- subset(obj_atac, subset = prediction.score.max > score_cut)
  # integrated = IntegrateData(transfer.anchors, dims = 1:n_pcs, k.weight = 50)


  ## ===========================================
  #        co-embedding (optional)
  #=============================================

  coembed <- merge(x = obj_rna, y = obj_atac, merge.data=TRUE)
  DefaultAssay(coembed) <- rna_assay_impute
  coembed$primer <- ifelse(is.na(coembed$primer), "ATAC", coembed$primer)


  coembed$RNAcluster = ifelse(is.na(coembed@meta.data[[name_trans]]),
                              coembed$predicted.id,
                              coembed@meta.data[[name_trans]])
  coembed$RNAcluster <- orderStr(coembed$RNAcluster)
  # table(coembed$RNAcluster)
  # levels(coembed$RNAcluster)

  coembed <- ScaleData(coembed, features = hvgs, split.by = "primer") # do.scale = FALSE
  coembed <- RunPCA(coembed, features = hvgs, verbose = FALSE)

  mindist = ifelse(stg == "E14", 0.5, 0.3)

  # using Harmony to correct the PC space
  # library("harmony")
  coembed <- RunHarmony(coembed, group.by.vars = "tech", dims.use = 1:n_pcs) # need pakage `Harmony`
  # embedding
  coembed <- RunUMAP(coembed, dim = 1:n_pcs, reduction = "harmony",
                     min.dist=mindist)


  ## ===========================================
  #        co-embedding - Visualization
  #=============================================


  # Vis - 1 (after Harmony)
  rd1 = "umap"
  plt1 <- DimPlot(coembed, reduction = rd1, group.by = "tech") +
    ggtitle(sprintf("Co-embedded data (stage: %s/%s)", sn_atac, stg))
  plt2 <- DimPlot(coembed, reduction = rd1, cols=COLORS,
                  group.by = "RNAcluster",
                  label = TRUE, repel = TRUE) +
    ggtitle("colored by clusters")
  ggsave(filename = sprintf("%s/coembed_%s_%s_%s.pdf", figdir, rd1, stg, sn_atac),
         plot = plot_grid(plt1, plt2),
         width=10, height = 5)

  # Vis - 2
  plt <- DimPlot(coembed,
                 reduction = rd1,
                 split.by = "tech",
                 group.by = "RNAcluster",
                 cols = COLORS,
                 label = TRUE,
                 repel = TRUE) +
    ggtitle(sprintf("Transfered labels from scRNA-seq (%s -> %s)", stg, sn_atac))#+ NoLegend()
  ggsave(filename = sprintf("%s/coembed_trans_split_%s_%s.pdf", figdir, stg, sn_atac),
         plot = plt,
         width=10, height = 5)


  #######################################
  # Vis - 000 (before Harmony; optional)
  coembed <- RunUMAP(coembed, dim = 1:n_pcs, reduction = "pca",
                     reduction.name = "umap0",
                     min.dist=mindist)
  rd = "umap0"
  plt1 <- DimPlot(coembed, reduction = rd, group.by = "tech") +
    ggtitle(sprintf("Co-embedded data (stage: %s/%s)", sn_atac, stg))

  plt2 <- DimPlot(coembed, reduction = rd, cols=COLORS,
                  group.by = "RNAcluster", label = TRUE, repel = TRUE) +
    ggtitle("colored by clusters")
  ggsave(filename = sprintf("%s/_coembed_%s_%s_%s.pdf", figdir, rd, stg, sn_atac),
         plot = plot_grid(plt1, plt2),
         width=10, height = 5)



  ## ===========================================
  #        Save results
  #=============================================

  saveRDS(coembed, file=sprintf("%s/coembed_%s_%s.rds", resdir, sn_rna, sn_atac))
  saveRDS(obj_atac, file=sprintf("%s/atac_%s.rds", resdir, sn_atac))

}




