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
datadir_ga = sprintf("%s/GA/%s", DATADIR, "0701-default")

# ========================================================

### parameters for gene activity matrix
ga_type = c("simple", "cicero")[2]
norm_mtd = c("norm", "log1f")[1]
upstr = 1500          # only used when `ga_tupe == "simple"`
recalculate_ga = T

mtd_trans = sprintf("knn_%s_%s%s", ga_type, norm_mtd, "")

### result directory
resdir = sprintf("%s/results/%s-%s", DATADIR, "0707", mtd_trans)
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)

# ========================================================
# merged metadata (with lineage labels)
metadata_all = read.csv(file.path(DATADIR_main, "afterQC_formal", "merged_metadata.csv"),
                        row.names = 1)
metadata_all$refined_group = str_split(metadata_all$stg_groups_new, "_", simplify = T)[,2]
str(metadata_all)
meta_columns = c("stage", "primer", "leiden_new", "refined_group", "lineage", "stage_name")
write.csv(metadata_all[, meta_columns],
          sprintf("%s/metadata_main.csv", resdir))

# ========================================================
# source("RFunxATAC.R")
NPCs_trans = list(
  E6 = 30, # B
  E7 = 30, # G3
  E8 = 30, # G4
  E9 = 30, # G5
  E10 = 30,# G6
  E11 = 40,# N1
  E12 = 40,# N2
  E13 = 40,# N3
  E14 = 40 # L0
)
stages = names(StageMap)
# stages = c("E9")
stg = "E9"


for (stg in stages){

  sn_rna = StageSname(stg)
  sn_atac = StageMap[[stg]]
  n_pcs = NPCs_trans[[stg]]

  if (stg == "E6"){sn_rna = paste0(sn_rna, "_sub1")}
  message(sprintf("processing stage %s & %s", stg, sn_atac))

  # ===============[ RNA-seq data ]==================
  obj_rna = readRDS(sprintf("%s/RNA/obj_rna_%s.rds", DATADIR, sn_rna))
  hvgs = VariableFeatures(obj_rna, assay = "RNA")


  ###=== adding `stage_name` and `linead` labels ===
  meta_rna = metadata_all[colnames(obj_rna), ]
  obj_rna$refined_group = orderStr(meta_rna$refined_group) # remember orderStr()
  obj_rna$lineage = as.character(meta_rna$lineage)
  obj_rna$stage_name = as.character(meta_rna$stage_name)
  # str(obj_rna@meta.data)

  ### group --> lineage mapping
  group_lin_df = obj_rna@meta.data[, c('refined_group', 'lineage')] %>% unique()
  rownames(group_lin_df) = group_lin_df$refined_group

  # ===============[ ATAC-seq data ]==================
  #### 0.2: load processed Seurat object =====
  obj_atac = readRDS(sprintf("%s/QC/processed_%s.rds", DATADIR, sn_atac))

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
  assay_name_ga = "RNA"
  obj_atac[[assay_name_ga]] = CreateAssayObject(counts = ga_raw)

  if (norm_mtd == "norm"){
    scale_fct = median(obj_atac@meta.data[[paste0("nCount_", assay_name_ga)]])
    # scale_fct = 1000
    message(sprintf("NOTE: using standard normalization, scale.factor = %.1f", scale_fct))
    obj_atac = NormalizeData(obj_atac, assay=assay_name_ga,
                             scale.factor = scale_fct)

  }else{
    obj_atac = NormalizeLog(obj_atac, assay_name_ga, s = NULL)
  }


  # ======================================================
  #     directly integrate using Harmony
  # ======================================================
  coembed <- merge(x = obj_rna, y = obj_atac, merge.data=TRUE)
  DefaultAssay(coembed) <- "RNA"
  coembed$nCount_RNA_log10 = log10(coembed$nCount_RNA)
  coembed$primer <- ifelse(is.na(coembed$primer), "ATAC", coembed$primer)
  coembed$RNAcluster <- ifelse(is.na(coembed$refined_group), "unlabeled", coembed$refined_group)
  # coembed$lineage <- ifelse(is.na(coembed$lineage), )

  # processing data
  coembed <- ScaleData(coembed, features = hvgs, split.by = "primer") # do.scale = FALSE
  coembed <- RunPCA(coembed, features = hvgs, verbose = FALSE)


  # n_pcs = 30
  use_pcs = 1:n_pcs
  ### using Harmony to correct the PC space ### library("harmony")
  coembed <- RunHarmony(coembed, group.by.vars = "tech", dims.use = use_pcs,
                        max.iter.harmony = 15,
                        sigima = 0.2)
  # embedding
  mindist = ifelse(stg == "E14", 0.5, 0.3)
  coembed <- RunUMAP(coembed, dim = use_pcs, reduction = "harmony",
                     min.dist=mindist)

  # ================================================================
  #     Label transfer using KNN on Harmony-corrected PC space
  # ================================================================

  name_trans = "refined_group"
  key_cluster = c("refined_group", "RNAcluster")[2]
  k = 19
  res_knn = labelByKNN(coembed, key_knn = "harmony",
                       col_trans=name_trans,
                       col_split = "tech", ref_groups = c("scRNA-seq"),
                       k=k)

  ### setting labels
  coembed$RNAcluster = orderStr(res_knn$all_labels)
  obj_atac$RNAcluster = coembed$RNAcluster
  obj_atac$RNAcluster_prob = res_knn$probs
  table(obj_atac$RNAcluster)

  ### cluster labels --> lineage labels
  coembed$lineage = group_lin_df[as.character(coembed$RNAcluster), "lineage"]
  obj_atac$lineage = group_lin_df[as.character(obj_atac$RNAcluster), "lineage"]

    # ============== prediction scores ============
  score_cut = 0.4
  histPredictScores(res_knn$probs,
                    score_cut = score_cut,
                    fname = sprintf("%s/predict_scores_%s_%s.pdf", figdir, stg, sn_atac),
                    main = sprintf("prediction scores (%s -> %s)", stg, sn_atac))

  # ============= save labels ===============
  # mtd_trans = sprintf("knn_%s_%s", ga_type, norm_mtd)
  dflb = data.frame(self_cluster = obj_atac$seurat_clusters)
  dflb[[mtd_trans]] = obj_atac$RNAcluster
  write.csv(dflb,
            sprintf("%s/transferred_labels_%s_%s.csv", resdir, stg, sn_atac))

  contin_mat = table(dflb)
  colnames(contin_mat) = paste(stg, colnames(contin_mat), sep="_")

  # source("RFunxATAC.R")
  plotContinMat(contin_mat,
                filename=sprintf("%s/contin_mat_%s_%s.pdf", figdir, stg, sn_atac))


  ## ===========================================
  #         Visualization
  #=============================================

  # ============== separated UMAP ================
  p1 <- DimPlot(obj_atac, reduction = "umap", cols = COLORS,
                group.by = key_cluster, label = TRUE, repel = TRUE) +
    ggtitle("scATAC-seq cells") +
    NoLegend()
  p2 <- DimPlot(obj_rna, group.by = name_trans,  cols = COLORS,
                label = TRUE, repel = TRUE) +
    ggtitle("scRNA-seq cells") +
    NoLegend()
  ggsave(filename = sprintf("%s/vis_sep_%s_%s.pdf", figdir, stg, sn_atac),
         plot = plot_grid(p1, p2),
         width=10, height = 5)


  # ============== co-embedded UMAP ================

  # Vis - 1 (after Harmony)
  rd1 = "umap"
  plt1 <- DimPlot(coembed, reduction = rd1, group.by = "tech") +
    ggtitle(sprintf("Co-embedded data (stage: %s/%s)", sn_atac, stg))
  plt2 <- DimPlot(coembed, reduction = rd1, cols=COLORS,
                  group.by = key_cluster,
                  label = TRUE, repel = TRUE) +
    ggtitle("colored by clusters")
  ggsave(filename = sprintf("%s/coembed_%s_%s_%s.pdf", figdir, rd1, stg, sn_atac),
         plot = plot_grid(plt1, plt2),
         width=10, height = 5)


  # Vis - 2
  plt <- DimPlot(coembed,
                 reduction = rd1,
                 split.by = "tech",
                 group.by = key_cluster,
                 cols = COLORS,
                 label = TRUE,
                 repel = TRUE) +
    ggtitle(sprintf("Transfered labels from scRNA-seq (%s -> %s)", stg, sn_atac))#+ NoLegend()
  ggsave(filename = sprintf("%s/coembed_trans_split_%s_%s.pdf", figdir, stg, sn_atac),
         plot = plt,
         width=10, height = 5)



  # Vis - 3: depth (log10_RNA_counts; optional)
  plt_depth = FeaturePlot(coembed, reduction = "umap",
                          "nCount_RNA_log10", split.by = "tech")
  ggsave(filename = sprintf("%s/coembed_depth_split_%s_%s.pdf", figdir, stg, sn_atac),
         plot = plt_depth,
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
                  group.by = key_cluster, label = TRUE, repel = TRUE) +
    ggtitle("colored by clusters")
  ggsave(filename = sprintf("%s/_coembed_%s_%s_%s.pdf", figdir, rd, stg, sn_atac),
         plot = plot_grid(plt1, plt2),
         width=10, height = 5)


  ## ===========================================
  #        Save results
  #=============================================


  saveRDS(coembed, file=sprintf("%s/coembed_%s_%s.rds", resdir, sn_rna, sn_atac))
  saveRDS(obj_atac, file=sprintf("%s/atac_%s.rds", resdir, sn_atac))

  ### with new labels added (lineage, modified cluster labels)
  saveRDS(obj_rna, file=sprintf("%s/obj_rna_%s.rds", resdir, sn_rna))

  print(resdir)

}













