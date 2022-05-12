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

DATADIR = "E:/Users/xyliu/data003/amph/ATAC"

### result directory
resdir = sprintf("%s/QC", DATADIR)
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)

# ========================================================

ATACStages
assay_name = "peaks"

snames_atac = c("blastula", "G3", "G6", "N1", "N3", "L0")
for (sn_atac in snames_atac){
{
### 0.1: raw peak counts ======
peak_cnt = Read10X_h5(sprintf("%s/data_peak/%s/filtered_peak_bc_matrix.h5", 
                              DATADIR, sn_atac))
meta_atac = read.csv(sprintf("%s/data_peak/%s/singlecell.csv", 
                             DATADIR, sn_atac), row.names = 1)

obj_atac0 <- CreateSeuratObject(
  counts = peak_cnt,
  assay = assay_name,
  project = 'ATAC',
  min.cells = 1,
  meta.data = meta_atac
)
obj_atac0$tech = "scATAC-seq"
obj_atac0
# str(obj_atac0@meta.data)
}

#==================[ peak matrix summary ]==================
{
summ0 = makeMatSummary(obj_atac0[[assay_name]]@counts, 
                      row_names = c("n_peaks_raw", "n_counts_raw"))
histMatSummary(summ0, fn_pdf = sprintf("%s/rawPmat_summary_%s.pdf", figdir, sn_atac))
print(summ0$summdf)
}

# =====================[ filtering ]======================
QCList_min_peaks = list(
  blastula = 600,
  G3 = 1500, # 2500 in v0
  YC = 1000,
  N1 = 5000,
  N15 = 1500,
  N3 = 1500,
  L2 = 1500
)
QCList_max_peaks = list(
  blastula = 3500,
  G3 = 1e4,
  YC = 1.5e4,
  N1 = 2e4,
  N15 = 1e4,
  N3 = 1e4,
  L2 = 1.5e4
)
{
obj_atac = subset(obj_atac0, 
                  nFeature_peaks > QCList_min_peaks[[sn_atac]] &
                    nFeature_peaks < QCList_max_peaks[[sn_atac]]
                  )
summ = makeMatSummary(obj_atac[[assay_name]]@counts, 
                      row_names = c("n_peaks", "n_counts"))
histMatSummary(summ, fn_pdf = sprintf("%s/pmat_summary_%s.pdf", figdir, sn_atac))
print(summ$summdf)

write.csv(rbind(summ0$summdf, summ$summdf), 
          sprintf("%s/pmat_summary_%s.csv", resdir, sn_atac))
obj_atac

}
   
# ======================[ process ]=========================
n_pcs = 50

{
# source("RFunxATAC.R")
obj_atac = processATAC(obj_atac, min.cutoff = "q0", 
                       n_pcs = n_pcs, use_pcs = 2:n_pcs,
                       cluster=F)

# DimPlot(obj_atac)
ggsave(sprintf("%s/depth_corr_%s.pdf", figdir, sn_atac),
       DepthCor(obj_atac),
       width = 6, height = 4)

obj_atac = clusterATAC(obj_atac, nneigh = 30, use_pcs = 2:n_pcs, reso=0.6)

### plot clusters and depth
obj_atac$nCount_peaks_log10 = log10(obj_atac$nCount_peaks)
n_cells = length(obj_atac$orig.ident)
cowplot::plot_grid(
  DimPlot(obj_atac, group.by = "ident") + 
    ggtitle(sprintf("%s (%d cells)", sn_atac, n_cells)),
  FeaturePlot(obj_atac, features = sprintf("nCount_%s_log10", assay_name))
) -> p1
p1
ggsave(sprintf("%s/umap_pre_%s.pdf", figdir, sn_atac), 
       p1, width = 9, height = 4)

fn_rds = sprintf("%s/processed_%s.rds", resdir, sn_atac)
saveRDS(obj_atac, fn_rds)
print(fn_rds)
}

saveRDS(obj_atac[["peaks"]]@counts,
        sprintf("%s/pmat_%s.rds", resdir, sn_atac))

}
#################################################################

for (sn_atac in ATACStages){
  fn_rds = sprintf("%s/processed_%s.rds", resdir, sn_atac)
  obj_atac = readRDS(fn_rds)
  
  pmat = obj_atac[["peaks"]]@counts
  str(pmat)
  fn_pmat_rds = sprintf("%s/pmat_%s.rds", resdir, sn_atac)
  saveRDS(pmat, fn_pmat_rds)
  print(fn_pmat_rds)
}


