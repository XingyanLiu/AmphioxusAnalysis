library("stringr")
library("dplyr")
library("ggplot2")
library("GenomicRanges")
library("cicero")
library("Seurat")

setwd("D:/Users/xyliu/003")
source("RFunxATAC.R")

DATADIR = "E:/Users/xyliu/data003/amph/ATAC"

peakdir = file.path(DATADIR, "QC")
gadir = sprintf("%s/GA", DATADIR)

resdir = file.path(gadir, "0701-default")
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = TRUE)


#=================================[ load data ]=================================
##### 0.1 load genome informations
size.genome = read.delim(file.path(DATADIR, "size.genome"), header = F, sep='\t')

##### 0.2 load processed annotations =================
gene_anno = readRDS(sprintf("%s/gene_anno_for_cicero.rds", DATADIR))
gene_anno_sub = prepareAnno4GA(gene_anno)

# ==============================================================================
##### 0.3 load peak matrix
for(sn_atac in ATACStages){
  # if(sn_atac == "L2"){next()}
  message(sprintf("processing ATAC stage: %s", sn_atac))

fn_pmat = sprintf("%s/pmat_%s.rds", peakdir, sn_atac)
fn_input_cds = sprintf("%s/CDS/input_cds_%s.rds", DATADIR, sn_atac)
fn_cicero_cds = sprintf("%s/CDS/cicero_cds_%s.rds", DATADIR, sn_atac)

if (! file.exists(fn_input_cds)){
  message("loading filtered peak matrix")
  peakmat = readRDS(fn_pmat)
  peak_meta = peakMetaFromNames(rownames(peakmat))
  rownames(peakmat) = rownames(peak_meta)
  
  # ======= Input CDS ========
  input_cds = new_cell_data_set(peakmat, gene_metadata = peak_meta)
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]
  
  # ======= preprocessing ========
  input_cds <- processCDS4ATAC(input_cds) # source("RFunxATAC.R")
  message(sprintf("saving `input_cds` object into:\n\t%s", fn_input_cds))
  saveRDS(input_cds, fn_input_cds)
  ggsave(sprintf("%s/umap_%s.pdf", figdir, sn_atac), plot_cells(input_cds),
         width = 5, height = 5)
  
}else{
  ### reload
  message(sprintf("loading existing `input_cds` object from:\n\t%s", fn_input_cds))
  input_cds = readRDS(fn_input_cds)

}


# ===== Cicero CDS ======
if (! file.exists(fn_cicero_cds)){
  umap_coords <- reducedDims(input_cds)$UMAP
  cicero_cds <- make_cicero_cds(input_cds, k=30,
                                reduced_coordinates = umap_coords)
  saveRDS(cicero_cds, fn_cicero_cds)
}else{
  message(sprintf("loading existing `cicero_cds` object from:\n\t%s", fn_input_cds))
  cicero_cds = readRDS(fn_cicero_cds)
}

cicero_cds

# ============================[ Run Cicero ]=============================
cicero_default = T
if(cicero_default){
  ##### way 1: one-step Cicero
  conns <- run_cicero(cicero_cds, size.genome, sample_num = 100)
}else{
  ##### way 2: step-by-step
  window_size = 350000 # 500000  by default
  distance_constraint = 150000 # 250000 by default
  s = 0.8 # 0.75 by default
  
  conns = runCustomedCicero(
    cicero_cds,
    genomic_coords = size.genome,
    window_size = window_size,
    distance_constraint = distance_constraint,
    s = s,
    sample_num = 100
  )
}
  

saveRDS(conns, sprintf("%s/peak_conns_%s.rds", resdir, sn_atac))

##### visualization of connections (unnecessary) #####
# plot_connections(conns, "Chr2", 9773451, 9848598,
#                  gene_model = gene_anno, 
#                  coaccess_cutoff = .25, 
#                  connection_width = .5, 
#                  collapseTranscripts = "longest" )

# =========[ Cicero gene activity scores ]===========

# gene_anno_sub = prepareAnno4GA(gene_anno)
input_cds <- annotate_cds_by_site(input_cds, gene_anno_sub)
tail(fData(input_cds))


# ======= Generate gene activity scores ============
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
str(unnorm_ga)

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
str(cicero_gene_activities)

saveRDS(unnorm_ga, sprintf("%s/cicero_ga_raw_%s.rds", resdir, sn_atac))
saveRDS(cicero_gene_activities, sprintf("%s/cicero_ga_%s.rds", resdir, sn_atac))

summ_ga_raw = makeMatSummary(unnorm_ga)
write.csv(summ_ga_raw$summdf, sprintf("%s/cicero_ga_summary_%s.csv", resdir, sn_atac))

xg = "active genes per cell"
xc = "active counts per cell"
colr = "lightblue"
fn_hist = sprintf("%s/activity_summary_%s.pdf", figdir, sn_atac)

histMatSummary(summ_ga_raw, fn_hist, xg, xc)

}




