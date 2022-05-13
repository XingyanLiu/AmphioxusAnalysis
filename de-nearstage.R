# differential expression analysis - current vs the previous stage
library("stringr")
library("dplyr")
library("Matrix")
library("ggplot2")
library("Seurat")

setwd("D:/Users/xyliu/003")
source("RFunx.R")
source("RFunxATAC.R")

# ========================================================
DATADIR_main = "E:/Users/xyliu/data003/amph"
datadir.m = file.path(DATADIR_main, "afterQC_formal")

resdir = sprintf("%s/res-srt/1116", DATADIR_main)
figdir = file.path(resdir, "figs")
dir.create(figdir, recursive = T)

### ==============[ loading Annos (gene annotations) ]=============
Annos = readRDS(sprintf("%s/Annos.rds", DATADIR_main))
str(Annos)

### ===================================
###     loading data (merged) etc.
### ===================================
.i = 2
sname = c("merged", "merged_E1-E15")[.i]
obj = readRDS(sprintf("%s/obj_%s.rds", datadir.m, sname))
levels(obj$stage_name) <-  c(
  "2cell", "4cell", "8cell", "32cell", "256cell", "B",
  "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0", "L2"
)
obj@commands

obj$stage_name %>% table
obj = NormalizeData(obj, scale.factor = 1000)
# saveRDS(obj, sprintf("%s/obj_%s_lognorm.rds", datadir.m, sname))
groupby_de = "stage_name"

tests = c('wilcox', 'MAST', 't')
test_use = tests[2]

for(test_use in tests[-1]){
  print(test_use)
  for(i in seq(StageOrd)){
    stg1 = StageOrd[i]
    if (stg1 == "G5"){break}
    if (stg1 == "L2"){break}
    stg2 = StageOrd[i + 1]
    message(sprintf('DE between stage %s and %s before', stg2, stg1))
  
    stg1.before = paste0(stg1, "_before")
  obj$tmp.ident = ifelse(obj$stage_name %in% StageOrd[1: i], 
                         yes=stg1.before,
                         no = as.character(obj$stage_name))
  Idents(obj) <- "tmp.ident"
  obj$tmp.ident %>% table()
  markers_all = FindMarkers(obj, 
                            ident.1 = stg2, 
                            ident.2 = stg1.before,
                            only.pos = T, min.pct = 0.025,
                            logfc.threshold = 0.15,# 2 ^ 0.25 = 1.189207
                            test.use = test_use)
  markers_all$gene = rownames(markers_all)
  marker_ids = gsub("-", "_", markers_all$gene)
  rownames(markers_all) = marker_ids
  length(marker_ids) %>% message()
  
  #================[ arrange ]================
  ### adding annotation columns (from Annos)
  markers_all$gene_short_name = Annos[marker_ids, "gene_short_name"]
  markers_all$tt = gsub("\n", "&", Annos[marker_ids, "tt"])
  str(markers_all)
  
  ### saving marker tables (mktb)
  ptag = sprintf("%s_%s(%s)", stg1, stg2, test_use)
  fn_mktb = sprintf("%s/DEGtable_%s.csv", resdir, ptag)
  write.csv(markers_all, fn_mktb)
  print(fn_mktb)
  
  #==============[ select genes for dotplot ]=============
  ### selecting top-significant markers for visualization
  ntop = 50
  markertb_plt = markers_all #%>% group_by(cluster) 
  if (test_use == 'roc'){
    markertb_plt = markertb_plt %>% top_n(ntop, power)
  }else{
    markertb_plt = markertb_plt %>% top_n(ntop, -p_val_adj) %>%
      top_n(ntop, avg_logFC)
  }
  markers_plt = unique(markertb_plt$gene)

  #==============[ dotplot ]=============
  source("RFunx.R")
  # xtlabs = Annos[marker_ids, "gene_short_name"]
  xtlabs0 = getGeneNames(markers_plt)
  xtlabs = shortMultiNames(xtlabs0)
  xtlabs = shortMultiNames(xtlabs, sep=' ')
  
  pdot = WrapperDotPlot(obj, markers_plt, groupby = groupby_de,
                        gene_labs = xtlabs,
                        transpose = T,
                        dir_fig = figdir, sname = ptag)
  pdot = WrapperDotPlot(obj, markers_plt, groupby = groupby_de,
                        gene_labs = xtlabs,
                        transpose = F,
                        dir_fig = figdir, sname = paste0(ptag, '-T'))
  
}

}








