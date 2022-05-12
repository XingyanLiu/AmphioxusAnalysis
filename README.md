# AmphioxusAnalysis

Code for analysis of snRNA-seq and scATAC-seq of amphioxus

Code for building the developmental tree can be found at: https://github.com/XingyanLiu/stagewiseNN

| file                           | description                                                                              |
|--------------------------------|------------------------------------------------------------------------------------------|
| rna_proc_stagewise.py          | preprocessing, clustering, and visualization of the snRNA-seq data                       |
| atac_preproc.R                 | preprocessing and visualization of the scATAC-seq data                                   |
| atac2ga_cicero.R               | compute the gene activities from the scATAC-seq data                                     |
| omics_integration-Harmony.R    | integration of snRNA_seq and scATAC-seq data using Harmony, and KNN-based label transfer |
| omics_integration-Seurat.R     | integration of snRNA_seq and scATAC-seq data using Seurat                                |
| lineage_dynamics-swnn-GA.ipynb | using StagewiseNN to integrately visualize the scATAC-seq data across six stages         |
| GRN-scenic.py                  | compute the TF-target pairs using pySCENIC                                               |
| GRN-intersection.py            | compute the species-conserved TF-target pairs                                            |
| src/*                          | general utils functions (Python)                                                         |
| RFunx.R                        | general utils functions (R)                                                              |
| RFunxATAC.R                    | utils functions for analysis related to scATAC-seq data                                  |
| RFunxTreePlot.R                | utils functions for plotting the tree                                                    |
| COLOR_SETS.R                   | sss                                                                                      |


