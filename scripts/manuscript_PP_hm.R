# Load
########
load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")
library(gplots)

# Select sites
################
sites<- c(
   'AEBP2_S167',
   'ATR_S435',
   'ATR_T1989',
   'ATRX_S978',
   'CASC3_S10',
   'DCK_S74',
   'E2F3_S166',
   'E2F3_S172',
   'EXO1_S746',
   'FANCD2_S1418',
   'FANCD2_T716',
   'FANCI_S730',
   'NBN_S615',
   'RPA2_S174',
   'SLBP_S110',
   'TOPBP1_S888',
   'UTP14A_S445',
   'WRN_S426'
  )

# Get PP counts
res_proc<- normalized_counts_PP[["4082"]]
res_proc<- res_proc[rownames(res_proc)%in%sort(sites),c(9:15)]
SQ_matrix<- SQ_sites[["4082"]]
isSQ<- rownames(res_proc)%in%rownames(SQ_matrix[SQ_matrix$isSQ,])

pdf("results/figs/PP_hm.pdf")
site_col<- rep("black", nrow(res_proc))
site_col[isSQ]<- "red"
heatmap.2(res_proc,Rowv = T, col=bluered(75), colRow = site_col, symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="row",density.info='none', cexRow=1.5, cexCol=1.75)
dev.off()
