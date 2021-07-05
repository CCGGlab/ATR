# Load 
######
load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")

# Get sequences around PP sites
####################################

# Select experiment & condition
exp<- "4082"
cond<- "Sync_BAY_DMSO"
res_proc<- res_diff_expr_PP[[exp]][[cond]]  
res_SQ<- SQ_sites[[exp]]
sites_DP<- get_DE(res_proc$logFC,res_proc$pvalue,rownames(res_proc),th_logFC=0.3,th_logP= -log10(0.05),curve=0.1)
sites_down<- sites_DP$down
sites_up<- sites_DP$up

# Get aa sequences around PP site
library(protr)
res_SQ$PP<- 0
res_SQ$PP[rownames(res_SQ)%in%sites_down]<- -1
res_SQ$PP[rownames(res_SQ)%in%sites_up]<- 1

res_SQ$seq5<- NA
for(i in 1:nrow(res_SQ)){
  cat(i," ")
  prot_seq_tmp<- getUniProt(res_SQ$uniProt[i])[[1]]
  res_SQ$seq5[i]<- substr(prot_seq_tmp,res_SQ$aa_pos[i]-5,res_SQ$aa_pos[i]+5)
}

# Save
saveRDS(res_SQ,"temp/PP_motifs.rds")

# Process: Add Xs to make sure positions & peptide lengths match
################################################################

PP_motifs <- readRDS("temp/PP_motifs.rds")

PP_motifs$seq5_corr<- PP_motifs$seq5
PP_motifs$aa_pos_corr<- PP_motifs$aa_pos
for(i in 1:nrow(PP_motifs)){
  l<- nchar(PP_motifs$seq5[i])
  p<- PP_motifs$aa_pos[i]
  if(p<6){
    PP_motifs$seq5_corr[i]<- paste0(paste0(rep("X",6-p),collapse = ""),PP_motifs$seq5_corr[i])
    PP_motifs$aa_pos_corr[i]<- 6
  }
  if(p>=6&l<11) PP_motifs$seq5_corr[i]<- paste0(PP_motifs$seq5_corr[i],paste0(rep("X",11-l),collapse = ""))
  if(substr(PP_motifs$seq5_corr[i],6,6)!=PP_motifs$aa[i]) PP_motifs$seq5_corr[i]<- NA # Remove motifs with wrong sites
}
PP_motifs<- PP_motifs[!is.na(PP_motifs$seq5_corr),]

# Create logo's for dePP & hyperPP
###################################

# Get dePP & hyperPP
PP_motifs_noChange<- PP_motifs[PP_motifs$PP==(0),"seq5_corr"]
PP_motifs_dePP<- PP_motifs[PP_motifs$PP==(-1),"seq5_corr"]
PP_motifs_hyperPP<- PP_motifs[PP_motifs$PP==(+1),"seq5_corr"]

# Use position-specific p-values calculated with kpLogo
# http://kplogo.wi.mit.edu/submit.html?
# Export files fr kpLogo
cat(PP_motifs_dePP,sep="\n",file="temp/kpLogo/ATR_dePP_motifs.txt")
cat(PP_motifs_hyperPP,sep="\n",file="temp/kpLogo/ATR_hyperPP_motifs.txt")
cat(PP_motifs_noChange,sep="\n",file="temp/kpLogo/ATR_noChange_motifs.txt") # Use this as background
