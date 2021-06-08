load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")
library(ggplot2)
library(gridExtra)
library(grid)
# library(ggrepel)
# library(fgsea)

# Load data
###########

# Select experiment & condition
exp<- "4082"
cond<- "Sync_BAY_DMSO"
res_proc<- res_diff_expr_PP[[exp]][[cond]]  

# SQ sites?
res_SQ<- SQ_sites[[exp]]

# General description
#####################

# Unique sites
length(unique(rownames(res_proc))) # 14509

# Unique proteins
length(unique(gsub("_.*","",rownames(res_proc)))) # 4085

# Type of sites
sites<- gsub(".*_","",rownames(res_proc))
table(substr(sites,1,1))
# S    T    Y 
# 12776  1676    57 

# How many unique peptides?
# length(unique(prot$`Annotated Sequence`)) # 13,517

# Volcano
#########

# Get DE
sites_DP<- get_DE(res_proc$logFC,res_proc$pvalue,rownames(res_proc),th_logFC=0.3,th_logP= -log10(0.05),curve=0.1)
length(sites_DP$DE) # 580 sites
length(unique(gsub("_.*","",sites_DP$DE))) # 399 proteins
length(sites_DP$up) # 528
length(unique(gsub("_.*","",sites_DP$up))) # 360 proteins
length(sites_DP$down) # 52
length(unique(gsub("_.*","",sites_DP$down))) # 43 proteins

# Sites in discussion and/or labelled
genes_SQ_down<- intersect(rownames(res_SQ)[res_SQ$isSQ], sites_DP$down)
genes_sel<- unique(c(genes_SQ_down, "ATR_S435","ATR_T1989","ATM_S2996","FANCD2_T716", "E2F3_S166", "E2F3_S172", "SLBP_S110", "TOPBP1_S888", "MCM3_S535", "TP53BP1_S1104", "FOXM1_S717", "FOXM1_S730", "DCK_S74", "LIG1_S91"))
res_proc[sort(genes_sel),]

# logFC       pvalue        FDR           d
# AEBP2_S167    -0.917154444 0.0001746158 0.00000000  10.8162042
# ATM_S2996      1.021270303 0.0076208216 0.07763975  -3.8931618
# ATR_S435      -0.425826571 0.0007575298 0.01570681   7.4632690
# ATR_T1989     -1.048618338 0.0041828865 0.05010438   4.5549616
# BCLAF1_S259   -0.412339298 0.0048855193 0.05566219   4.3792193
# CASC3_S10     -0.333837614 0.0040308774 0.05010438   4.5965098
# DCK_S74       -0.516238440 0.0032819629 0.04400978   4.8874125
# E2F3_S166     -0.470451951 0.0001458405 0.00000000  11.5564346
# E2F3_S172     -0.470451951 0.0001458405 0.00000000  11.5564346
# FANCD2_T716   -1.230070355 0.0035768144 0.04895105   4.7545393
# FOXM1_S717     0.445993844 0.0135294989 0.12183908  -3.3308007
# FOXM1_S730     0.445993844 0.0135294989 0.12183908  -3.3308007
# LIG1_S91      -0.005602599 0.9356539734 0.99875450   0.1169554
# MCM3_S535      2.249209128 0.0008015714 0.01762115  -7.3690540
# RECQL4_S27    -0.541334811 0.0289177752 0.21389539   2.6780122
# RPA2_S174     -0.700078082 0.0079621614 0.08153846   3.8478629
# SLBP_S110     -0.622422564 0.0048748363 0.05566219   4.3821897
# TOPBP1_S888   -0.372125921 0.0069176718 0.07491857   3.9983971
# TP53BP1_S1104  0.955042675 0.0000381832 0.00000000 -18.1053553
# UTP14A_S445   -0.583256455 0.0005418016 0.01092896   8.1333318

p_volc<- plot_volcano(DE_results = res_proc, gene=genes_sel, isProt=T, useHyperbolicTH = T, p_cu = 0.05, logFC_cu = 0.3, curve = 0.1, labelAll = T, labelCol="black", plotTH = F)
p_volc<- p_volc +
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7)
  ) +
  scale_x_continuous(name = "log2(Fold Change)", limits = c(-1.5,2.5)) +
  scale_y_continuous(name = "-log10(P)")

# Label ATM/ATR in red
idx_sel<- grep("ATR_|ATM_",rownames(res_proc))
df_sel<- cbind(res_proc[idx_sel,],gene_id=rownames(res_proc[idx_sel,]))
colnames(df_sel)[1]<- "log2FoldChange"
p_volc<- p_volc + geom_point(data=df_sel, colour="red", size=2) # this adds a red point

# Label diff PP SQ sites as well?
idx_sel<- which(rownames(res_proc)%in%intersect(rownames(res_SQ)[res_SQ$isSQ], sites_DP$DE))
df_sel<- cbind(res_proc[idx_sel,],gene_id=rownames(res_proc[idx_sel,]))
colnames(df_sel)[1]<- "log2FoldChange"
p_volc<- p_volc + geom_point(data=df_sel, colour="red", size=3, shape=21) # this adds a red point

# GSEA
######

# Select genes
genes_up<- unique(gsub("_.*","",sites_DP$up))
genes_down<- unique(gsub("_.*","",sites_DP$down))
genes_all<- unique(gsub("_.*","",rownames(res_proc)))

# GSEA
GSEA_up<- do_GSEA2(genes_retrieved = genes_up, genes_all = genes_all, GSEA_db = geneset_ls$CP_ls, min_genes = 4, isList = T)
GSEA_down<- do_GSEA2(genes_retrieved = genes_down, genes_all = genes_all, GSEA_db = geneset_ls$CP_ls, min_genes = 4, isList = T)

# Focus on selection
pws_to_plot<-c(
  "REACTOME_DNA_REPAIR",
  "REACTOME_FANCONI_ANEMIA_PATHWAY",
  "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
  "REACTOME_CELL_CYCLE",
  "PID_ATR_PATHWAY",
  "PID_ATM_PATHWAY",
  "WP_ATR_SIGNALING",
  "WP_ATM_SIGNALING_PATHWAY",
  "WP_DNA_IRDOUBLE_STRAND_BREAKS_DSBS_AND_CELLULAR_RESPONSE_VIA_ATM",
  "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR"
)
GSEA_df<- data.frame(
  logq= c(-log10(as.numeric(GSEA_down[pws_to_plot,"q"])),-log10(as.numeric(GSEA_up[pws_to_plot,"q"]))),
  n= c(GSEA_down[pws_to_plot,"n_genes_pw_pos"],GSEA_up[pws_to_plot,"n_genes_pw_pos"]),
  prop= c(GSEA_down[pws_to_plot,"prop_pos"],GSEA_up[pws_to_plot,"prop_pos"]),
  OR= c(GSEA_down[pws_to_plot,"OR"],GSEA_up[pws_to_plot,"OR"]),
  genes= c(GSEA_down[pws_to_plot,"genes"],GSEA_up[pws_to_plot,"genes"]),
  cond=rep(c("hypo","hyper"),each=length(pws_to_plot)),
  pw= factor(rep(pws_to_plot,2),levels=rev(rownames(GSEA_down[rownames(GSEA_down)%in%pws_to_plot,])))
)

# Plot
p_GSEA<- ggplot(data=GSEA_df, aes(y=logq, x=pw, fill=cond)) +
  # geom_bar(stat="identity", fill="blue") +
  geom_bar(stat="identity", position = position_dodge()) +
  coord_flip(ylim=c(0, 5)) + 
  # ylim(0,5) + # Don't do this, bars get removed!, previous option
  xlab("") +
  ylab("-log10(Padj)") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text.x = element_text(size=7),  
    axis.text.y = element_text(size=7),
    axis.title = element_text(size=7),
    legend.position="top",
    legend.title = element_blank()
  )

# Kninjnebyrg set of annotated DDR genes?
knij <- readxl::read_excel("downloads/pub/knijnenburg_2018/NIHMS962902-supplement-2.xlsx", skip = 3)
knij_genes<- knij$`Gene Symbol`

intersect(knij_genes, genes_down) # 9/43 = 20.9%
# [1] "RPA2"    "ATR"     "DCLRE1A" "FANCD2"  "RBBP8"   "RECQL4"  "TOPBP1"  "UNG"     "WRN"    length(genes_sel) # 43
isDo<- genes_all%in%genes_down
isKnij<- genes_all%in%knij_genes
f_t<- table(isDo, isKnij)
f_prop<- prop.table(f_t,1)[,"TRUE"]
f_prop
# FALSE       TRUE 
# 0.03092528 0.20930233 
f_prop["TRUE"]/f_prop["FALSE"] # 6.768
fisher.test(f_t) # 7.35e-06; OR=8.28

# Enrichment SQ after ATRi?
############################

# Create table
PP_all<- rownames(res_proc)
DQ_df<- data.frame(isDP=PP_all%in%sites_DP$DE, isDown=PP_all%in%sites_DP$down, isUp=PP_all%in%sites_DP$up, isSQ=res_SQ[PP_all,"isSQ"])

# n DQ for DP sites?
DQ_DP_t<- tapply(DQ_df$isSQ,DQ_df$isDP,"mean")
DQ_DP_t
# FALSE       TRUE 
# 0.03273745 0.10172414 
DQ_DP_t["TRUE"]/DQ_DP_t["FALSE"]
# 3.11

# Pie charts
for(i in 1:3){
  
  # condition
  c<- c("isDown","isUp","Other")[i]
  
  # Table
  if(c=="Other") DQ_t<- table(!DQ_df$isDP,DQ_df$isSQ)
  else DQ_t<- table(DQ_df[,c],DQ_df$isSQ)
  data<- data.frame(n=DQ_t["TRUE",c("TRUE","FALSE")])
  data$prop<- data$n/sum(data$n)
  data$group<- c("ST/Q", "Other")
  data$group_label<- paste0(data$group, " (n=", data$n, "; ",signif(100*data$prop,3), "%)")
  data$ypos<- cumsum(data$prop)- 0.5*data$prop
  
  # p value?
  ft<- fisher.test(DQ_t)
  cat(i, " ", ft$p.value,"\n")
  # 1   0.002283675 
  # 2   2.225388e-11 
  # 3   2.199581e-13 
  
  # Numbers
  n_DP<- data["TRUE","n"]
  n_tot<- sum(data$n)
  prop<- 100*round(n_DP/n_tot,3)
  
  # Plot piechart
  p<- ggplot(data, aes(x="", y=prop, fill=group_label)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(
      title = c("HypoPP", "HyperPP", "No change")[i],
      caption = paste0(prop,"% (",n_DP,"/",n_tot,")")
      ) + 
    theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5, size=8),
      plot.caption = element_text(hjust= 0.5, size=7),
      legend.title = element_blank(),
      legend.position = "none"
      ) + 
    scale_fill_manual(values=c("#999999", "#D55E00"))
  assign(paste0("p",i), p)
}

p_pie<- grid.arrange(p1,p2,p3,nrow=1)

# Plot
###########

p<- grid.arrange(p_volc, p_pie, p_GSEA, widths=c(1,1), heights=c(1,1,3,2,2),layout_matrix=rbind(
 c(NA,NA),
 c(NA,2),
 c(1,NA),
 c(3,3),
 c(NA,NA)
))
ggsave("results/figs/manuscript_PP_fig.pdf",  p, width = 178, height = 265, units = "mm")

# Save GSEA results in excel
############################
WriteXLS::WriteXLS(c("GSEA_down","GSEA_up"),ExcelFileName = "results/tables/manuscript_GSEA_table.xlsx",SheetNames = c(paste0("HypoPP proteins (",length(genes_down),")"),paste0("HyperPP proteins (",length(genes_up),")")),row.names = T, col.names = T)


