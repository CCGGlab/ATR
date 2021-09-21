# Load
########
load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")
library(ggplot2)
library(gridExtra)
library(grid)
library(visNetwork)
library(geomnet)
library(igraph)
library(dplyr)
library(gplots)

# Get PP data
#############

# Select data
res_proc_PP<- res_diff_expr_PP[["BAR"]]

# Get DP genes
genes_DP_ls<- get_DE(res_proc_PP$logFC,res_proc_PP$pvalue,rownames(res_proc_PP),th_logFC=0.3,th_logP= -log10(0.05),curve=0.1)
sites_DP<- genes_DP_ls$DE
genes_DP<- unique(gsub("_.*","",sites_DP))
genes_DP_up<- unique(gsub("_.*","",genes_DP_ls$up))
genes_DP_down<- unique(gsub("_.*","",genes_DP_ls$down))

# Create graph
###############

pdf("results/figs/PP_network_new_cu_onlyDown_withP2.pdf")

# SHow both general network and zoom on hub
for(i in 1:2){
  
  # Select genes for network
  genes<- genes_DP_down
  # genes<- c(genes_DP_down,"BRCA1","FANCI","EXO1","WRN","BLM","ATRX","XRCC1")
  # if(i==2){
  #   genes<- setdiff(genes,c("EIF3B","EIF4B","CASC3","SRRM2")) # Seperate clusters
  #   genes<- setdiff(genes,c("SMARCA5","GSK3B","CEBPB","MED1","TSC2")) # Focus on hub only
  # } 
  
  # Create graph
  ppi_db_sel<- ppi_db[ppi_db$protein1%in%genes&ppi_db$protein2%in%genes,]
  if(i==2){
    res_proc_PP_sel<- res_proc_PP[gsub("_.*","",rownames(res_proc_PP))%in%genes,]
    res_proc_PP_sel<- res_proc_PP_sel[rownames(res_proc_PP_sel)%in%sites_DP,]
    res_proc_PP_sel<- data.frame(protein1=gsub("_.*","",rownames(res_proc_PP_sel)),protein2=rownames(res_proc_PP_sel),combined_score=1)
    ppi_db_sel<- rbind(ppi_db_sel,res_proc_PP_sel[res_proc_PP_sel$protein1%in%union(ppi_db_sel$protein1,ppi_db_sel$protein2),])
  }
  g_ppi <- graph_from_data_frame(ppi_db_sel, directed=FALSE)
  
  # Get nodes
  nodes<- as.character(V(g_ppi)$name)
  
  # Node colors
  nodes_sites<- nodes[grep("_",nodes)]
  logFC_norm<- round(50+50*res_proc_PP[nodes_sites,"logFC"]/2)
  names(logFC_norm)<- nodes_sites
  node_col<- bluered(100)[logFC_norm[nodes]]
  # Proteins
  node_col[is.na(node_col)]<- "grey"

  # Node shapes
  node_shape<- rep("circle",length(nodes))
  if(i==2) node_shape[nodes%in%genes]<- "sphere"

  # Node Labels
  node_label<- nodes
  node_label[node_label%in%nodes_sites]<- NA # Don't label sites
  
  # Node border colors: SQ sites
  PP_SQ<- rownames(SQ_sites[["BAR"]])[SQ_sites[["BAR"]][,"isSQ"]]
  node_border_col<- rep("black",length(nodes))
  node_border_col[nodes%in%PP_SQ]<- "red"
  # node_border_col[nodes%in%nodes_genes]<- "black"
  
  # Edge line types
  E(g_ppi)$lty<- 1
  E(g_ppi)$lty[ppi_db_sel$protein2%in%nodes_sites]<- 2
  
  plot.igraph(
    g_ppi,
    vertex.frame.color=node_border_col,
    vertex.label=node_label,
    vertex.color=node_col,
    vertex.label.color="black",
    vertex.label.cex=0.5,
    vertex.label.font=3,
    vertex.size=(3+2*log2(degree(g_ppi,mode = "in"))),
    vertex.shape=node_shape, 
    # edge.curved=T,
    layout=layout.kamada.kawai
    # layout=layout_with_fr,
  )
}
dev.off()

save(g_ppi, node_border_col, node_label, node_shape, node_col, file="results/data/PP_network.RData")

# tkplot locally

# Some data
#############

# N nodes?
length(nodes) #41
length(unique(gsub("_.*","",nodes))) #17

# Knij genes?
knij <- readxl::read_excel("downloads/pub/knijnenburg_2018/NIHMS962902-supplement-2.xlsx", skip = 3)
knij_genes<- knij$`Gene Symbol`
intersect(nodes, knij_genes)
# "TOPBP1"  "FANCD2"  "WRN"     "ATR"     "RPA2"    "DCLRE1A" "RBBP8"   "RECQL4" 
setdiff(nodes, knij_genes) # Only DTL & SLBP not!

# Check sites
res_proc_PP[nodes,]

# Col key
pdf("results/figs/PP_network_key.pdf")
plot(1:100,rep(0,100),col=bluered(100)[1:100],pch="|",cex=5,axes=F,xlab="log2FC",ylab="")
axis(1,at = c(0,100), labels = c(-2,2))
dev.off()
