load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")
library(ggplot2)
library(gridExtra)
library(grid)
library(DEGreport)
library(ggrepel)
library(fgsea)
library(gplots)

# Volcano
############

# Same genes as CLs
genes_to_label<- sort(c("BRCA1", "BRCA2", "E2F1", "E2F2", "E2F8", "TOP2A", "MCM4", "MCM5", "MCM6", "MCM10", "MKI67", "CCNB1", "CCNB2", "CCNA2", "FOXM1", "BRIP1", "BARD1", "FANCD2", "FANCI", "FANCB", "CDKN1A","FAS", "TP53INP1", "GDF15", "NOTCH1"))
genes_to_label_MGI<- MGI_to_HGNC[MGI_to_HGNC$HGNC.symbol%in%genes_to_label,"MGI.symbol"]

# ALK & ALKAL
genes_DE<- list()
for(i in 1:2){
  
  model<- c("ALK","ALKAL")[i]
  res_proc<- as.data.frame(res_diff_expr_mice[[paste0("DE_",model,"_ATRi_Ctrl")]])
  p_tmp<- plot_volcano(DE_results = res_proc, gene = genes_to_label_MGI, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5, plotTH=F)
  p_tmp<- p_tmp + 
    ggtitle(label = paste0(c("Alk-F1178S;Th-MYCN","Rosa26_Alkal2;Th-MYCN")[i])) +
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
    scale_x_continuous(name = "log2(Fold Change)", limits = c(-15,15)) +
    scale_y_continuous(name = "-log10(Padj)")
  assign(paste0("p_volc_",model),p_tmp)
  # n DE?
  res_proc<- na.omit(res_proc)
  genes_DE[[i]]<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$MGI, th_logFC = 1.5, th_logP = 2, curve = 0.5)
}

sapply(genes_DE, function(x) sapply(x,"length"))
# up    569 4195
# down  603 2220
# DE   1172 6415

# fGSEA + RS plots main pathways
##################################

pws<- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_P53_PATHWAY"
)

for(model in c("ALK","ALKAL")){
  
  res_proc<- as.data.frame(res_diff_expr_mice[[paste0("DE_",model,"_ATRi_Ctrl")]])
  GSEA_res<- do_fGSEA(geneset_ls[["Ha_ls"]], get_GSEA_stat(res_proc, isMice = T, MGI_to_HGNC))
  GSEA_res$isLabel<- NA
  GSEA_res$isLabel[GSEA_res$padj< 1e-20]<- T
  assign(paste0("fGSEA_res_",model),GSEA_res)
  
  p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway)) +
    ggtitle(label = paste0(c("Alk-F1178S;Th-MYCN","Rosa26_Alkal2;Th-MYCN")[i])) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Normalized Enrichment Score") + 
    ylab("-log10(Padj)") + 
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
    ) +
    geom_point(size=2)
  
  # Label
  p<- p + 
    geom_text_repel(data = subset(GSEA_res, isLabel),
                    aes(label = pathway, fontface=3),
                    size = 2,
                    colour="black",
                    box.padding   = 0.5, 
                    point.padding = 0.5,
                    segment.color = 'grey50')
  
  assign(paste0("p_fGSEA_",model),p)
  
  # RS
  stat<- get_GSEA_stat(res_proc, isMice = T, MGI_to_HGNC = MGI_to_HGNC)
  for(j in 1:length(pws)){
    p_pw<- format(signif(GSEA_res[GSEA_res$pathway==pws[j],"padj"],3),scientific=T)
    p<- plotEnrichment(geneset_ls[["Ha_ls"]][[pws[j]]],stat) +
      ggtitle(paste0(pws[j],"\n(Padj=",p_pw,")")) + 
      geom_line(size=1.5, col="green") +
      theme(
        plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),  
        axis.title = element_text(size=7),
      ) +
      scale_x_continuous(name = "Rank") +
      scale_y_continuous(name = "Enrichment Score", limits = c(-1,1))
    assign(paste0("p_RS_",model,j), p)
  }
}

# Save fGSEA plot for ALKAL
p_tmp<- grid.arrange(p_volc_ALK, p_volc_ALKAL, p_fGSEA_ALKAL, ncol=3)
ggsave(filename = "results/figs/manuscript_mice.pdf", plot = p_tmp, width = 178, height = 60, units = "mm")

# Keep fGSEA results for supplementare table with all GSEA results
save(fGSEA_res_ALK, fGSEA_res_ALKAL, file = "results/data/fGSEA_Ha_mice.RData")

# Correlate RNA to mice
#######################

# Get CL data
res_tmp<- res_diff_expr$GE$'48'
cl<- res_tmp$log2FoldChange
names(cl)<- rownames(res_tmp)

# HGNC-MGI conversion table
MGI_to_HGNC<- MGI_to_HGNC[!duplicated(MGI_to_HGNC$MGI.symbol),]
rownames(MGI_to_HGNC)<- MGI_to_HGNC$MGI.symbol

# Correlate with mice data
for(model in c("ALK","ALKAL")){

  # Get Mice data
  res_tmp<- as.data.frame(res_diff_expr_mice[[paste0("DE_",model,"_ATRi_Ctrl")]])
  # Convert to HGNC
  res_tmp$HGNC<- MGI_to_HGNC[rownames(res_tmp),"HGNC.symbol"]
  res_tmp<- res_tmp[!is.na(res_tmp$HGN),]
  # Get logFC
  mm<- res_tmp$log2FoldChange
  names(mm)<- res_tmp$HGNC

  # Common genes CL-mice
  common_genes<- intersect(names(mm),names(cl))
  
  # Correlate
  for(j in 1:length(pws)){
    pw<- pws[j]
    genes_sel<- geneset_ls$Ha_ls[[pw]]
    cor_df<- data.frame(mice=mm[genes_sel], CL=cl[genes_sel])
    # Plot
    p <- ggplot(cor_df, aes(mice, CL)) +
      geom_point() + 
      geom_smooth(method = lm) + # Add regression line
      xlab("Log2(Fold Change) Mice NB") +
      ylab("Log2(Fold Change) CLB-GE (48h)") +
      geom_cor(method = "pearson",cex=3) +
      theme(
        plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),  
        axis.title = element_text(size=7),
      )
    assign(paste0("p_cor_",model,j), p)
  }
}

# Suppl figure with RS & correlations
p_ALK<- grid.arrange(p_RS_ALK1, p_RS_ALK2, p_RS_ALK3, p_cor_ALK1, p_cor_ALK2, p_cor_ALK3, ncol=3, nrow=2, top="Alk-F1178S;Th-MYCN")
p_ALKAL<- grid.arrange(p_RS_ALKAL1, p_RS_ALKAL2, p_RS_ALKAL3, p_cor_ALKAL1, p_cor_ALKAL2, p_cor_ALKAL3, ncol=3, nrow=2, top="Rosa26_Alkal2;Th-MYCN")
p_tmp<- grid.arrange(p_ALK,p_ALKAL,nrow=2)
ggsave(filename = "results/figs/manuscript_mice_GSEA_suppl.pdf", device = cairo_pdf, plot = p_tmp, width = 178, height = 265, units = "mm")

# Heatmap G2M & E2F
#########################

# Identify pathways
# for(pw in c("HALLMARK_E2F_TARGETS","HALLMARK_INFLAMMATORY_RESPONSE")){
for(pw in c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")){
  genes_pw<- geneset_ls$Ha_ls[[pw]]
  genes_pw_MGI<- MGI_to_HGNC[MGI_to_HGNC$HGNC.symbol%in%genes_pw,"MGI.symbol"]
  
  cat(pw, " ",length(genes_pw_MGI),"\n")
  # HALLMARK_E2F_TARGETS   196 
  # HALLMARK_G2M_CHECKPOINT   195 

  # Take log(x+1), tried rlogs as well, nothing different
  data<- data.matrix(na.omit(log2(normalized_counts_mice[genes_pw_MGI,]+1)))
  
  # Samples ALK & ALKAL, Atri vs control
  samples_ALK<- c(rownames(sample_info_mice[!is.na(sample_info_mice$protID)&sample_info_mice$condition=="ALK",]),rownames(sample_info_mice[sample_info_mice$condition=="ALK_ATRi",]))
  samples_ALKAL<- c(rownames(sample_info_mice[!is.na(sample_info_mice$protID)&sample_info_mice$condition=="clay",]),rownames(sample_info_mice[sample_info_mice$condition=="clay_ATRi",]))
  
  # DE 
  res_proc_ALK<- res_diff_expr_mice$DE_ALK_ATRi_Ctrl[genes_pw_MGI,]
  res_proc_ALKAL<- res_diff_expr_mice$DE_ALKAL_ATRi_Ctrl[genes_pw_MGI,]
  
  # Order counts
  data_ALK<- data[order(res_proc_ALK[rownames(data),"stat"]),samples_ALK]
  data_ALKAL<- data[order(res_proc_ALKAL[rownames(data),"stat"]),samples_ALKAL]
  
  # Colors
  cond<- c("ALK_ATRi", "clay_ATRi", "clay", "ALK")
  cond_cols<- rainbow(length(cond))
  names(cond_cols)<- cond
  
  # Plot
  svglite::svglite(paste0("results/figs/mice_hm_",pw,".svg"),height=28)
  heatmap.2(data_ALKAL,Rowv = F, dendrogram = "col", col=bluered(75),ColSideColors=  cond_cols[sample_info_mice[colnames(data_ALKAL),"condition"]],symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="row",density.info='none', cexRow=0.5, cexCol=0.8,key.title = "logFC")
  dev.off()
}

