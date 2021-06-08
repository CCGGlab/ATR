load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")
library(ggplot2)
library(gridExtra)
library(grid)

# Basic description: Volcano plot
####################################

# N proteins?
sapply(res_diff_expr_prot, function(x) sapply(x, function(y) nrow(y)))
# $`3336`
# 24   48 
# 6635 6635 
# 
# $`3367`
# 24   48 
# 5189 5290 
# 
# $`3403`
# 24   48 
# 6792 6791 

# BAR & GE 24-48h
# genes_to_label<- sort(c("CDKN1A","TOP2A","CCNA2","RRM2","MCM4","MKI67","MCM10","FANCI","CCNB1","MCM5","FOXM1","MCM6","E2F1","E2F2","E2F8","TP53INP1","BRIP1","FANCD2", "NOTCH1", "BRCA1", "BRCA2", "FAS"))
genes_to_label<- sort(c("BRCA1", "BRCA2", "E2F1", "E2F2", "E2F8", "TOP2A", "MCM4", "MCM5", "MCM6", "MCM10", "MKI67", "CCNB1", "CCNB2", "CCNA2", "FOXM1", "BRIP1", "BARD1", "FANCD2", "FANCI", "FANCB", "CDKN1A","FAS", "TP53INP1", "GDF15", "NOTCH1"))
genes_DE<- list()
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
    if(CL=="BAR") exp<- "3403"
    if(CL=="GE") exp<- "3367"
    res_proc<- as.data.frame(res_diff_expr_prot[[exp]][[t]])
    p_tmp<- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 0.3,curve = 0.1, plotTH=F, isProt = T, p_cu = 0.05)
      p_tmp<- p_tmp +
        ggtitle(label = paste0("CLB-",CL,": ",t,"h")) +
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
        scale_x_continuous(name = "log2(Fold Change)", limits = c(-3,3)) +
        scale_y_continuous(name = "-log10(Padj)")
      
      assign(paste("p",CL,t,sep="_"),p_tmp)
      # n DE?
      res_proc<- na.omit(res_proc)
      genes_DE[[paste(CL,t,sep="_")]]<- get_DE(logFC = res_proc$logFC, P = res_proc$pvalue, genes = rownames(res_proc), th_logFC = 0.5, th_logP = 0.05, curve = 0.1)
  }
}

# p_GE_48

# n DE genes: 
############
n_DE_BAR<- sapply(genes_DE[1:2], function(x) sapply(x, length))
n_DE_GE<- sapply(genes_DE[3:4], function(x) sapply(x, length))

n_DE_BAR
# BAR_24 BAR_48
# up       12     33
# down     14     58
# DE       26     91

n_DE_GE
# GE_24 GE_48
# up      70   259
# down   219   291
# DE     289   550

g_E2F<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_E2F_TARGETS)
g_G2M<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_G2M_CHECKPOINT)
g_TP53<- intersect(genes_DE$GE_48$up,geneset_ls$Ha_ls$HALLMARK_P53_PATHWAY)
g_MTORC1<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_MTORC1_SIGNALING)

length(intersect(genes_DE$BAR_48$down,genes_DE$BAR_24$down)) # 13/58
length(intersect(genes_DE$GE_48$down,genes_DE$GE_24$down)) # 92/291
length(intersect(genes_DE$BAR_48$up,genes_DE$BAR_24$up)) # 7/33
length(intersect(genes_DE$GE_48$up,genes_DE$GE_24$up)) # 35/259
length(intersect(genes_DE$BAR_48$up,genes_DE$GE_48$up)) # 7/259
length(intersect(genes_DE$BAR_48$down,genes_DE$GE_48$down)) # 18/291

# GSEA 48h --> E2F & G2M + RS
###############################

# BAR & GE 24-48h?
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
    if(CL=="BAR") exp<- "3403"
    if(CL=="GE") exp<- "3367"
    res_proc<- as.data.frame(res_diff_expr_prot[[exp]][[t]])
    GSEA_res<- do_fGSEA(geneset_ls[["Ha_ls"]], get_GSEA_stat(res_proc, isProt = T))
    GSEA_res$isLabel<- NA
    GSEA_res$isLabel[GSEA_res$padj< 0.05]<- T
    
    p<- ggplot(GSEA_res, aes(x = NES, y = -log10(padj), key=pathway)) +
      ggtitle(paste0("CLB-",CL," ", t,"h")) +
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
    GSEA_res$pathway_short<- gsub("HALLMARK_"," ",GSEA_res$pathway)
    p<- p + 
      geom_text_repel(data = subset(GSEA_res, isLabel),
                    aes(label = pathway_short, fontface=3),
                    size = 2,
                    colour="black",
                    box.padding   = 0.5, 
                    point.padding = 0.5,
                    segment.color = 'grey50')
    
    assign(paste("p_GSEA",CL,t,sep="_"),p)
    assign(paste("fGSEA",CL,t,sep="_"),GSEA_res)
  }
}

# Combine with RS
pws<- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_P53_PATHWAY"
)

stat<- get_GSEA_stat(res_proc, isProt = T)
for(i in 1:length(pws)){
  p_pw<- format(signif(GSEA_res[GSEA_res$pathway==pws[i],"padj"],3),scientific=T)
  pw_short<- gsub("HALLMARK_"," ",pws[i])
  p<- plotEnrichment(geneset_ls[["Ha_ls"]][[pws[i]]],stat) +
    ggtitle(paste0(pw_short)) + 
    # ggtitle(paste0(pws[i],"\n(Padj=",p_pw,")")) + 
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
    scale_y_continuous(name = "Enrichment Score")
  assign(paste0("p_RS",i), p)
}

# Create Plots
##############

# Group GSEA & RS plots
p_GSEA<- grid.arrange(p_GSEA_GE_48 + ggtitle("Preranked GSEA (Hallmark genesets)"), p_RS1, p_RS2, p_RS3,
                    widths=c(2,1),layout_matrix=rbind(c(1,2),c(1,3),c(1,4)))

ggsave("results/figs/manuscript_Prot_fGSEA.pdf",  p_GSEA, width = 120, height = 100, units = "mm")

# Volcano
p_volc<- grid.arrange(p_BAR_24, p_BAR_48, p_GE_24, p_GE_48,nrow=4,ncol=2)
ggsave("results/figs/manuscript_Prot_volc.pdf", p_volc, width = 178, height = 265, units = "mm")

