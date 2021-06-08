load("data/ATR_data.RData")
source("scripts/functions/ATR_functions.R")

# Basic description: Volcano plot
####################################

# BAR & GE 24-48h
genes_to_label<- sort(c("BRCA1", "BRCA2", "E2F1", "E2F2", "E2F8", "TOP2A", "MCM4", "MCM5", "MCM6", "MCM10", "MKI67", "CCNB1", "CCNB2", "CCNA2", "FOXM1", "BRIP1", "BARD1", "FANCD2", "FANCI", "FANCB", "CDKN1A","FAS", "TP53INP1", "GDF15", "NOTCH1", "RAD51"))
genes_DE<- list()
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
      res_proc<- as.data.frame(res_diff_expr[[CL]][[t]])
      p_tmp<- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5,curve = 0.1, plotTH=F)
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
        scale_x_continuous(name = "log2(Fold Change)", limits = c(-7,7)) +
        scale_y_continuous(name = "-log10(Padj)")
      
      assign(paste("p",CL,t,sep="_"),p_tmp)
      # n DE?
      res_proc<- na.omit(res_proc)
      genes_DE[[paste(CL,t,sep="_")]]<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1.5, th_logP = 2, curve = 0.1)
  }
}
# p_GE_48

# n DE genes
############
n_DE_BAR<- sapply(genes_DE[1:2], function(x) sapply(x, length))
n_DE_GE<- sapply(genes_DE[3:4], function(x) sapply(x, length))

n_DE_BAR
#         BAR_24 BAR_48
# up       29     89
# down      0      7
# DE       29     96

n_DE_GE
#        GE_24 GE_48
# up      70   198
# down     3   370
# DE      73   568

g_E2F<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_E2F_TARGETS)
g_G2M<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_G2M_CHECKPOINT)
g_TP53<- intersect(genes_DE$GE_48$up,geneset_ls$Ha_ls$HALLMARK_P53_PATHWAY)
g_MTORC1<- intersect(genes_DE$GE_48$down,geneset_ls$Ha_ls$HALLMARK_MTORC1_SIGNALING)

length(intersect(genes_DE$BAR_48$up,genes_DE$BAR_24$up)) # 29(/29) from 89
length(intersect(genes_DE$GE_48$up,genes_DE$GE_24$up)) # 63(/70) from 198
length(intersect(genes_DE$BAR_48$up,genes_DE$GE_48$up)) # 45
length(intersect(genes_DE$BAR_48$down,genes_DE$GE_48$down)) # 0 (none at th in BAR)

# Dynamics for examples (including Prot)
#######################################
p_ls<- list()
for(gene in genes_to_label){
  res_proc<- as.data.frame(sapply(res_diff_expr,function(x) sapply(x,function(y) y[gene,"log2FoldChange"])))
  res_proc_P<- as.data.frame(sapply(res_diff_expr_prot[c("3367","3403")],function(x) sapply(x,function(y) y[gene,"logFC"])))
  colnames(res_proc_P)<- c("GE","BAR")
  res_gene<- rbind(
    data.frame(exp="RNA-Seq", CL="BAR", time=as.numeric(c(0,rownames(res_proc))), counts=c(0,res_proc$BAR)),
    data.frame(exp="RNA-Seq", CL="GE", time=as.numeric(c(0,rownames(res_proc))), counts=c(0,res_proc$GE)),
    data.frame(exp="Proteomics", CL="BAR", time=as.numeric(c(0,rownames(res_proc_P))), counts=c(0,res_proc_P$BAR)),
    data.frame(exp="Proteomics", CL="GE", time=as.numeric(c(0,rownames(res_proc_P))), counts=c(0,res_proc_P$GE))
  ) 
  res_gene$exp<- factor(res_gene$exp,levels=c("RNA-Seq", "Proteomics")) 
  p<- ggplot(res_gene, aes(x=time, y=counts, group=interaction(CL,exp), color=CL)) +
    geom_line(size=1, aes(linetype = exp)) +
    ggtitle(label = gene) +
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face='italic'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7)
    ) +
    scale_x_continuous(name = "log2(Fold Change)", breaks = seq(0,48,24)) +
    scale_y_continuous(name = "-log10(Padj)")
  
    p_ls[[gene]]<- p + theme(legend.position = "none")
}

# GSEA 48h --> E2F & G2M + RS
###############################

# BAR & GE 24-48h?
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
    res_proc<- as.data.frame(res_diff_expr[[CL]][[t]])
    GSEA_res<- do_fGSEA(geneset_ls[["Ha_ls"]], get_GSEA_stat(res_proc))
    GSEA_res$isLabel<- NA
    GSEA_res$isLabel[GSEA_res$padj< 1e-10]<- T
    
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

# p_GSEA_GE_48

# Combine with RS
pws<- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_P53_PATHWAY"
)

stat<- get_GSEA_stat(res_proc)
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

# TF enrichment (barplots)
###########################

# BAR & GE 24-48h?
CL<- "GE"
t<- "48"
res_proc<- as.data.frame(res_diff_expr[[CL]][[t]])

genes_sel<- genes_DE$GE_48$DE
genes_all<- rownames(res_proc)
GSEA<- do_GSEA2(genes_retrieved = genes_sel, genes_all = genes_all, GSEA_db = geneset_ls$TFT_RN, min_genes = 4, isList = T)
GSEA$pathway<- factor(rownames(GSEA),levels=rev(rownames(GSEA)))
GSEA$q<- as.numeric(GSEA$q)
GSEA_sel<- GSEA[GSEA$q<0.05,]
    
p_TFT<- ggplot(GSEA_sel, aes(x = pathway, y = -log10(q))) +
  geom_bar(stat="identity", fill="darkblue") +
  ggtitle(paste0("CLB-",CL," ", t)) +
  coord_flip(expand = F) +
  ylab("-log10(Padj)") +
  xlab("") +
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

# p_TFT

# Create Plots
##############

# Group GSEA & RS plots
p_RS<- grid.arrange(p_GSEA_GE_48, p_RS1, p_RS2, p_RS3,
                    widths=c(2,1),layout_matrix=rbind(c(1,2),c(1,3),c(1,4)))

# Group dynamics & TF GSEA
p_dyn<- grid.arrange(p_ls$E2F1,p_ls$BRCA1,
                      p_ls$FANCD2,p_ls$FANCI,
                      p_ls$CDKN1A,p_ls$TP53INP1,
                      p_ls$GDF15,
                      p_TFT,
                      nrow=5, heights=c(1,1,1,1,3), layout_matrix=rbind(
                        c(1,2),
                        c(3,4),
                        c(5,6),
                        c(7,NA),
                        c(8,8)
                      ))

# Combine with volcano
p<- grid.arrange(p_GE_48,p_RS, p_dyn, widths=c(2,1), heights=c(1,3,3,3),layout_matrix=rbind(
  c(NA,3),
  c(1,3),
  c(2,3),
  c(NA,NA)
))

# Save main plot
ggsave("results/figs/manuscript_RNA_fig.pdf",  p, width = 178, height = 265, units = "mm")

# Suppl fig with all Volcano
p<- grid.arrange(p_BAR_24, p_BAR_48, p_GE_24, p_GE_48,nrow=4,ncol=2)
ggsave("results/figs/manuscript_RNA_volcan_suppl.pdf", p, width = 178, height = 265, units = "mm")

# Save GSEA results
GSEA_GE_48_TFT<- GSEA
fGSEA_GE_48_Ha<- fGSEA_GE_48
save(fGSEA_GE_48_Ha, GSEA_GE_48_TFT, file="results/data/GSEA_RNA.RData")

