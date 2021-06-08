# Summary table
load("data/ATR_data.RData")

# RNA
###########
RNA_df<- list()
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
    res_tmp<- as.data.frame(res_diff_expr[[CL]][[t]])[,c(1,2,5,6)]
    colnames(res_tmp)<- paste(CL,t,colnames(res_tmp),sep="_")
    RNA_df<- c(RNA_df,list(res_tmp))
  }
  }
RNA_df<- cbind(RNA_df[[1]],RNA_df[[2]],RNA_df[[3]],RNA_df[[4]])

# Proteomics
###########
P_df<- NULL
for(CL in c("BAR","GE")){
  for(t in c("24","48")){
    if(CL=="GE") exp<- "3367"
    if(CL=="BAR") exp<- "3403"
    res_tmp<- as.data.frame(res_diff_expr_prot[[exp]][[t]])[,c(1,2,3)]
    colnames(res_tmp)<- paste(CL,t,colnames(res_tmp),sep="_")
    res_tmp$id<- rownames(res_tmp)
    if(is.null(P_df)) P_df<- res_tmp
    else P_df<- merge(P_df, res_tmp, all=T,by = "id")
  }
}
rownames(P_df)<- P_df$id
P_df$id<- NULL

# Phosphoprot
###############

PP_df<- as.data.frame(res_diff_expr_PP[["4082"]][["Sync_BAY_DMSO"]])[,c(1,2,3)]
colnames(PP_df)<- paste("GE_6",colnames(PP_df),sep="_")

# RNA-Seq mice
################
res_ALK<- as.data.frame(res_diff_expr_mice$DE_ALK_ATRi_Ctrl)[,c(1,2,5,6)]
colnames(res_ALK)<- paste("ALK",colnames(res_ALK),sep="_")
res_ALKAL<- as.data.frame(res_diff_expr_mice$DE_ALKAL_ATRi_Ctrl)[,c(1,2,5,6)]
colnames(res_ALKAL)<- paste("ALKAL",colnames(res_ALKAL),sep="_")
RNA_mice_df<- cbind(res_ALK,res_ALKAL)

# Save
WriteXLS::WriteXLS(c("RNA_df","P_df","PP_df","RNA_mice_df"),"results/tables/manuscript_summary_table.xlsx",row.names = T, SheetNames = c("RNA-Seq","Proteomics","Phosphoproteomics","RNA-Seq mice"))
