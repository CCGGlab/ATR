install.packages("shinyBS")
install.packages("corrplot")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("plotly") 
install.packages("htmlwidgets") # Need to update for plotly, issues with DT have bene reported!
install.packages("qgraph") 
# visnetwork
install.packages("visNetwork")
install.packages("geomnet")

# install.packages("dqshiny")  # Currently not available in CRAN, use following workaround for now
install.packages("remotes")
remotes::install_github("daqana/dqshiny")

install.packages("UniprotR") # TO get protein sequence from uniprot (tried to derive SQ sites); not working (dependencies)
install.packages("protr") 

install.packages("survminer")

BiocManager::install("DEGreport") # Geom_cor
BiocManager::install("DEP") 

install.packages("ggseqlogo")
