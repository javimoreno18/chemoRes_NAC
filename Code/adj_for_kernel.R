# Read the network in SIF format. Produce adjacency matrices for kernel calculations


# Load required libraries
source("../../../../../bin/0_loadLibraries.R")

loadpkg("dplyr") 
loadpkg("igraph")

library(dplyr)
library(igraph)

# set working directory
workingDir <- "/Users/aurelio/OneDrive - Universidad de Málaga/p.BCGE_NAC/results/exp2019-04-16/chemoRes_NAC/Results/ct_kernels"
subnetDir <- "/Users/aurelio/OneDrive - Universidad de Málaga/p.BCGE_NAC/results/exp2019-04-16/chemoRes_NAC/Results/expanded_modules_networks"
resultsDir = workingDir
setwd(workingDir)



gene.name2gene.id <- function(gene_name) {
        require(clusterProfiler)
        gid <-tryCatch(bitr(gene_name, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db"),error=function(e){return(bitr(gene_name, fromType="ALIAS", toType="ENTREZID", annoDb="org.Hs.eg.db"))})
        #gid <- bitr(gene_name, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
        gid <- gid$ENTREZID
        return(gid[1])
}




get_adjacency <- function(sif){
  sif_name = deparse(substitute(sif))
  cat("processing ", sif_name, " subnetwork", "\n")
  grafo = graph_from_data_frame(d = sif,directed = F)
  grafo = igraph::simplify(graph=grafo, remove.multiple = TRUE, remove.loops = TRUE,edge.attr.comb ='max')
  kk<- clusters(grafo)
  conectados <- names(kk$membership[kk$membership == 1])
  connected<- induced.subgraph(graph = grafo,vids = conectados)
  ## Write out the connected unweighted graph
  df.connected <- get.data.frame(x = connected)
  write.table(df.connected,file =  file.path(resultsDir, paste(sif_name,"_subnetwork_connected.tsv", sep = "")) , sep = '\t', quote = F, row.names = F)
  write.graph(connected, file = file.path(resultsDir, paste(sif_name,"_subnetwork_connected.gml", sep = "")), format = 'gml')
  ## Adjacency matrix for kernel calculation
  cat("computing adjacency matrix", "\n")
  adj <- as_adj(connected,type = "both",sparse = F, names = T)
  k <- adj[1,]
  adj_names <- names(k)
  cat("writing adjacency matrix to file", "\n")
  write.table(adj,file = file.path(resultsDir, paste(sif_name,"_adj.txt", sep = "")), sep = "\t",quote = F,row.names = F,col.names = F)
  zip(zipfile = file.path(resultsDir, paste(sif_name,"_adj.gz", sep = "")), files =  file.path(resultsDir, paste(sif_name,"_adj.txt", sep = "")))
  file.remove(file.path(resultsDir, paste(sif_name,"_adj.txt", sep = "")))
  write(adj_names, file = file.path(resultsDir, paste(sif_name,"_names.txt", sep = "")), ncolumns = 1)
  cat("done", "\n")
  
}


##### Load Data

ghiassian100grey.df <- read.csv(file.path(subnetDir, "Ghiassian_interactome_n=100_expanded_grey_module_network.csv"), skip = 1, header = F) 
ghiassian100grey <- ghiassian100grey.df[,2:3]
ghiassian100turq.df <- read.csv(file.path(subnetDir, "Ghiassian_interactome_n=100_expanded_turquoise_module_network.csv"), skip = 1, header = F) 
ghiassian100turq <- ghiassian100turq.df[,2:3]
ghiassian100yellow.df <- read.csv(file.path(subnetDir, "Ghiassian_interactome_n=100_expanded_yellow_module_network.csv"), skip = 1, header = F) 
ghiassian100yellow <- ghiassian100yellow.df[,2:3]

iref100grey.df <- read.csv(file.path(subnetDir, "iRef_interactome_n=100_expanded_grey_module_network.csv"), skip = 1, header = F) 
iref100grey <- iref100grey.df[,2:3]
iref100turq.df <- read.csv(file.path(subnetDir, "iRef_interactome_n=100_expanded_turquoise_module_network.csv"), skip = 1, header = F) 
iref100turq <- iref100turq.df[,2:3]
iref100yellow.df <- read.csv(file.path(subnetDir, "iRef_interactome_n=100_expanded_yellow_module_network.csv"), skip = 1, header = F) 
iref100yellow <- iref100yellow.df[,2:3]

string100grey.df <- read.csv(file.path(subnetDir, "STRINGdb_interactome_n=100_expanded_grey_module_network.csv"), skip = 1, header = F) 
string100grey <- string100grey.df[,2:3]
string100turq.df <- read.csv(file.path(subnetDir, "STRINGdb_interactome_n=100_expanded_turquoise_module_network.csv"), skip = 1, header = F) 
string100turq <- string100turq.df[,2:3]
string100yellow.df <- read.csv(file.path(subnetDir, "STRINGdb_interactome_n=100_expanded_yellow_module_network.csv"), skip = 1, header = F) 
string100yellow <- string100yellow.df[,2:3]


metabol100grey.df <- read.csv(file.path(subnetDir, "metabolism-centered_interactome_n=100_expanded_grey_module_network.csv"), skip = 1, header = F) 
metabol100grey <- metabol100grey.df[,2:3]
metabol100turq.df <- read.csv(file.path(subnetDir, "metabolism-centered_interactome_n=100_expanded_turquoise_module_network.csv"), skip = 1, header = F) 
metabol100turq <- metabol100turq.df[,2:3]
metabol100yellow.df <- read.csv(file.path(subnetDir, "metabolism-centered_interactome_n=100_expanded_yellow_module_network.csv"), skip = 1, header = F) 
metabol100yellow <- metabol100yellow.df[,2:3]

sifs = c(ghiassian100grey,ghiassian100turq,ghiassian100yellow,iref100grey,iref100turq,iref100yellow,metabol100grey,metabol100turq,metabol100yellow,string100grey,string100turq,string100yellow)

for(i in sifs){
        names(i)
}

get_adjacency(string)
get_adjacency(iref)
get_adjacency(metabol)
get_adjacency(string)
