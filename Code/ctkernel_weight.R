



#!/usr/bin/env Rscript

# Get the ctkernel-weighted network of each interactome.
# This script uses the interactomes and their corresponding normalised CT-kernels which were computed in picasso.


args = commandArgs(trailingOnly=TRUE)

# test if the number of arguments is right: if not, return an error
if (length(args)< 3) {
  stop("3 arguments must be supplied: graph file, names file and kernel file (input file).n", call.=FALSE)
} 


# set working directory
workingDir <- "/Users/aurelio/OneDrive - Universidad de Málaga/p.BCGE_NAC/results/exp2019-04-23/NAC-BL/Results/expanded_modules_networks"
setwd(workingDir)

options(stringsAsFactors = FALSE)

# Load required libraries
source("../../../../../bin/0_loadLibraries.R")

loadpkg("dplyr") 
loadpkg("igraph")
loadpkg("reshape2")



assignKernelprob <- function(grafo, kernel.df){
  # Get the kernel probabilities from the kernel data frame
  # assign kernel probabilities to a network derived from the Kernel
  df <-get.data.frame(grafo)
  mer <- merge(df,kernel.df, by = c('from', 'to'))
  gr <- graph.data.frame(d = mer,directed = T)
  # return(subgraph.edges(graph = gr, eids= E(gr)[which(E(gr)$weight > 0.1)]))
  return(gr)
}


subnetwork_kernel = function(graphfile,namesfile,kernelfile){
  graphname = strsplit(graphfile, "_")
  graphname = graphname[[1]][1]
  cat(paste("Load the ", graphname, " CT kernel and graph \n", sep = ""))
  grafo <- read.graph(file = graphfile, format ='gml')
  nodesnames<-read.table(namesfile)
  nodenames <- nodesnames$V1
  ctdf <-read.table(kernelfile, sep='\t')
  names(ctdf)<-nodenames
  row.names(ctdf) <- nodenames
  
  cat("convert the ct kernel to matrix \n")
  ct <- as.matrix(ctdf)
  ct.mat.df <- melt(ct)
  
  cat("add  weights to the ct-kernel dataframe \n")
  ct.mat.df <- mutate(ct.mat.df, weightlog = -log(value))
  ct.mat.df$weightlog[is.na(ct.mat.df$weightlog)] <- max(ct.mat.df$weightlog, na.rm = T) ####REplace NaN with the highest distance
  ct.mat.df$weightlog[is.infinite(ct.mat.df$weightlog)] <- max(ct.mat.df$weightlog[is.finite(ct.mat.df$weightlog)], na.rm = T) ####REplace Inf with the highest distance
  ct.mat.df <- dplyr::select(ct.mat.df, Var1,Var2,weightlog) 
  names(ct.mat.df)<-c('from','to','weight')
  
  cat(paste("Save the ", graphname, " weighted interactome with CT kernel probabilities\n", sep = ""))
  interactome_CT <- assignKernelprob(grafo, kernel.df = ct.mat.df)
  write.graph(interactome_CT, file=paste(graphname,'_ct.gml', sep = ""), format = 'gml')
  
  interactome_CT.df <- get.data.frame(interactome_CT)
  write.csv(interactome_CT.df, file =paste(graphname,'_ct.csv', sep = "") ,quote = F, row.names = F)
}

graphfile = args[1]
namesfile = args[2]
kernelfile = args[3]

subnetwork_kernel(graphfile, namesfile, kernelfile)



