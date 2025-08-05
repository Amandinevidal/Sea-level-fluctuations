#### ---------------------------Informations--------------------------####
#### Date : 2025-01-14
#### Author : Amandine Vidal-Hosteng
#### Encoding : UTF-8
#### Email : amandine.vidal-hosteng@univ-tlse3.fr
#### 
#### - This script calculates the statistical gamma and the sackin index 
#### for phylogenies generated with the archipelago model, the values
#### are saved in two files "sim1.irun_trees_statistic_gamma.txt" and 
#### "sim1.irun_trees_statistic_sackin.txt" (irun is the replicate number)
#### - Then it computes the average over the replicates of simulations per 
#### time step, the values are saved in the file "mean_statistics.txt"
#### - These calculations are parallelized, and 3 cpus are used (see line 87) !!!!!!
#### --------------------------------------------------------------------#

source("parameters.r")
library(foreach)
library(doParallel)
library(ape)
library(apTreeshape)
library(phytools)

#### 1. Statistics ####

# Function Sackin index
parallelization_sackin <- function(table_extant_phylo, cluster) {
  tryCatch({
    doParallel::registerDoParallel(cluster)
  }, error = function(e) {
    stop(sprintf("ERROR: Cluster register failed for sackin index. Message: %s\n", e$message))
  })
  
  res <- foreach::foreach(t = seq(1, max(table_extant_phylo$time), by = 100), 
                          .combine = 'rbind', 
                          .packages = c('ape', 'apTreeshape')) %dopar% {
                            tryCatch({
                              phylo_extant_w <- ape::read.tree(text = paste(table_extant_phylo[table_extant_phylo$time == t, 2], sep = ""))
                              phylo_list <- list(phylo_extant_w)
                              vec <- data.frame(time = t, sackin = NA)
                              
                              for (i in 1:length(phylo_list)) {
                                if (!is.null(phylo_list[[i]]) && phylo_list[[i]]$Nnode > 1) {
                                  vec$sackin <- apTreeshape::sackin(apTreeshape::as.treeshape(phylo_list[[i]]), norm = "yule")
                                }
                              }
                              return(vec)
                            }, error = function(e) {
                              cat(sprintf("ERROR: Sackin calcul fails t=%d. Message: %s\n", t, e$message))
                              return(data.frame(time = t, sackin = NA))
                            })
                          }
  return(res)
}

# Function gamma statistic
parallelization_gamma <- function(table_extant_phylo, cluster) {
  tryCatch({
    doParallel::registerDoParallel(cluster)
  }, error = function(e) {
    stop(sprintf("ERROR: Cluster register failed for gamma statistic. Message: %s\n", e$message))
  })
  
  res <- foreach::foreach(t = seq(1, max(table_extant_phylo$time), by = 100), 
                          .combine = 'rbind', 
                          .packages = c('ape', 'phytools')) %dopar% {
                            tryCatch({
                              phylo_extant_w <- ape::read.tree(text = paste(table_extant_phylo[table_extant_phylo$time == t, 2], sep = ""))
                              phylo_list <- list(phylo_extant_w)
                              vec <- data.frame(time = t, gamma = NA)
                              
                              for (i in 1:length(phylo_list)) {
                                if (!is.null(phylo_list[[i]])) {
                                  vec$gamma <- phytools::ltt(phylo_list[[i]], plot = FALSE)$gamma
                                }
                              }
                              return(vec)
                            }, error = function(e) {
                              cat(sprintf("ERROR: Gamma statistic calcul fails t=%d. Message: %s\n", t, e$message))
                              return(data.frame(time = t, gamma = NA))
                            })
                          }
  return(res)
}

# Create cluster
ncpus <- 3
cluster <- tryCatch(
  parallel::makeCluster(ncpus),
  error = function(e) {
    stop(sprintf("ERROR: Unable to create cluster. Message: %s\n", e$message))
  }
)

# Verify package 
if (!requireNamespace("doParallel", quietly = TRUE)) {
  stop("Missing 'doParallel' package.")
}

# Test cluster
# tryCatch({
#   test <- parallel::clusterCall(cluster, function() Sys.info())
#   print(test)
# }, error = function(e) {
#   stop(sprintf("ERROR: Problem with cluster nodes. Message: %s\n", e$message))
# })

for (irun in 0:(nrun - 1)) {
  
  # Create results files
  file.create(paste("results/sim1.", irun, "_trees_statistic_gamma.txt", sep = ""))
  file.create(paste("results/sim1.", irun, "_trees_statistic_sackin.txt", sep = ""))
  
  # Load data
  extant_phylo <- read.table(paste("results/sim1.", irun, "_insular_tree.txt", sep = ""))
  colnames(extant_phylo) <- c("time", "extant_spc_withoutconti")
  
  # Simulation time
  time <- max(extant_phylo[, 1])
  
  # Stat
  gamma_res <- parallelization_gamma(extant_phylo, cluster)
  sackin_res <- parallelization_sackin(extant_phylo, cluster)
  
  # Save results
  write.table(gamma_res, file = paste("results/sim1.", irun, "_trees_statistic_gamma.txt", sep = ""), 
              col.names = F, row.names = F, append = T)
  write.table(sackin_res, file = paste("results/sim1.", irun, "_trees_statistic_sackin.txt", sep = ""), 
              col.names = F, row.names = F, append = T)
  
  rm(extant_phylo, gamma_res, sackin_res)
  gc()
}

parallel::stopCluster(cluster)

#### 2. Mean statistics ####

#### Mean/Confidence-interval function ####
mean_ci <- function(vec){ # calcul mean and confidence interval for each line (time)
  result <- c(NA,NA,NA)
  if(!is.null(vec)){
    mean <- mean(vec,na.rm=T) # mean
    sd <- sd(vec,na.rm=T) # standard deviation
    se <- sd/(sqrt(length(vec))-1) # standard error
    ci <- 1.96*se # confidence interval
    result[1] <- mean # save for each time
    result[2] <- mean+ci # save upper values ci
    result[3] <- mean-ci # save lower values ci
  }
  return(result) # return a table of mean and ci values for each time
}

# create results file
file.create(paste("results/mean_statistics.txt",sep=""))

# list results file
list_files_gamma <- list.files(path=paste0("results/"),pattern=paste0("_gamma.txt"))
list_files_sackin <- list.files(path=paste0("results/"),pattern=paste0("_sackin.txt"))

# list of df
list_df_gamma <- vector("list", length(list_files_gamma))
list_df_sackin <- vector("list",length(list_files_sackin))
for (i in seq_along(list_files_gamma)) {
  file <- list_files_gamma[i]
  data <- read.table(paste0("results/",file),col.names = c("t","phylo"))
  list_df_gamma[[i]] <- data
  file <- list_files_sackin[i]
  data <- read.table(paste0("results/",file),col.names = c("t","phylo"))
  list_df_sackin[[i]] <- data
}

# find simulation time of the shortest replicate
max_t <- min(sapply(list_df_gamma, function(df) max(df$t, na.rm = TRUE)))

for(i in seq(1,min(max_t),by=100)) {
  vec_gamma <- sapply(list_df_gamma, function(df) {df[df$t==i,]$phylo})
  vec_sackin <- sapply(list_df_sackin, function(df) {df[df$t==i,]$phylo})
  res_gamma <- mean_ci(vec_gamma)
  res_sackin <- mean_ci(vec_sackin)
  res <- data.frame(t=i,
                    mg=res_gamma[1],
                    cug=res_gamma[2],
                    clg=res_gamma[3],
                    ms=res_sackin[1],
                    cus=res_sackin[2],
                    cls=res_sackin[3])
  write.table(res, file=paste("results/mean_statistics.txt", sep=""),col.names=F,row.names=F,append=T)
  
}







