# library
library(ape)
library(TreeSim)
library(paleotree)
library(ggplot2)
library(dplyr)
library(cowplot)

# Functions ####

stat <- function(phylo) {
  if(!is.null(phylo) && Ntip(phylo) > 2) {
    sackin <- apTreeshape::sackin(apTreeshape::as.treeshape(phylo), norm = "yule")
    gamma <- phytools::ltt(phylo, plot = FALSE)$gamma
  } else {
    sackin <- NA
    gamma <- NA
  }
  return(c(gamma,sackin))
}

search_phylo <- function(t,data){
  tree <- data[data$t==t,]$tree
  if(tree != "NULL") {
    phylo <- read.tree(text=tree)
    depth <- node.depth.edgelength(phylo)
    phylo$root.time <- max(depth)
  } else {
    phylo <- NULL
  }
  return(phylo)
}

inf_round <- function(x) {
  return(floor(x / 100) * 100 +1)
}

# Code: past statistics with present phylogenies ####

ifelse(!dir.exists("results/plotstatpast"), dir.create("results/plotstatpast"), FALSE)

for (irun in 0:(nrun-1)) {
  
  # data
  file.create(paste0("results/sim1.", irun, "_trees_statistics_past_estimates.txt", sep = ""))
  data <- read.table(paste0("results/sim1.",irun,"_insular_tree.txt"))
  colnames(data) <- c("t","tree")
  
  for (tf in seq(1,max(data$t),by=1000)) { # loop over simulation time
    
    phylo <- search_phylo(tf,data)
    stat_tf <- stat(phylo)
    
    if(!is.null(phylo) && phylo$root.time>1000) {
      
      vector_t <- seq(1,phylo$root.time - 1000,by=1000)
      
      for (t in vector_t) { # loop over phylogeny lifetime
        phylo_cut <- timeSliceTree(
          phylo,
          sliceTime = t,
          plot = FALSE,
          drop.extinct = FALSE
        )
        stat_est <- stat(phylo_cut)
        t_obs <- inf_round(phylo$root.time-t) # arrondir au nombre rond inférieur parce que sinon dans data ça marche pas
        phylo_obs <- search_phylo(t_obs,data)
        stat_obs <- stat(phylo_obs)
        save <- data.frame(tf = tf, t = t, gamma_est = stat_est[1], gamma_obs = stat_obs[1], sackin_est = stat_est[2], sackin_obs = stat_obs[2], error_gamma = round(stat_obs[1] - stat_est[1],2), error_sackin = round(stat_obs[2] - stat_est[2],2), irun = irun)
        
        write.table(save, 
                    file = paste("results/sim1.", irun, "_trees_statistics_past_estimates.txt", sep = ""), 
                    col.names = F, row.names = F, append = T)
        # # plot
        # layout(1:3)
        # if(!is.null(phylo)) {
        #   plot(ladderize(phylo), show.tip.label = FALSE)
        #   axisPhylo()
        #   title(paste("Phylo t",tf,"- Gamma:",round(stat_tf[1],2),"- Sackin:",round(stat_tf[2],2)))
        # } else {
        #   plot(1,1)
        # }
        # if(!is.null(phylo_cut)) {
        #   plot(ladderize(phylo_cut), show.tip.label = FALSE)
        #   axisPhylo()
        #   title(paste("Phylo cut",t,"- Gamma:",round(stat_est[1],2),"- Sackin:",round(stat_est[2],2)))
        # } else {
        #   plot(1,1)
        # }
        # if(!is.null(phylo_obs)) {
        #   plot(ladderize(phylo_obs), show.tip.label = FALSE)
        #   axisPhylo()
        #   title(paste("Phylo osb",phylo$root.time-t+1,"- Gamma:",round(stat_obs[1],2),"- Sackin:",round(stat_obs[2],2)))
        # } else {
        #   plot(1,1)
        # }
        # layout(1)
        
      }
      
    }
    
  }
  
  cat("irun",irun,"done\n")
}

# Code: replicate average ####

# create results file
file.create(paste("results/mean_statistics_past_estimates.txt",sep=""))

# list results file
list_files <- list.files(path=paste0("results/"),pattern=paste0("_trees_statistics_past_estimates.txt"))

# list of df
df <- c()
for (i in seq_along(list_files)) {
  file <- list_files[i]
  data <- read.table(paste0("results/",file),col.names = c("tf","t","gamma_est","gamma_obs","sackin_est","sackin_obs","error_gamma","error_sackin","irun"))
  df <- rbind(df,data)
}

df$t_cor <- df$tf - df$t

data_est <- df %>%
  filter(!is.na(gamma_est)) %>%
  filter(!is.na(sackin_est)) %>%
  group_by(tf,t_cor) %>%
  summarize(
    mg = mean(gamma_est, na.rm = TRUE),
    sd_gamma_est = sd(gamma_est, na.rm = TRUE),
    n_gamma_est = n(),
    ms = mean(sackin_est, na.rm = TRUE),
    sd_sackin_est = sd(sackin_est, na.rm = TRUE),
    n_sackin_est = n()
  ) %>%
  mutate(se_gamma_est = sd_gamma_est / sqrt(n_gamma_est),
         clg = mg - qt(1 - (0.05 / 2), n_gamma_est - 1) * se_gamma_est,
         cug = mg + qt(1 - (0.05 / 2), n_gamma_est - 1) * se_gamma_est,
         se_sackin_est = sd_sackin_est / sqrt(n_sackin_est),
         cls = ms - qt(1 - (0.05 / 2), n_sackin_est - 1) * se_sackin_est,
         cus = ms + qt(1 - (0.05 / 2), n_sackin_est - 1) * se_sackin_est)

write.table(data_est, 
            file = paste("results/mean_statistics_past_estimates.txt", sep = ""), 
            col.names = F, row.names = F, append = T)

data_obs <- read.table("results/mean_statistics.txt",col.names = c("t","mg","cug","clg","ms","cus","cls"))

# Plots: average data ####

# for (time in unique(data_est$tf)) {
#   data_est_t <- data_est[data_est$tf == time,]
#   data_obs_t <- data_obs[data_obs$t <= time,]
#   if(nrow(data_est_t[which(!is.na(data_est_t$mg)),]) < 2){
#     plot_gamma <- ggplot(data_obs_t) +
#       geom_ribbon(data = data_obs_t, aes(x = t, ymax = cug, ymin = clg), fill = "black", alpha = 0.3)+
#       geom_ribbon(data = data_est_t, aes(x = t_cor, ymax = cug, ymin = clg), fill = "red", alpha = 0.3)+
#       geom_line(data = data_obs_t, aes(x=t,y = mg), col = "black") +
#       geom_point(data = data_est_t, aes(x=t_cor,y = mg), col ="red") +
#       labs(x = "Times (in generations)", y = "Gamma statistic",
#            title = paste("Time",time))
#   } else {
#     plot_gamma <- ggplot(data_obs_t) +
#       geom_ribbon(data = data_obs_t, aes(x = t, ymax = cug, ymin = clg), fill = "black", alpha = 0.3)+
#       geom_ribbon(data = data_est_t, aes(x = t_cor, ymax = cug, ymin = clg), fill = "red", alpha = 0.3)+
#       geom_line(data = data_obs_t, aes(x=t,y = mg), col = "black") +
#       geom_line(data = data_est_t, aes(x=t_cor,y = mg), col ="red") +
#       labs(x = "Times (in generations)", y = "Gamma statistic",
#            title = paste("Time",time))
#   }
#   if(nrow(data_est_t[which(!is.na(data_est_t$ms)),]) < 2){
#     plot_sackin <- ggplot(data_obs_t) +
#       geom_ribbon(data = data_obs_t, aes(x = t, ymax = cus, ymin = cls), fill = "black", alpha = 0.3)+
#       geom_ribbon(data = data_est_t, aes(x = t_cor, ymax = cus, ymin = cls), fill = "red", alpha = 0.3)+
#       geom_line(data = data_obs_t, aes(x=t,y = ms), col = "black") +
#       geom_point(data = data_est_t, aes(x=t_cor,y = ms), col ="red") +
#       labs(x = "Times (in generations)", y = "Sackin Index",
#            title = paste("Time",time))
#   } else {
#     plot_sackin <- ggplot(data_obs_t) +
#       geom_ribbon(data = data_obs_t, aes(x = t, ymax = cus, ymin = cls), fill = "black", alpha = 0.3)+
#       geom_ribbon(data = data_est_t, aes(x = t_cor, ymax = cus, ymin = cls), fill = "red", alpha = 0.3)+
#       geom_line(data = data_obs_t, aes(x=t,y = ms), col = "black") +
#       geom_line(data = data_est_t, aes(x=t_cor,y = ms), col ="red") +
#       labs(x = "Times (in generations)", y = "Sackin Index",
#            title = paste("Time",time))
#   }
#   ggsave(plot_grid(plot_gamma,plot_sackin,nrow=1,ncol=2),path = "results/plotstatpast", file=paste0("stat_past_t",time,".png"), width=12, height=8)
# }

data_est_filt <- data_est %>%
  filter(tf %in% c(seq(1001,15001,by=1000),seq(20001,300001,by=20000)))

plot_all_gamma <- ggplot(data_obs) +
  geom_ribbon(data = data_obs, aes(x = t, ymax = cug, ymin = clg), fill = "black", alpha = 0.3)+
  geom_ribbon(data = data_est_filt, aes(x = t_cor, ymax = cug, ymin = clg), fill = "red", alpha = 0.3)+
  geom_line(data = data_obs, aes(x=t,y = mg), col = "black") +
  geom_line(data = data_est_filt, aes(x=t_cor,y = mg, group=tf), col ="red") +
  labs(x = "Times (in generations)", y = "Gamma statistic")+
  theme(legend.position = "up")
plot_all_sackin <- ggplot(data_obs) +
  geom_ribbon(data = data_obs, aes(x = t, ymax = cus, ymin = cls), fill = "black", alpha = 0.3)+
  geom_ribbon(data = data_est_filt, aes(x = t_cor, ymax = cus, ymin = cls), fill = "red", alpha = 0.3)+
  geom_line(data = data_obs, aes(x=t,y = ms), col = "black") +
  geom_line(data = data_est_filt, aes(x=t_cor,y = ms, group=tf), col ="red") +
  labs(x = "Times (in generations)", y = "Sackin index")+
  theme(legend.position = "up")

ggsave(plot_all_gamma,path = "results/plotstatpast", file=paste0("gamma_past_t.png"), width=12, height=8)
ggsave(plot_all_sackin,path = "results/plotstatpast", file=paste0("sackin_past_t.png"), width=12, height=8)

# Plots: brut replicates ####

# list_files_obs <- list.files(path=paste0("../sim/dynamic/slf/",model,"/results/"),pattern=paste0("_trees_statistic_gamma.txt"))
# list_files_est <- list.files(path=paste0("../sim/dynamic/slf/",model,"/results/"),pattern=paste0("_trees_statistics_past_estimates.txt"))
# 
# file_obs <- list_files_obs[[1]]
# file_est <- list_files_est[[1]]
# 
# data_obs <- read.table(paste0("../sim/dynamic/slf/",model,"/results/",file_obs),col.names = c("t","gamma"))
# data_est <- read.table(paste0("../sim/dynamic/slf/",model,"/results/",file_est),col.names = c("tf","t","gamma_est","gamma_obs","sackin_est","sackin_obs","error_gamma","error_sackin","irun"))
# data_est$t_corrected <- data_est$tf-data_est$t
# 
# ifelse(!dir.exists("../sim/dynamic/slf/baseline/results/plotstatpast"), dir.create("../sim/dynamic/slf/baseline/results/plotstatpast"), FALSE)
# 
# for (time in unique(data_est$tf)) {
#   data_est_t <- data_est[data_est$tf == time,]
#   data_obs_t <- data_obs[data_obs$t <= time,]
#   
#   plot <- ggplot(data_obs) +
#     geom_line(data = data_obs_t, aes(x=t,y = gamma), col = "black") +
#     geom_line(data = data_est_t, aes(x=t_corrected,y = gamma_est), col ="red") +
#     labs(x = "Times (in generations)", y = "Gamma statistic",
#          title = paste("Time",time)) 
#   
#   ggsave(plot,path = "../sim/dynamic/slf/baseline/results/plotstatpast", file=paste0("gamma_past_t",time,".png"), width=12, height=8)
# }
