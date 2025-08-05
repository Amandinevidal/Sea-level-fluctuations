#### ---------------------------Informations--------------------------####
#### Date : 2025-01-31
#### Author : Amandine Vidal-Hosteng
#### Encoding : UTF-8
#### Email : amandine.vidal-hosteng@univ-tlse3.fr
#### ------------------------------------------------------------------ #

rm(list = ls())

# Library ####
library(plotly)
library(dplyr)
library(htmlwidgets)
source("parameters.r")

# plot file
ifelse(!dir.exists("results/plottopo"), dir.create("results/plottopo"), FALSE)

#### PLOTS TOPO PARAMETERS ####
data <- read.table("out_arch")
colnames(data) <- c("t",paste0("K",1:nbi),paste0("dc",1:nbi),paste0("di", rep(1:nbi, each = nbi), rep(1:nbi, times = nbi)),"sl",paste0("Kadj",1:nbi),paste0("phi",1:nbi),paste0("h",1:nbi),paste0("hz",1:nbi),paste0("r",1:nbi),paste0("rz",1:nbi),"rc")
brut_data <- read.table("out_arch")
colnames(brut_data) <- c("t",paste0("K",1:nbi),paste0("dc",1:nbi),paste0("di", rep(1:nbi, each = nbi), rep(1:nbi, times = nbi)),"sl",paste0("Kadj",1:nbi),paste0("phi",1:nbi),paste0("h",1:nbi),paste0("hz",1:nbi),paste0("r",1:nbi),paste0("rz",1:nbi),"rc")

main <- ""

png(filename = paste0("results/plottopo/K.png"), width = 800, height = 600, res = 120)
plot(brut_data$t,brut_data$K1,type='l',lwd = 5, xlab = "Time", ylab = "Carrying capacty",main = main, xlim=c(min(brut_data$t),max(brut_data$t)))
lines(brut_data$t,brut_data$K2,lwd = 5)
lines(brut_data$t,brut_data$K3,lwd = 5)
lines(brut_data$t,brut_data$K4,lwd = 5)
lines(brut_data$t,brut_data$Kadj1, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$Kadj2, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$Kadj3, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$Kadj4, col="lightblue", lty = 1, lwd = 1)
legend(x = "topright",          # Position
       legend = c("K", "Kadj"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("black","lightblue"),           # Line colors
       lwd = 2) 
dev.off()  # Termine la sauvegarde du fichier

png(filename = paste0("results/plottopo/radius.png"), width = 800, height = 600, res = 120)
plot(brut_data$t,brut_data$r1,type='l',lwd = 5, xlab = "Time", ylab = "Island's radius",main = main,xlim=c(min(brut_data$t),max(brut_data$t)))
lines(brut_data$t,brut_data$r2,lwd = 5)
lines(brut_data$t,brut_data$r3,lwd = 5)
lines(brut_data$t,brut_data$r4,lwd = 5)
lines(brut_data$t,brut_data$rz1, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$rz2, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$rz3, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$rz4, col="lightblue", lty = 1, lwd = 1)
legend(x = "topright",          # Position
       legend = c("r", "rz"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("black","lightblue"),           # Line colors
       lwd = 2) 
dev.off()

png(filename = paste0("results/plottopo/height.png"), width = 800, height = 600, res = 120)
plot(brut_data$t,brut_data$h1,type='l',lwd = 5, xlab = "Time", ylab = "Island's height",main = main,xlim=c(min(brut_data$t),max(brut_data$t)))
lines(brut_data$t,brut_data$h2, lwd = 5)
lines(brut_data$t,brut_data$h3, lwd = 5)
lines(brut_data$t,brut_data$h4, lwd = 5)
lines(brut_data$t,brut_data$hz1, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$hz2, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$hz3, col="lightblue", lty = 1, lwd = 1)
lines(brut_data$t,brut_data$hz4, col="lightblue", lty = 1, lwd = 1)
legend(x = "topright",          # Position
       legend = c("h", "hz"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("black","lightblue"),           # Line colors
       lwd = 2) 
dev.off()

png(filename = paste0("results/plottopo/theta.png"), width = 800, height = 600, res = 120)
plot(brut_data$t,brut_data$phi1,type='l',lwd = 5, xlab = "Time", ylab = "Island's apex angle",main = main,xlim=c(min(brut_data$t),max(brut_data$t)))
lines(brut_data$t,brut_data$phi2, lwd = 5)
lines(brut_data$t,brut_data$phi3, lwd = 5)
lines(brut_data$t,brut_data$phi4, lwd = 5)
legend(x = "topright",          # Position
       legend = c("theta"),  # Legend texts
       lty = c(1),           # Line types
       col = c("black"),           # Line colors
       lwd = 2) 
dev.off()

png(filename = paste0("results/plottopo/dc.png"), width = 800, height = 600, res = 120)
plot(brut_data$t, brut_data$dc1, type = 'l', lwd = 1, 
     xlab = "Time", ylab = "Mainland distance", 
     main = main, xlim = c(min(brut_data$t), max(brut_data$t)), 
     ylim = c(min(brut_data$dc1, brut_data$dc2), max(brut_data$dc1, brut_data$dc2)))
lines(brut_data$t, brut_data$dc2, col = "black", lwd = 1)
lines(brut_data$t, brut_data$dc3, col = "black", lwd = 1)
lines(brut_data$t, brut_data$dc4, col = "black", lwd = 1)
par(new = TRUE)
plot(brut_data$t, brut_data$K1/1000,type='l', col = 'grey', lwd = 2, 
     axes = FALSE, xlab = "", ylab = "", 
     ylim = c(0, max(brut_data$K1/1000))
)
lines(brut_data$t, brut_data$K2/1000, col = 'cornsilk4', lwd = 2)
axis(4)  # Ajouter un axe y à droite
dev.off()

png(filename = paste0("results/plottopo/di.png"), width = 800, height = 600, res = 120)
plot(brut_data$t, brut_data$di12, type = 'l', lwd = 1, 
     xlab = "Time", ylab = "Inter-island distance", 
     main = main, xlim = c(min(brut_data$t), max(brut_data$t)), 
     ylim = c(min(brut_data$di12, brut_data$di21), max(brut_data$di12, brut_data$di21)))
lines(brut_data$t, brut_data$di13, col = "black", lwd = 1)
lines(brut_data$t, brut_data$di14, col = "black", lwd = 1)
lines(brut_data$t, brut_data$di23, col = "red", lwd = 1)
lines(brut_data$t, brut_data$di24, col = "red", lwd = 1)
lines(brut_data$t, brut_data$di34, col = "green", lwd = 1)
par(new = TRUE)
plot(brut_data$t, brut_data$K1/1000,type='l', col = 'grey', lwd = 2, 
     axes = FALSE, xlab = "", ylab = "", 
     ylim = c(0, max(brut_data$K1/1000))
)
lines(brut_data$t, brut_data$K2/1000, col = 'cornsilk4', lwd = 2)
axis(4)  # Ajouter un axe y à droite
dev.off()

png(filename = paste0("results/plottopo/slv.png"), width = 800, height = 600, res = 120)
plot(brut_data$t,brut_data$sl,type='l',lwd = 1, xlab = "Time", ylab = "Sea level",main = main)
legend(x = "topright",          # Position
       legend = c("slv"),  # Legend texts
       lty = c(1),           # Line types
       col = c("black"),           # Line colors
       lwd = 2) 
dev.off()

#### PLOTS RATES PARAMETERS ####
# data <- read.table("results/sim1.0_rates")
data <- read.table("out_rates")
data_baseline <- read.table("/home/avidalhosteng/Bureau/slv/sim/sim_test_baseline_comparison/two_isl/two_isl_baseline/results/sim1.0_rates")
colnames(data) <- c("t",paste0("n_mig_conti",1:nbi),paste0("n_mig_isl", rep(1:nbi, each = nbi), rep(1:nbi, times = nbi)),paste0("n_ana_isl",1:nbi),paste0("n_ana_conti",1:nbi),paste0("n_clado",1:nbi),paste0("n_spe_event",1:nbi),paste0("n_spe_event_ori_mig_conti",1:nbi),paste0("n_spe_event_ori_mig_isl",1:nbi),paste0("n_spe_event_ori_clado",1:nbi),paste0("n_pop_ext",1:nbi),paste0("n_spe_ext",1:nbi),"n_spe_ext_arch")
colnames(data_baseline) <-  c("t",paste0("n_mig_conti",1:nbi),paste0("n_mig_isl", rep(1:nbi, each = nbi), rep(1:nbi, times = nbi)),paste0("n_ana_isl",1:nbi),paste0("n_ana_conti",1:nbi),paste0("n_clado",1:nbi),paste0("n_spe_event",1:nbi),paste0("n_spe_event_ori_mig_conti",1:nbi),paste0("n_spe_event_ori_mig_isl",1:nbi),paste0("n_spe_event_ori_clado",1:nbi),paste0("n_pop_ext",1:nbi),paste0("n_spe_ext",1:nbi),"n_spe_ext_arch")

main <- ""
png(filename = paste0("results/plottopo/mig_isl.png"), width = 800, height = 600, res = 120)
plot(data$t, data$n_mig_isl12, type = 'l')
lines(data_baseline$t, data_baseline$n_mig_isl12, col = "red", lwd = 5)
dev.off()

#### 3D MODELING ####

# Data ####
source("parameters.r")
data <- read.table("results/sim1.0_arch")
colnames(data) <- c("t",paste0("K",1:nbi),paste0("dc",1:nbi),paste0("di", rep(1:nbi, each = nbi), rep(1:nbi, times = nbi)),"sl",paste0("Kadj",1:nbi),paste0("phi",1:nbi),paste0("h",1:nbi),paste0("hz",1:nbi),paste0("r",1:nbi),paste0("rz",1:nbi),"rc")

# time sequence for the animation (number of image you want)
time_seq <- 1000
specified_time <- seq(min(data$t), max(data$t), by = time_seq)

# filter data set
data_filtered <- data %>% filter(t %in% specified_time)

# space limits
xlim_fixed <- c(-max(data_filtered$r1, data_filtered$rz1) * 1.5, max(data_filtered$r1, data_filtered$rz1) * 3.5) # width
ylim_fixed <- c(min(data_filtered$sl) * 1.2, max(data_filtered$h1) * 1.2) # height

# Offset for each island (based on inter island centroid distance, first one on zero)
x_offset <- c(0,rep(di[1],nbi-1))

# Functions ####
create_cone <- function(base_radius, height, z_offset, x_offset = 0) {
  angle <- seq(0, 2 * pi, length.out = 100)
  x_base <- base_radius * cos(angle) + x_offset
  y_base <- base_radius * sin(angle)
  
  x <- c(x_offset, x_base)
  y <- c(0, y_base)
  z <- c(height, rep(z_offset, length(angle)))
  
  return(list(x = x, y = y, z = z))
}

# Code ####

# create list of frames that contains 3d mesh of each cones
frames <- vector(mode = "list", length = length(specified_time))    
names(frames) <- paste("Temps:", specified_time) 

for (current_time in 1:length(specified_time)) {
  cone_data <- data_filtered %>% filter(t == specified_time[current_time])
  cone_list <- list()
  
  for (i in 1:nbi) { # for each island
    r <- cone_data[, paste0("r", i)]
    rz <- cone_data[, paste0("rz", i)]
    h <- cone_data[, paste0("h", i)]
    z <- 0
    zz <- cone_data[, "sl"]
    rmax <- (h - (-amp)) * tan(cone_data[, paste0("phi", i)])
    cone_list[[paste0("full_cone", i)]] <- create_cone(r, h, z, x_offset[i])
    cone_list[[paste0("new_cone", i)]] <- create_cone(rz, h, zz, x_offset[i])
    cone_list[[paste0("sea_cone", i)]] <- create_cone(rmax, h, -amp, x_offset[i])
  }
  
  frames[[current_time]] <- c(cone_list, list(time = specified_time[current_time]))
}

fig_animation <- plot_ly() # object for animation

for (i in 1:length(frames)) { # for each frame (time seq)
  # figure init
  fig_animation <- fig_animation
  
  for (j in 1:nbi) { # for each island
    fig_animation <- fig_animation %>%
      add_trace(x = frames[[i]][[paste0("sea_cone", j)]]$x,
                y = frames[[i]][[paste0("sea_cone", j)]]$y,
                z = frames[[i]][[paste0("sea_cone", j)]]$z,
                type = "mesh3d", intensity = 1, 
                opacity = 0.9, colorscale = list(c(0, 'lightblue'), c(1, 'lightblue')), 
                frame = i, showscale = FALSE) %>%
      add_trace(x = frames[[i]][[paste0("full_cone", j)]]$x,
                y = frames[[i]][[paste0("full_cone", j)]]$y,
                z = frames[[i]][[paste0("full_cone", j)]]$z,
                type = "mesh3d", intensity = 1,
                opacity = 0.9, colorscale = list(c(0, 'lightblue'), c(1, 'lightblue')),
                frame = i, showscale = FALSE) %>%
      add_trace(x = frames[[i]][[paste0("new_cone", j)]]$x,
                y = frames[[i]][[paste0("new_cone", j)]]$y,
                z = frames[[i]][[paste0("new_cone", j)]]$z,
                type = "mesh3d", intensity = 1,
                opacity = 0.9, colorscale = list(c(0, 'orange'), c(1, 'orange')),
                frame = i, showscale = FALSE)
  }
  
  # animation parameters
  fig_animation <- fig_animation %>%
    layout(scene = list(aspectmode = "cube",
                        xaxis = list(title = "X Axis", range = xlim_fixed, autorange = FALSE),
                        yaxis = list(title = "Y Axis", range = xlim_fixed, autorange = FALSE),
                        zaxis = list(title = "Z Axis", range = ylim_fixed, autorange = FALSE)))
}

fig_animation <- fig_animation %>% animation_opts(frame = 100, transition = 0, redraw = TRUE) %>% animation_slider(currentvalue = list(prefix = "New Cone t = ", font = list(size = 14)))

# see live
# fig_animation

saveWidget(fig_animation, "results/plottopo/animation_topo.html")

