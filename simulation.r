## include R functions and libraries
source("Rfunctions.r")

## include parameter values
source("parameters.r")
for (i in 1:nbi) di[i,i] <- 0  ## distance between an island and itself must be zero
if (min(ti) == 0) ti <- ti - 1 ## the first island must be emerged when the simulation starts
tmax <- max(ti+tl)             ## duration of simulations

## check for obvious errors in parameter values
checkparameters(simulnumber, beta, gamma, sigma, alpha, tau, mui, ci, muc, cc, Nc, theta, nbi, islarch, Kmax, ti, pmax, tl, dc, di, tmax, tiss, wtup, nrun, per, rngseed, static, teq, tint)

## write the constants, compile the C program and load library
writeconstants(nbi, islarch, static)
if (file.exists("Cfunctions.o")) file.remove("Cfunctions.o")
if (file.exists("Cfunctions.so")) file.remove("Cfunctions.so")
# dyn.unload("Cfunctions.so")
## next line is valid on linux OS. To be adapted for other OS. 
system("R CMD SHLIB -lm -lgsl -lgslcblas Cfunctions.c")
dyn.load("Cfunctions.so")


## run the simulations and plot the results
if (static) { ## static archipelago
    vtini <- seq(0,tmax,tint)
    for (simulnumber in 1:length(vtini)) { # for each possible static state of the archipelago    
        tini <- vtini[simulnumber]
        output <- runsimul(simulnumber, beta, gamma, sigma, alpha, tau, mui, ci, muc, cc, Nc, theta, nbi, Kmax, ti, pmax, tl, dc, di, tmax, tiss, wtup, nrun, rngseed, teq, tini, phi_max, d, amp, period, phi_min, psi, h_break_point)
   #   source("plotresults.r")
    }
    source("plotequilibrium.r")
} else { ## dynamic archipelago
    tini <- 0 
    output <- runsimul(simulnumber, beta, gamma, sigma, alpha, tau, mui, ci, muc, cc, Nc, theta, nbi, Kmax, ti, pmax, tl, dc, di, tmax, tiss, wtup, nrun, rngseed, teq, tini, phi_max, d, amp, period, phi_min, psi, h_break_point)
    source("build_phylo.r")
    source("plotresults.r")
    # source("stat_phylo.r")
}


