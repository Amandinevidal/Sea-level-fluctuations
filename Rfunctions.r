## THIS FUNCTION CHECKS THAT THERE IS NO OBVIOUS ERROR IN PARAMETER
## VALUE. IF THERE IS ONE, STOP SIMULATION, ELSE SAVE PARAMETERS
checkparameters <- function(simulnumber, beta, gamma, sigma, alpha, tau, mui, ci, muc, cc, Nc, theta, nbi, islarch, Kmax, ti, pmax, tl, dc, di, tmax, tiss, wtup, nrun, per, rngseed, static, teq, tint) {
    
    if (!is.numeric(simulnumber)) {
        cat("\nERROR: simulnumber must be an integer\n\n")
        stop()
    }
    if (simulnumber != floor(simulnumber)) {
        cat("\nERROR: simulnumber must be an integer\n\n")
        stop()
    }
    if (beta <= 0) {
        cat("\nERROR: beta must be > 0\n\n")
        stop()
    }
    if (gamma <= 0) {
        cat("\nERROR: gamma must be > 0\n\n")
        stop()
    }
    if (sigma < 0) {
        cat("\nERROR: sigma must be >= 0\n\n")
        stop()
    }
    if (alpha < 0) {
        cat("\nERROR: alpha must be >= 0\n\n")
        stop()
    }
    if (tau < 0) {
        cat("\nERROR: tau must be >= 0\n\n")
        stop()
    }
    if (tau > tmax) {
        cat("\nERROR: tau must be <= tmax\n\n")
        stop()
    }
    if (mui < 0) {
        cat("\nERROR: mui must be >= 0\n\n")
        stop()
    }
    if (ci <= 0) {
        cat("\nERROR: ci must be > 0\n\n")
        stop()
    }
    if (muc < 0) {
        cat("\nERROR: muc must be >= 0\n\n")
        stop()
    }
    if (cc <= 0) {
        cat("\nERROR: cc must be > 0\n\n")
        stop()
    }
    if (Nc < 1) {
        cat("\nERROR: Nc must be >= 1\n\n")
        stop()
    }
    if (theta < 0) {
        cat("\nERROR: theta must be >= 0\n\n")
        stop()
    }
    if (nbi != floor(nbi)) {
        cat("\nERROR: nbi must be an integer\n\n")
        stop()
    }    
    if (nbi < 1) {
        cat("\nERROR: nbi must be >= 1\n\n")
        stop()
    }
    if (islarch != 0 && islarch != 1 && islarch != 2) {
        cat("\nERROR: islarch must be 0, 1 or 2\n\n")
        stop()
    }
    if (sum(Kmax <= 0) != 0) {
        cat("\nERROR: all Kmax must be > 0\n\n")
        stop()
    }
    if (length(Kmax) != nbi) {
        cat("\nERROR: Kmax must have nbi element(s)\n\n")
        stop()
    }
    if (length(ti) != nbi) {
        cat("\nERROR: ti must have nbi element(s)\n\n")
        stop()
    }
    if (ti[1] >= 0) {
        cat("\nERROR: first ti must be < 0\n\n")
        stop()
    }
    if (length(pmax) != nbi) {
        cat("\nERROR: pmax must have nbi element(s)\n\n")
        stop()
    }
    if (sum(pmax <= 0)) {
        cat("\nERROR: all pmax must be > 0\n\n")
        stop()
    }
    if (sum(pmax >= 1)) {
        cat("\nERROR: all pmax must be < 1\n\n")
        stop()
    }
    if (length(tl) != nbi) {
        cat("\nERROR: tl must have nbi element(s)\n\n")
        stop()
    }
    if (sum(tl <= 0) != 0) {
        cat("\nERROR: all tl must be > 0\n\n")
        stop()
    }    
    if (sum(dc <= 0) != 0) {
        cat("\nERROR: all dc must be > 0\n\n")
        stop()
    }
    if (length(dc) != nbi) {
        cat("\nERROR: dc must have nbi element(s)\n\n")
        stop()
    }
    if (!is.matrix(di)) {
        cat("\nERROR: di must be a matrix\n\n")
        stop()
    }
    if (dim(di)[1] !=  dim(di)[2]) {
        cat("\nERROR: di must be a square matrix\n\n")
        stop()
    }
    if (dim(di)[1] != nbi) {
        cat("\nERROR: the dimension of di must be nbi\n\n")
        stop()
    }
    if (sum(di < 0) != 0) {
        cat("\nERROR: all di must be >= 0\n\n")
        stop()
    }
    for (i in 1:nbi) {
        if (di[i,i] != 0) {
            cat("\nERROR: all diagonal elements of di must be 0\n\n")
            stop()
        } 
    }
    if (nbi > 1) {
        for (i in 1:(nbi-1)) {
            for (j in (i+1):nbi) {
                if (di[i,j] != di[j,i]) {
                    cat("\nERROR: di must be symmetrical\n\n")
                    stop()
                }
                if (di[i,j] == 0) {
                    cat("\nERROR: non-diagonal elements of di must be > 0\n\n")
                    stop()
                }
            }
        }   
    }
    if (tmax <= 0) {
        cat("\nERROR: tmax must be > 0\n\n")
        stop()
    }
    archexists <- rep(0,floor(tmax))
    for (i in 1:nbi) {
        for (t in 1:floor(tmax)) {            
            if (t >= ti[i] && t <= ti[i]+tl[i])
                archexists[t] <- 1
        }
    }    
    if (sum(archexists) < length(archexists)) {
         cat("\nERROR: there is at least one moment before tmax with no archipelago\n\n")
         stop()
    }
    if (tiss <= 0) {
        cat("\nERROR: tiss must be > 0\n\n")
        stop()
    }
    if (tiss > tmax) {
        cat("\nERROR: tiss must be <= tmax\n\n");
        stop()
    }
    if (wtup <= 0) {
        cat("\nERROR: wtup must be > 0\n\n")
        stop()
    }
    if (wtup > tmax) {
        cat("\nERROR: wtup must be <= tmax\n\n");
        stop()
    }
    if (wtup >= min(tl)) {
        cat("\nERROR: wtup must be < min(tl)\n\n");
        stop()
    }
    if (wtup > tiss) {
        cat("\nERROR: wtup must be <= tiss\n\n");
        stop()
    }
    if (nrun != floor(nrun)) {
        cat("\nERROR: nrun must be an integer\n\n")
        stop()
    }    
    if (nrun < 1) {
        cat("\nERROR: nrun must be >= 1\n\n")
        stop()
    }
    if (per != 0 && per != 1) {
        cat("\nERROR: per must be 0 or 1\n\n")
        stop()
    }
    if (rngseed != floor(rngseed)) {
        cat("\nERROR: rngseed must be an integer\n\n")
        stop()
    }    
    if (rngseed < 0) {
        cat("\nERROR: rngseed must be >= 0\n\n")
        stop()
    }    
    if (rngseed != 0 && nrun != 1) {
        cat("\nERROR: rngseed = 0 but nrun != 1\n\n")
        stop()
    }
    if (static != 0 && static != 1) {
        cat("\nERROR: static must be 0 or 1\n\n")
        stop()
    }
    if (static == 1) {
        if (teq <= 0) {
            cat("\nERROR: teq must be > 0\n\n")
            stop()
        }    
        if (tint <= 0) {
            cat("\nERROR: tint must be > 0\n\n")
            stop()
        }    
    }
    
    ## save the parameter file
    if (!file.exists("results")) dir.create("results")
    if (static == 1)
        file.copy("parameters.r", "results/simall.parameters.r", overwrite=TRUE)
    else
        file.copy("parameters.r", paste("results/sim", simulnumber, ".parameters.r", sep=""), overwrite=TRUE)
}

## THIS FUNCTION WRITES THE CONSTANTS FOR THE C PROGRAM
writeconstants <- function(nbi, islarch, static) {
    write("#define NB_POP_MAX 5000", file="Cconstants.c", append=FALSE)
    write("#define NB_SP_CONTI_MAX 2500", file="Cconstants.c", append=TRUE)
    write(paste0("#define NB_ISLAND_MAX ", nbi), file="Cconstants.c", append=TRUE)
    write(paste0("#define ISL_ARCH ", islarch), file="Cconstants.c", append=TRUE)
    write(paste0("#define STATIC ", static), file="Cconstants.c", append=TRUE)   
    write("#define TPI 6.283185   // 2*PI", file="Cconstants.c", append=TRUE)
    write("#define PI_TH 5.568328 // PI^(3/2)", file="Cconstants.c", append=TRUE)
    write("#define PI M_PI // pi", file="Cconstants.c",append=TRUE)
    write("#define NB_SP_PHYLO_MAX 50000", file="Cconstants.c", append=TRUE)
    write("#define rad_coef 0.011834057 // normalization radius", file="Cconstants.c", append=TRUE)
}

## THIS FUNCTION IS A CALL TO THE C PROGRAM
runsimul <- function(simulnumber, beta, gamma, sigma, alpha, tau, mui, ci, muc, cc, Nc, theta, nbi, Kmax, ti, pmax, tl, dc, di, tmax, tiss, wtup, nrun, rngseed, teq, tini, phi_max, d, amp, period, phi_min, psi, h_break_point) {

  .C("simul", as.integer(simulnumber), as.double(beta), as.double(gamma), as.double(sigma), as.double(alpha), as.double(tau), as.double(mui), as.double(ci), as.double(muc), as.double(cc), as.integer(Nc), as.double(theta), as.integer(nbi), as.double(Kmax), as.double(ti), as.double(pmax), as.double(tl), as.double(dc), as.double(di), as.double(tmax), as.double(tiss), as.double(wtup), as.integer(nrun), as.integer(rngseed), as.double(teq), as.double(tini),as.double(phi_max),as.double(d),as.double(amp),as.double(period),as.double(phi_min),as.double(psi),as.double(h_breakpoint))

}

