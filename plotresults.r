plotmat <- function(time, mat, vectwindow, adresse1, adresse2, graphename, nbi, tiss) {
    if (sum(!is.na(mat)) != 0) { ## if there are data to plot
        nlig <- dim(mat)[1]
        ncol <- dim(mat)[2]
        coul <- rainbow(ncol+1) ## +1 allow to avoid the flash yellow which is invisible

        ## window = nb of time steps over which means are computed
        nwindow <- length(vectwindow)
        png(paste0(adresse2, "_", graphename, ".png"), width=500*nwindow, height=500)
        if (nwindow > 1)
            layout(matrix(1:nwindow, 1, nwindow))
        
        for (window in vectwindow) {
            
            if (window != 1) {
                ## build the mean matrix; for the last time steps, the mean is computed
                ## over the available time steps
                meanmat <- matrix(NA, nlig, ncol)
                for (ilig in 1:nlig) 
                    for (icol in 1:ncol)
                        meanmat[ilig,icol] <- mean(mat[ilig:(min(c(ilig+window,nlig))),icol],na.rm=TRUE)
                mattoplot <- meanmat
            } else {
                mattoplot <- mat
            }

            ## plot the mean matrix
            par(cex.axis=1+0.5*nwindow, cex.lab=1+0.5*nwindow, cex.main=1+0.25*nwindow, mar=c(5,5.5,4,2))
            plot(time, mattoplot[,1], type="l", col=coul[1], lwd=2, ylim=c(0,max(mattoplot, na.rm=TRUE)), ylab=graphename, main=paste0(graphename, " (mean over ", tiss*window, " gen.)"))
            if (ncol > 1) {
                for (i in 2:ncol)
                    lines(time, mattoplot[,i], col=coul[i], lwd=2)
            }

            ## superimpose (scaled) K        
            coul <- rainbow(nbi+1)
            d <- read.table(paste0(adresse1, "_arch"))
            K <- d[, 2:(2+nbi-1)]
            K <- max(mattoplot, na.rm=TRUE) * K / max(K)
            for (i in 1:nbi)
                lines(time, as.matrix(K)[,i], col=coul[i], lty=2)
        }
        dev.off()
    } else {
        nwindow <- 1
        png(paste0(adresse2, "_", graphename, ".png"), width=500*nwindow, height=500)
        par(cex.axis=1+0.5*nwindow, cex.lab=1+0.5*nwindow, cex.main=1+0.25*nwindow, mar=c(5,5.5,4,2))
        plot(0, 0, type="n", ylab=graphename, main=graphename, xlab="time")
        dev.off()
    }
}

##source("parameters.r")
##for (i in 1:nbi) di[i,i] <- 0
##if (min(ti) == 0) ti <- ti - 1
##tmax <- max(ti+tl)


##instead, if entered manually, needed parameters are:
##simulnumber <- 1
##nrun <- 1
##nbi <- 4 
##islarch <- 0
##tiss <- 25
##tmax <- 3.5e5-1
##static <- 0
##per <- 0


## read time vector
d <- read.table(paste0("results/sim", simulnumber, ".0_arch"))
time <- d[,1]
nt <- length(time)

## means over simulation replicates
Mdivisl <- matrix(0, nt, nbi)
Mendem <- matrix(0, nt, nbi)
Mpopsizeisl <- matrix(0, nt, nbi)
Mpresdurisl <- matrix(0, nt, nbi)
Mspecdurisl <- matrix(0, nt, nbi)
Mspsizeisl <- matrix(0, nt, nbi)
Mrateanagconti <- matrix(0, nt, nbi)
Mratecladog <- matrix(0, nt, nbi)
Mratemigconti <- matrix(0, nt, nbi)
Mratepopext <- matrix(0, nt, nbi)
Mratespecallisl <- matrix(0, nt, nbi)
Mratespeccladog <- matrix(0, nt, nbi)
Mratespecmigconti <- matrix(0, nt, nbi)
Mratespextisl <- matrix(0, nt, nbi)
if (nbi > 1 && islarch == 0) {       ## some outputs are meaningful only if
    Mpopsizearch <- matrix(0, nt, 2) ## we are in a real archipelago
    Mdivarch <- matrix(0, nt, 1)
    MendemSIE <- matrix(0, nt, nbi)
    MendemMIE <- matrix(0, nt, nbi)
    Mspsizearch <- matrix(0, nt, 1)
    Mpresdurarch <- matrix(0, nt, 1)
    Mspecdurarch <- matrix(0, nt, 1)
    Mratemigisl <- matrix(0, nt, nbi)
    Mrateanagisl <- matrix(0, nt, nbi)
    Mratespecmigisl <- matrix(0, nt, nbi)
    Mratespextarch <- matrix(0, nt, 1)
    Mratespecallarch <- matrix(0, nt, 1)
}

for (irun in 0:(nrun-1)) {

    adresse <- paste0("results/sim", simulnumber, ".", irun)

    
    ## PLOT K ################################################################

    d <- read.table(paste0(adresse, "_arch"))
    K <- d[, 2:(2+nbi-1)]
    if (static == 0 && per == 1) {
        plotmat(time, as.matrix(K), 1, adresse, adresse, "K", nbi, tiss)
    }
##    Kmax <- max(K)
    
    ## PLOT CONTINENT ########################################################

    dconti <- read.table(paste0(adresse, "_conti"))
    spidcontimax <- dconti[dim(dconti)[1],1]

    ## species-rank abundance curve on the continent
    spcontifreq <- c(dconti$V2[1], diff(dconti$V2)) ## freq of each species
    if (static == 0 && per == 1) {
        png(paste0(adresse, "_conti.png"), width=500, height=500)
        par(cex.axis=1.5, cex.lab=1.5, cex.main=1.25, mar=c(5,4.5,4,2))    
        plot(sort(spcontifreq, decreasing=TRUE), xlab="Abundance rank", ylab="Relative abundance", pch=20, type="b", main="Species-rank abundance on the continent")
        dev.off()
    }
    
    ## PLOT POP SIZES AND DIVERSITIES ########################################

    d <- read.table(paste0(adresse, "_metapop"))
    div <- matrix(0, nt, nbi)  ## nb of different species on each island
    divarch <- matrix(0, nt, 1)         ## nb of different species on archipelago
    N <- matrix(0, nt, nbi)    ## nb of indiv on each island
    SIE <- matrix(0, nt, nbi)  ## nb of SIE
    MIE <- matrix(0, nt, nbi)  ## nb of MIE
    E <- matrix(0, nt, nbi)    ## nb of endemic species
    
    for (i in 1:nt) {
 	subd <- d[d$V1==time[i],]   ## table of this time step
 	if (!is.na(subd$V2[1])) {   ##  if there are populations existing
            for (j in 1:nbi) {
                subsubd <- subd[subd$V2==(j-1),]   ## on this island
                if (dim(subsubd)[1] > 0) {
                    div[i,j] <- sum(!duplicated(subsubd$V3))
                    N[i,j] <- sum(subsubd$V5)
                    ## for each species, check its endemic status
                    spid <- subsubd$V3[!duplicated(subsubd$V3)] ## list of sp id on this isl
                    for (k in 1:length(spid)) {
                        if (spid[k] > spidcontimax) { ## this species is endemic
                            ## check if it exists on other islands
                            ## select pops on other islands, of the same species
                            check <- subd[subd$V2!=(j-1),]
                            if (dim(check)[1] > 0) { ## if there are species on other islands
                                checkcheck <- check[check$V3==spid[k],]
                                if (dim(checkcheck)[1] > 0) ## if we find indiv of this species
                                    MIE[i,j] <- MIE[i,j] + 1
                                else
                                    SIE[i,j] <- SIE[i,j] + 1
                            } else {
                                SIE[i,j] <- SIE[i,j] + 1 ## SIE because only 1 island has species
                            }
                        }
                    }
                }
            }
            divarch[i,1] <- sum(!duplicated(subd$V3))
 	}
    }
    E <- SIE + MIE
    ## endemic status as a proportion of diversity
    ## the proportion of SIE et MIE are not computed "knowing that it is endemic"
    ## therefore pSIE + pMIE != 1 but pSIE + pMIE = pE. 
    divNA <- div
    divNA[divNA==0] <- NA
    pSIE <- SIE / divNA
    pMIE <- MIE / divNA
    pE <- E / divNA

    if (per) {
        if (static)
            plotmat(time, div, c(1,10,100), adresse, adresse, "div_isl", nbi, tiss)
        else
            plotmat(time, div, c(200,400,800), adresse, adresse, "div_isl", nbi, tiss)
        if (!static) {
            plotmat(time, N, c(4,40), adresse, adresse, "pop-size_isl", nbi, tiss)
            plotmat(time, pE, c(1,80,200), adresse, adresse, "endem_all", nbi, tiss)
        }
    }
    Mdivisl <- Mdivisl + div ## save results of this run
    Mpopsizeisl <- Mpopsizeisl + N
    Mendem <- Mendem + pE

    if (nbi > 1 && islarch == 0) {
        Ntot <- as.matrix(apply(N,1,sum))
        Ktot <- as.matrix(apply(K,1,sum))
        if (per) {
            if (static)
                plotmat(time, divarch, c(1,10,100), adresse, adresse, "div_arch", nbi, tiss)
            else
                plotmat(time, divarch, c(200,400,800), adresse, adresse, "div_arch", nbi, tiss)
            if (!static) {
                plotmat(time, cbind(Ntot,Ktot), c(4,40), adresse, adresse, "pop-size_arch", nbi, tiss)
                plotmat(time, pSIE, c(1,80,200), adresse, adresse, "endem_SIE", nbi, tiss)
                plotmat(time, pMIE, c(1,80,200), adresse, adresse, "endem_MIE", nbi, tiss)
            }
        }
        Mpopsizearch <- Mpopsizearch + cbind(Ntot,Ktot)
        Mdivarch <- Mdivarch + divarch
        MendemSIE <- MendemSIE + pSIE
        MendemMIE <- MendemMIE + pMIE        
    }


    ## PLOT SPECIES SIZE #################################################

    ## species size distribution/average on each island
    spsizeisl <- matrix(0, nt, nbi)
##    divmax <- max(div)
##    pdf(paste0(adresse, "_sppopsize_isl.pdf"), height=7, width=7*nbi)
##    par(cex.axis=1+0.2*nbi, cex.lab=1+0.2*nbi, cex.main=1+0.15*nbi, mar=c(5,5.5,4,2))
    for (i in 1:nt) {
 	subd <- d[d$V1==time[i],]   ## table of this time step
 	if (!is.na(subd$V2[1])) {   ## if there are populations existing
##            layout(matrix(1:nbi,1,nbi))
            for (j in 1:nbi) {
                subsubd <- subd[subd$V2==(j-1),]   ## on this island
                if (dim(subsubd)[1] > 0) { ## if there are populations existing
                    spid <- subsubd$V3[!duplicated(subsubd$V3)] ## list of sp id on this isl
                    sppopsize <- matrix(0,length(spid),2)
                    sppopsize[,1] <- spid
                    for (k in 1:dim(subsubd)[1]) { ## sum pop size of all pops of the same species
                        sppopsize[sppopsize[,1]==subsubd$V3[k],2] <- sppopsize[sppopsize[,1]==subsubd$V3[k],2] + subsubd$V5[k]
                    }
                    spsizeisl[i,j] <- mean(sppopsize[,2])
##                    sppopsize[,2] <- sppopsize[,2] / sum(sppopsize[,2]) ## pop size as frequency
##                    plot(sort(sppopsize[,2], decreasing=TRUE), xlab="Abundance rank", ylab="Relative abundance", pch=20, type="b", main=paste0("Island ", j, " - time=", time[i]), xlim=c(0,divmax), ylim=c(0,1))
                }## else {
##                    plot(0, 0, xlab="", ylab="", type="n", main=paste0("Island ", j, " - time=", time[i]), xaxt="n", yaxt="n")
##                }
            }
        }
    }
##    dev.off()
    if (static == 0 && per == 1) {
        plotmat(time, spsizeisl, c(1,50,200), adresse, adresse, "sp-size_isl", nbi, tiss)
    }
    Mspsizeisl <- Mspsizeisl + spsizeisl
   
    ## species size distribution/average on the archipelago
    if (nbi > 1 && islarch == 0) {
        spsizearch <- matrix(0, nt, 1)
    ##    divmax <- max(divarch)
    ##    pdf(paste0(adresse, "_sppopsize_arch.pdf"), height=7, width=7)
    ##    par(cex.axis=1.25, cex.lab=1.25, cex.main=1.2, mar=c(5,5,4,2))    
        for (i in 1:nt) {
     	subd <- d[d$V1==time[i],]   ## table of this time step
     	if (!is.na(subd$V2[1])) {   ## if there are populations existing
                spid <- subd$V3[!duplicated(subd$V3)] ## list of sp id
                sppopsize <- matrix(0,length(spid),2)
                sppopsize[,1] <- spid
                for (k in 1:dim(subd)[1]) { ## sum pop size of all pops of the same species
                    sppopsize[sppopsize[,1]==subd$V3[k],2] <- sppopsize[sppopsize[,1]==subd$V3[k],2] + subd$V5[k]
                }
                spsizearch[i,1] <- mean(sppopsize[,2])            
    ##            sppopsize[,2] <- sppopsize[,2] / sum(sppopsize[,2]) ## pop size as frequency
    ##            plot(sort(sppopsize[,2], decreasing=TRUE), xlab="Abundance rank", ylab="Relative abundance", pch=20, type="b", main=paste0("Whole archipelago - time=", time[i]), xlim=c(0,divmax), ylim=c(0,1))
            }
        }
    ##    dev.off()
        if (static == 0 && per == 1) {
            plotmat(time, spsizearch, c(1,50,200), adresse, adresse, "sp-size_arch", nbi, tiss)
        }
        Mspsizearch <- Mspsizearch + spsizearch
    }

    
    ## PLOT SPECIES PRESENCE DURATION ###########################################
        
    ## presence duration on each island (continuous presence; if the species was present, went extinct
    ## and is again present, only the last event is taken into account; extinction erase memory)
    ## presence duration is average NOT weighted by pop size
    presdurisl <- matrix(0, nt, nbi)
##    divmax <- max(div)
##    pdf(paste0(adresse, "_presence_isl.pdf"), height=7, width=7*nbi)
##    par(cex.axis=1+0.2*nbi, cex.lab=1+0.2*nbi, cex.main=1+0.15*nbi, mar=c(5,5.5,4,2))
    for (i in 1:nt) {
 	subd <- d[d$V1==time[i],]   ## table of this time step
 	if (!is.na(subd$V2[1])) {   ## if there are populations existing
##            layout(matrix(1:nbi,1,nbi))
            for (j in 1:nbi) {
                subsubd <- subd[subd$V2==(j-1),]   ## on this island
                if (dim(subsubd)[1] > 0) { ## if there are populations existing
                    spid <- subsubd$V3[!duplicated(subsubd$V3)] ## list of sp id on this isl
                    nspid <- length(spid)
                    sppres <- numeric(nspid)
                    for (k in 1:nspid)  ## presence duration computed from min date of origin
                        sppres[k] <- time[i] - min(subsubd$V7[subsubd$V3==spid[k]])
                    presdurisl[i,j] <- mean(sppres)
##                    plot(sort(sppres, decreasing=TRUE), xlab="Rank", ylab="Presence duration", pch=20, type="b", main=paste0("Island ", j, " - time=", time[i]), xlim=c(0,divmax), ylim=c(0,Kmax))
                }## else {
##                    plot(0, 0, xlab="", ylab="", type="n", main=paste0("Island ", j, " - time=", time[i]), xaxt="n", yaxt="n")
##                }
            }
        }
    }
##    dev.off()
    if (static == 0 & per == 1) {
        plotmat(time, presdurisl, c(1,50,200), adresse, adresse, "pres-dur_isl", nbi, tiss)
    }
    Mpresdurisl <- Mpresdurisl + presdurisl
    
    ## presence duration on the archipelago (continuous prensence; see above)
    if (nbi > 1 && islarch == 0) {
        presdurarch <- matrix(0, nt, 1)
    ##    divmax <- max(divarch)
    ##    pdf(paste0(adresse, "_presence_arch.pdf"), height=7, width=7)
    ##    par(cex.axis=1.25, cex.lab=1.25, cex.main=1.2, mar=c(5,5,4,2))    
        for (i in 1:nt) {
     	subd <- d[d$V1==time[i],]   ## table of this time step
     	if (!is.na(subd$V2[1])) {   ## if there are populations existing
                spid <- subd$V3[!duplicated(subd$V3)] ## list of sp id
                nspid <- length(spid)
                sppres <- numeric(nspid)
                for (k in 1:nspid)  ## presence duration computed from min date of origin
                    sppres[k] <- time[i] - min(subd$V7[subd$V3==spid[k]])
                presdurarch[i,1] <- mean(sppres)
    ##            plot(sort(sppres, decreasing=TRUE), xlab="Rank", ylab="Presence duration", pch=20, type="b", main=paste0("Whole archipelago - time=", time[i]), xlim=c(0,divmax), ylim=c(0,Kmax))
            }
        }
        ##    dev.off()
        if (static == 0 && per == 1) {
            plotmat(time, presdurarch, c(1,50,200), adresse, adresse, "pres-dur_arch", nbi, tiss)
        }
        Mpresdurarch <- Mpresdurarch + presdurarch
    }

    ## PLOT SPECIATION DURATION #################################################

    specdurisl <- matrix(NA, nt, nbi) ## mean speciation duration of existing species on each island
    specdurarch <- matrix(NA, nt, 1)  ## same on the archipelago
        
    for (i in 1:nt) {
 	subd <- d[d$V1==time[i],]   ## table of this time step
 	if (!is.na(subd$V2[1])) {   ##  if there are populations existing
            for (j in 1:nbi) {
                subsubd <- subd[subd$V2==(j-1),]   ## on this island
                if (dim(subsubd)[1] > 0) {
                    subsubsubd <- subsubd[subsubd$V4>tmax,] ## species which have speciated
                    if (dim(subsubsubd)[1] > 0) {
                        specdurisl[i,j] <- mean(subsubsubd$V8)
                    }
                }
            }
            specdurarch[i,1] <- mean(specdurisl[i,], na.rm=TRUE)
        }
    }
    if (static == 0 && per == 1) {
        plotmat(time, specdurisl, c(1,20,50), adresse, adresse, "spec-dur_isl", nbi, tiss)
    }
    Mspecdurisl <- Mspecdurisl + specdurisl
    
    if (nbi > 1 && islarch == 0) {
        if (static == 0 && per == 1) {
            plotmat(time, specdurarch, c(1,20,50), adresse, adresse, "spec-dur_arch", nbi, tiss)
        }
        Mspecdurarch <- Mspecdurarch + specdurarch
    }

    
    ## PLOT RATES ################################################################
    
    d <- read.table(paste0(adresse, "_rates"))
    ## outputs are counts: divide them by tiss to convert them as rates
    ncol <- dim(d)[2]
    d[,2:ncol] <- d[,2:ncol] / tiss
    
    c <- 2 ## column number
    migconti <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, migconti, c(800,2000,3200), adresse, adresse, "rate_mig_conti", nbi, tiss)
    }
    Mratemigconti <- Mratemigconti + migconti
    c <- c + nbi

    if (nbi > 1 && islarch == 0) {
        migisl <- d[, c:(c+nbi*nbi-1)]
        ## total migration (from any island) to each island    
        totmigisl <- matrix(0,nt,nbi)
        for (i in 1:nbi)
            totmigisl[,i] <- apply(migisl[,seq(i,i+(nbi-1)*nbi,nbi)], 1, sum)
        if (static == 0 && per == 1) {
            plotmat(time, totmigisl, c(200,400,800), adresse, adresse, "rate_mig_isl", nbi, tiss)
        }
        Mratemigisl <- Mratemigisl + totmigisl
    }
    c <- c + nbi*nbi

    if (nbi > 1 && islarch == 0) {
        anaisl <- d[, c:(c+nbi-1)]
        if (static == 0 && per == 1) {
            plotmat(time, anaisl, c(200,400,800), adresse, adresse, "rate_anag_isl", nbi, tiss)
        }
        Mrateanagisl <- Mrateanagisl + anaisl
    }
    c <- c + nbi

    anaconti <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, anaconti, c(800,2000,3200), adresse, adresse, "rate_anag_conti", nbi, tiss)
    }
    Mrateanagconti <- Mrateanagconti + anaconti
    c <- c + nbi

    clado <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, clado, c(4,200), adresse, adresse, "rate_cladog", nbi, tiss)
    }
    Mratecladog <- Mratecladog + clado
    c <- c + nbi

    ratespecallisl <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, ratespecallisl, c(100,200,400), adresse, adresse, "rate_spec_all_isl", nbi, tiss)
    }
    Mratespecallisl <- Mratespecallisl + ratespecallisl
    c <- c + nbi
    spe <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, spe, c(100,200,400), adresse, adresse, "rate_spec_mig-conti", nbi, tiss)
    }
    Mratespecmigconti <- Mratespecmigconti + spe
    c <- c + nbi
    if (nbi > 1 && islarch == 0) {
        spe <- d[, c:(c+nbi-1)]
        if (static == 0 && per == 1) {
            plotmat(time, spe, c(100,200,400), adresse, adresse, "rate_spec_mig-isl", nbi, tiss)
        }
        Mratespecmigisl <- Mratespecmigisl + spe
    }
    c <- c + nbi
    spe <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, spe, c(100,200,400), adresse, adresse, "rate_spec_cladog", nbi, tiss)
    }
    Mratespeccladog <- Mratespeccladog + spe
    c <- c + nbi

    popext <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, popext, c(4,200), adresse, adresse, "rate_pop-ext", nbi, tiss)
    }
    Mratepopext <- Mratepopext + popext
    c <- c + nbi

    speextisl <- as.matrix(d[, c:(c+nbi-1)])
    if (static == 0 && per == 1) {
        plotmat(time, speextisl, c(200, 400), adresse, adresse, "rate_sp-ext_isl", nbi, tiss)
    }
    Mratespextisl <- Mratespextisl + speextisl
    c <- c + nbi

    if (nbi > 1 && islarch == 0) {
        speextarch <- matrix(d[, c], length(time), 1)
        if (static == 0 && per == 1) {
            plotmat(time, speextarch, c(200, 400), adresse, adresse, "rate_sp-ext_arch", nbi, tiss)
        }
        Mratespextarch <- Mratespextarch + speextarch
    }
    
    if (nbi > 1 && islarch == 0) {
        ratespecallarch <- as.matrix(apply(ratespecallisl, 1, sum))
        if (static == 0 && per == 1) {
            plotmat(time, ratespecallarch, c(100, 200, 400), adresse, adresse, "rate_spec_all_arch", nbi, tiss)
        }
        Mratespecallarch <- Mratespecallarch + ratespecallarch
    }
    
}


## plot means over simulations replicates
adresse2 <- paste0("results/sim", simulnumber, ".all")
if (nrun > 1) {
    Mdivisl <- Mdivisl / nrun
    Mendem <- Mendem / nrun
    Mpopsizeisl <- Mpopsizeisl / nrun
    Mpresdurisl <- Mpresdurisl / nrun
    Mspecdurisl <- Mspecdurisl / nrun
    Mspsizeisl <- Mspsizeisl / nrun
    Mrateanagconti <- Mrateanagconti / nrun
    Mratecladog <- Mratecladog / nrun
    Mratemigconti <- Mratemigconti / nrun
    Mratepopext <- Mratepopext / nrun
    Mratespecallisl <- Mratespecallisl / nrun
    Mratespeccladog <- Mratespeccladog / nrun
    Mratespecmigconti <- Mratespecmigconti / nrun
    Mratespextisl <- Mratespextisl / nrun


    plotmat(time, Mdivisl, c(1,10,100), adresse, adresse2, "div_isl", nbi, tiss)
    if (!static) {
        plotmat(time, Mendem, c(1,10,100), adresse, adresse2, "endem_all", nbi, tiss)
        plotmat(time, Mpopsizeisl, c(1,10,100), adresse, adresse2, "pop-size_isl", nbi, tiss)
        plotmat(time, Mpresdurisl, c(1,10,100), adresse, adresse2, "pres-dur_isl", nbi, tiss)
        plotmat(time, Mspecdurisl, c(1,10,100), adresse, adresse2, "spec-dur_isl", nbi, tiss)
        plotmat(time, Mspsizeisl, c(1,10,100), adresse, adresse2, "sp-size_isl", nbi, tiss)
        plotmat(time, Mrateanagconti, c(100,200,500), adresse, adresse2, "rate_anag_conti", nbi, tiss)
        plotmat(time, Mratecladog, c(1,10,100), adresse, adresse2, "rate_cladog", nbi, tiss)
        plotmat(time, Mratemigconti, c(100,200,500), adresse, adresse2, "rate_mig_conti", nbi, tiss)
        plotmat(time, Mratepopext, c(1,10,100), adresse, adresse2, "rate_pop-ext", nbi, tiss)
        plotmat(time, Mratespecallisl, c(100,200,500), adresse, adresse2, "rate_spec_all_isl", nbi, tiss)
        plotmat(time, Mratespeccladog, c(100,200,500), adresse, adresse2, "rate_spec_cladog", nbi, tiss)
        plotmat(time, Mratespecmigconti, c(1,10,100), adresse, adresse2, "rate_spec_mig-conti", nbi, tiss)
        plotmat(time, Mratespextisl, c(100,200,500), adresse, adresse2, "rate_sp-ext_isl", nbi, tiss)

    }

    if (nbi > 1 && islarch == 0) {
        Mpopsizearch <- Mpopsizearch / nrun
        Mdivarch <- Mdivarch / nrun
        MendemSIE <- MendemSIE / nrun
        MendemMIE <- MendemMIE / nrun
        Mspsizearch <- Mspsizearch / nrun
        Mpresdurarch <- Mpresdurarch / nrun
        Mspecdurarch <- Mspecdurarch / nrun
        Mratemigisl <- Mratemigisl / nrun
        Mrateanagisl <- Mrateanagisl / nrun
        Mratespecmigisl <- Mratespecmigisl / nrun
        Mratespextarch <- Mratespextarch / nrun
        Mratespecallarch <- Mratespecallarch / nrun
        
        plotmat(time, Mdivarch, c(1,10,100), adresse, adresse2, "div_arch", nbi, tiss)
        if (!static) {
            plotmat(time, Mpopsizearch, c(1,10,100), adresse, adresse2, "pop-size_arch", nbi, tiss)
            plotmat(time, MendemSIE, c(1,10,100), adresse, adresse2, "endem_SIE", nbi, tiss)
            plotmat(time, MendemMIE, c(1,10,100), adresse, adresse2, "endem_MIE", nbi, tiss)
            plotmat(time, Mspsizearch, c(1,10,100), adresse, adresse2, "sp-size_arch", nbi, tiss)
            plotmat(time, Mpresdurarch, c(1,10,100), adresse, adresse2, "pres-dur_arch", nbi, tiss)
            plotmat(time, Mspecdurarch, c(1,10,100), adresse, adresse2, "spec-dur_arch", nbi, tiss)
            plotmat(time, Mratemigisl, c(50,100,250), adresse, adresse2, "rate_mig_isl", nbi, tiss)
            plotmat(time, Mrateanagisl, c(50,100,250), adresse, adresse2, "rate_anag_isl", nbi, tiss)
            plotmat(time, Mratespecmigisl, c(100,200,500), adresse, adresse2, "rate_spec_mig-isl", nbi, tiss)
            plotmat(time, Mratespextarch, c(100,200,500), adresse, adresse2, "rate_sp-ext_arch", nbi, tiss)
            plotmat(time, Mratespecallarch, c(100,200,500), adresse, adresse2, "rate_spec_all_arch", nbi, tiss)
        }
    }
}

## save means over simulations replicates
Mtosave <- list()
Mtosave[[1]] <- Mdivisl
Mtosave[[2]] <- Mendem
Mtosave[[3]] <- Mpopsizeisl
Mtosave[[4]] <- Mpresdurisl
Mtosave[[5]] <- Mspecdurisl
Mtosave[[6]] <- Mspsizeisl
Mtosave[[7]] <- Mrateanagconti
Mtosave[[8]] <- Mratecladog
Mtosave[[9]] <- Mratemigconti
Mtosave[[10]] <- Mratepopext
Mtosave[[11]] <- Mratespecallisl
Mtosave[[12]] <- Mratespeccladog
Mtosave[[13]] <- Mratespecmigconti
Mtosave[[14]] <- Mratespextisl
if (nbi > 1 && islarch == 0) {
    Mtosave[[15]] <- Mpopsizearch
    Mtosave[[16]] <- Mdivarch
    Mtosave[[17]] <- MendemSIE
    Mtosave[[18]] <- MendemMIE
    Mtosave[[19]] <- Mspsizearch
    Mtosave[[20]] <- Mpresdurarch
    Mtosave[[21]] <- Mspecdurarch
    Mtosave[[22]] <- Mratemigisl
    Mtosave[[23]] <- Mrateanagisl
    Mtosave[[24]] <- Mratespecmigisl
    Mtosave[[25]] <- Mratespextarch
    Mtosave[[26]] <- Mratespecallarch
}
dput(Mtosave, file=paste0(adresse2, "_means"))

## for the null model of a static archipelago,
## compute the equilibrium values and save them
## I assume that the equilibrium is reached at around
## the middle of simulations (to be checked by hand!)
## and the equilibrium value is then computed as the
## temporal average of the considered output (already
## averaged time step by time step over replicates) over
## the second half of the simulation
if (static) {
    mid <- floor(nt/2) ## index of the middle of the time of simulation
    Eqdivisl <- apply(as.matrix(Mdivisl[mid:nt,]), 2, mean)
    Eqendem <- apply(as.matrix(Mendem[mid:nt,]), 2, function(x){mean(x,na.rm=TRUE)})
    Eqpopsizeisl <- apply(as.matrix(Mpopsizeisl[mid:nt,]), 2, mean)
    Eqpresdurisl <- apply(as.matrix(Mpresdurisl[mid:nt,]), 2, mean)
    Eqspecdurisl <- apply(as.matrix(Mspecdurisl[mid:nt,]), 2, function(x){mean(x,na.rm=TRUE)})
    Eqspsizeisl <- apply(as.matrix(Mspsizeisl[mid:nt,]), 2, mean)
    Eqrateanagconti <- apply(as.matrix(Mrateanagconti[mid:nt,]), 2, mean)
    Eqratecladog <- apply(as.matrix(Mratecladog[mid:nt,]), 2, mean)
    Eqratemigconti <- apply(as.matrix(Mratemigconti[mid:nt,]), 2, mean)
    Eqratepopext <- apply(as.matrix(Mratepopext[mid:nt,]), 2, mean)
    Eqratespecallisl <- apply(as.matrix(Mratespecallisl[mid:nt,]), 2, mean)
    Eqratespeccladog <- apply(as.matrix(Mratespeccladog[mid:nt,]), 2, mean)
    Eqratespecmigconti <- apply(as.matrix(Mratespecmigconti[mid:nt,]), 2, mean)
    Eqratespextisl <- apply(as.matrix(Mratespextisl[mid:nt,]), 2, mean)
    if (nbi > 1 && islarch == 0) {
        Eqpopsizearch <- apply(as.matrix(Mpopsizearch[mid:nt,]), 2, mean)
        Eqdivarch <- apply(as.matrix(Mdivarch[mid:nt,]), 2, mean)
        EqendemSIE <- apply(as.matrix(MendemSIE[mid:nt,]), 2, function(x){mean(x,na.rm=TRUE)})
        EqendemMIE <- apply(as.matrix(MendemMIE[mid:nt,]), 2, function(x){mean(x,na.rm=TRUE)})
        Eqspsizearch <- apply(as.matrix(Mspsizearch[mid:nt,]), 2, mean)
        Eqpresdurarch <- apply(as.matrix(Mpresdurarch[mid:nt,]), 2, mean)
        Eqspecdurarch <- apply(as.matrix(Mspecdurarch[mid:nt,]), 2, function(x){mean(x,na.rm=TRUE)})
        Eqratemigisl <- apply(as.matrix(Mratemigisl[mid:nt,]), 2, mean)
        Eqrateanagisl <- apply(as.matrix(Mrateanagisl[mid:nt,]), 2, mean)
        Eqratespecmigisl <- apply(as.matrix(Mratespecmigisl[mid:nt,]), 2, mean)
        Eqratespextarch <- apply(as.matrix(Mratespextarch[mid:nt,]), 2, mean)
        Eqratespecallarch <- apply(as.matrix(Mratespecallarch[mid:nt,]), 2, mean)
    }
    Eqtosave <- list()
    Eqtosave[[1]] <- Eqdivisl
    Eqtosave[[2]] <- Eqendem
    Eqtosave[[3]] <- Eqpopsizeisl
    Eqtosave[[4]] <- Eqpresdurisl
    Eqtosave[[5]] <- Eqspecdurisl
    Eqtosave[[6]] <- Eqspsizeisl
    Eqtosave[[7]] <- Eqrateanagconti
    Eqtosave[[8]] <- Eqratecladog
    Eqtosave[[9]] <- Eqratemigconti
    Eqtosave[[10]] <- Eqratepopext
    Eqtosave[[11]] <- Eqratespecallisl
    Eqtosave[[12]] <- Eqratespeccladog
    Eqtosave[[13]] <- Eqratespecmigconti
    Eqtosave[[14]] <- Eqratespextisl
    if (nbi > 1 && islarch == 0) {
        Eqtosave[[15]] <- Eqpopsizearch
        Eqtosave[[16]] <- Eqdivarch
        Eqtosave[[17]] <- EqendemSIE
        Eqtosave[[18]] <- EqendemMIE
        Eqtosave[[19]] <- Eqspsizearch
        Eqtosave[[20]] <- Eqpresdurarch
        Eqtosave[[21]] <- Eqspecdurarch
        Eqtosave[[22]] <- Eqratemigisl
        Eqtosave[[23]] <- Eqrateanagisl
        Eqtosave[[24]] <- Eqratespecmigisl
        Eqtosave[[25]] <- Eqratespextarch
        Eqtosave[[26]] <- Eqratespecallarch
    }
    dput(Eqtosave, file=paste0(adresse2, "_equilibrium"))
}

