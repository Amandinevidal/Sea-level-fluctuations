simulnumber <- 1          ## integer; simulation named "simX"; ignored for static models
beta <- 1                 ## per capita birth rate
gamma <- 0.1              ## strength of density dependent mortality
sigma <- 5e-5             ## per capita rate of cladogenesis
alpha <- 100               ## magnitude of the effect of gene flow on speciation
tau <- 500                 ## min delay for speciation
mui <- 1e-1                ## per capita migration rate between islands
ci <- 10                   ## mean dispersal distance for individuals between islands
muc <- 1                   ## per capita migration rate from the continent
cc <- 25                   ## mean dispersal distance for individuals from the continent
Nc <-  1e6                 ## number of individuals on the continent
theta <- 20                ## fundamental biodiversity number
nbi <- 4                   ## number of islands in the archipelago
islarch <- 0               ## 0=archipel; 1/2=single island archipel equivalent (null models)
Kmax <- rep(40e3, nbi)     ## vector; maximum carrying capacity for each island
ti <- 5e4*(0:(nbi-1))      ## vector; time when each island emerges
pmax <- rep(0.2, nbi)      ## vector; fraction of island life time at which size is max
tl <- rep(20e4, nbi)       ## vector; life time of each island
dc <- rep(150, nbi)        ## vector; initial distance of each island from the continent
di <- matrix(22,nbi,nbi)   ## matrix; initial distances between islands
tiss <- 25                 ## time interval between each saved state
wtup <- 1                  ## waiting time between updates of the landscapes and rates
nrun <- 15                 ## number of simulation replicates
per <- 0                   ## 1=plot graphs for each run; 0=don't
rngseed <- 0               ## seed used for random number; 0 means seed=time
static <- 0                ## 0=dynamic model; 1=static model, teq & tint only for this case
teq <- 10e4                ## duration until equilibrium is estimated to be achieved
tint <- 1e4	            	 ## time interval between two equilibria estimates
phi_max <- 79*3.141592/180 ## maximum angle for apex of the cone (80 degres) 3.141592/4
phi_min <- 10*3.141592/180 ## minimum apex angle of the island's cone  
d <- -log(0.2)/(0.8*tl[1]) ## increase coefficient for theta function, with 80% (1-0.9) of its maximum opening at 80% of island's lifetimex 
amp <- 0.15              ## amplitude of sea level variations (sinusoidal function)
period <- 5000             ## period of sea level variations (sinusoidal function) moitié moins
psi <-  5*(3.141592/180)   ## mainland slope (here 5 degrees converted in radians)
h_breakpoint <- 4          ## mainland height (here > max elevation of islands ~ 3)