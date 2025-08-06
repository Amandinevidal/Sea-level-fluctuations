
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Archipelago with sea-level fluctuations IBM

This repository implements a stochastic, individual-based, and neutral
model inspired by the Aguilée et al. (2021) archipelago framework. The
model simulates a dynamic archipelago composed of multiple islands that
successively emerge and submerge along with sinusoidal sea-level
fluctuations. Each island follows a hump-shaped ontogeny (as described
in the General Dynamic Model, Whittaker et al. 2008): it grows in area
after emergence, then gradually erodes. This erosion implies that
sea-level fluctuations have a greater impact in later island stages due
to the gentler slopes formed over time. We simulate neutral ecological
dynamics at the individual level, tracking all individuals to
reconstruct complete phylogenies from which we extract extant
phylogenies. From these, we compute the gamma statistic (speed of
lineage accumulation) and the Sackin index (tree imbalance).

## Contents

- [:file_folder: Simulation](R)
  - [Cfunctions.C](C) - Gathers the functions constituting the model and
    written in C
  - [Rfunctions.r](R) - Gathers the functions called by the main file
    and written in R
  - [build_phylo.r](R) - Builds the phylogenies of existing species on
    the archipelago at 100 generation intervals
  - [build_phylo_functions.r](R) - Gathers the functions called in
    build_phylo.r
  - [parameters.r](R) - Gathers the parameters values
  - [plotresults.r](R) - Computes and plot the outputs characterizing
    the state and evolution of the population
  - [plottopo.r](R) - Computes and plot the outputs characterizing the
    landscape parameters
  - [simulation.r](R) - Is the main program file. This is the one to run
  - [stat_phylo.r](R) - Computes the phylogenetic metrics :warning:
    *Caution: This R script is **parallelized**. Make sure your system
    supports parallel execution and adjust the script according to your
    system (check ‘ncpus’).*

## References

Aguilée, R., Pellerin, F., Soubeyrand, M., Choin, J., & Thébaud, C.
(2021). Biogeographic drivers of community assembly on oceanic islands:
The importance of archipelago structure and history. Journal of
Biogeography, 48 (10), 2616–2628. <https://doi.org/10.1111/jbi.14228>  
Whittaker, R. J., Triantis, K. A., & Ladle, R. J. (2008). A general
dynamic theory of oceanic island biogeography. Journal of Biogeography,
35 (6), 977–994. <https://doi.org/10.1111/j.1365-2699.2008.01892.x>

### Licenses

**Code :** See the [DESCRIPTION](DESCRIPTION) file
