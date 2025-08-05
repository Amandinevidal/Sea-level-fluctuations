#### ---------------------------Informations--------------------------####
#### Date : 2025-01-14
#### Author : Amandine Vidal-Hosteng
#### Encoding : UTF-8
#### Email : amandine.vidal-hosteng@univ-tlse3.fr
####
#### this file contains all functions needed in the build_phylo scripts
#### ------------------------------------------------------------------ #

# reorganize and rename phylo data
if (static == 0) {
  phylo_table <- function(irun) {
    file_path <- paste("results/sim1.", irun, "_phylo", sep = "")
    phylo <- read.table(file_path)
    colnames(phylo) <- c("lineage", "date_origin", "mother", "birth", "death")
    phylo <- phylo[,-2]
    phylo[,3:4] <- round(phylo[,3:4])
    return(phylo)
  }
} else {
  phylo_table <- function(jrun, irun) {
    file_path <- paste("results/sim",jrun,".", irun, "_phylo", sep = "")
    phylo <- read.table(file_path)
    colnames(phylo) <- c("lineage", "date_origin", "mother", "birth", "death")
    phylo <- phylo[,-2]
    phylo[,3:4] <- round(phylo[,3:4])
    return(phylo)
  }
}

# build continental phylogeny
conti_tree <- function(first_endemic_species,t){
  continental_species = seq(0,first_endemic_species-1,by=1) 
  conti <- matrix(continental_species,ncol=4,nrow=length(continental_species)) # create a matrix same structure as phylo
  colnames(conti) <- c("lineage","mother","birth","death")
  arbre_continental <- ape::rtree(nrow(conti),tip.label = conti[,1])
  node_table = tidytree::as_tibble(arbre_continental)
  mil_time = mil_time = 1/1000 # put all length to 25
  for (n in 1:length(node_table$parent)) {
    if (is.na(node_table$label[n])) {node_table$branch.length[n]=mil_time}
  }
  arbre_continental = ape::as.phylo(node_table)
  return(arbre_continental)
}

# add continental species into phylo data
create_phylo_with_continental_species <- function(metapop, phylo, first_endemic_species) {
  continental_species <- seq(0, first_endemic_species - 1, by = 1) # continental species seq
  conti <- matrix(continental_species, ncol = 4, nrow = length(continental_species))
  colnames(conti) <- c("lineage", "mother", "birth", "death")
  remove <- c() # unobserved spc vector
  for (i in 1:length(continental_species)) { # save birth and death time of each continental species
    if (length(which(metapop[, 3] == continental_species[i])) > 0) { # if observed
      conti[i, 3] <- 0
      conti[i, 4] <- max(metapop[which(metapop[, 3] == continental_species[i]), 1])
    } else { # if not
      remove <- c(remove, i)
    }
  }
  if (length(remove) > 0) {
    conti <- conti[-remove, ]
  }
  phylo <- rbind(conti,phylo) # add continental species in phylo
  return(phylo)
}

# filter only species that existed on the archipelago into phylo data
update_species_table <- function(phylo, t, times) {
  table <- phylo
  if (t != times) {                         # if the experience time isn't over
    if (length(which(table[,3]>t))>0) {     # species which birth after times[t]
      table = table[-(which(table[,3]>t)),] # remove these species
    }
  }
  table[which(table[,4]>t),4] = t # set death time of species which die after times[t] = times[t]
  return(table)
}

# filter only species that are extant at t
filter_living_species <- function(table, t) {
  # filter
  living_spc <- table[table[, 3] <= t & table[, 4] >= t, ]
  # new table
  table_extant <- data.frame(lineage = living_spc[, 1],
                             mother = living_spc[, 2],
                             birth = living_spc[, 3],
                             death = living_spc[, 4])
  return(table_extant)
}

# filter continental species that migrated on the archipelago
identify_continental_species <- function(table_extant, metapop, t, first_endemic_species) {
  continental <- table_extant[which(table_extant[,1]<first_endemic_species),1] # who is continental
  # remove continental species that are not present in the archipelago at t
  time_obs_1 <- t-1
  time_obs_2 <- t+24
  metapop_obs <- rbind(metapop[which(metapop[,1]==time_obs_1),],metapop[which(metapop[,1]==time_obs_2),])
  for (i in 1:length(continental)){
    lines <- which(metapop_obs[,3]==continental[i])
    line <- which(table_extant[,1]==continental[i])
    if (length(lines)==0) {
      table_extant <- table_extant[-line,]
    }
  }
  return(table_extant)
}

# transform a single species lineage into Newick format
parse_single_species_tree <- function(lineage, t) {
  if (startsWith(as.character(lineage), "(")) {
    return(paste0(lineage))
  } else {
    return(paste0(lineage, ":",t))
  }
}

# add coma at the end of the phylogeny for Newick format
add_coma <- function(newick) {
  if (substr(newick, nchar(newick), nchar(newick)) != ";") {  # Vérifie si le point-virgule est absent à la fin
    newick <- paste0(newick, ";")  # Ajoute un point-virgule à la fin
  }
  return(newick)
}

# combine trees
combine_trees <- function(tree1_newick, tree2_newick) {
  tree1 <- read.tree(text = tree1_newick)
  tree2 <- read.tree(text = tree2_newick)
  combined_tree <- bind.tree(tree1, tree2, where = "root")
  return(combined_tree)
}

# combine insular trees from the same continental species
assemble_trees <- function(mom) {
  subset_mother <- table_trees_extant_spc_withoutconti %>%
    filter(mother == mom) %>%
    arrange(birth_time) %>%
    mutate(lineage = mapply(parse_single_species_tree, lineage, t))
  first_tree <- subset_mother$lineage[1]
  for (i in 2:nrow(subset_mother)) {
    second_tree <- subset_mother$lineage[i]
    ifelse(substr(first_tree, nchar(first_tree), nchar(first_tree)) == ")", type1 <- "tree", type1 <- "single")
    ifelse(substr(second_tree, nchar(second_tree), nchar(second_tree)) == ")", type2 <- "tree", type2 <- "single")
    if (type1 == "tree" && type2 == "tree") { 
      first_phylo <- ape::read.tree(text=paste0(first_tree,";"))
      first_t <- t-max(ape::branching.times(first_phylo))
      second_phylo <- ape::read.tree(text=paste0(subset_mother$lineage[i],";"))
      second_t <- t-max(ape::branching.times(second_phylo))
      first_tree <-paste0("(",first_tree,":",first_t,",",second_tree,":",second_t,")")
    } 
    else if (type1 == "tree" && type2 == "single") {
      first_phylo <- ape::read.tree(text=paste0(first_tree,";"))
      first_t <- t-max(ape::branching.times(first_phylo))
      first_tree <- paste0("(",first_tree,":",first_t,",",second_tree,")")
    }
    else if (type1 =="single" && type2 == "tree") {
      second_phylo <- ape::read.tree(text=paste0(subset_mother$lineage[i],";"))
      second_t <- t-max(ape::branching.times(second_phylo))
      first_tree <- paste0("(",second_tree,":",second_t,",",first_tree,")")
    }
    else {
      first_tree <- paste0("(",first_tree,",",second_tree,")")
    }
  }
  new_tree <- first_tree
  return(new_tree)
}

# combine insular trees with continental tree
build_tree <- function(table_trees_extant_spc_withoutconti,continental_tree){
  branches_selectionnees <- c(table_trees_extant_spc_withoutconti$mother)
  branches_selectionnees <- as.character(branches_selectionnees)
  arbre_selectionne <- ape::keep.tip(continental_tree, tip = branches_selectionnees)
  bp_continental_tree_t <<- arbre_selectionne
  node_table = tidytree::as_tibble(arbre_selectionne)
  
  # attach each row to its corresponding branch
  for (i in table_trees_extant_spc_withoutconti$mother) {
    branch_to_connect <- table_trees_extant_spc_withoutconti[which(table_trees_extant_spc_withoutconti[,3]==i),]
    tip_to_connect <- which(node_table$label==i)
    if (branch_to_connect[,1]==branch_to_connect[,2]|branch_to_connect[,1]==branch_to_connect[,6]) { # 1 spc to attach
      node_table[tip_to_connect,]$branch.length = as.numeric(branch_to_connect[,5]) # false spc life time but need to be ultrametric bcs of coalescence
      node_table[tip_to_connect,]$label <- as.character(table_trees_extant_spc_withoutconti[which(table_trees_extant_spc_withoutconti[,3]==i),1])
      new_tree = ape::as.phylo(node_table)
    } else { # subtree to attach
      node_table[tip_to_connect,]$branch.length = as.numeric(branch_to_connect[,5])-as.numeric(branch_to_connect[,7])# spc life time
      random_tree = ape::as.phylo(node_table)
      lineages <- branch_to_connect[,1]
      tree_newick <- add_coma(lineages)
      tree <- ape::read.tree(text=tree_newick)
      new_tree = ape::bind.tree(random_tree,tree,where=node_table[tip_to_connect,]$node)
    }
    random_tree = new_tree
    node_table = tidytree::as_tibble(random_tree)
  }
  
  return(random_tree)
}

# check if errors in reconstruction
if (static == 0) {
  tree_check <- function(table_extant,first_endemic_species,phylogeny_extant_spc_withoutconti) {
    # create table_extant without continental species
    lines = which(table_extant[,1]<first_endemic_species)
    if (length(lines)!=0) {
      table_extant_withoutconti = table_extant[-lines,]
    } else {
      table_extant_withoutconti = table_extant
    }
    if (!is.null(phylogeny_extant_spc_withoutconti)) {
      # Dichotomous phylogeny
      binary <- ape::is.binary(phylogeny_extant_spc_withoutconti)[1]
      if (binary == FALSE) {
        trees_error <- print(paste("Error : t = ",t,"tree",i,": Not dichotomous"))
        write.table(trees_error, file=paste("results/sim1.",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
      }
      # Number of species
      tree_nb_spc = ape::Ntip(phylogeny_extant_spc_withoutconti)
      data_nb_spc = length(table_extant_withoutconti[,1])
      data = table_extant_withoutconti[,1]
    }
    # Data and trees taxa matching
    checking <- setdiff(data,as.numeric(phylogeny_extant_spc_withoutconti$tip.label))
    if (tree_nb_spc != data_nb_spc) {
      trees_error <- print(paste("Error t = ",t,"tree",i," : Not all species recorded"))
      write.table(trees_error, file=paste("results/sim1.",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
    }
    if (length(checking)>0) {
      trees_error <- print(paste("Error : t = ",t,"tree",i," : Tree branch species are not matching data"))
      write.table(trees_error, file=paste("results/sim1.",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
    }
  }
} else {
  tree_check <- function(table_extant,first_endemic_species,phylogeny_extant_spc_withoutconti) {
    # create table_extant without continental species
    lines = which(table_extant[,1]<first_endemic_species)
    if (length(lines)!=0) {
      table_extant_withoutconti = table_extant[-lines,]
    } else {
      table_extant_withoutconti = table_extant
    }
    if (!is.null(phylogeny_extant_spc_withoutconti)) {
      # Dichotomous phylogeny
      binary <- ape::is.binary(phylogeny_extant_spc_withoutconti)[1]
      if (binary == FALSE) {
        trees_error <- print(paste("Error : t = ",t,"tree",i,": Not dichotomous"))
        write.table(trees_error, file=paste("results/sim",jrun,".",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
      }
      # Number of species
      tree_nb_spc = ape::Ntip(phylogeny_extant_spc_withoutconti)
      data_nb_spc = length(table_extant_withoutconti[,1])
      data = table_extant_withoutconti[,1]
    }
    # Data and trees taxa matching
    checking <- setdiff(data,as.numeric(phylogeny_extant_spc_withoutconti$tip.label))
    if (tree_nb_spc != data_nb_spc) {
      trees_error <- print(paste("Error t = ",t,"tree",i," : Not all species recorded"))
      write.table(trees_error, file=paste("results/sim",jrun,".",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
    }
    if (length(checking)>0) {
      trees_error <- print(paste("Error : t = ",t,"tree",i," : Tree branch species are not matching data"))
      write.table(trees_error, file=paste("results/sim",jrun,".",irun,"_trees_error.txt", sep=""),col.names=F,row.names=F,append=T)
    }
  }
}

# is the tree perfectly dichotomic ?
check_binarity <- function(tree) {
  internal_nodes <- unique(tree$edge[,1])
  non_binary_nodes <- list()
  for (node in internal_nodes) {
    descendants <- tree$edge[tree$edge[,1] == node, 2]
    if (length(descendants) != 2) {
      non_binary_nodes <- c(non_binary_nodes, list(node))
    }
  }
  return(non_binary_nodes)
}