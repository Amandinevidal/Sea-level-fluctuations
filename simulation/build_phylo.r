#### ---------------------------Informations--------------------------####
#### Date : 2025-01-14
#### Author : Amandine Vidal-Hosteng
#### Encoding : UTF-8
#### Email : amandine.vidalhosteng@gmail.com
####
#### this script builds the phylogenies of existing species on the 
#### archipelago at 100-generation intervals, which are stored in the 
#### file ending with "sim1.irun_insular_tree.txt" (irun being the simulation
#### replicate number); the subphylogenies of continental species whose 
#### descendants are on the archipelago are also recorded in the file 
#### "sim1.irun_continental_tree.txt"
#### ------------------------------------------------------------------ #

library(dplyr)
source("parameters.r")
source("build_phylo_functions.r")

for (irun in 0:(nrun-1)) { #### 1. Irun loop ####
  
  # phylo and metapop files 
  phylo <- phylo_table(irun)
  first_endemic_species <- phylo[1, 1] # save the id of the first endemic species saved 
  times = max(phylo[,3])
  metapop = read.table(paste("results/sim1.",irun,"_metapop", sep=""))
  
  # continental tree
  continental_tree <- conti_tree(first_endemic_species)
  
  # remove continental spc not observed on the archipelago from phylo file
  phylo <- create_phylo_with_continental_species(metapop,phylo,first_endemic_species) 
  
  # save phylo
  write.table(phylo, paste("results/sim1.",irun,"_phylo.txt", sep=""), col.names=F, row.names=F) # save phylo file
  
  # create files
  file.create(paste("results/sim1.",irun,"_trees_error.txt",sep=""))
  file.create(paste("results/sim1.",irun,"_insular_tree.txt",sep=""))
  file.create(paste("results/sim1.",irun,"_continental_tree.txt", sep=""))
  
  for (t in seq(1,times,by=100)) { #### 2. Time loop  ####
    
    # remove non-born species
    table <- update_species_table(phylo,t,times)
    # keep only extant species
    table_extant <- filter_living_species(table, t)
    # who are the continental species that are present on the arch in table extant
    table_extant <- identify_continental_species(table_extant,metapop,t,first_endemic_species)
    
    # data spc
    table_extant_spc <- data.frame(
      lineage=table_extant[,1],
      ancestor_of_lineage=table_extant[,1],
      mother=table_extant[,2],
      birth_time=table_extant[,3],
      death_time=table_extant[,4],
      last_daughter_registered=rep(0,length(table_extant[,1])),
      last_daughter_life_time=rep(0,length(table_extant[,1]))
    )
    
    # data spc without continental species
    table_extant_spc_withoutconti <- data.frame(
      table_extant_spc[which(as.numeric(table_extant_spc[,2])<first_endemic_species),]
    )
    
    # data trees without continental species
    table_trees_extant_spc_withoutconti <- table_trees_extant_spc <- data.frame(
      lineage = character(),
      ancestor_of_lineage = numeric(),
      mother=numeric(),
      birth_time=numeric(),
      death_time=numeric(),
      last_daughter_registered=numeric(),
      last_daughter_life_time=numeric()
    )
    
    if (nrow(table_extant[table_extant$lineage>=first_endemic_species,])>2) {
      
      while (nrow(table_extant_spc)>0) { # for each endemic species in table
        
        # daughter info 
        daughter_line = which(table_extant_spc[,4]==max(table_extant_spc[,4])) # who is the last species that appeared = daughter
        if (length(daughter_line)>1) { # if some species appeared at the same time
          daughter_line=daughter_line[length(daughter_line)] # take the first species of them
        }
        
        if (as.numeric(table_extant_spc[daughter_line,2]) < first_endemic_species) { # daughter species is continental
          
          daughter_line_withoutconti = which(table_extant_spc_withoutconti[,2]==table_extant_spc[daughter_line,2])
          
          if (as.character(table_extant_spc[daughter_line,1])==as.character(table_extant_spc[daughter_line,2])) {
            line <- data.frame(table_extant_spc[daughter_line,])
            names(line) <- colnames(table_trees_extant_spc)
            table_trees_extant_spc <- rbind(table_trees_extant_spc,line)
          } else {
            line <- data.frame(table_extant_spc[daughter_line,])
            names(line) <- colnames(table_trees_extant_spc)
            table_trees_extant_spc <- rbind(table_trees_extant_spc,line)
            line <- data.frame(table_extant_spc_withoutconti[daughter_line_withoutconti,])
            names(line) <- colnames(table_trees_extant_spc_withoutconti)
            table_trees_extant_spc_withoutconti <- rbind(table_trees_extant_spc_withoutconti,line)
          }
          table_extant_spc <- table_extant_spc[-daughter_line,] # rem daughter species in table_all_spc
          table_extant_spc_withoutconti <- table_extant_spc_withoutconti[-daughter_line_withoutconti,] # rem daughter species in table_all_spc
          
        } else { # daughter species is endemic
          
          if (as.numeric(table_extant_spc[daughter_line,3]) >= first_endemic_species) { # mother is endemic
            
            if (length(which(table_extant_spc[,2]==table_extant_spc[daughter_line,3]))==1) { # mother is extant
              
              mother_line=which(table_extant_spc[,2]==table_extant_spc[daughter_line,3]) # mother line
              
              if (as.character(table_extant_spc[daughter_line,1])==as.character(table_extant_spc[daughter_line,2])) { # no daughter subtree
                
                if (as.character(table_extant_spc[mother_line,1])==as.character((table_extant_spc[mother_line,2]))) { # no mother subtree
                  table_extant_spc[mother_line,1]= paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                } else  { # mother subtree
                  table_extant_spc[mother_line,1]= paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[mother_line,7]),")",sep="")
                }
                
              } else { # daughter subtree
                
                if (as.character(table_extant_spc[mother_line,1])==as.character((table_extant_spc[mother_line,2]))) { # no mother subtree
                  table_extant_spc[mother_line,1]= paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                } else  { # mother subtree
                  table_extant_spc[mother_line,1]= paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[mother_line,7]),")",sep="")
                }
                
              } # end type of insertion
              
              table_extant_spc[mother_line,6]=table_extant_spc[daughter_line,2]
              table_extant_spc[mother_line,7]=table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]
              table_extant_spc <- table_extant_spc[-daughter_line,] # rem daughter species in table_extant_spc
              
            } else { # mother is extinct
              
              # mother species and its descents are extinct, replace the mother species by its own ancestor
              table_extant_spc[daughter_line,4] <- phylo[which(phylo[,1]==table_extant_spc[daughter_line,3]),3] # replace daughter' birth by mother's birth
              table_extant_spc[daughter_line,3] <- phylo[which(phylo[,1]==table_extant_spc[daughter_line,3]),2] # replace daughter's mother by grand mother
              
            } # end extinct endemic mother
            
          } else { # mother is continental
            
            if (length(which(table_extant_spc[,2]==table_extant_spc[daughter_line,3]))==1) { # continental mother is extant
              
              mother_line=which(table_extant_spc[,2]==table_extant_spc[daughter_line,3]) # mother line
              mother_line_withoutconti=which(table_extant_spc_withoutconti[,2]==table_extant_spc[daughter_line,3])
              
              if (as.character(table_extant_spc[daughter_line,1])==as.character(table_extant_spc[daughter_line,2])) { # no daughter subtree
                
                if (as.character(table_extant_spc[mother_line,1])==as.character((table_extant_spc[mother_line,2]))) { # no mother subtree
                  
                  table_extant_spc[mother_line,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                  
                  table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste(table_extant_spc[daughter_line,1],sep="")
                  
                } else { # mother subtree
                  
                  table_extant_spc[mother_line,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[mother_line,7]),")",sep="")
                  
                  if (as.character(table_extant_spc_withoutconti[mother_line_withoutconti,1])==as.character(table_extant_spc_withoutconti[mother_line_withoutconti,6])) { # mother subtree = 1 species
                    table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc_withoutconti[mother_line_withoutconti,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                  } else { # mother subtree = >1 species
                    table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),",",table_extant_spc_withoutconti[mother_line_withoutconti,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc_withoutconti[mother_line_withoutconti,7]),")",sep="")
                  }
                  
                }
                
                table_extant_spc_withoutconti[mother_line_withoutconti,6]=table_extant_spc[daughter_line,2]
                table_extant_spc_withoutconti[mother_line_withoutconti,7]=table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]
                
              } else { # daughter subtree
                
                if (as.character(table_extant_spc[mother_line,1])==as.character((table_extant_spc[mother_line,2]))) { # no mother subtree
                  
                  table_extant_spc[mother_line,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                  table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste(table_extant_spc[daughter_line,1],sep="")
                  table_extant_spc_withoutconti[mother_line_withoutconti,6:7]=table_extant_spc[daughter_line,6:7]
                  
                } else { # mother subtree
                  
                  table_extant_spc[mother_line,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc[mother_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[mother_line,7]),")",sep="")
                  
                  if (as.character(table_extant_spc_withoutconti[mother_line_withoutconti,1])==as.character(table_extant_spc_withoutconti[mother_line_withoutconti,6])) { # mother subtree = 1 species
                    table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc_withoutconti[mother_line_withoutconti,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]),")",sep="")
                  } else {
                    table_extant_spc_withoutconti[mother_line_withoutconti,1]=paste("(",table_extant_spc[daughter_line,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc[daughter_line,7]),",",table_extant_spc_withoutconti[mother_line_withoutconti,1],":",as.character(table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]-table_extant_spc_withoutconti[mother_line_withoutconti,7]),")",sep="")
                  }
                  
                  table_extant_spc_withoutconti[mother_line_withoutconti,6]=table_extant_spc[daughter_line,2]
                  table_extant_spc_withoutconti[mother_line_withoutconti,7]=table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]
                  
                }
                
              }
              
              table_extant_spc[mother_line,6]=table_extant_spc[daughter_line,2]
              table_extant_spc[mother_line,7]=table_extant_spc[daughter_line,5]-table_extant_spc[daughter_line,4]
              table_extant_spc <- table_extant_spc[-daughter_line,] # rem daughter species in table_extant_spc
              
            } else { # continental mother is extinct
              
              if (length(which(table_extant_spc[,3]==table_extant_spc[daughter_line,3]))>1) {
                lines <- table_extant_spc[which(table_extant_spc[,3]==table_extant_spc[daughter_line,3]),]
                new_mother_line <- lines[which(lines[,4]==min(lines[,4])),]
                if (dim(new_mother_line)[1]>1) {
                  new_mother_line = new_mother_line[1,]
                }
                mother_line = which(table_extant_spc[,2]==new_mother_line[,2])
                table_extant_spc[mother_line,6]=table_extant_spc[mother_line,2]
                # table_extant_spc[mother_line,7] <- table_extant_spc[mother_line,5]-table_extant_spc[mother_line,4] modif 1
                table_extant_spc[mother_line,2]=table_extant_spc[daughter_line,3]
                table_extant_spc[mother_line,4] <- 0 # MODIF 
                table_extant_spc_withoutconti<- rbind(table_extant_spc_withoutconti,table_extant_spc[mother_line,])
              } else {
                line <- table_extant_spc[daughter_line,]
                names(line) <- colnames(table_trees_extant_spc)
                table_trees_extant_spc <- rbind(table_trees_extant_spc,line)
                table_trees_extant_spc_withoutconti <- rbind(table_trees_extant_spc_withoutconti,line)
                table_extant_spc <- table_extant_spc[-daughter_line,] # rem daughter species in table_all_spc
              }
              
            } # end continental mother is extant
            
          } # end mother is continental
          
        } # end daughter spc is endemic
        
      } # end while (nrow(table_extant_spc)>0)
      
      # is there continental species that have colonie multple time ?
      duplicate_mother <- table_trees_extant_spc_withoutconti %>%
        filter(duplicated(mother)) %>%
        pull(mother) %>%
        unique()
      if(length(duplicate_mother)[1]>0){
        for(d in duplicate_mother){
          new_tree_newick <- assemble_trees(d) # combine mother trees
          new_tree_phylo <- ape::read.tree(text=paste0(new_tree_newick,";")) # convert into phylo object
          tree_life_time <- max(ape::branching.times(new_tree_phylo)) # branch total length
          new_line <- data.frame(lineage=new_tree_newick, ancestor_of_lineage=d, mother=d, birth_time=0, death_time=t, last_daughter_registered=0,last_daughter_life_time=tree_life_time)
          table_trees_extant_spc_withoutconti <- table_trees_extant_spc_withoutconti[which(table_trees_extant_spc_withoutconti$mother!=d),]
          table_trees_extant_spc_withoutconti <- rbind(table_trees_extant_spc_withoutconti,new_line)
        }
      }
      
      branches = nrow(table_trees_extant_spc_withoutconti)
      if (branches > 1) {
        # build phylogeny
        bp_continental_tree_t <<- NULL
        phylogeny_extant_spc_withoutconti <- build_tree(table_trees_extant_spc_withoutconti,continental_tree)
        tree_extant_spc_withoutconti = ape::write.tree(phylogeny_extant_spc_withoutconti)
        bp_continental_tree_t <- data.frame(t=t,tree=ape::write.tree(bp_continental_tree_t))
        write.table(bp_continental_tree_t,file=paste("results/sim1.",irun,"_continental_tree.txt", sep=""),col.names=F,row.names=F,append=T)
        
        # check if there are errors in the tree
        tree_check(table_extant,first_endemic_species,phylogeny_extant_spc_withoutconti)
        
      } else {
        bp_continental_tree_t <- data.frame(t=t,tree="NULL")
        write.table(bp_continental_tree_t,file=paste("results/sim1.",irun,"_continental_tree.txt", sep=""),col.names=F,row.names=F,append=T)
        tree_extant_spc_withoutconti <- paste0(table_trees_extant_spc_withoutconti$lineage[1],";")
        phylogeny_extant_spc_withoutconti <- ape::read.tree(text=tree_extant_spc_withoutconti)
        tree_check(table_extant,first_endemic_species,phylogeny_extant_spc_withoutconti)
      }
      
    } else { # if no rows of species, or only one species, no tree
      
      tree_extant_spc_withoutconti=c("NULL")
      bp_continental_tree_t <- data.frame(t=t,tree="NULL")
      write.table(bp_continental_tree_t,file=paste("results/sim1.",irun,"_continental_tree.txt", sep=""),col.names=F,row.names=F,append=T)
      
    } # end if (nrow(table_extant_spc)>2)
    
    # save phylogeny
    save_tree_extant_species=data.frame(t=t,phylo=tree_extant_spc_withoutconti)
    write.table(save_tree_extant_species,file=paste("results/sim1.",irun,"_insular_tree.txt", sep=""),col.names=F,row.names=F,append=T)
    
  } # end of time loop
  
}

