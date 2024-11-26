##### ------------- TEST ----------------- #####

#Occurence_data = Sampled_Pre_Abs_Data[[x]][[y]]
#Phylo_tree = List_Phylo[[x]]
#PW = PW_Obs[PW_Obs$tree_NM == x,-2]
#Vector <- Occurence_data[1,]
# 
#TEST <- Ex_quad_entropy(Occurence_data = Occurence_data, Phylo_tree = Phylo_tree, PW = PW, FUN = sum)

# Occurence_data = Sampled_Pre_Abs_Data[[x]][[y]][,-(1:2)]
# Phylo_tree = List_Phylo[[x]]
# PW = Obs_PW[[x]][,-(1:2)]
# FUN = sum
# Vector <- Occurence_data[1,]

##### ------------ TEST ------------------ ######

# Install/load tidyverse
suppressPackageStartupMessages(if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)})

Ex_quad_entropy <- function(Occurence_data, Phylo_tree, PW ,FUN = sum){
  
  ### Create a function to apply a function (Sum here) on a subset of all the Patristic Pairwise distances
  F_ex_quad_ent <- function(Vector,                # A vector of species abundance or occurrence within a plot
                            PW,                    # A global matrix of PW distance
                            FUN2 = FUN){           # A function to apply on the results
  
  # Find the species present in the plot
  Sp.names <- names(Vector[which(Vector > 0)])
  # Sp.names <- Sp.names[- which(Sp.names == "rep" | Sp.names == "sample" )] # Remove eventually a element that is not a species name.
  
  # Extract the PW distance from the global dataset
  PW_distances <- PW[rownames(PW) %in% Sp.names,colnames(PW) %in% Sp.names]
  
  # Get only a triangular matrix and transform it into a vector (We keep the diagonal)
  PW_distances[upper.tri(PW_distances,diag = FALSE)] <- NA # Make the matrix triangular
  PW_distances <- na.omit(unlist(PW_distances))            # Make a vector from only numeric values
  
  # Apply the function "FUN" entered
  result <- FUN2(PW_distances)
   
  # Return the value to the outer function
  return(result)
  }
  
  # Apply the inner function to the community dataset
  out <- apply(Occurence_data, MARGIN = 1, F_ex_quad_ent, PW = PW )
  
  # Return the results
  return(out)
}
