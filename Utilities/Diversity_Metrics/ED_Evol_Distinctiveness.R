##### ------------- TEST ----------------- #####


# Occurence_data <- PreAbs_1[[i]]
# Vector <- Occurence_data[1,] 
# # # Get the Global ED 
# Evol_distinct <- ED_NM[ED_NM$rep == i,c(1,3)]
# # 
# # Give it a try
# TEST <- ED_Tucker(Occurence_data = Occurence_data, Evol_distinct = Evol_distinct)


##### ------------ TEST ------------------ ######


ED_Tucker <- function(Occurence_data , Evol_distinct, FUN = sum){
  
  # Verify that the needed is present and loaded
  if (!require("dplyr")) install.packages("dplyr")
  library(dplyr)

### Create a function to apply a function a function on a subset of all the Evolutionary Distinctiveness
ED_metric <- function(Vector,               # A vector of species abundance or occurrence within a plot
                      EV = Evol_distinct,    # A global vector of species Evolutionary distinctiveness
                      FUN2 = FUN){           # A function to apply on the results
  
  
  # Find the species present in the plot
  Sp.names <- names(Vector[which(Vector > 0)])
  # Extract the ED from the global dataset
  ED <- Evol_distinct[which(Evol_distinct$Species %in% Sp.names),2]
  
  # Apply the function "FUN" entered
  result <- FUN2(ED)
  
  # Return the value to the outer function
  return(result)
}

# Apply the inner function to the community dataset
out <- apply(Occurence_data, MARGIN = 1, ED_metric)

# Return the results
return(out)
}
