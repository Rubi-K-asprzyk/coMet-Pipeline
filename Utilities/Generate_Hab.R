Generate_Hab <- function(
    L,          # Length of the grid
    W,          # Width of the grid
    sq_x,       # sq_x (Number of squares along the x axis)
    sq_y,       # sq_y (Number of squares along the y axis)
    sq_size,    # sq_size 
    Nhabitat  # Number of habitat
    # HabitatSize = NULL # HabitatSize only used if type = "random" instead of "two_sym"
) {   
  
  # Transformation of sq_x and sq_y into vectors
  sq_x <- 1:sq_x
  sq_y <- 1:sq_y
  
  # Creation of the habitat structure
  if (Nhabitat == 2) {
    
    # Tests -----
    if (max(sq_x) * sq_size + sq_size * 2 > W - sq_size) {
      stop("The number of and size of nested habitat expand beyong the width of the grid")
    }
    
    if (max(sq_y) * sq_size + sq_size * 2 > L - sq_size) {
      stop("The number of and size of nested habitat expand beyong the height of the grid")
    }
    # Tests end -----
    
    # Transformation of Nhabitat into a vector 
    Nhabitat <- 1:Nhabitat
    
    # Creation of the matrix of habitat structure filled with the value of Nhabitat[1] = 1 
    hab_mat <- matrix(nrow = L-sq_size, ncol = W-sq_size, Nhabitat[1])
    
    # Modification of the matrix to fill the wanted places with Nhabitat[2] = 2
    for (i in 1:length(sq_x)) {
      for (j in 1:length(sq_y)) {
        hab_mat[(((sq_x[i] * sq_size - sq_size + sq_x[i] * sq_size)+1)-(sq_size/2)):((sq_x[i] * sq_size + sq_x[i] * sq_size)-(sq_size/2)),
                (((sq_y[j] * sq_size - sq_size + sq_y[j] * sq_size)+1)-(sq_size/2)):((sq_y[j] * sq_size + sq_y[j] * sq_size)-(sq_size/2))] <- Nhabitat[2]
      }
    }
  }
  
  # Return the habitat structure
  return(hab_mat)
}
