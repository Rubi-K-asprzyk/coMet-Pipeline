### ----------------TEST------------------- ###

### ----------------TEST------------------- ###

mpd_V2 <- function(Sample, PW_distances, abundance.weighted=FALSE , intra = F)
{
  N <- nrow(Sample)  # Find the size of the sample (the number of communities)
  mpd <- numeric(N)  # Initialize the MPD vector
  
  if (identical(rownames(PW_distances),colnames(PW_distances)) == FALSE) {stop("ERROR: Colnames and rownames of the pairwise distance matrix are not identical.")} # Be sure that the row names are the same as the col names
  
  for (i in 1:N) {
    # Look for the species present in the samples
    sppInSample <- names(Sample[i,which(Sample[i,] > 0)])
    # If they are species 
    if (length(sppInSample) > 1) {
      # Extract the values of the pairwise comparisons for the species present in the community
      sample.dis.intra <- as.matrix(PW_distances[sppInSample, sppInSample])
      # Create the inter versions (without the diagonals)
      # sample.dis.inter <- as.matrix(sample.dis.intra)
      # diag(sample.dis.inter) <- NA
      
      # Abundance.weighted version
      if (abundance.weighted) {
        stop("WARNING: Abundance weighted version not implemented yet.")
        # Create the samples weights
        sample.weights.inter <- t(as.data.frame(sample.dis.inter[i,sppInSample,drop=FALSE])) %*% as.matrix(sample.dis.inter[i,sppInSample,drop=FALSE]) %>% as.data.frame()
        # test <- as.matrix(sample.dis.inter,na.omit = T) %*% as.matrix(sample.dis.inter)
        sample.weights.intra <- t(as.matrix(sample.dis.intra[i,sppInSample,drop=FALSE])) %*% as.matrix(sample.dis.intra[i,sppInSample,drop=FALSE]) %>% as.data.frame()
        # COmpute the weighted mean. 
        if (intra == F) {mpd[i] <- stats::weighted.mean(as.vector(sample.dis.inter),as.vector(sample.weights.inter),na.rm = TRUE)}
        if (intra == T) {mpd[i] <- weighted.mean(sample.dis.intra,sample.weights.intra)}
      }
      else {
        # OLD
        # Initial portion of the script
        #  mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
        
        # NEW
        # Just create a parameter to keep or not the diagonal. 
        
        if (intra == F) {mpd[i] <- mean(sample.dis.intra[lower.tri(sample.dis.intra)])} # We do not want the diagonal. 
        if (intra == T) {mpd[i] <- mean(sample.dis.intra[lower.tri(sample.dis.intra,diag = T)])} # We want the diagonal
        
      }
    }
    else{
      mpd[i] <- NA
    }
  }
  names(mpd) <- rownames(Sample)
  
  return(mpd)
}
