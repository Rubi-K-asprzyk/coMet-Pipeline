# ----- TESTS ---- # 
PW_Dis <- Obs_PW_Inter[[x]][,-(1:2)]
Sample <- Sampled_Abundance_Data[[x]][[y]][,-c(1:2)]
Sample2 <- Sample[,colnames(Sample) %in% rownames(PW_Dis)] 
samp = Sample ; dis = PW_Dis; abundance.weighted = F

# UNtriangle the matrix.




 
TEST <- mntd(samp = Sample, dis = PW_Dis, abundance.weighted = F) 

# ----- TESTS ---- # 

mntd <- function(samp, dis, abundance.weighted=FALSE)
{
  N <- dim(samp)[1]
  mntd <- numeric(N)
  for (i in 1:N) {
    sppInSample <- names(samp[i,samp[i,] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample,]
      diag(sample.dis) <- NA
      if (abundance.weighted)
      {
        mntds <- apply(sample.dis,2,min,na.rm=TRUE)
        sample.weights <- samp[i,sppInSample]
        mntd[i] <- weighted.mean(mntds, sample.weights)
      }
      else
      {
        mntd[i] <- mean(apply(sample.dis,1,min,na.rm=TRUE))
      }
    }
    else {
      mntd[i] <- NA
    }
  }
  mntd
}
