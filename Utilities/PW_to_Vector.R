##### ----- TEST ----- #####
# PW_matrix <- Metric_List[[1]] ; Type <- "Intra" # Intra
# PW_matrix <- Metric_List[[6]] ; Type <- "Inter" # Inter
# 
# Test <- PW_to_Vector(Intra, Colname = "Test")

##### ----- TEST ----- ##### 

PW_to_Vector <- function(PW_matrix, Colname, Diag = T, Type = "Intra"){
  
  if (Type == "Intra"){  
    
    if (Diag == T) {Compound_new <- upper.tri(PW_matrix , diag = TRUE)}
    if (Diag == F) {Compound_new <- upper.tri(PW_matrix , diag = FALSE)}
    i2 <- which(Compound_new, arr.ind=TRUE)
    Compound_final <- data.frame(Sample = paste0(rownames(PW_matrix)[i2[,1]],"-",colnames(PW_matrix)[i2[,2]]),
                                 Sample_A = rownames(PW_matrix)[i2[,1]], 
                                 Sample_B = colnames(PW_matrix)[i2[,2]], 
                                 Colname = PW_matrix[Compound_new]) %>%
      setNames(nm = c("Sample","Sample_A","Sample_B",Colname))
    
    return(Compound_final) } # End of type "Intra"
  
  if (Type == "Inter"){
    Compound_new <- matrix(data = TRUE, nrow = nrow(PW_matrix), ncol = ncol(PW_matrix))
    i2 <- which(Compound_new, arr.ind=TRUE)
    Compound_final <- data.frame(Sample = paste0(rownames(PW_matrix)[i2[,1]],"-",colnames(PW_matrix)[i2[,2]]),
                                 Sample_A = rownames(PW_matrix)[i2[,1]], 
                                 Sample_B = colnames(PW_matrix)[i2[,2]], 
                                 Name = PW_matrix[Compound_new]) %>%
      setNames(nm = c("Sample","Sample_A","Sample_B",Colname))
    
    return(Compound_final)  } # End of type "Inter"
  
}
