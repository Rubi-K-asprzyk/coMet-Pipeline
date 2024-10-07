
##### ------------- TEST ----------------- #####
# x <- 1
# y <-"Hab1Core"  #  "Hab2Core"    "Hab1Ecotone" "Hab2Ecotone"
#weighted.faith( my.phylo = List_Phylo[[x]] , my.sample = Sampled_Abundance_Data[[x]][[y]][,-(1:2)], type = "V" )
# my.sub.sample <- Sampled_Abundance_Data[[x]][[y]][1,-(1:2)]
##### ------------- TEST ----------------- #####

weighted.faith <-
  function(my.phylo, my.sample, type){
    
    weighted.faith.function = function(my.sub.sample){

      ## Modify the sub.samples to only keep present species
      my.sub.sample <- my.sub.sample[which(my.sub.sample > 0)]
      
      ## Add a condition where we need at least two tips to draw a tree
      if (length(my.sub.sample) < 2){
        # Mectric can't be computed, therefore <- NA
        out <- NA
        return(out)
      }
      
      ## extract the names of species in a community with an abundance greater than zero and use that information to make a pruned phylogeny for that community.
      tmp.tree = treedata(my.phylo,my.sub.sample,warnings=F)$phy
      
      ## Create empty branches matrix
      branches = matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
      
      ## Fill first two columns of the matrix with node numbers defining each edge.
      branches[,1:2] = tmp.tree$edge
      
      ## Fill the third column with the length of each branch
      branches[,3] = tmp.tree$edge.length
      
      get.leaves<-function(x){
        leaves.node<-tips(tmp.tree,x[2])
      }
      
      ## Apply the get.leaves() function to each row in the branches object. This will retrieve species names subtended by each branch (i.e. ## row) in the branches matrix
      leaves = apply(branches, MARGIN = 1, get.leaves)
      
      ## Now quickly loop through each set of leaves to ## calculate the mean abundance (Ai) for those species.
      for(i in 1:length(leaves)){
        branches[i,4] = mean(my.sub.sample[leaves[[i]]], na.rm = T) 
        
      }
      
      ## Lastly calculated the Weighted Faith’s Index sensu VELEND
      if(type == "V") {
        
        out <- nrow(tmp.tree$edge) * ((sum(branches[,3] * branches[,4])) / sum(branches[,4]))}
      
      # ------- ADDED PART TO COMPUTE Weighted Faith’s Index sensu BARKER -------- # 
      ## Lastly calculated the Weighted Faith’s Index sensu BARKER
      if(type == "B") {

        out <- sum(branches[,3] * branches[,4])}
      
      # return the last value
      return(out)
    }
    
    if(type == "V") { Name <- "PDab" } else {Name <- "DnPD"}

    outt = apply(my.sample, MARGIN = 1, weighted.faith.function)
    outt <- as.data.frame(outt)
    colnames(outt) <- Name
    return(outt)
    
  }
