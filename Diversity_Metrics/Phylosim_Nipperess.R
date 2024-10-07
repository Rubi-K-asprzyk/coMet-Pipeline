## R function for the calculation of pairwise similarity in phylogenetic diversity across multiple samples.
## by David Nipperess, Macquarie University, Australia (dnippere@bio.mq.edu.au)
## Last edited 3rd February 2009.

##  Copyright (C) 2009  David Nipperess

##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  any later version.

##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License <http://www.gnu.org/licenses/> for more details.

## You must download and install the ape package (available from http://cran.r-project.org/) to use this function.

## ---------------

## Usage: phylosim (x, phy, incidence = T, method = "sorensen")
## Arguments:
##   "x" is a community data table (as in the vegan package) with species/OTUs as columns and samples/sites as rows. Columns are labelled with the names of the species/OTUs. Rows are labelled with the names of the samples/sites. Data can be either abundance or incidence (0/1). 
##   "phy" is a rooted phylogenetic tree with branch lengths stored as a phylo object (as in the ape package) with terminal nodes labelled with names matching those of the community data table. Note that the function trims away any terminal taxa not present in the community data table, so it is not necessary to do this beforehand. 
##   "incidence" is a logical indicating whether the data are to be treated as incidence (binary presence-absence) or abundance. 
##   "method" indicates the particular form of the similarity index you wish to use. Current options are: "sorensen" (default - 2a/a+b+c), "jaccard" (a/a+b+c), "simpson" (a/a+min{b,c}) and "faith" (a+0.5d/a+b+c+d). 
## Details: phylosim takes a community data table and a rooted phylogenetic tree (with branch lengths) and calculates the similarity in Phylogenetic Diversity (PD-similarity) of all pairwise combinations of samples/sites. The principles for calculating PD-similarity on incidence data are discussed by Ferrier et al. (2007). I have extended this approach to include abundance data (manuscript in preparation). 
## Value: phylosim returns a dist object giving the PD-similarity of all pairwise combinations of sample/sites in x. 
## References: Ferrier S, Manion G, Elith J & Richardson K. 2007. Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. Diversity & Distributions 13: 252-264. 

## ---------------

phylosim <- function (x, phy, incidence = T, method = "sorensen") {
  
  METHODS <- c("sorensen", "jaccard", "simpson", "faith")
  method <- pmatch(method, METHODS)
  
  ### step 1: trimming the tree to match the community data table thus creating a "community tree" (sensu Cam Webb).
  
  if (length(phy$tip.label) > length(x[1,])) {
    phy <- drop.tip (phy, which(!phy$tip.label %in% colnames(x))) }
  
  # script is modified from that of Paradis (2006) "Analysis of Phylogenetics and Evolution with R", Springer.
  
  ### step 2: converting a community tree into a MRP matrix
  
  # A MRP matrix, used in supertree methods, is where the membership of an OTU in a clade spanned by a branch is indicated by a 0 (no) or 1 (yes). Unlike supertree MRP matrices, our matrix includes terminal branches.
  # the new phylo object model for ape 2.0 broke the original code. The following replaces that code.
  
  phylomatrix <- matrix (0, length(phy$tip.label), length(phy$edge.length))
  for (i in 1:length(phy$tip.label)) {
    lineage <- which (phy$edge[,2] == i)
    node <- phy$edge[lineage,1]
    while (node > length(phy$tip.label)+1) {
      branch <- which (phy$edge[,2] == node)
      lineage <- c (lineage, branch)
      node <- phy$edge[branch,1]
    }
    phylomatrix[i,lineage] = 1
  }
  
  # this script fills a matrix of 0's with 1's corresponding to the incidence of an OTU on a particular branch.
  # the code is pretty slow on large trees.
  
  ### step 3: re-ordering the OTUs of the occurrence table and MRP matrix to match.
  
  phylomatrix <- phylomatrix[sort.list(phy$tip.label), ]
  x <- x[ ,sort.list(colnames(x))]
  
  # data are sorted to a common ordering standard, that is alphabetic order, so that OTUs match up.
  
  ### step 4: creating a community phylogeny matrix from a MRP matrix and an occurrence matrix
  
  x <- as.matrix (x)
  commphylo <- x %*% phylomatrix
  
  # the above code performs simple matrix multiplication to produce a composite branch by sample matrix (the community phylogeny matrix) where each cell now contains either OTU richness or total abundance for each branch in each sample.
  
  ### step 5 (optional): converting a community phylogeny matrix into an incidence (0/1) form
  
  if (incidence == T) {
    commphylo <- ifelse (commphylo > 0, 1, 0) }
  
  ### step 6: calculating shared diversity (a)
  
  a <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
  for (i in 1:length(x[ ,1])) {
    for (j in 1:length(x[ ,1])) {
      a[i,j] = sum (phy$edge.length * pmin (commphylo[i, ], commphylo[j, ])) }
  }
  
  ### step 7: calculating unshared diversity (b)
  
  b <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
  for (i in 1:length(x[ ,1])) {
    for (j in 1:length(x[ ,1])) {
      b[i,j] = sum (phy$edge.length * (pmax(commphylo[i, ], commphylo[j, ]) - commphylo[j, ])) }
  }
  
  ### step 8: calculating unshared diversity (c)
  
  c <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
  for (i in 1:length(x[ ,1])) {
    for (j in 1:length(x[ ,1])) {
      c[i,j] = sum (phy$edge.length * (pmax(commphylo[i, ], commphylo[j, ]) - commphylo[i, ])) }
  }
  
  ### step 9: calculating absent diversity (d)
  
  d <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
  maxabundance <- vector (mode = "numeric", length = length(commphylo[1, ]))
  for (i in 1:length(commphylo[1, ])) {
    maxabundance[i] = max(commphylo[,i])
  }
  for (i in 1:length(x[ ,1])) {
    for (j in 1:length(x[ ,1])) {
      d[i,j] = sum (phy$edge.length * (maxabundance - pmax(commphylo[i, ], commphylo[j,]))) }
  }
  
  ### step 10: calculating the pairwise similarity Index of samples
  
  if (method == 1) {
    pdsim <- (2 * a)/((2 * a) + b + c)
  }
  else if (method == 2) {
    pdsim <- a/(a + b + c)
  }
  else if (method == 3) {
    pdsim <- a/(a + pmin(b,c))
  }
  else if (method == 4) {
    pdsim <- (a + (0.5 * d))/(a + b + c + d)
  }
  
  rownames (pdsim) <- rownames (x)
  colnames (pdsim) <- rownames (x)
  pdsim <- as.dist (pdsim)
  
  return (pdsim)
  
}