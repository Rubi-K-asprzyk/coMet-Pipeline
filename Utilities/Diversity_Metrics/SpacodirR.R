#### ___________ TEST _______________ ####
# OK <- t(Abund_1[[3]][,-1])
# NOTOK <- t(Abund_1[[2]][,-1])
# phy <- List_Phylo[[2]]
# 
# check.distmat(OK)
# 
# THE ERROR WAS IN CHECK.DISTMAT THAT SOMETIMES THINK THAT MY COMMUNITY MATRIX WAS A DISTANCE MATRIX
# Metric <- spacodi.calc(
#   sp.plot = t(Sampled_Abund_Data_Inter[[x]][,-c(1:2)]),   # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
#   phy = List_Phylo[[x]],            # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
#   all.together = TRUE,              # whether to treat all traits together or separately
#   prune = TRUE,
#   pairwise = TRUE)          

#### ------------------------------- ####

# Loading of the internal functions
match.spacodi.data<-function(sp.plot, phy=phy, sp.traits=NULL, prune=TRUE, verbose=TRUE) {
  # major error checking
  sp.plot=as.matrix(sp.plot)
  if(is.null(row.names(sp.plot))) stop("Check that sp.plot has row names as species.") 
  if(is.null(colnames(sp.plot))) {
    warning("sp.plot does not appear to have plots as column names.")
    names(sp.plot)=paste("plot",1:ncol(sp.plot),sep="")
  }
  
  if(!missing(phy)) {
    if(class(phy)=="phylo") { 
      if(length(unique(phy$tip.label))!=length(phy$tip.label)) stop("Redundant taxa were found in tree.")
      if(!is.null(phy$node.label)) phy$node.label=NULL
    }
  }
  
  if(length(unique(names(sp.plot)))!=length(names(sp.plot))) stop("Redundant plots were found in sp.plot.")
  
  # check poor values in sp.plot
  if(any(!is.finite(sp.plot))) {
    poor=which(!is.finite(sp.plot))
    poor%%nrow(sp.plot)->poor.rows
    ceiling(poor/nrow(sp.plot))->poor.cols
    poor.data=as.list(paste(rownames(sp.plot)[poor.rows],colnames(sp.plot)[poor.cols]))
    stop("Poor data values found in supplied sp.plot (NA, NaN, or Inf):\n\n\t",toString(poor.data),"\n")
  }
  
  # find undersampled plots
  prune.sp=function(sp.plot, verbose=FALSE){
    rr=rownames(sp.plot)
    l.spp=nrow(sp.plot)
    drop.plots=vector()
    for(sp in 1:ncol(sp.plot)) {
      l.nulls=length(which(sp.plot[,sp]==0))
      if((l.spp-l.nulls)<2) {
        drop.plots=cbind(drop.plots, sp)
      }
    }
    
    plot.names.orig=colnames(sp.plot)
    dropped.plots=plot.names.orig[drop.plots]
    if(length(drop.plots)!=0) {
      if(verbose)message({cat("\nThe following plots were dropped from sp.plot:\n\t");cat(dropped.plots, sep=" "); cat("\n")})
      sp.plot=as.matrix(sp.plot[,-as.numeric(drop.plots[!is.na(drop.plots)])])
      rownames(sp.plot)=rr
      colnames(sp.plot)=plot.names.orig[which(!plot.names.orig%in%plot.names.orig[drop.plots])]
      if(is.null(ncol(sp.plot))) sp.plot=NULL
    }
    return(sp.plot)
  }
  
  if(prune) sp.plot=prune.sp(sp.plot, verbose=verbose)
  
  # match ordering and size of all data objects	
  usable.species=rownames(sp.plot)
  
  # phy only
  if(missing(sp.traits) && !missing(phy)) {   
    if(class(phy)=="phylo") {
      usable.species.tmp=intersect(phy$tip.label, usable.species)
      phy=reorderspacodiobj(phy,usable.species.tmp)             
      usable.species=usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
      if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
      r.out=list(sp.plot=sp.plot,sp.tree=phy)
    } else if(check.distmat(phy)) {
      usable.species.tmp=intersect(rownames(phy), usable.species)
      phy=reorderspacodiobj(phy,usable.species.tmp)
      usable.species=usable.species.tmp[match(rownames(phy), usable.species.tmp)]
      if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
      r.out=list(sp.plot=sp.plot,sp.tree=phy)			
    }
  }
  
  # sp.traits only
  if(!missing(sp.traits) & missing(phy)) {
    usable.species=intersect(usable.species, rownames(sp.traits))
    if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
    r.out=list(sp.plot=sp.plot,sp.traits=reorderspacodiobj(sp.traits, usable.species))			
  }
  
  # both sp.traits and phy
  if(!missing(sp.traits) & !missing(phy)) {
    if(class(phy)=="phylo") {			
      usable.species.tmp=intersect(phy$tip.label, intersect(rownames(sp.traits), usable.species))
      phy=reorderspacodiobj(phy,usable.species.tmp)
      usable.species=usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
      if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
      r.out=list(sp.plot=sp.plot,sp.tree=reorderspacodiobj(phy, usable.species), sp.traits=reorderspacodiobj(sp.traits, usable.species))
    } else if(check.distmat(phy)) {
      usable.species.tmp=intersect(rownames(phy), intersect(rownames(sp.traits), usable.species))
      phy=reorderspacodiobj(phy,usable.species.tmp)
      usable.species=usable.species.tmp[match(rownames(phy), usable.species.tmp)]
      if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
      r.out=list(sp.plot=sp.plot,sp.tree=reorderspacodiobj(phy, usable.species), sp.traits=reorderspacodiobj(sp.traits, usable.species))			
    }
  }
  
  # neither sp.traits nor phy
  if(missing(sp.traits) & missing(phy)) {
    if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
    r.out=list(sp.plot=sp.plot)
  }
  
  return(r.out)
}	

reorderspacodiobj = function(obj, names) {
  if (any(grepl("phylo", class(obj)))) {
    obj.labels = obj$tip.label
    if (any(!obj.labels %in% names))
      obj <- ape::drop.tip(obj, obj.labels[!obj.labels %in% names])
    else
      obj <- obj
  } else if (check.distmat(obj)) {
    obj.labels <- rownames(obj)
    obj <- as.matrix(obj[match(names, obj.labels), match(names, obj.labels)])
    rownames(obj) <- colnames(obj) <- names
  } else {
    nn <- colnames(obj)
    obj.labels <- rownames(obj)
    obj <- as.data.frame(obj[match(names, obj.labels), ])
    names(obj) <- nn
  }
  if (!all(names %in% obj.labels))
    warning(paste(
      paste(names[!names %in% obj.labels], sep = " ", collapse = " "),
      "were not found in the supplied object",
      sep = " "
    ))
  return(obj)
}



##### ERROR HERE #####
  # The function check.distmat return a True when the obj is a community matrix, not a distance matrix #
  # Therefore, it was replaced by another version because we don't need this internal function. 
# check.distmat <- function(obj, tol = 1e-9) {
#   if (any(grepl("matrix", class(obj)))) {
#     if (all(diag(obj) - 0 <= tol)) {
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   } else {
#     return(FALSE)
#   }
# }


check.distmat <- function(obj) {return(FALSE)}
cophen <- function(phy){
  n=Ntip(phy)
  out <- .Call("cophen", tree=list(
    ROOT = as.integer(n+1),
    MAXNODE = as.integer(max(phy$edge[,1])),
    ENDOFCLADE = as.integer(dim(phy$edge)[1]),
    ANC = as.integer(phy$edge[,1]),
    DES = as.integer(phy$edge[,2]),
    EDGES = as.double(c(phy$edge.length,0)),
    COPHEN = as.double(array(matrix(0, n, n)))),
    PACKAGE = "spacodiR")
  cc=matrix(out$COPHEN,nrow=n,byrow=FALSE)
  rownames(cc)<-colnames(cc)<-phy$tip.label
  return(cc)
}

# Definition of the function
spacodi.calc <-function(sp.plot, phy = NULL, sp.traits = NULL, all.together=TRUE, prune=TRUE, pairwise=TRUE,...){
  if(!missing(phy) && !missing(sp.traits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
  
  rematrix=function(arr, names){
    n=sqrt(length(arr))
    counter=1
    M=matrix(0,n,n)
    for(i in 1:n){
      j=i+1
      while(j<=n) {
        M[i,j]=arr[counter]
        counter=counter+1
        j=j+1
      }
    }
    M[lower.tri(M)]=t(M)[lower.tri(t(M))]
    dimnames(M)[[1]]<-dimnames(M)[[2]]<-names
    return(M)
  }
  
  # determine type of abundance
  stripped=unname(unlist(c(sp.plot)))
  if(all(!is.na(match(stripped, c(0,1))))) {
    abundt = 0		# presence|absence
  } else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
    abundt = 1		# relative abundance
  } else {
    abundt = 2		# abundance (n.individuals)
  }
  
  # iter is 1 unless more than a single trait is being used for separate analyses
  iter=1
  traits.tmp<-distmat<-NULL
  
  # INTERPRET whether to compute trait diversity or phylogenetic diversity & prepare data
  if(!missing(phy)) {
    # distmat is phylogenetic: Pst
    sp.data=match.spacodi.data(sp.plot=sp.plot, phy=phy, prune=prune)
    sp.plot=sp.data$sp.plot
    phy=sp.data$sp.tree
    if(check.distmat(phy)) distmat <- phy else distmat <- cophen(phy)
    
  } else if(missing(sp.traits) & missing(phy)) {
    # distmat is null: Ist
    distmat <- matrix(1, ncol = nrow(sp.plot), nrow = nrow(sp.plot))
    
  } else if(!missing(sp.traits)){
    # distmat is trait-based: Tst
    if(class(sp.traits)=="dist") sp.traits=as.matrix(sp.traits)
    if(ncol(sp.traits)==1) all.together=TRUE
    if(all(is.null(names(sp.traits)))) names(sp.traits)=paste("trt",seq(1:ncol(sp.traits)),sep="")	
    if(all.together==TRUE | check.distmat(sp.traits)) {
      sp.data=match.spacodi.data(sp.plot=sp.plot, sp.traits=sp.traits, prune=prune, ...)
      sp.plot=sp.data$sp.plot
      sp.traits=sp.data$sp.traits
      if(check.distmat(sp.traits)) distmat=sp.traits else distmat=as.matrix(dist(sp.traits))
    } else {	
      iter=ncol(sp.traits)
      traits.tmp=lapply(1:ncol(sp.traits), function(x) {
        trait.tt<-data.frame(sp.traits[,x])
        row.names(trait.tt)<-row.names(sp.traits)
        sp.data<-match.spacodi.data(sp.plot=sp.plot, sp.traits=trait.tt, prune=prune, ...)
        distmat <- as.matrix(dist(sp.data$sp.traits))
        return(list(distmat=distmat, sp.plot=sp.data$sp.plot))
      }
      )
    }
  } else if(is.null(distmat)){ 
    stop("Cannot decipher input object(s).")
  }
  
  if(is.null(names(sp.plot))) pnames=paste("plt",seq(1,ncol(sp.plot)),sep=".") else pnames=names(sp.plot)
  
  
  # PREPARE output
  gen.out=list()
  prw.out=list()
  
  for(tt in 1:iter) {	
    if(is.null(traits.tmp)) {
      sp.plot <-	as.matrix(sp.plot)
    } else {
      sp.plot <-	as.matrix(traits.tmp[[tt]]$sp.plot)
      distmat <-	as.matrix(traits.tmp[[tt]]$distmat)
    }
    
    diag(distmat) <- 0
    np <- ncol(sp.plot)
    out <- NA
    out <- .C("spacodi", 
              np = as.integer(np),
              ns = as.integer(nrow(sp.plot)),
              sp.plot = as.double(as.vector(sp.plot)),
              distmat = as.double(as.vector(as.matrix(distmat))),
              abundtype = as.integer(abundt),
              Ndclass = as.integer(0),
              dclass = as.double(c(0,0)),
              Ist = 0,
              Pst = 0,
              Bst = 0,
              PIst = 0,
              pairwiseIst = as.double(as.vector(matrix(0,np,np))),
              pairwisePst = as.double(as.vector(matrix(0,np,np))),
              pairwiseBst = as.double(as.vector(matrix(0,np,np))),
              pairwisePIst = as.double(as.vector(matrix(0,np,np))),
              PACKAGE="spacodiR"
    )
    
    # compile results for phylogenetic turnover
    if(missing(phy) & missing(sp.traits)) {
      r.out=as.numeric(out[8])
      names(r.out)="Ist"
      gen.out[[tt]]=as.data.frame(t(r.out))
      if(pairwise) {
        prw.out[[tt]]=list(rematrix(unlist(out[12]),pnames))
        names(prw.out[[tt]])<-paste("pairwise","Ist",sep=".")
      }
    } else if(!missing(phy)){
      r.out=as.numeric(c(out[8:11]))
      names(r.out)<-c("Ist","Pst","Bst","PIst")
      gen.out[[tt]]=as.data.frame(t(r.out))
      if(pairwise) {
        prw.out[[tt]]=lapply(out[12:15], function(x) rematrix(x,pnames))
        names(prw.out[[tt]])<-paste("pairwise",c("Ist","Pst","Bst","PIst"),sep=".")
      }
    } else if(!missing(sp.traits)) {
      r.out=as.numeric(c(out[8:11]))
      names(r.out)=c("Ist","Tst","Ust","TAUst")
      gen.out[[tt]]=as.data.frame(t(r.out))
      if(pairwise) {
        prw.out[[tt]]=lapply(out[12:15], function(x) rematrix(x,pnames))
        names(prw.out[[tt]])<-paste("pairwise",c("Ist","Tst","Ust","TAUst"),sep=".")
      }
    } 
  } 
  if(iter>1 & !all.together) {
    RES=lapply(1:iter, function(x) {
      if(pairwise) return(c(gen.out[[x]], prw.out[[x]])) else return(c(gen.out[[x]]))
    }
    )
    names(RES)<-names(sp.traits)
  } else {
    prw.out=unlist(prw.out,recursive=FALSE)
    gen.out=unlist(gen.out,recursive=FALSE)
    if(pairwise) RES=c(gen.out,prw.out) else RES=gen.out
  }
  
  return(unlist(list(RES,sp.plot=list(sp.plot),sp.tree=list(phy),sp.traits=list(sp.traits)),recursive=FALSE)) 
  
} 



