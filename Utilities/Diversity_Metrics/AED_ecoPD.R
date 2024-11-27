# internal function to compute distances from one node to all others
# NOTE: assumes node IDs are always 1:n (currently true of phylo4 objs)
dijkstra <- function(phy, node) {
  edge <- edges(phy, drop.root=TRUE)
  elen <- edgeLength(phy)
  nodes <- nodeId(phy, "all")
  n <- length(nodes)
  d <- rep(Inf, length=n)
  names(d) <- nodes
  d[node] <- 0
  while (length(nodes)>0) {
    u <- nodes[which.min(d[nodes])[[1]]]
    if (is.infinite(d[u])) break
    nodes <- nodes[-match(u, nodes)]
    anc <- edge[edge[,2]==u,1]
    for (parent in anc[anc %in% nodes]) {
      alt <- d[u] + elen[paste(parent, u, sep="-")]
      if (alt < d[parent]) {
        d[parent] <- alt
      }
    }
    des <- edge[edge[,1]==u,2]
    for (child in des[des %in% nodes]) {
      alt <- d[u] + elen[paste(u, child, sep="-")]
      if (alt < d[child]) {
        d[child] <- alt
      }
    }
  }
  d
}

# tip length extractor
tipLength <- function(phy, from=c("parent", "root")) {
  from <- match.arg(from)
  tips <- nodeId(phy, type="tip")
  if (from=="parent") {
    len <- edgeLength(phy, tips)
  } else if (from=="root") {
    len <- dijkstra(phy, rootNode(phy))
    len <- len[match(tips, names(len))]
  }
  names(len) <- tipLabels(phy)
  return(len)
}

# abundance extractor
abundance <- function(phy, comm, tip, na.zero=FALSE) {
  communities <- names(phy@metadata$comms)
  if (missing(comm)) {
    comm <- communities
  }
  doNotExist <- !comm %in% communities
  if (any(doNotExist)) {
    stop("one or more communities not found in phy: ",
         paste(comm[doNotExist], collapse=", "))
  }
  if (missing(tip)) {
    tip <- tipLabels(phy)
  } else {
    tip <- getNode(phy, tip, type="tip", missing="warn")
    tip <- names(tip)[!is.na(tip)]
  }
  N <- tipData(phy)[tip, comm, drop=FALSE]
  if (na.zero) N[is.na(N)] <- 0
  return(N)
}

# abundance assignment function
`abundance<-` <- function(phy, comm, tip, value) {
  if (!is.atomic(comm) || length(comm)!=1) {
    stop("when replacing, comm must be a vector of length 1")
  } else if (!comm %in% names(phy@metadata$comms)) {
    stop(paste("community", comm, "not found in phy", sep=" "))
  }
  if (missing(tip)) {
    tip <- tipLabels(phy)
  } else {
    tip <- names(getNode(phy, tip, type="tip", missing="fail"))
  }
  tipData(phy)[tip, comm] <- value
  return(phy)
}

presence <- function(phy, comm, tip, na.zero=FALSE) {
  N <- abundance(phy, comm, tip, na.zero=na.zero)
  N[N > 0] <- 1
  N[N <= 0] <- 0
  N
}

richness <- function(phy, comm, na.zero=FALSE) {
  P <- presence(phy, comm, na.zero=na.zero)
  colSums(P)
}

# minTL extractor
minTL <- function(phy) {
  minTL <- tipData(phy)$minTL
  if (!is.null(minTL)) {
    names(minTL) <- row.names(tipData(phy))
  }
  return(minTL)
}

# Note: function assumes values are in nodeId(phy, "tip") order
# minTL assignment function
`minTL<-` <- function(phy, value) {
  if (!is.numeric(value)) {
    stop("minTL values must be a numeric vector")
  } else if (length(value)!=nTips(phy)) {
    stop("number of minTL values must equal number of species")
  }
  tipData(phy)[["minTL"]] <- value
  return(phy)
}

# community labels extractor
communities <- function(x) {
  names(x@metadata$comms)
}

## this works as implementation of dist.nodes for phylo4 objects, albeit
## about 1.5x slower than dist.nodes
pairdist <- function(phy, type=c("all", "tip"), use.labels=FALSE) {
  
  if (hasPoly(phy) || hasRetic(phy)) {
    stop("pairdist can't currently handle polytomies or reticulation")
  }
  
  type <- match.arg(type)
  edge <- edges(phy, drop.root=TRUE)
  elen <- edgeLength(phy)
  elen <- elen[!names(elen) %in% edgeId(phy, "root")]
  nodes <- nodeId(phy, "all")
  n <- length(nodes)
  d <- matrix(NA_real_, nrow=n, ncol=n)
  diag(d) <- 0
  d[edge] <- d[edge[,2:1]] <- elen
  ntips <- nTips(phy)
  tips <- nodeId(phy, "tip")
  tip.parents <- sapply(tips, function(n) edge[edge[,2]==n, 1])
  nodes.todo <- tip.parents[duplicated(tip.parents)]
  done <- tips
  root <- if (isRooted(phy)) rootNode(phy) else ntips+1
  descendants <- matrix(edge[order(edge[,1]),2], nrow=2)
  ancestors <- edge[match(nodes, edge[,2]), 1]
  
  while(length(nodes.todo)>0) {
    nod <- nodes.todo[1]
    if (nod==root && length(nodes.todo)>1) {
      nod <- nodes.todo[2]
    }
    des1 <- descendants[1, nod-ntips]
    des2 <- descendants[2, nod-ntips]
    if (des1>ntips) des1 <- which(!is.na(d[des1, ]))
    if (des2>ntips) des2 <- which(!is.na(d[des2, ]))
    for (y in des1) {
      d[y, des2] <- d[des2, y] <- d[nod, y] + d[nod, des2]
    }
    if (nod!=root) {
      anc <- ancestors[nod]
      l <- elen[paste(anc, nod, sep="-")]
      d[des2, anc] <- d[anc, des2] <- d[nod, des2] + l
      d[des1, anc] <- d[anc, des1] <- d[nod, des1] + l
      done <- c(done, nod)
      if (all(descendants[, anc-ntips] %in% done)) {
        nodes.todo <- c(nodes.todo, anc)
      }
    }
    nodes.todo <- nodes.todo[nod!=nodes.todo]
  }
  
  if (type=="tip") {
    d <- d[1:ntips, 1:ntips]
  }
  
  if (use.labels) {
    rownames(d) <- colnames(d) <- unname(labels(phy, type=type))
  } else {
    rownames(d) <- colnames(d) <- nodeId(phy, type=type)
  }
  
  return(d)
}



setClass("phylo4com",
         representation(subtrees="list"),
         prototype = list(subtrees = list()),
         validity = checkPhylo4,
         contains="phylo4d")

## create single tree with all data in it, plus a list of subtrees (sans
## data) for each unique species composition
setGeneric("phylo4com", function(x, n, ...) {
  standardGeneric("phylo4com")
})

setMethod("phylo4com", c("phylo", "ANY"),
          function(x, n, ..., check.node.labels="keep") {
            x <- phylo4(x, check.node.labels)
            phylo4com(x, n, ...)
          })

setMethod("phylo4com", c("phylo4", "numeric"),
          function(x, n, communityID, species,
                   missing=c("warn", "OK", "fail")) {
            
            ## create site-by-species abundance matrix
            comm <- factor(communityID, levels=unique(communityID))
            taxa <- factor(species)
            dat <- matrix(0, nrow=nlevels(taxa), ncol=nlevels(comm),
                          dimnames=list(levels(taxa), unique(comm)))
            dat[cbind(taxa, comm)] <- n
            
            ## hand off to the phylo4com matrix method
            phylo4com(x, dat, missing)
            
          })

setMethod("phylo4com", c("phylo4", "matrix"),
          function(x, n, missing=c("warn", "OK", "fail")) {
            
            taxa <- rownames(n)
            if (any(duplicated(taxa))) {
              stop("duplicated taxa are not permitted")
            }
            comm <- colnames(n)
            if (any(duplicated(comm))) {
              stop("duplicated community IDs are not permitted")
            }
            
            phy <- subset(x, tips.include=taxa)
            phyd <- addData(phy, tip.data=n, extra.data="OK")
            phylo4com(phyd, cols=comm)
            
          })

setMethod("phylo4com", c("phylo4", "data.frame"),
          function(x, n, missing=c("warn", "OK", "fail")) {
            n <- as.matrix(n)
            phylo4com(x, n, missing)
          })

setMethod("phylo4com", c("phylo4d", "missing"),
          function(x, n, cols) {
            
            if (missing(cols)) {
              cols <- names(tipData(x))
            }
            if (is.null(cols)) {
              res <- as(x, "phylo4com")
              res@metadata$comms <- NULL
              return(res)
            }
            x@metadata$comms <- setNames(rep(NA, length(cols)),
                                         make.names(cols))
            
            # create trees for each unique community wrt composition only
            P <- presence(x)
            phy <- extractTree(x)
            ids <- as.character(as.numeric(factor(sapply(P, paste,
                                                         collapse=""))))
            subtrees <- lapply(P[!duplicated(ids)],
                               function(n) subset(phy, rownames(P)[n %in% 1]))
            names(subtrees) <- ids[!duplicated(ids)]
            
            res <- as(x, "phylo4com")
            res@subtrees <- subtrees
            res@metadata$comms[] <- ids
            
            return(res)
          })

##
## ED and related methods
##

setGeneric("ed", function(x, na.rm=TRUE) {
  standardGeneric("ed")
})

setMethod("ed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  ed(phyc, na.rm=na.rm)
})

# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
setMethod("ed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  
  comms <- x@metadata$comms
  if (is.null(comms)) {
    stop("no community data specified in phylo4com object")
  }
  subtrees <- x@subtrees[unique(as.character(comms))]
  .edi <- function(tree) {
    # set length of root edge to zero
    edgeLength(tree)[edgeId(tree, "root")] <- 0
    
    all.nodes <- nodeId(tree, type = "all")
    des <- descendants(tree, all.nodes, type="tips")
    nv <- edgeLength(tree, all.nodes) / sapply(des, length)
    names(nv) <- all.nodes
    
    tip.nodes <- nodeId(tree, "tip")
    anc <- ancestors(tree, tip.nodes, "ALL")
    
    res <- sapply(anc, function(n) sum(nv[as.character(n)], na.rm=TRUE))
    names(res) <- tipLabels(tree)
    res
  }
  res <- lapply(subtrees, .edi)[as.character(comms)]
  names(res) <- names(comms)
  return(res)
  
})

setGeneric("hed", function(x, na.rm=TRUE) {
  standardGeneric("hed")
})

setMethod("hed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  hed(phyc, na.rm=na.rm)
})

setMethod("hed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  ED <- ed(x)
  PD <- pd(x)
  res <- sapply(seq_len(length(PD)), function(i) {
    scaledED <- ED[[i]] / PD[[i]]
    -sum(scaledED * log(scaledED))
  })
  names(res) <- names(PD)
  return(res)
})

setGeneric("eed", function(x, na.rm=TRUE) {
  standardGeneric("eed")
})

setMethod("eed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  eed(phyc, na.rm=na.rm)
})

setMethod("eed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  comms <- x@metadata$comms
  if (is.null(comms)) {
    stop("no community data specified in phylo4com object")
  }
  subtrees <- x@subtrees[unique(as.character(comms))]
  hed(x) / log(sapply(subtrees, nTips)[as.character(comms)])
})

setGeneric("aed", function(x, na.rm=TRUE) {
  standardGeneric("aed")
})

setMethod("aed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  aed(phyc, na.rm=na.rm)
})

# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
setMethod("aed", signature(x="phylo4com"),
          function(x, na.rm=TRUE) {
            
            comms <- x@metadata$comms
            if (is.null(comms)) {
              stop("no community data specified in phylo4com object")
            }
            subtrees <- x@subtrees[unique(as.character(comms))]
            .isD <- function(tree) {
              # get all node IDs, but excluding root node
              nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
              # Create logical matrix indicating which tips (in columns) are
              # descendants of each node (in rows), self-inclusive
              t(sapply(ancestors(tree, tipLabels(tree), "ALL"),
                       function(n) nonroot.nodes %in% n))
            }
            
            .elen <- function(tree) {
              nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
              edgeLength(tree, nonroot.nodes)
            }
            isDescendant <- lapply(subtrees, .isD)[as.character(comms)]
            edge.length <- lapply(subtrees, .elen)[as.character(comms)]
            N <- abundance(x)
            
            res <- lapply(seq_along(N), function(i) {
              spp <- row.names(isDescendant[[i]])        
              dAbund <- N[spp, i] * isDescendant[[i]]
              
              # Calculate individual-based AED of each species
              AED <- colSums(edge.length[[i]] * t(prop.table(dAbund, margin=2)))
              AED/N[spp, i]
            })
            names(res) <- names(comms)
            return(res)
            
          })

#    .isD <- function(tree, template) {
#        # get all node IDs, but excluding root node
#        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
#        # Create logical matrix indicating which tips (in columns) are
#        # descendants of each node (in rows), self-inclusive
#        res <- t(sapply(ancestors(tree, tipLabels(tree), "ALL"),
#          function(n) nonroot.nodes %in% n))
#        template[match(rownames(res), rownames(template)), ] <- res
#        template
#    }
#
#    .el <- function(tree) {
#        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
#        edgeLength(tree, nonroot.nodes)
#    }
#    tmp <- setNames(rep(NA, nTips(x)), tipLabels(x))
#    tmp <- matrix(NA, nrow=nTips(x), ncol=nNodes(x)+nTips(x)-1)
#    rownames(tmp) <- tipLabels(x)
#    isDescendant <- lapply(subtrees, .isD, tmp)
#    isDescendant <- array(unlist(isDescendant), dim=c(nTips(x),
#        nTips(x)+nNodes(x)-1, length(subtrees)))
#
#    # Create vector of ancestral edge lengths
#    edge.length <- sapply(subtrees, .elen)
#
#    # Create matrix containing number of individuals of each species
#    # descending from each interior node
#    N <- as.matrix(abundance(x))
#    dAbund <- sweep(isDescendant, c(1,3), N, "*")
#
#    # Calculate individual-based AED of each species
#    pt <- prop.table(dAbund, margin=c(2,3))
#    AED <- sweep(pt, c(2,3), edge.length, "*")
#    AED <- apply(AED, 3, rowSums) / N
#
#    return(AED)

setGeneric("haed", function(x, na.rm=TRUE) {
  standardGeneric("haed")
})

setMethod("haed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  haed(phyc, na.rm=na.rm)
})

setMethod("haed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  # Recast AED in terms of individuals
  AED <- aed(x)
  PD <- pd(x)
  N <- abundance(x)
  scaled.AED <- lapply(seq_along(N), function(i) {
    spp <- names(AED[[i]])        
    rep(unname(AED[[i]]), N[spp, i]) / PD[[i]]
  })
  res <- sapply(scaled.AED, function(x) -sum(x * log(x)))
  names(res) <- names(AED)
  return(res)
})

setGeneric("eaed", function(x, na.rm=TRUE) {
  standardGeneric("eaed")
})

setMethod("eaed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
  phyc <- phylo4com(x)
  eaed(phyc, na.rm=na.rm)
})

setMethod("eaed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  haed(x) / log(colSums(abundance(x)))
})

setGeneric("value",
           function(x, na.rm=TRUE) {
             standardGeneric("value")
           })

setMethod("value", signature(x="phylo4d"),
          function(x, na.rm=TRUE) {
            phyc <- phylo4com(x)
            value(phyc, na.rm=na.rm)
          })

setMethod("value", signature(x="phylo4com"),
          function(x, na.rm=TRUE) {
            AED <- aed(x)
            lapply(AED, function(x) x/sum(x))
          })