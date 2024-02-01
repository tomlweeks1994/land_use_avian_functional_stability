rao_diversity_gaw <- function (comm, forage_true = c("yes", "no"), traits = NULL, weights, dist_matrix = NULL, phylodist = NULL, checkdata = TRUE,
                               ord = "metric", put.together = NULL, standardize = TRUE,
                               ...) 
{
  diver.internal <- function(community, distance) {
    if (any(is.na(distance))) {
      distance.na <- ifelse(is.na(distance), 0, 1)
      inter.na <- community %*% distance.na
      adjustment <- rowSums(sweep(community, 1, inter.na, 
                                  "*", check.margin = FALSE))
      distance[is.na(distance)] <- 0
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
      res <- ifelse(adjustment > 0, res/adjustment, res)
    }
    else {
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
    }
    return(res)
  }
  res <- list(call = match.call())
  if (inherits(comm, "metacommunity.data")) {
    if (!is.null(traits) | !is.null(phylodist) | !is.null(put.together)) {
      stop("\n When you use an object of class metacommunity.data the arguments traits, phylodist and put.together must be null. \n")
    }
    traits <- comm$traits
    phylodist <- comm$phylodist
    put.together <- comm$put.together
    comm <- comm$community
  }
  list.warning <- list()
  if (checkdata) {
    organize.temp <- organize.syncsa(comm, traits = traits, 
                                     phylodist = phylodist, check.comm = TRUE)
    if (!is.null(organize.temp$stop)) {
      organize.temp$call <- match.call()
      return(organize.temp)
    }
    list.warning <- organize.temp$list.warning
    comm <- organize.temp$community
    traits <- organize.temp$traits
    phylodist <- organize.temp$phylodist
  }
  if (length(list.warning) > 0) {
    res$list.warning <- list.warning
  }
  if (any(is.na(comm))) {
    stop("\n community data with NA\n")
  }
  comm <- as.matrix(comm)
  N <- nrow(comm)
  S <- ncol(comm)
  ## we have a 1s and 0s distace matrix
  dist.1 <- 1 - diag(x = rep(1, S))
  if (!is.null(traits)) {
    traits <- as.data.frame(traits)
    m <- ncol(traits)
    weights <- rep(1, m)
    make.names <- is.null(colnames(traits))
    colnames(traits) <- colnames(traits, do.NULL = FALSE, 
                                 prefix = "T")
    names(weights) <- colnames(traits)
    if (!is.null(put.together)) {
      if (!inherits(put.together, "list")) {
        stop("\n put.together must be a object of class list\n")
      }
      if (make.names) {
        for (k in 1:length(put.together)) {
          put.together[[k]] <- paste("T", put.together[[k]], 
                                     sep = "")
        }
      }
      if (max(table(unlist(put.together))) > 1) {
        stop("\n The same trait appears more than once in put.together\n")
      }
      if (length(setdiff(unlist(put.together), colnames(traits))) > 
          0) {
        stop("\n Check traits names in put.together\n")
      }
      for (k in 1:length(put.together)) {
        weights[put.together[[k]]] <- 1/length(put.together[[k]])
      }
    }
    
    if(forage_true == "yes") {
  
    dist.functional <- dist_matrix
      
    } else {
    dist.functional <- as.matrix(cailliez(gawdis::gawdis(x = traits, w.type = "optimized", groups = c(1,2,3), opti.maxiter = 1000)))
    }
    
    if (checkdata) {
      if (any(is.na(dist.functional))) {
        warning("Warning: NA in distance between species", 
                call. = FALSE)
      }
    }
  }
  if (!is.null(phylodist)) {
    dist.phylogenetic <- as.matrix(phylodist)
    if (checkdata) {
      if (any(is.na(dist.phylogenetic))) {
        warning("Warning: NA in phylodist", call. = FALSE)
      }
    }
    if (standardize) {
      dist.phylogenetic <- dist.phylogenetic/max(dist.phylogenetic, 
                                                 na.rm = TRUE)
    }
  }
  comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/") ## relative abundance matrix
  SD <- diver.internal(comm, dist.1)
  res$Simpson <- SD
  if (!is.null(traits)) {
    FD <- diver.internal(comm, dist.functional)
    res$FunRao <- FD
    res$FunRedundancy <- SD - FD
  }
  if (!is.null(phylodist)) {
    PD <- diver.internal(comm, dist.phylogenetic)
    res$PhyRao <- PD
    res$PhyRedundancy <- SD - PD
  }
  return(res)
}
