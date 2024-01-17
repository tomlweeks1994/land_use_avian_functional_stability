TPDRich <- function (TPDc = NULL, TPDs = NULL) 
{
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm") {
    stop("The class of one object do not match the expectations,\n         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp") {
    stop("The class of one object do not match the expectations,\n         Please, specify if your object is a TPDc or a TPDs")
  }
  results <- list()
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  if (!is.null(TPDc)) {
    results$communities <- list()
    message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)

  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" | 
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    }
    else {
      results$populations <- list()
      message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness\n              and Divergence are meaningless!!")
    }
  }
  return(results)
}


TPDRed <- function (x = NULL) 
  {
  if (class(x) != "TPDcomm") {
    stop("TPDc must be an object of class TPDcomm generated with the TPDc\n\t\t    function")
  }
  results <- list()
  results$redundancy <- results$richness <- numeric()
  for (i in 1:length(x$TPDc$TPDc)) {
    
    TPDc_aux <- x$TPDc$TPDc[[i]]
    RTPDs_aux <- x$TPDc$RTPDs[[i]]
    RTPDs_aux[RTPDs_aux > 0] <- 1
    M <- rowSums(RTPDs_aux)
    results$redundancy[i] <- sum(M * TPDc_aux) - 1
    
  } 
  names(results$redundancy) <- names(x$TPDc$TPDc)
  return(results)
}

