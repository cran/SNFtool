estimateNumberOfClustersGivenGraph <- function(W, NUMC=2:5) {
    # Estimates the best number of clusters from a vector of choices, using 
    #   the eigen-gap & rotation cost heuristics.
    #
    # Args:
    #   W: Affinity matrix (usually result from SNF)
    #   NUMC: A vector of integers specifying which cluster numbers to check
    #
    # Returns:
    #   A vector of the top two suggested number of clusters using
    #       the eigen-gap and rotation cost heuristics. 
    #

    #Put this check after the length(NUMC) check?
    if (min(NUMC) == 1) {
        warning('Note that we always assume there are more than one cluster.')
        NUMC <- NUMC[NUMC > 1]
    }

    #Why is this performed here?
    W <- (W + t(W))/2
    diag(W) <- 0

    #NUMC validity check
    if (length(NUMC) <= 0) {
        warning(paste("Invalid NUMC provided, must be an integer vector",
             "with atleast one other number than 1.",
              "Using default NUMC=c(2,3,4,5)",sep=""))
        NUMC <- 2:5 
    }

    # compute unnormalized Laplacian
    degs <- rowSums(W)
    degs[degs == 0] <- .Machine$double.eps    
    D <- diag(degs)    
    L <- D - W
    Di <- diag(1 / sqrt(degs))
    L <- Di %*% L %*% Di
    #print(dim(L))

    # compute the eigenvectors corresponding to the k smallest
    eigs <- eigen(L)
    eigs_order <- sort(eigs$values, index.return=T)$ix
    eigs$values <- eigs$values[eigs_order]
    eigs$vectors <- eigs$vectors[, eigs_order]
    eigengap <- abs(diff(eigs$values))
#    eigengap <- eigengap * (1 - eigs$values[1:length(eigs$values) - 1]
#        ) / (1 - eigs$values[2:length(eigs$values)])

    quality <- list()
    for (c_index in 1:length(NUMC)) {
        ck <- NUMC[c_index]
        UU <- eigs$vectors[, 1:ck]
        EigenvectorsDiscrete <- .discretisation(UU)[[1]]
        EigenVectors <- EigenvectorsDiscrete^2
      
        #MATLAB: sort(EigenVectors,2, 'descend');
        temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors),
             function(i) EigenVectors[, i])), ]
        temp1 <- t(apply(temp1, 1, sort, TRUE))  
  
        quality[[c_index]] <- (1 - eigs$values[ck + 1]) / 
            (1 - eigs$values[ck]) * 
            sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*%
            temp1[, 1:max(2, ck-1)] ))
    }
    #Eigen-gap best two clusters
    t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)$ix
    K1 <- NUMC[t1[1]]
    K12 <- NUMC[t1[2]]

    #Rotation cost best two clusters
    t2 <- sort(unlist(quality), index.return=TRUE)$ix
    K2 <- NUMC[t2[1]]
    K22 <- NUMC[t2[2]]    
  
    output <- list("Eigen-gap best"=K1, "Eigen-gap 2nd best"=K12,
        "Rotation cost best"=K2, "Rotation cost 2nd best"=K22)
    return (output)
}
