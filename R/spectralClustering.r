spectralClustering <- function(affinity, K, type=3) {
    # Implements spectral clustering on given affinity matrix into K clusters.
    #
    # Args:
    #   affinity: Affinity matrix (size NxN) to perform clustering on
    #   K: Number of clusters 
    #   type (default 3): Used to speciy the type of spectral clustering
    #
    # Returns:
    #   labels: A vector of length N assigning a label 1:K to each sample

    d <- rowSums(affinity)
    d[d == 0] <- .Machine$double.eps
    D <- diag(d)
    L <- D - affinity

    if (type == 1) {
        NL <- L

    } else if (type == 2) {
        Di <- diag(1 / d)
        NL <- Di %*% L

    } else if(type == 3) {
        Di <- diag(1 / sqrt(d))
        NL <- Di %*% L %*% Di
    }

    eig <- eigen(NL)
    res <- sort(abs(eig$values),index.return = TRUE)
    U <- eig$vectors[,res$ix[1:K]]
    normalize <- function(x) x / sqrt(sum(x^2))

    if (type == 3) {
        U <- t(apply(U,1,normalize))
    }

    eigDiscrete <- .discretisation(U)
    eigDiscrete <- eigDiscrete$discrete
    labels <- apply(eigDiscrete,1,which.max)
  
  
 
  return(labels)
}
