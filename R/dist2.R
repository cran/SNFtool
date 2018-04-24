dist2 <- function(X,C) {
    # Calculates the squared euclidean distance between two matrices where rows 
    #   represent a single data point or patient.
    #
    # Args:
    #   X: Matrix with each row representing a single data point (or patient)
    #   C: Matrix with each row representing a single data point (or patient)
    #
    # Returns:
    #   res: A NxM matrix where nrow(X) == N and nrow(C) == M. Element [n,m] 
    #       is the squared euclidean distance between rows N[n,] and C[m,].

    ndata <- nrow(X)
    ncentres <- nrow(C)
    
    sumsqX <- rowSums(X^2)
    sumsqC <- rowSums(C^2)
      
    XC <- 2 * (X %*% t(C))
    res <- matrix(rep(sumsqX, times=ncentres), ndata, ncentres) + 
        t(matrix(rep(sumsqC, times=ndata), ncentres, ndata)) - XC
    res[res < 0] <- 0

    return(res)
}
