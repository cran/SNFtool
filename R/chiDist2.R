chiDist2 <- function(A){
    # Calculates the pairwise chi-square distance between all rows in a given matrix.
    # Uses chi2Dist from 'ExPosition' package
    #
    # Args:
    #   A: Matrix with rows representing samples
    # 
    # Returns:
    #   D: NxN matrix where N is the number of rows in A. Element i,j in 
    #       the returned matrix is the chi-square distance between A[i,]
    #        and A[j,].

    return(chi2Dist(A)$D)
}
