standardNormalization = function(x) {
    # Normalizes each column of x to have mean of 0 and standarad deviation of 1
    #
    # Args:
    #   x: Matrix
    #
    # Returns:
    #   xNorm: Normalized matrix
 
    x <- as.matrix(x)
    mean <- apply(x, 2, mean)
    sd <- apply(x, 2, sd)
    sd[sd==0] <- 1
    xNorm <- t((t(x) - mean) / sd)

    return(xNorm)
}
