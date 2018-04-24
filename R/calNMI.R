calNMI <- function(x, y) {
    # Calculates normalized mutual information between two vectors
    #
    # Args: 
    #   x: a vector
    #   y: a vector
    #
    # Returns:
    #   The normalized mutual information between vectors x and y.
	
    x <- as.vector(x)
	y <- as.vector(y)

    mutual.info <- (.mutualInformation(x, y)/
        sqrt(.entropy(x) * .entropy(y)))

    return(max(0, mutual.info, na.rm=TRUE))
}


