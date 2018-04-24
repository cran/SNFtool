concordanceNetworkNMI <- function(Wall,C) {
    # Calculates all pairwise NMIs between matrices in Wall and the fusion of 
    #   these matrices.
    #
    # Args:
    #   Wall: List of affinity matrices
    #   C: Number of clusters
    # 
    # Returns:
    #   A nxn matrix containing NMIs between cluster assignments made by spectral
    #        clustering for all n matrices in Wall and the fusion of those matrices.

    # Get the cluster labels for each of the networks
    labels <- lapply(Wall, function(x) spectralClustering(x, C))

    # Calculate the NMI between each pair clusters
    LW <- length(Wall)
    NMIs <- matrix(NA, LW, LW)
    for (i in 1:LW) {
        for (j in 1:LW) {
            NMIs[i, j] <- calNMI(labels[[i]], labels[[j]])
        }
    }
  
  return(NMIs)
}
