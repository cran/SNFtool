rankFeaturesByNMI <- function(data, W){  
    # Calculates the normalized mutual information (NMI) score between each 
    # features clustering and the clustering of the fused matrix W. Each feature
    # is ranked based on how similar it is to the clustering of the fused matrix.
    #
    # Args:
    #   data: A list of matrices 
    #   W: Fused affinity matrix from all data types in data
    #
    # Returns:
    #   A list containing NMI score for each feature from all data types
    #   and their NMI score ranks.
    #

    stopifnot(class(data) == "list")
    
    NUM.OF.DATA.TYPES <- length(data)
    NMI.scores <- vector(mode="list", length=NUM.OF.DATA.TYPES)
    NMI.ranks <- vector(mode="list", length=NUM.OF.DATA.TYPES)
    num.of.clusters.fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
    clustering.fused <- spectralClustering(W, num.of.clusters.fused)
    
    for (data.type.ind in 1:NUM.OF.DATA.TYPES){
        NUM.OF.FEATURES <- dim(data[[data.type.ind]])[2] 
        NMI.scores[[data.type.ind]] <- vector(mode="numeric", 
            length=NUM.OF.FEATURES)
      
        for (feature.ind in 1:NUM.OF.FEATURES){
            affinity.matrix <- affinityMatrix(
                dist2(as.matrix(data[[data.type.ind]][, feature.ind]), 
                as.matrix(data[[data.type.ind]][, feature.ind])))      

            clustering.single.feature <- spectralClustering(affinity.matrix, 
                num.of.clusters.fused)

            NMI.scores[[data.type.ind]][feature.ind] <- calNMI(clustering.fused, 
                clustering.single.feature)      
        }    

        NMI.ranks[[data.type.ind]] <- rank(-NMI.scores[[data.type.ind]],
             ties.method="first")
    }
    
    return(list(NMI.scores, NMI.ranks))  
}
