displayClustersWithHeatmap <- function (W, group, ColSideColors=NULL,...) {
    # Visualize the clusters present in given similarity matrix with sample 
    #   information displayed by given colors.
    #
    # Args:
    #   W: Affinity matrix
    #   group: labels of cluster groups 
    #   ColSideColors: Character vector of length(group) containing color names 
    #       for horizontal side bar to annotate columns of W, OR a character 
    #       matrix with number of rows matching number of rows in W.
    #
    # Returns:
    #   NULL

    normalize <- function(X) X/rowSums(X)
    ind <- sort(as.vector(group), index.return = TRUE)
    ind <- ind$ix

    diag(W) <- median(as.vector(W))
    W <- normalize(W)
    W <- W + t(W)

    if(is.null(ColSideColors)){
        heatmap(W[ind, ind],scale="none",Rowv=NA,Colv=NA,...)
    }
    else{
        if(is.vector(ColSideColors)){
            heatmap(W[ind, ind],scale="none",Rowv=NA,Colv=NA,
                ColSideColors=ColSideColors[ind],...)
        }

        else{
            heatmapPlus(W[ind, ind],scale="none",Rowv=NA,Colv=NA,
                ColSideColors=ColSideColors[ind,],...)
        }
    }
    return()
}
