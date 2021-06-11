plotAlluvial <- function(W, clust.range, color.vect="gray") {
    # Plots alluvial of patient clusterings for clustering into clust.range clusters.
    #
    # Args:
    #   W: A similarity matrix to be used in spectral clustering
    #   clust.range: An integer vector specifying the number of clusters to be
    #	    chosen for spectral clustering
    #   color.vect: A colour vector of length dim(W)[[1]] specifying the 
    #       colorings for patients in the alluvial (default all gray)
    #
    # Returns:
    #    NULL alluvial plot is output to k
    
    #Wrapper of spectralClustering to account for when clust.num=1 
    get_spectral_labels <- function(aff, clust.num){
        if(clust.num == 1){
            return(rep(1,dim(aff)[[1]]))
        }
        else{
            return(spectralClustering(aff, clust.num))
        }
    }
    
    #Error check for any value in clust.range ==1  or greater than N (dim(W)[[1]])
    if(any(clust.range) < 1){
        stop('All numbers in clust.range must be greater than or equal to 1.')
    }
   
    #Initialize clust.map matrix specifying sample-groups for clust.range clusterings 
    n.samples <- dim(W)[[1]]
    n.clusterings <- length(clust.range)
    clust.map <- matrix(rep(NA, n.samples*n.clusterings), n.samples, n.clusterings)
    
    clust.range <- unique(clust.range)
    

    #Generates (patient X num.clusters) mapping for clustering 
    for(i in c(1:length(clust.range))){
        clust.map[,i] <- get_spectral_labels(W, clust.range[[i]])
    }
    colnames(clust.map) <- paste("",clust.range,sep="")
    alluvial(clust.map, freq=rep(1,n.samples), col=color.vect)
}
