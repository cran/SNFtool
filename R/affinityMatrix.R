affinityMatrix <- function(diff,K=20,sigma=0.5) {
    # Computes the affinity matrix for a given distance matrix
    # 
    # Args:
    #   diff: Distance matrix 
    #   K: Number of nearest neighbours to sparsify similarity
    #   sigma: Variance for local model
    #
    # Returns:
    #   Affinity matrix using exponential similarity kernel scaled by k nearest
    #       neighbour average similarity
    #

    N <- nrow(diff)
    
    diff <- (diff + t(diff)) / 2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff,2,sort)))

    finiteMean <- function(x){
        return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
    
    avg <- function(x,y){
        return((x+y)/2)
    }
    Sig <- outer(means,means,avg)/3*2 + diff/3 + .Machine$double.eps;
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(diff, 0, sigma*Sig, log = FALSE)
    
    W <- (densities + t(densities)) / 2
    return(W)
}
