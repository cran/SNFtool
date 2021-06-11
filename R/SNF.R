SNF <- function(Wall, K=20, t=20) {
    # Similarity Network Fusion takes multiple views of a network (Wall) and
    # fuses them together to create a overall affinity matrix.
    #
    # Args:
    #   Wall: List of matrices, each element is a square symmetric affinity 
    #       matrix.
    #   K: Number of neighbors used in the K-nearest neighbours step,??? more details???
    #   t: Number of iterations for the diffusion process
    #
    # Returns:  
    #   W: Unified similarity graph of all data types in Wall. 

    check_wall_names <- function(Wall){
        # Checks if dimnames are consistant across all matrices in Wall
        #   #Move to internal functions?
        # Args:
        #   Wall: List of matrices
        # Returns:
        #   logical: True/False indicator of dimnames equivalence
        name_match <- function(names_A, names_B){
            return(identical(dimnames(names_A), dimnames(names_B)))
        }

        return(all(unlist(lapply(Wall, FUN=name_match, Wall[[1]]))))
    }

    #Check if Wall names are consistant across all matrices in Wall
    wall.name.check <- check_wall_names(Wall)
    wall.names <- dimnames(Wall[[1]])
    if(!wall.name.check){
        warning("Dim names not consistent across all matrices in Wall.
            Returned matrix will have no dim names.")
    }
 
    LW <- length(Wall)

    #Normalization method for affinity matrices
    normalize <- function(X){
        row.sum.mdiag <- rowSums(X) - diag(X) 
        #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
        row.sum.mdiag[row.sum.mdiag == 0] <- 1   
        X <- X/(2*(row.sum.mdiag))
        diag(X) <- 0.5
        return(X)
    }
    
    #Normalize different networks to avoid scale problems.
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    for(i in 1:LW){
      Wall[[i]] <- normalize(Wall[[i]])
      Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
    }
    
    ### Calculate the local transition matrix. (KNN step?)
    for(i in 1:LW){
      newW[[i]] <- (.dominateset(Wall[[i]], K))
    }
    
    #Perform the diffusion for t iterations
    for (i in 1:t) {
        for(j in 1:LW){
            sumWJ <- matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
            for(k in 1:LW){
                if(k != j) {
                    sumWJ <- sumWJ + Wall[[k]]
                }
            }
            nextW[[j]] <- newW[[j]] %*% (sumWJ/(LW-1)) %*% t(newW[[j]])
        }

        #Normalize each new obtained networks.
        for(j in 1 : LW){
          Wall[[j]] <- normalize(nextW[[j]])
          Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2;
        }
    }
    
    # Construct the combined affinity matrix by summing diffused matrices
    W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:LW){
        W <- W + Wall[[i]]
    }

    W <- W/LW
    W <- normalize(W)
    W <- (W + t(W)) / 2

    if(wall.name.check){
        dimnames(W) <- wall.names
    } 

    return(W)  
}
