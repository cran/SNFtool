groupPredict <- function(train, test, groups, K=20, alpha=0.5, t=20, method=1){
    # Predicts subtype of new patients from labeled training set of patients 
    #   using label propigation or local and global consistency.
    #
    # Args:
    #   train: List affinity matrices for samples with known labels
    #   test: List affinity matrices for samples with unknown labels.
    #       Length of test must match length of train (and order?)
    #   groups: Labels specifying the groups in train
    #   K: SNF parameter for number of neighbours in KNN step
    #   alpha: SNF Hyperparameter 
    #   t: SNF varaible - number of iterations
    #   method: 0/1 specifies method used (1) Label propagation or
    #       (0) Local & global consistency.
    #
    # Returns: 
    #   Vector of new labels assigned to the test samples 
    Wi <- vector("list", length=length(train))
    
    for (i in 1:length(train)){
        view <- standardNormalization(rbind(train[[i]],test[[i]]))
        Dist1 <- dist2(view, view)
        Wi[[i]] <- affinityMatrix(Dist1, K, alpha)
    }
    
    W <- SNF(Wi,K,t)
    Y0 <- matrix(0,nrow(view), max(groups))
    for (i in 1:length(groups)){
        Y0[i,groups[i]] <- 1
    }

    Y <- .csPrediction(W,Y0,method)
    newgroups <- rep(0,nrow(view))
    for (i in 1:nrow(Y)){
        newgroups[i] <- which(Y[i,] == max(Y[i,]))
    }
    
    return (newgroups)
}
