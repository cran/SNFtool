getColorsForGroups <- function(group,
    colors=c("red","blue","green","purple","grey","cyan","brown","pink")){
    # Constructs a vector of colours given a numeric vector of cluster labels.
    #   If the number of groups > 8 a vector with length equal to the number of 
    #   groups  must be provided.
    # 
    # Args:
    #   group: A numeric vector of group labels
    #   colors: Provided to override default colour options, must be provided if
    #       number of groups exceeds eight. Length of provided colours must be 
    #        >= number of groups.
    #
    # Returns:
    #   A vector of characters specifying colours for each group
    #

    cluster.colors <- group

    if(max(group) <= length(colors)){
        for(i in 1:max(group)){
            cluster.colors[which(group==i)] <- colors[i]
        }
    return(cluster.colors)

    } else {
        warning(paste("ERROR: Not enough colors using the default color argument",
            "for the different groups, PLEASE inform the colors argument",
             sep=""))
        return(NULL)
  }
}
