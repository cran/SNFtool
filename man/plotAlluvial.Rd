\name{plotAlluvial}
\alias{plotAlluvial}
\title{
Plot Alluvial 
}
\description{
This function plots an alluvial (Parallel coordinate plot) of sample clusterings for a specified number of clusters. Samples can be coloured by providing a vector of colours, allowing for the visualization of sample properties over a range of clustering number choices.

*This is a wrapper function calling the \href{https://CRAN.R-project.org/package=alluvial}{Alluvial Package} (Bojanowski M. & Edwards R)

}
\usage{
plotAlluvial(W, clust.range, color.vect)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
  Affinity matrix of dimension n.samples by n.samples
}
  \item{clust.range}{
  Integer vector specifying the number of clusters for each clustering
}
  \item{color.vect}{
  A vector of color's of length n.samples to colour the samples
}
}
\value{
Plots an alluvial plot for range of clustering choices. 
}
\author{
Daniel Cole
}
\seealso{More information on \href{https://CRAN.R-project.org/package=alluvial}{Alluvial Package}} 
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
K <- 20
alpha <- 0.5
iter <- 20

data(Data1)
data(Data2)

dist1 <- (dist2(as.matrix(Data1), as.matrix(Data1)))^(1/2)
dist2 <- (dist2(as.matrix(Data2), as.matrix(Data2)))^(1/2)

W1 <- affinityMatrix(dist1, K, alpha)
W2 <- affinityMatrix(dist2, K, alpha)

W <- SNF(list(W1, W2), K, iter)

#Plots the alluvial with no colouring
plotAlluvial(W, 2:5)

#Change the colour of all samples a single colour
plotAlluvial(W, 2:5, col="red")

colour.breaks <- 30
#This will assign each sample to one of colour.breaks colour bins between green and red.
colFunc <- colorRampPalette(c("green", "red"))
colours <- colFunc(colour.breaks)[as.numeric(cut(Data1[,1],breaks=colour.breaks))]
plotAlluvial(W, 2:5, col=colours)
}
