\name{displayClustersWithHeatmap}
\alias{displayClustersWithHeatmap}
\title{
Display the similarity matrix by clusters with some sample information
}
\description{
Visualize the clusters present in the given similarity matrix as well as some sample information.
}
\usage{
displayClustersWithHeatmap(W, group, ColSideColors=NULL, ...)
}
\arguments{
  \item{W}{
Similarity matrix
}
  \item{group}{
A numeric vector containing the groups information for each sample in W such as the result of the spectralClustering function. The order should correspond to the sample order in W.
}
  \item{ColSideColors}{
(optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x, used by the heatmap function, OR a character matrix with number of rows matching number of rows in x.  Each column is plotted as a row similar to heatmap()'s ColSideColors by the heatmap.plus function.
}
  \item{...}{
other paramater that can be pass on to the heatmap (if ColSideColor is a NULL or a vector) or  heatmap.plus function (if ColSideColors is matrix)  
}
}
\details{
Using the heatmap or heatmap.plus function to display the similarity matrix
For representation purpose, the similarity matrix diagonal is set to the median value of W, the matrix is normalised and W = W + t(W) is applied
In this presentation no clustering method is ran the samples are ordered in function of their group label present in the group arguments.
}
\value{
Plots the similarity matrix using the heatmap function. Samples are ordered by the clusters provided by the argument groups with sample information displayed with a color bar if the ColSideColors argument is informed.
}
\author{
Florence Cavalli
}
\examples{
## First, set all the parameters:
K = 20;    # number of neighbors, usually (10~30)
alpha = 0.5;    # hyperparameter, usually (0.3~0.8)
T = 20;   # Number of Iterations, usually (10~20)

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

## Here, the simulation data (SNFdata) has two data types. They are complementary to each other. 
## And two data types have the same number of points. 
## The first half data belongs to the first cluster; the rest belongs to the second cluster.
truelabel = c(matrix(1,100,1),matrix(2,100,1)); ## the ground truth of the simulated data

## Calculate distance matrices
## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)

## If the data are all continuous values, we recommend the users to perform 
## standard normalization before using SNF, 
## though it is optional depending on the data the users want to use.  
# Data1 = standardNormalization(Data1);
# Data2 = standardNormalization(Data2);

## Calculate the pair-wise distance; 
## If the data is continuous, we recommend to use the function "dist2" as follows 
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2), K, T)

## With this unified graph W of size n x n, 
## you can do either spectral clustering or Kernel NMF. 
## If you need help with further clustering, please let us know. 

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2   							# number of clusters
group = spectralClustering(W,C); 	# the final subtypes information

## Get a matrix containing the group information 
## for the samples such as the SpectralClustering result and the True label
M_label=cbind(group,truelabel)
colnames(M_label)=c("spectralClustering","TrueLabel")

## ****
## Comments
## rownames(M_label)=names(spectralClustering) To add if the spectralClustering function 
## pass the sample ID as names.
## or rownames(M_label)=rownames(W) Having W with rownames and colmanes 
## with smaple ID would help as well.
## ***

## Use the getColorsForGroups function to assign a color to each group
## NB is more than 8 groups, you will have to input a vector 
## of colors into the getColorsForGroups function
M_label_colors=t(apply(M_label,1,getColorsForGroups))
## or choose you own colors for each label, for example:
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],
colors=c("blue","green")),"TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],
colors=c("orange","cyan")))

## Visualize the clusters present in the given similarity matrix 
## as well as some sample information
## In this presentation no clustering method is ran the samples 
## are ordered in function of their group label present in the group arguments
displayClustersWithHeatmap(W, group, M_label_colors[,"spectralClustering"]) 
displayClustersWithHeatmap(W, group, M_label_colors)
}