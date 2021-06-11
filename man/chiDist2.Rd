\name{chiDist2}
\alias{chiDist2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Chi-squared distances
}
\description{
Wrapper function chi2Dist imported from 'ExPosition' package.
Computes the Chi-squared distances between all pairs of data point given
}
\usage{
chiDist2(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
A data matrix where each row is a different data point
}
}
\value{
Returns an N x N matrix where N is the number of rows in X. element (i,j) is the squared Chi-squared distance between ith data point in X and jth data point in X.
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir

}
\examples{

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

## Calculate distance matrices(here we calculate Euclidean Distance, 
## you can use other distance, e.g. correlation)
Dist1 = chiDist2(as.matrix(Data1))
Dist2 = chiDist2(as.matrix(Data2))

}
