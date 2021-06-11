\name{heatmapPlus}
\alias{heatmapPlus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display heatmap for clusters 
}
\description{
Visualize clusters with heatmap
}
\usage{
heatmapPlus(x, 
Rowv = NULL, 
Colv = if (symm) "Rowv" else NULL, 
distfun = dist, 
hclustfun = hclust, 
reorderfun = function(d, w) reorder(d, w), 
add.expr, 
symm = FALSE, 
revC = identical(Colv, "Rowv"), 
scale = c("row", "column", "none"), 
na.rm = TRUE, 
margins = c(5, 5), 
ColSideColors, 
RowSideColors, 
cexRow = 0.2 + 1/log10(nr), 
cexCol = 0.2 + 1/log10(nc), 
labRow = NULL, 
labCol = NULL, 
main = NULL, 
xlab = NULL, 
ylab = NULL, 
keep.dendro = FALSE, 
verbose = getOption("verbose"), ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric matrix of the values to be plotted.
}
 \item{Rowv}{
determines if and how the row dendrogram should be computed and reordered. Either a dendrogram or a vector of values used to reorder the row dendrogram or NA to suppress any row dendrogram (and reordering) or by default, NULL.
}
 \item{Colv}{
determines if and how the column dendrogram should be reordered. Has the same options as the Rowv argument above and additionally when x is a square matrix, Colv = "Rowv" means that columns should be treated identically to the rows.
}
 \item{distfun}{
function used to compute the distance (dissimilarity) between both rows and columns. Defaults to dist.
}
 \item{hclustfun}{
function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to hclust.
}
 \item{reorderfun}{
function(d,w) of dendrogram and weights for reordering the row and column dendrograms. The default uses reorder.dendrogram.
}
 \item{add.expr}{
expression that will be evaluated after the call to image. Can be used to add components to the plot.
}
 \item{symm}{
logical indicating if x should be treated symmetrically; can only be true when x is a square matrix.
}
 \item{revC}{
logical indicating if the column order should be reversed for plotting, such that e.g., for the symmetric case, the symmetry axis is as usual.
}
 \item{scale}{
character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "row" if symm false, and "none" otherwise. 
}
 \item{na.rm}{
logical indicating whether NA's should be removed.

}
 \item{margins}{
numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively.
}
 \item{ColSideColors}{
(optional) character matrix with number of rows matching number of rows in x. Each column is plotted as a row similar to heatmap()'s ColSideColors. colnames() will be used for labels if present.
}
 \item{RowSideColors}{
(optional)character matrix with number of rows matching number of columns in x. Each column is plotted as a row similar to heatmap()'s RowSideColors. colnames() will be used for labels if present.
}
 \item{cexRow}{
positive numbers, used as cex.axis in for the row axis labeling. The defaults currently only use number of rows.
}
 \item{cexCol}{
positive numbers, used as cex.axis in for the column axis labeling. The defaults currently only use number of columns.

}
 \item{labRow}{
character vectors with row labels to use; these default to rownames(x).
}
 \item{labCol}{
character vectors with column labels to use; these default to colnames(x). 
}
 \item{main}{
y-axis titles; defaults to none.
}
 \item{xlab}{
x-axis titles; defaults to none.
}
 \item{ylab}{
main titles; defaults to none.
}
 \item{keep.dendro}{
logical indicating if the dendrogram(s) should be kept as part of the result (when Rowv and/or Colv are not NA).
}
 \item{verbose}{
logical indicating if information should be printed.
}
 \item{...}{
additional arguments passed on to image, e.g., col specifying the colors
}


}
\value{
Invisibly, a list with components

rowInd
row index permutation vector as returned by order.dendrogram.

colInd
column index permutation vector.

Rowv
the row dendrogram; only if input Rowv was not NA and keep.dendro is true.

Colv
the column dendrogram; only if input Colv was not NA and keep.dendro is true.
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
  z = matrix(rnorm(30),nrow=5,ncol=6);
  rlab = matrix(as.character(c(1:5,2:6,3:7,4:8)),nrow=5,ncol=4);
  clab = matrix(as.character(c(1:6,6:1)),nrow=6,ncol=2);
  colnames(rlab) = LETTERS[1:dim(rlab)[2]];
  colnames(clab) = 1:dim(clab)[2];
  heatmapPlus(z,ColSideColors=clab,RowSideColors=rlab);
}
