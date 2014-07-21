\name{low_rank}
\alias{low_rank}
\title{
Implementation of HSSVD framework for biclusering with heterogeneous variance.
}
\description{
HSSVD is a recently developed data mining tool for discovering subgroups of 
patients and genes which simultaneously display unusual levels of variability 
compared to other genes and patients. Previous biclustering methods were 
restricted to mean level detection, while the new method can detect both 
mean and variance biclusters.}
\usage{
low_rank(x, ranks = c(rep(NA,3)), sparse = rep(FALSE, 3), 
est = rep("wold", 3), add=c(1,0,0),...)
}
\arguments{
  \item{x}{
Data matrix to bicluster, if data range is [0; 1] then logit transform is
recommended. If [0;inf] then log transform is recommended, otherwise no transform or
standardization is needed.}
  \item{ranks}{
Integer vector (dim=3) of rank input for: rescale, mean approximation and
variance approximation steps. Default is c(NA,NA,NA), which means that
Wold or Gabriel-style cross validation will be used for rank estimation.}
  \item{sparse}{
Vector of logical indicators (dim=3). TRUE=sparse input of left and right singular vector
is needed. Default is c(FALSE,FALSE,FALSE), see Details for insight.}
  \item{est}{
Character vector (dim=3) giving estimation methods for each step. Only "GB" or "Wold"
method are available. Default is "Wold"}
  \item{add}{
Integer vector (dim=3), where 0 <= add <= 2, which is added to the rank estimate
at each step for dense matrices.}
  \item{...}{Additional input variables for FIT.SSVD. See Details section
for more information.}
}
\details{
FIT-SSVD method (3) assumes that over half of the rows are approximate null rows
with little or no signal and similar for the columns, their preselecting row and column
steps for initial values and rank estimation will under-select informative rows/columns if
all rows/columns are informative or the dimension of rows/columns are not large enough,
both of which can happen in biological research. So in practice, the initial value would
be the vanilla singular value decomposition's result. For the rank approximation, the
Gabriel-style "block" holdout tend to give a larger rank estimation than the Wold-style
"speckled" holdout. For our experience Wold method may be better for the
large data set. The rank estimation steps are slow, especially when using Wold-style
 method, so we suggest to run in batch mode or run on clustering if rank estimation
is needed. 

The FIT.SSVD has several input variables that can be adjusted by the
user if desired. However, these are currently limited to \cr
dothres : 'hard' or 'soft' thresholding; default='hard' \cr
n.step : maximum number of iterations in Algorithm 1; default=100 \cr
n.err : number of Bootstrap samples in Algorithm 3; default=100 \cr
}
\value{
\item{call}{record of call to low_rank}
\item{rescale}{Sparse SVD approximation result of rescale, typically this type
object contain singular value (d), left (u) and right (v) singular vectors.
(Step 2)}
\item{result_mean}{Sparse SVD approximation result of Y . (Step 3)}
\item{result_var}{Sparse SVD approximation result of Z. (Step 4)}
\item{ranks}{vector of estimated ranks for (U, Y_tilde, Z_tilde).}
\item{bgmean}{Mean for null cluster. (Step 5)}
\item{bgstd}{Standard deviation for null cluster. (Step 5)}
\item{back}{Logical matrix where TRUE=in null cluster.}
\item{mean_app}{Mean approximation of X. (Step 6)}
\item{std_app}{Standard deviation approximation of X. (Step 6)}
}
\references{
Chen G, Sullivan PF, Kosorok MR (2013) Biclustering with heterogeneous variance.
Proceedings of the National Academy of Sciences, doi:10.1073/pnas.1304376110.

Owen AB, Perry PO (2009) Bi-cross-validation of the SVD and the nonnegative matrix
factorization. The Annals of Applied Statistics 3:564-594.

Yang, D., Ma, Z., and Buja, A. (2011) A sparse SVD method for high-dimensional
data arXiv:1112.2433.
}
\author{
Authors: Guanhua Chen and Michael R. Kosorok \cr
Contributor: Shannon T. Holloway sthollow@ncsu.edu \cr
Maintainer: Guanhua Chen <guanhuac@live.unc.edu> \cr
}


\examples{
data(Methylation)
beta_name <- colnames(Methylation)[grep("AVG_Beta",colnames(Methylation))]
ds <- as.matrix(Methylation[beta_name],ncol=length(beta_name))
info <- t(ds)
## The methylation data takes value between [0,1]; 
## therefore, we do a logit tranformation ##
info <- log(info/(1-info))

result <- low_rank(info,ranks=c(8,4,4))
}
