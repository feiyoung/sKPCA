\name{recerrfun}
\alias{recerrfun}
\title{Evaluate the empirical reconstruction error}
\description{
  Evaluate the reconstruction error based on object with 'sKPCA' class or 'oKPCA' class.
}
\usage{
  recerrfun(obj)
}
\arguments{
  \item{obj}{an object with 'sKPCA' class or 'oKPCA' class.}
}
\value{
  return a numeric scalar.
}
\references{
  nothing.
}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  nothing
}
\examples{
  library(kernlab)
  data(iris)
  test <- sample(1:150,20)
  # Sketched Kernel PCA
  d <- 2
  fkpca <- sKPCA(iris[-test,-5], d=d, seed = 8)
  recerrfun(fkpca)
}
\keyword{Reconstruction error}
