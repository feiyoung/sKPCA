\name{sigmaselect}
\alias{sigmaselect}
\title{Select the kernel width}
\description{
  Select the kernel width sigma based on the average interpoint distance
}
\usage{
  sigmaselect(X)
}
\arguments{
  \item{X}{a n*p matrix, the data matrix.}
}

\value{
  return a positive scalar.
}
\examples{
  library(kernlab)
  data(iris)
  test <- sample(1:150,20)
  # Sketched Kernel PCA
  sigma <- sigmaselect(iris[-test,-5])
  kern <- rbfdot(sigma)
  d <- 2
  fkpca <- sKPCA(iris[-test,-5], d=d, kern=kern, seed = 8)
}
\keyword{Width slection}
