\name{sigmaslect}
\alias{sigmaslect}
\title{Select the kernel width}
\description{
  Select the kernel width sigma based on the average interpoint distance
}
\usage{
  sigmaslect(X)
}
\arguments{
  \item{X}{a n*p matrix, the data matrix.}
}
\details{
   sigma = sqrt(frac(1)(n^2)sumisumj |X_i -X_j|_2^2)
}
\value{
  return a positive scalar.
}
\references{
  Wang, S., Gittens, A., and Mahoney, M. W. (2019). Scalable kernel k-means clustering with nystrom approximation: relative-error bounds. The Journal of Machine Learning Research, 20(1):431– 262.
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
  sigma <- sigmaslect(iris[-test,-5])
  kern <- rbfdot(sigma)
  d <- 2
  fkpca <- sKPCA(iris[-test,-5], d=d, kern=kern, seed = 8)
}
\keyword{Width slection}
