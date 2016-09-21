\name{nlsw88}
\alias{nlsw88}
\docType{data}
\title{
NLSW, 1988 extract
}
\description{
National Longitudinal Surveys, Women sample
}
\usage{data(nlsw88)}
\format{
   A data frame with 2246 observations on the following 17 variables:
  \describe{
    \item{\code{idcode}}{a numeric vector, NLS id}
    \item{\code{age}}{a numeric vector, age at current year}
    \item{\code{race}}{a numeric vector, race}
    \item{\code{married}}{a numeric vector}
    \item{\code{never_married}}{a numeric vector}
    \item{\code{grade}}{a numeric vector, current grade completed}
    \item{\code{collgrad}}{a numeric vector, college graduate}
    \item{\code{south}}{a numeric vector, lives in south}
    \item{\code{smsa}}{a numeric vector, lives in SMSA}
    \item{\code{c_city}}{a numeric vector, lives in central city}
    \item{\code{industry}}{a numeric vector}
    \item{\code{occupation}}{a numeric vector}
    \item{\code{union}}{a numeric vector, union worker}
    \item{\code{wage}}{a numeric vector, hourly wage}
    \item{\code{hours}}{a numeric vector, usual hours worked}
    \item{\code{ttl_exp}}{a numeric vector, total work experience}
    \item{\code{tenure}}{a numeric vector, job tenure(years)}
  }
}
\details{
The NLSW88 data contains data of a group of women in their 30s and early 40s to study labor force patterns.
}
\source{
Stata website
}
\references{
Stata website: http://www.stata-press.com/data/r10/g.html
}
\examples{
data(nlsw88)
attach(nlsw88) 
plot(wage,tenure)
}
\keyword{datasets}
