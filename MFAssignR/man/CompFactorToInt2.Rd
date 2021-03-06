% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities_Halogen.R
\name{CompFactorToInt2}
\alias{CompFactorToInt2}
\title{Component factor to integer}
\usage{
CompFactorToInt2(x)
}
\arguments{
\item{x}{component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})}
}
\value{
numeric; integer representation for passed in factor component
}
\description{
Given a component as a factor returns the correponding integer representation.
For example, "C" -> 1, "H" -> 2, ..., "Z" -> 9
}
\details{
This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
 environment.
}
\examples{
CompFactorToInt("C")

CompFactorToInt('C')

}
