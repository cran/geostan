% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geary.R
\name{lg}
\alias{lg}
\title{Local Geary}
\source{
Anselin, Luc. "Local indicators of spatial association—LISA." Geographical analysis 27, no. 2 (1995): 93-115.

Chun, Yongwan, and Daniel A. Griffith. Spatial Statistics and Geostatistics: Theory and Applications for Geographic Information Science and Technology. Sage, 2013.
}
\usage{
lg(x, w, digits = 3, scale = TRUE, na.rm = FALSE, warn = TRUE)
}
\arguments{
\item{x}{Numeric vector of length \code{n}. By default, this will be standardized using the \code{scale} function.}

\item{w}{An \verb{n x n} spatial connectivity matrix. See \link[geostan]{shape2mat}.}

\item{digits}{Number of digits to round results to.}

\item{scale}{If \code{TRUE}, then \code{x} will automatically be standardized using \code{scale(x, center = TRUE, scale = TRUE)}.}

\item{na.rm}{If \code{na.rm = TRUE}, observations with \code{NA} values will be dropped from both \code{x} and \code{w}.}

\item{warn}{If \code{FALSE}, no warning will be printed to inform you when observations with \code{NA} values have been dropped, or if any observations without neighbors have been found.}
}
\value{
The function returns a vector of numeric values, each value being a local indicator of spatial association (or dissimilarity), ordered as \code{x}.
}
\description{
A local indicator of spatial association based on the Geary Ratio (Geary's C) for exploratory spatial data analysis. Large values of this statistic highlight local outliers, that is, values that are not like their neighbors.
}
\details{
Local Geary's C is found in the numerator of the Geary Ratio (GR). For the \eqn{i^{th}} observation, the local Geary statistic is
\deqn{C_i = \sum_j w_{i,j} * (x_i - x_j)^2}
Hence, local Geary values will be largest for those observations that are most unlike their neighboring values. If a binary connectivity matrix is used (rather than row-standardized), then having many neighbors can also increase the value of the local Geary statistic. For most purposes, the row-standardized spatial weights matrix may be the more appropriate choice.
}
\examples{
library(ggplot2)
data(georgia)
x <- log(georgia$income)
w <- shape2mat(georgia, "W")
lisa <- lg(x, w)
hist(lisa)
ggplot(georgia) +
  geom_sf(aes(fill = lisa)) +
scale_fill_gradient(high = "navy",  
                   low = "white") 
 ## or try: scale_fill_viridis() 

}
