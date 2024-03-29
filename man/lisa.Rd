% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moran.R
\name{lisa}
\alias{lisa}
\title{Local Moran's I}
\source{
Anselin, Luc. "Local indicators of spatial association—LISA." Geographical Analysis 27, no. 2 (1995): 93-115.
}
\usage{
lisa(x, w, type = TRUE, scale = TRUE, digits = 3)
}
\arguments{
\item{x}{Numeric vector of length \code{n}.}

\item{w}{An \verb{n x n} spatial connectivity matrix. See \link[geostan]{shape2mat}. If \code{w} is not row standardized (\code{all(Matrix::rowSums(w) == 1)}), it will automatically be row-standardized.}

\item{type}{Return the type of association also (High-High, Low-Low, High-Low, and Low-High)? Defaults to \code{FALSE}.}

\item{scale}{If \code{TRUE}, then \code{x} will automatically be standardized using \code{scale(x, center = TRUE, scale = TRUE)}. If \code{FALSE}, then the variate will be centered but not scaled, using \code{scale(x, center = TRUE, scale = FALSE)}.}

\item{digits}{Number of digits to round results to.}
}
\value{
If \code{type = FALSE} a numeric vector of lisa values for exploratory analysis of local spatial autocorrelation. If \code{type = TRUE}, a \code{data.frame} with columns \code{Li} (the lisa value) and \code{type}.
}
\description{
A local indicator of spatial association (LISA) based on Moran's I (the Moran coefficient) for exploratory data analysis.
}
\details{
The values of \code{x} will automatically be centered first with \code{z = scale(x, center = TRUE, scale = scale)} (with user control over the \code{scale} argument). The LISA values are the product of each \code{z} value with the weighted sum of their respective surrounding value: \deqn{I_i = z_i \sum_j w_{ij} z_j} (or in R code: \code{lisa = z * (w \%*\% z)}). These are for exploratory analysis and model diagnostics.

An above-average value (i.e. positive z-value) with positive mean spatial lag indicates local positive spatial autocorrelation and is designated type "High-High"; a low value surrounded by high values indicates negative spatial autocorrelation and is designated type "Low-High", and so on.

This function uses Equation 7 from Anselin (1995). Note that the \code{spdep} package uses Formula 12, which divides the same value by a constant term \eqn{\sum_i z_i^2/n}. So the \code{geostan} version can be made equal to the \code{spdep} version by dividing by that value.
}
\examples{
library(ggplot2)
library(sf)
data(georgia)
w <- shape2mat(georgia, "W")
x <- georgia$ICE
li = lisa(x, w)
head(li)
ggplot(georgia, aes(fill = li$Li)) +
  geom_sf() +
  scale_fill_gradient2()
}
\seealso{
\code{\link[geostan]{moran_plot}}, \code{\link[geostan]{mc}}, \code{\link[geostan]{aple}}, \code{\link[geostan]{lg}}, \code{\link[geostan]{gr}}
}
