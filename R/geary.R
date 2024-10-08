

GRfun <- function(i, var, weights) matrix(weights[i,], nrow=1) %*% ((var[i] - var)^2)

#' Local Geary 
#'
#' @export
#' @description A local indicator of spatial association based on the Geary Ratio (Geary's C) for exploratory spatial data analysis. Large values of this statistic highlight local outliers, that is, values that are not like their neighbors. 
#'
#' @param x Numeric vector of length `n`. By default, this will be standardized using the `scale` function.
#' @param w An `n x n` spatial connectivity matrix. See \link[geostan]{shape2mat}.
#' @param digits Number of digits to round results to.
#' @param scale If `TRUE`, then `x` will automatically be standardized using `scale(x, center = TRUE, scale = TRUE)`.
#' @param warn If `FALSE`, no warning will be printed to inform you when observations with `NA` values have been dropped, or if any observations without neighbors have been found.
#' @param na.rm If `na.rm = TRUE`, observations with `NA` values will be dropped from both `x` and `w`.
#'
#' @details
#' 
#' Local Geary's C is found in the numerator of the Geary Ratio (GR). For the \eqn{i^{th}} observation, the local Geary statistic is
#' \deqn{C_i = \sum_j w_{i,j} * (x_i - x_j)^2}
#' Hence, local Geary values will be largest for those observations that are most unlike their neighboring values. If a binary connectivity matrix is used (rather than row-standardized), then having many neighbors can also increase the value of the local Geary statistic. For most purposes, the row-standardized spatial weights matrix may be the more appropriate choice.
#'
#' @return The function returns a vector of numeric values, each value being a local indicator of spatial association (or dissimilarity), ordered as `x`.
#'
#' @source
#'
#' Anselin, Luc. "Local indicators of spatial association&mdash;LISA." Geographical analysis 27, no. 2 (1995): 93-115.
#'
#' Chun, Yongwan, and Daniel A. Griffith. Spatial Statistics and Geostatistics: Theory and Applications for Geographic Information Science and Technology. Sage, 2013.
#'
#' @examples
#' library(ggplot2)
#' data(georgia)
#' x <- log(georgia$income)
#' w <- shape2mat(georgia, "W")
#' lisa <- lg(x, w)
#' hist(lisa)
#' ggplot(georgia) +
#'   geom_sf(aes(fill = lisa)) +
#' scale_fill_gradient(high = "navy",  
#'                    low = "white") 
#'  ## or try: scale_fill_viridis() 
#' 
lg <- function(x, w, digits = 3, scale = TRUE, na.rm = FALSE, warn = TRUE) {
    check_sa_data(x, w)
    na_idx <- which(is.na(x))
    if (na.rm == TRUE && length(na_idx) > 0) {   
        if (warn) message(length(na_idx), " NA values found in x. They will be removed from the data before calculating the Local Geary statistic. If matrix w was row-standardized, it may no longer be. To address this, you can use a binary connectivity matrix, using style = 'B' in shape2mat.")
        x <- x[-na_idx]
        w <- w[-na_idx, -na_idx]
    }
    z <- scale(x, center = scale, scale = scale)
    n <- length(z)    
    LISD <- sapply(1:n, FUN = GRfun, var = z, weights = w)
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        LISD[zero.idx] <- NA        
        if (warn) message(length(zero.idx), " observations with no neighbors found. Their local C values are being set to NA.")
    }             
    return(round(LISD, digits = digits))
}


#' The Geary Ratio 
#'
#' @export
#' @description An index for spatial autocorrelation. Complete spatial randomness (lack of spatial pattern) is indicated by a Geary Ratio (GR) of 1; positive autocorrelation moves the index towards zero, while negative autocorrelation will push the index towards 2. 
#'
#' @param x Numeric vector of length `n`. By default, this will be standardized using the `scale` function.
#' @param w An `n x n` spatial connectivity matrix. See \link[geostan]{shape2mat}.
#' @param digits Number of digits to round results to.
#' @param warn If `FALSE`, no warning will be printed to inform you when observations with `NA` values have been dropped, or if any observations without neighbors have been found.
#' @param na.rm If `na.rm = TRUE`, observations with `NA` values will be dropped from both `x` and `w`.
#'
#' @return Returns the Geary ratio (a single numeric value).
#' 
#' @details
#'
#' The Geary Ratio is an index of spatial autocorrelation. The
#' numerator contains a series of sums of squared deviations, which
#' will be smaller when each observation is similar to its
#' neighbors. This term makes the index sensitive to local outliers,
#' which is advantageous for detecting such outliers and for measuring
#' negative autocorrelation. The denominator contains the total sum of
#' squared deviations from the mean value. Hence, under strong
#' positive autocorrelation, the GR approaches zero. Zero spatial
#' autocorrelation is represented by a GR of 1. Negative
#' autocorrelation pushes the GR above 1, towards 2.  \deqn{GR =
#' \frac{n-1}{2K} \frac{M}{D}} \deqn{M = \sum_i \sum_j w_{i,j} (x_i -
#' x_j)^2 } \deqn{D = \sum_i (x_i - \overline{x})^2 }
#'
#' Observations
#' with no neighbors are removed before calculating the GR. (The
#' alternative would be for those observations to contribute zero to the
#' numerator---but zero is not a neutral value, it represents strong
#' positive autocorrelation.)
#' 
#' @source
#' 
#' Chun, Yongwan, and Daniel A. Griffith. Spatial Statistics and Geostatistics: Theory and Applications for Geographic Information Science and Technology. Sage, 2013.
#' 
#' Geary, R. C. "The contiguity ratio and statistical mapping." The Incorporated Statistician 5, no. 3 (1954): 115-127_129-146.
#'
#' Unwin, Antony. "Geary's Contiguity Ratio." The Economic and Social Review 27, no. 2 (1996): 145-159.
#'
#' @examples
#' data(georgia)
#' x <- log(georgia$income)
#' w <- shape2mat(georgia, "W")
#' gr(x, w)
#'
#' 
gr <- function(x, w, digits = 3, na.rm = FALSE, warn = TRUE) {
    check_sa_data(x, w) 
    na_idx <- which(is.na(x))   
    if (na.rm == TRUE && length(na_idx) > 0) {   
        if (warn) message(length(na_idx), " NA values found in x. They will be removed from the data before calculating the Geary Ratio. If matrix w was row-standardized, it may not longer be. To address this, you can use a binary connectivity matrix, using style = 'B' in shape2mat.")
        x <- x[-na_idx]
        w <- w[-na_idx, -na_idx]
    }
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        if (warn) message(length(zero.idx), " observations with no neighbors found. They will be removed from the data before calculating the Geary Ratio.")
        x <- x[-zero.idx]
        w <- w[-zero.idx, -zero.idx]
    }
    n <- length(x)    
    top <- (n-1) * sum(sapply(1:n, FUN = GRfun, var = x, weights = w))
    A <- sum(Matrix::rowSums(w))    
    bottom <- 2 * A * sum(scale(x, center = TRUE, scale = FALSE)^2)
    GR <- top / bottom
    return(round(GR, digits = digits))
}

