#' Data projectction
#'
#' Transforms the data by \eqn{- \boldsymbol{x}} and projects it on the
#' direction \eqn{\boldsymbol{v}}
#'
#' @param data Data points given as a matrix.
#' @param x Location as a vector.
#' @param v Direction as a vector.
#' @return The transformed and projected points. Either as a vector or dataframe
#'   if multiple points or directions are considered.
#' @examples
#' data <- SimCirlce(200)
#' x    <- c(-1, 0)
#' v    <- c( 1, 0)
#' z <- ProjectData(data, x, v)
#' ggplot(dataframe(z = t(z)), aes(z)) + geom_density()
#' @export
ProjectData <- function(data, x, v){
  if(is.matrix(x) | is.matrix(v)){
    stop("Only 1 point and 1 direction is implementend so far!")
  }
  # Normalizing the direction v
  v <- v / sqrt(sum(v^2))
  z <- colSums(v * (data - x))
  return(as.vector(z))
}

#' Weights
#'
#' Calculats the weights of the data for a given point \code{x}
#'
#' @param data Data in a matrix.
#' @param x Location as a vector.
#' @param h Bandwidth parameter.
#' @param type Type of the used kernel function.
#' @return Vector of the probability weights corresponding to the data.
#' @export
Weights <- function(data, x, h, type = "Gaussian"){
  w <- K((data - x) / h, type = type)
  w_norm <- w / sum(w)
  return(w_norm)
}

#' Weighted Variance
#'
#' Calculates the variance of a finite distribution with probability weigth \code{w} at \code{z}
#'
#' @param z Points with positive probability as vector.
#' @param w Probability weights as vector.
#' @return Variance of the finite distribution
#' @export
Weighted_Variance2 <- function(z, w){
  return(sum(w * z^2) - sum(w * z)^2)
}


#' Find Direction of Smallest Weighted Projected data variance
#'
#' Finding the direction minimizing the variance of the weighted projected data.
#'
#' @param data Data as matrix [d,n].
#' @param x Location as vector [d].
#' @param w Weights as vector [n].
#' @return List of \describe{ \item{v}{Vector [d] which minimizes the
#'   variance of the weighted projected data.}
#'   \item{directional_variance}{weighted projected data variance
#'   in direction of \code{v}} }
#' @export
FindDirection <- function(data, x, w){
  loc_data <- data - x
  d <- dim(loc_data)[1]
  Q <- loc_data %*% (diag(w) - w %*% t(w)) %*% t(loc_data)
  eigQ <- eigen(Q, symmetric = TRUE)
  var <- eigQ$values[d]
  v <- eigQ$vectors[ , d]
  return(list(v = v, directional_variance = var))
}

#' Find Reweighted Mode
#'
#' Finding the univariate mode of the weighted projected data, which is
#' reweighted by the kernel.
#'
#' @param z Projected data [n]
#' @param w Weights of the data points [n]
#' @param h Bandwidth
#' @return A list with the values: \item{m}{Mode}
#'   \item{logcon}{The fitted
#'   log-concave density of the function call
#'   \code{\link[logcondes]{activeSetLogCon}}}
#'   \item{reweighted_density}{The
#'   values of the reweighted density at the points of \code{logcon$x}}
#' @export
FindReweightedMode <- function(z, w, h){
  logcon <- activeSetLogCon(z, w = w)
  reweighted_density <- exp(logcon$phi + logcon$x^2 / h^2)
  # Calculation of the closest mode of the reweighted univariate log-concave
  # density to 0.
  knots <- logcon$x[as.logical(logcon$IsKnot)]
  phi   <- logcon$phi[as.logical(logcon$IsKnot)]
  nn    <- length(knots)
  # Calculate left and right derivative of reweighted phi at the knots
  slope <- (phi[-1] - phi[-nn]) / (knots[-1] - knots[-nn])
  left_derivative  <- c(Inf, slope + 2 * knots[-1] / h^2)
  right_derivative <- c(slope + 2 * knots[-nn] / h^2, -Inf)
  # Calculating the modes and closest mode
  modes_index <- which(left_derivative > 0 & right_derivative < 0)
  dist <- abs(knots[modes_index])
  m <- knots[modes_index[which.min(dist)]]
  return(list(m = m, logcon = logcon, reweighted_density = reweighted_density))
}

#' Find Smoothed Mode
#'
#' Finding the mode of the smoothed log-concave density with a bisection search.
#'
#' @param x Vector [m] of support points of the data.
#' @param phi Vector [m] with the values of the log-concave density at the
#'   support points.
#' @param gam Bandwidth parameter of the Gaussian smoothing.
#' @param isKnot Logical vector [m] with the information if the corresponding
#'   support point is a knot of the log-concave density. By default all support
#'   points are treated as knots.
#' @param F.smoothed Vector [m] with the values of the smoothe log-concave
#'   density at the support points caluclated with bandwidth \code{gam}; see
#'   also \pkg{logcondens}. If not \code{NULL} the starting interval of the
#'   bisection search will be chosen accordingly.
#' @param delta Stopping criteria. (default = 1e-10)
#' @return Mode of the smoothed log-concave density.
#' @export
FindSmoothedMode <- function(x, phi, gam, isKnot = TRUE, F.smoothed = NULL, delta = 1e-10){
  isKnot <- as.logical(isKnot)
  m <- length(x)
  if (!is.null(F.smoothed)){
    max <- which.max(F.smoothed)
    l <- ifelse(max > 1, x[max - 1], x[max])
    r <- ifelse(max < m, x[max + 1], x[max])
  } else{
    l <- x[1]
    r <- x[m]
  }
  knots <- x[isKnot]
  phi   <- phi[isKnot]
  while (r - l > delta){
    c     <- (l + r) / 2
    fc    <- SmoothedDensity(c, knots, phi, gam)$f1
    if (fc > 0){
      l <- c
    } else{
      r <- c
    }
  }
  return((l + r) / 2)
}

#' Derivatves of the smoothed log-concave density
#'
#' @param z Where the function is evaluated.
#' @param x Vector [m] of knots.
#' @param phi Vector [m] of values of the log-concave density corresponding to
#'   the knots \code{x}.
#' @param gam Bandwidth paramter of the Gaussian smoothing.
#' @return A list with the value, 1st  and second order derivative of the smoothed
#'   log-concave density at \code{z}.
#' @export
SmoothedDensity <- function(z, x, phi, gam){
  m    <- length(x)
  I1   <- 1:(m-1)
  I2   <- 2:m
  s    <- (phi[I2] - phi[I1]) / (x[I2] - x[I1])
  f    <- exp(phi[I1])
  tmp  <- exp(s * (z - x[I1]) + s^2 * gam^2 / 2)
  q0   <- tmp * (pnorm(z - x[I1] + s * gam^2, 0, gam) - pnorm(z - x[I2] + s * gam^2, 0, gam))
  tmp2 <- tmp * (dnorm(z - x[I1] + s * gam^2, 0, gam) - dnorm(z - x[I2] + s * gam^2, 0, gam))
  q1   <- s * q0 + tmp2
  q2   <- s * q1 + s * tmp2 + tmp * ((z - x[I2] + s * gam^2) / gam^2 *
          dnorm(z - x[I2] + s * gam^2, 0, gam) - (z - x[I1] + s * gam^2) / gam^2 *
          dnorm(z - x[I1] + s * gam^2, 0, gam))
  f0 <- sum(f * q0)
  f1 <- sum(f * q1)
  f2 <- sum(f * q2)
  return(list(f0 = f0, f1 = f1, f2 = f2))
}

#' Weight cutoff
#'
#' Cutoff weights with very small influance
#'
#' @param w Vector of weights [n]
#' @param alpha Percentage of weighted data to keep
#'
#'
#' @return A list with the values:
#' \item{w}{Vector of new weights [n]}
#' \item{isPositive}{Logical vector wheter or not the entry in \code{w} is 0 or not. [n]}
#' @export
NewWeights <- function(w, alpha = 0.99){
  ord <- order(w, decreasing = TRUE)
  pos_n <- sum(cumsum(w[ord]) < alpha) + 1
  isPositive <- rep(FALSE, length(w))
  isPositive[ord[1:pos_n]] <- TRUE
  new_w <- w * isPositive
  new_w <- new_w / sum(new_w)
  return(list(w = new_w, isPositive = isPositive))
}

##################################################
# Functions for the Subsapce Constraint Mean Shift
##################################################

#' Mean Shift
#'
#' Find the mean shift at the current point.
#'
#' @param data data as matrix [d, n].
#' @param x current point [d].
#' @param h Bandwidth of the Gaussian kernel
#' @return A list with the following values:
#' \item{m}{Mean shift vector [d].}
#' \item{w}{Weights for each datapoint [n].}
#' \item{loc_data}{Localized data by \code{x} and \code{h} [d, n].}
#' @export
GetMeanShift <- function(data, x, h){
  d <- dim(data)[1]
  loc_data <- (data - x) / h
  w <- exp(-colSums(loc_data^2) / 2)
  m <- rowSums(data * rep(w, each = d)) / sum(w) - x
  return(list(m = m, w = w, loc_data = loc_data))
}

#' Subsapce Projection
#'
#' Finds the subspace and projects the given vector on it
#'
#' @param data data as matrix [d, n].
#' @param m vector to project [d].
#' @param h Bandwidth
#' @return Projected vector [d].
#' @export
GetProjectedMeanShift <- function(data, x, h){
  d <- dim(data)[1]
  tmp <- GetMeanShift(data, x, h)
  loc_data <- tmp$loc_data
  w <- tmp$w
  m <- tmp$m
  theta <- LogDensity(x, h, method = 'KDE', data = data)
  # Estimated Hessian of the log-density at x.
  f <- sum(w)
  H <- matrix(-theta[-(1:d)] / f + kronecker(theta[1:d], theta[1:d]) / f^2, d, d)
  eigenH <- eigen(H, symmetric = TRUE)
  v <- eigenH$vectors[ , 1]
  return(sum(v * m) * v)
}

###################################################################
# Auxiliary Functions for the Newton step for the smoothed density
###################################################################

#' Derivatives of smooth (re-weighted) density
#'
#' Calculation of the smooth (re-weighted) density and the first and second order
#' derivative. To be precise, the density is estimated without normalization and
#' the derivatives thereof ar calculated. However, calculation of the modes are
#' not effected.
#'
#'
#' @param y Where the function should be evaluated.
#' @param logcon Fitted log-concave density of class \code{dlc} from package
#'   \code{logcondens}.
#' @param h Bandwith used for localization.
#' @return Vector of length3 containing \eqn{\hat{f}^*_h(z), \hat{f}^*_h(z)',
#'   \hat{f}^*(z)''}.
#' @export
GetDerivativeSRLC <- function(y, logcon, h){
  names(y) <- NULL
  IsKnot <- as.logical(logcon$IsKnot)
  x   <- logcon$x[IsKnot]
  phi <- logcon$phi[IsKnot]
  m   <- length(x)
  f   <- exp(phi)[-m]
  s   <- (phi[-1] - phi[-m]) / (x[-1] - x[-m])
  gam <- GetSmoothingBandwidth(logcon)
  # Ingriedients for the calculation of q0, q1, q2
  a <- s * (y - x[-m]) + s^2 * gam^2 / 2
  exp_a <- exp(a)
  u <- (x[-m] - y - s * gam^2) / gam
  v <- (x[-1] - y - s * gam^2) / gam
  # Auxiliary Functions
  q0 <- ifelse(is.infinite(exp_a), Inf, exp_a * (pnorm(v) - pnorm(u)))
  # Replace values wiht theoretical bounds if necessary.
  tmp_m     <- (u + v) / 2
  tmp_delta <- (v - u) / 2
  tmp_exp   <- exp(a - tmp_m^2 / 2)
  Phi_delta <- pnorm(tmp_delta) - pnorm(-tmp_delta)
  q0 <- pmax(tmp_exp * Phi_delta, pmin(tmp_exp * cosh(tmp_m * tmp_delta) *
                                         Phi_delta, q0))
  q1 <- s * q0 + (exp(a  - u^2 / 2) - exp(a - v^2 / 2)) / gam / sqrt(2 * pi)
  # q1 <- s * q0 + exp_a / gam * (dnorm(u) - dnorm(v))
  q2 <- 2 * s * q1 - s^2 * q0 + (u * exp(a - u^2 / 2) - v * exp(a - v^2 / 2)) / gam^2 / sqrt(2 * pi)
  # q2 <- 2 * s * q1 - s^2 * q0 + exp_a / gam^2 * (u * dnorm(u) - v * dnorm(v))

  # Derivatives of density smooth density
  f0 <- sum(f * q0)
  f1 <- sum(f * q1)
  f2 <- sum(f * q2)
  # Derivatives of the smooth re-weighted density (without) normalization
  phi_h <- dnorm(y, sd = h)
  f0_h <- f0 / phi_h
  f1_h <- (f1 + y / h^2 * f0) / phi_h
  f2_h <- (f2 + 2 * y / h^2 * f1 + (1 / h^2 + y^2 / h^4) * f0) / phi_h
  return(c(f0 = f0, f1 = f1, f2 = f2, f0_h = f0_h, f1_h = f1_h, f2_h = f2_h))
}

#' Get Smoothing Bandwidth
#'
#' Returns the bandwith for the smoothed density, such that that its variance
#' co-incides with the variance of the empirical distribution.
#'
#' @param logcon Object of class `dlc` from package \code{link[logcondens]{LogConDens}}.
#' @return Smoothing Bandwidth
#' @export
GetSmoothingBandwidth <- function(dlc){
  VarFn <- logcondens::LocalVariance(x = dlc$x, w = dlc$w, phi = dlc$phi)
  gam <- sqrt(dlc$sig ^ 2 - VarFn)
  return(gam)
}

#' Newton step for mode of smoothed (re-weighted) density
#'
#' Newton step with step-size correction, if \code{f2 < 0}. Ohterwise, goes to
#' direction of the next knot.
#'
#' @param y Starting point
#' @param lcd Fitted log-concave density of class \code{dlc} from package
#'   \code{logcondens}.
#' @param h Re-weighting bandwidth
#' @param reweighted Which function to use. The re-weighted one or not.
#' @param knotes Knotes of the log-concave density estimator
#' @export
GetNewtonStep <- function(y, lcd, h, reweighted = TRUE, knots, delta = 1e-10 / 2){
  if (reweighted){
    f0 <- 'f0_h'
    f1 <- 'f1_h'
    f2 <- 'f2_h'
  } else{
    f0 <- 'f0'
    f1 <- 'f1'
    f2 <- 'f2'
  }
  f <- GetDerivativeSRLC(y, lcd, h)
  if (f[f2] < 0){
    # Newton step
    y_new <- y - f[f1] / f[f2]
    # parameter for step size correction
    lambda <- abs(f[f1]^2 / f[f2]) / 3
  } else{
    if (y <= min(knots) | y >= max(knots)){
      warning(print('Current point not between knots with positive 2nd derivative.\nThis error occured in GetNewtonStep'))
      return(NA)
    }
    knot_tmp <- sort(c(knots[abs(knots - y) > delta], y))
    index <- which(knot_tmp == y)
    if (f[f1] < 0){
      y_new <- knot_tmp[index - 1]
    } else{
      y_new <- knot_tmp[index + 1]
    }
    lambda <- f[f1] * (y_new - y)
  }
  f0_new <- GetDerivativeSRLC(y_new, lcd, h)[f0]
  cond <- f0_new - f[f0] < lambda & abs(y - y_new) > delta
  if (is.na(cond)) return(NA)
  while (cond){
    lambda <- lambda / 2
    y_new  <- (y + y_new) / 2
    f0_new <- GetDerivativeSRLC(y_new, lcd, h)[f0]
    cond <- f0_new - f[f0] < lambda & abs(y - y_new) > delta
    if (is.na(cond)) return(NA)
  }
  names(y_new) <- NULL
  return(y_new)
}

#' Newton Method finding closest mode
#'
#' A Newton method finding the closest mode to zero of the smooothed
#' (re-weighted) density.
#'
#' @param y Starting point (default: the mode of the non-smooth estimator)
#' @param lcd Fitted log-concave denstiy of class \code{dlc} from package
#'   \code{logcondens}.
#' @param h Re-weighting bandwidth
#' @param delta accuracy in the Newton method.
#' @return mode of the smooth (re-weigted) density.
#' @export
FindModeNewton <- function(y = NULL, lcd, h, reweighted = TRUE, delta = 1e-10, max_iter = 20){
  if (is.null(y)){
    if (reweighted){
      y <- FindReweightedMode(lcd$x, lcd$w, h)$m
    } else{
    y <- lcd$mode
    }
  }
  y_new <- GetNewtonStep(y, lcd, h, reweighted, lcd$knots)
  j <- 1
  cond <- abs(y - y_new) > delta
  if (is.na(cond)) return(NA)
  while (cond & j <= max_iter){
    y <- y_new
    y_new <- GetNewtonStep(y, lcd, h, reweighted, lcd$knots)
    j <- j + 1
    cond <- abs(y - y_new) > delta
    if (is.na(cond)) return(NA)
  }
  if (j >= max_iter){
    warnings(print(paste('Newton stopped after', max_iter, 'steps at', y_new)))
  }
  return(y_new)
}

#' Get Confidence Interval
#'
#' Calculates the confidence interval of the mode for the weighted projected
#' data in direction of minimal variance.
#'
#' The points in \code{x} should correspond to the estimated ridgepoints of the
#' loc-concave method. Furthermore, \code{h} and \code{alpha} should be the
#' bandwidth which has been used for estimating \code{x}. The returned values
#' give the (\code{1 - level})-confidence interval for each point in \code{x}.
#'
#' Whenever \code{x} is not a vector, a for-loop is used. To speed up things,
#' one may use vectors for \code{x} with some multi-threat method outside the
#' function.
#'
#' @param data Data as dataframe/tibble [d, n].
#' @param x Either a vector [d] or a tibble [m, d] with column-names \code{X1,
#'   X2} corresponding to m different points
#' @param h Bandwidth (sould correspond to the one \code{x} was estimated with).
#' @param level The confidence level of each confidence interval is \code{1 -
#'   level}.
#' @param alpha Cutoff value of the weight vector in percent.
#' @return A named vector or a tibble with the columns \code{X1, X2, X1_left,
#'   X2_left, X1_right, X2_right}.
#' @export
GetCI <- function(data, x, h, level = 0.05, alpha = 1){
  GetCIone <- function(data, x, h, level, alpha){
    w <- Weights(data, x, h)
    if (alpha < 1){
      tmp <- NewWeights(w, alpha)
      w <- tmp$w[tmp$isPositive]
      data <- data[ , tmp$isPositive]
    }
    v <- FindDirection(data, x, w)$v
    z <- ProjectData(data, x, v)
    index <- order(z)
    CI_univar <- LCLRCImode(z[index], w = w[index], alpha = level)
    return(c(X1 = x[1], X2 = x[2],
           X1_left = x[1] + CI_univar[1] * v[1],
           X2_left = x[2] + CI_univar[1] * v[2],
           X1_right = x[1] + CI_univar[2] * v[1],
           X2_right = x[2] + CI_univar[2] * v[2]))
  }
  if (is.vector(x)){
    return(GetCIone(data, x, h, level, alpha))
  } else {
    m <- nrow(x)
    res <- tibble(X1 = rep(NA, m), X2 = rep(NA, m),
                  X1_left = rep(NA, m), X2_left = rep(NA, m),
                  X1_right = rep(NA, m), X2_right = rep(NA, m))
    for (i in 1:m){
      res[i, ] <- GetCIone(data, x[i, c('X1', 'X2')] %>% unlist(), h, level, alpha)
    }
    return(res)
  }
}

#' Post Smoothing
#'
#' After calculating the "lc"-algorithm, it performce a mode smoothing after
#' convergence.
#'
#' @param x Matrix of the points to consider [d, n]
#' @export
PostSmoothing <- function(data, x, h, alpha = 1){
  x <- as.matrix(x, nrow = 2)
  n <- dim(x)[2]
  smoothMode <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n){
    w <- Weights(data, x[ , i], h)
    tmp <- NewWeights(w, alpha)
    w <- tmp$w
    isPos <- tmp$isPositive
    v <- FindDirection(data, x[ , i], w)$v
    z <- ProjectData(data, x[ , i], v)
    dlc <- activeSetLogCon(z[isPos], w = w[isPos])
    gam <- GetSmoothingBandwidth(dlc)
    mode <- FindSmoothedMode(dlc$x, dlc$phi, gam)
    smoothMode[i, ] <- x[ , i] + mode * v
  }
  return(data.frame(X1 = smoothMode[ , 1],
                    X2 = smoothMode[ , 2]))
}

#' Post Interval
#'
#' Calculates for each ridge point calculated by "lc" an interval, where
#' the mode could also lie, based on a threshold.
#'
#' @param x Matrix of the points to consider [d, n]
#' @param thres percent value according to the maximum of the density
#' @export
PostInterval <- function(data, x, h, thres = 0.9, alpha = 1, print = FALSE){
  x <- as.matrix(x, nrow = 2)
  n <- dim(x)[2]
  interval <- matrix(NA, nrow = n, ncol = 4)
  for (i in 1:n){
    w <- Weights(data, x[ , i], h)
    tmp <- NewWeights(w, alpha)
    w <- tmp$w
    isPos <- tmp$isPositive
    v <- FindDirection(data, x[ , i], w)$v
    z <- ProjectData(data, x[ , i], v)
    dlc <- activeSetLogCon(z[isPos], w = w[isPos])
    IsKnot <- as.logical(dlc$IsKnot)
    knots <- dlc$x[IsKnot]
    phi <- dlc$phi[IsKnot]
    max <- max(phi)
    arg_max <- which.max(phi)
    y <- max + log(thres)
    ind1 <- ifelse(phi[1] < max + log(thres),
                   max(which(phi[1:arg_max] < max + log(thres))),
                   0)
    low <- ifelse(ind1 > 0,
                  (y - phi[ind1]) * (knots[ind1 + 1] - knots[ind1]) /
                    (phi[ind1 + 1] - phi[ind1]) + knots[ind1],
                  knots[1])
    ind2 <- ifelse(phi[length(phi)] < max + log(thres),
                   max(which(phi > max + log(thres))),
                   0)
    up  <- ifelse(ind2 > 0,
                  (y - phi[ind2]) * (knots[ind2 + 1] - knots[ind2]) /
                    (phi[ind2 + 1] - phi[ind2]) + knots[ind2],
                  knots[length(knots)])
    interval[i, 1:2] <- x[ , i] + v * low
    interval[i, 3:4] <- x[ , i] + v * up
    if (print){
      p <- ggplot(tibble(knots = knots, phi = phi)) +
        geom_line(aes(knots, phi)) +
        geom_vline(xintercept = c(low, up), col = "red") +
        geom_hline(yintercept = max + log(thres))
      print(p)
    }
  }
  return(data.frame(X1_start = interval[ , 1],
                    X2_start = interval[ , 2],
                    X1_end   = interval[ , 3],
                    X2_end   = interval[ , 4]))
}

#' Get Grid
#'
#' Calculates a starting point grid with some given maximal distance to the data.
#'
#' @param data data as matrix [d, n]
#' @param d distance between grid points
#' @param r maximal distance of each grid point from data. (default r = d)
#' @return matrix [d, n] containing one grid point per column.
#' @export
GetGrid <- function(data, d, r = NULL){
  if (is.null(r)) r <- d
  xx <- seq(min(data[1, ]) - r, max(data[1, ]) + r, by = d)
  yy <- seq(min(data[2, ]) - r, max(data[2, ]) + r, by = d)
  grid <- rbind(rep(xx, length(yy)), rep(yy, each = length(xx)))
  accept <- apply(grid, 2, function(y) any(colSums((data - y)^2) < r^2))
  return(grid[ , accept])
}
