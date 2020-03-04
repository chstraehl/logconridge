#' Log-Concave Ridge search
#'
#' Finds the closest Ridge Point using log-concave variante.
#'
#' @param data Data as  matrix [d, n].
#' @param x Starting location as vector [d].
#' @param h Bandwidth for the weight function
#' @param smoothed Logical, if \code{TRUE} the mode of the smoothed log-concave
#'   density estimator will be used. (default: \code{FALSE})
#' @param delta Stopping criteria (default: \code{delta = 1e-6}).
#' @param print If TRUE prints the values at each step (default \code{print =
#'   FALSE}).
#' @param save If TRUE returns a list, which contains for each step a list with
#' startpoint, endpoint, direction, log-concave density.
#' @param max_iter The maximal amount of iterations. If the maximal iterations
#'   are reached a warning will be printed.
#' @param alpha Cutoff value for the weight vector in percent.
#' @return Vector [d] the algorithm converges to on the ridge.
#' @examples
#' data <- SimCircle(200)
#' x <- c(-0.5, 0.1)
#' h <- 0.5
#' LogConRidge(data, x, h, print = TRUE)
#' @export
LogConRidge <- function(data, x, h, delta = 1e-6, print = FALSE, save = FALSE, max_iter = 20, reweight = TRUE, smooth = FALSE, alpha = 1){
  x_start = x
  j <- 1
  path <- x
  steps <- list()
  cond <- TRUE
  while(cond & (j <= max_iter)){
    w <- Weights(data, x, h)
    if (alpha < 1){
      tmp <- NewWeights(w, alpha)
      w   <- tmp$w
      isPositive <- tmp$isPositive
    } else{
      isPositive <- rep(TRUE, length(w))
    }
    v <- FindDirection(data[ , isPositive], x , w[isPositive])$v
    z <- ProjectData(data[ , isPositive], x, v)
    if (reweight){
      tmp <- FindReweightedMode(z, w[isPositive], h)
      m <- tmp$m
      if (smooth){
        m <- FindModeNewton(m, tmp$logcon, h)
      }
    } else{
      res <- activeSetLogCon(z, w = w[isPositive])
      m <- res$mode
      if (smooth){
        m <- FindModeNewton(m, res, h, reweighted = FALSE)
      }
    }
    x_old <- x
    x <- x + m * v
    if (print){
      path <- cbind(path, x)
      # Use the vector z from above
      #z <- colSums(v * (data - x_old))
      print(paste("Number of iterations:", j))
      print(paste("Direction:", v[1], v[2]))
      print(paste("Mode:", m))
      print(paste("New x:", x[1], x[2]))
      res <- activeSetLogCon(z, w = w[isPositive])
      res2 <- data.frame(x = res$x, y = exp(res$phi), w = w[isPositive])
      if (reweight){
        if (smooth){
          p <- PlotDensity(z, w[isPositive], h, 4) +
            geom_vline(xintercept = m)
        } else{
          p <- PlotReweightedDensity(z, w[isPositive], h)
        }

      } else{
        if (smooth){
          p <- PlotDensity(z, w[isPositive], h, 2) +
            geom_vline(xintercept = m)
        } else{
          p <- PlotDensityBoth(z, w[isPositive], x_old, v)
        }
      }
      print(p)
      readline(prompt = "PRESS ENTER!")
      if (j == 1){
        q <- PlotData(data) +
          geom_point(aes(x, y), data = data.frame(x = x_old[1], y = x_old[2]), color = "red")
      }
      q <- q +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "red",
                     data = data.frame(x = x_old[1], y = x_old[2], xend = x[1], yend = x[2])) +
        geom_point(aes(x, y), color = "red", data = data.frame(x = x[1], y = x[2]))
      print(q)
      readline(prompt = "PRESS ENTER!")
    }
    if (save){
      # Saves all steps into a list for printing
      lcd <- activeSetLogCon(z, w = w[isPositive])
      steps[[j]] <- list(x_old = x_old, x_new = x, v = v, lcd = lcd)
    }
  j <- j + 1
  cond <- abs(m) > delta
  if (is.na(cond)){
    warning(sprintf('NA is returned for startingpoint (%.3f, %.3f)', x_start[1], x_start[2]))
    return(c(NA, NA))
  }
  }
  if (j > max_iter) warning(paste("Stopped after", max_iter,
                                  "iterations for the starting point", x_start[1], x_start[2]))
  if (save){
    return(steps)
  } else{
    return(x)
  }
}

#' Subspace Constraint Mean Shift
#'
#' Implementation of the subspace constraint mean shift algorithm
#'
#' @param data Data as matrix [d, n]
#' @param x Starting point [d]
#' @param h Bandwidth of the kernel density estimator
#' @param delta Achivent accuracy
#' @param print If \code{TRUE}, the path of the algorithm is plotted.
#' @return Point on the ridge closest to \code{x} as vector
#' @export
SCMS <- function(data, x, h, delta = 1e-6, print = FALSE, iter_max = 100){
  if (print) path <- x
  x_start <- x
  m <- GetProjectedMeanShift(data, x, h)
  x <- x + m
  j <- 1
  while( sqrt(sum(m^2)) > delta & (j <= iter_max)){
    if (print) path <- cbind(path, x)
    m <- GetProjectedMeanShift(data, x, h)
    x <- x + m
    j <- j + 1
  }
  if (print){
    p <- PlotData(data) +
      geom_path(aes(x, y), data.frame(x = path[1, ], y = path[2, ]), color = 'red') +
      geom_point(aes(x, y), data.frame(x = path[1, ], y = path[2, ]), color = 'red')
    print(p)
    print(sprintf('Number of steps: %i', j))
  }
  if (j > iter_max) warnings(paste('Stopped after', iter_max, 'iterations at starting point',
                                   x_start[1], x_start[2]))
  return(x)
}


#' Get Bandwidth
#'
#' Calculation of the bandwidth based on minimum spanning trees
#'
#' @param data Datapoints
#' @export
GetBandwidth <- function(data){
  tmp <- ComputeMST(t(data))
  return((sum(tmp$distance) / dim(data)[2])^(1 / (4 + dim(data)[1])))
}
