#' Weighted Density Plot
#'
#' Plots of the estimated densities for different methods.
#'
#' @param z Univariate data [n].
#' @param w Weights of the data [n].
#' @param h Bandwith needed for the smoothed densities. (Default = 1)
#' @param which A (sub)vector of \code{1:5}. Each number
#'   corresponding to a certain estimator: 1: log-concave, 2: smoothed
#'   log-concave, 3: reweigted log-concave, 4: reweigted smoothed log-concave, 5: KDE.
#' @return ggplot2 object showing the desired densities.
#' @export
PlotDensity <- function(z, w, h = 1, which = c(1, 5)){
  bound <- 0.1 * (max(z) - min(z))
  xx <- seq(min(z) - bound, max(z) + bound, length.out = 500)
  x_weight <- (max(xx) - min(xx)) / 500
  h_kde <- bw.nrd0(z)
  res <- activeSetLogCon(z, w = w)
  res2 <- evaluateLogConDens(xx, res)
  rlc_raw <- exp(res2[ , 'log-density'] + xx^2 / (2 * h^2))
  srlc_raw <- res2[ , 'smooth.density'] * exp(xx^2 / (2 * h^2))
  data_raw <- tibble(
    xx   = xx,
    lc   = res2[ , 'density'],
    slc  = res2[ , 'smooth.density'],
    rlc  = rlc_raw / ( sum(rlc_raw) * x_weight),
    # rslc = exp(log(res2[ , 'smooth.density'] + xx^2 / h^2)),
    srlc = srlc_raw / (sum(srlc_raw) * x_weight),
    kde  = rowMeans(exp(-outer(xx, z, '-')^2 / h_kde^2 / 2)) / sqrt(2 * pi) / h_kde
  )
  data <- data_raw %>%
    select(1, which + 1) %>%
    gather('Estimator', 'y', -xx)
  p <- ggplot(data) + geom_line(aes(xx, y, color = Estimator))
  return(p)
}

#' Weighted and Log-Concave Density Plot
#'
#' @param z Projected Data [n].
#' @param w Weights of the projected data [n].
#' @param x Location vector [d]. If \code{x} and \code{v} are \code{NULL} the
#'   information will be plotted in the subtitle.
#' @param v Direction vector [d]. If \code{x} and \code{v} are \code{NULL} the
#'   information will be plotted in the subtitle.
#' @export
PlotDensityBoth <- function(z, w, x = NULL, v = NULL){
  res <- activeSetLogCon(z, w = w)
  res2 <- data.frame(x = res$x, y = exp(res$phi), w = w)
  p <- ggplot() +
    geom_density(aes(z, weight = w), data = data.frame(z = z, w = w)) +
    geom_line(aes(x, y), data = res2, color = "red") +
    geom_vline(xintercept = res$mode, color = "red")
  if (!is.null(v) & !is.null(x)){
    p <- p + ggtitle("Weighted Kernel and Log-Concave Density Estimation",
                     paste0("at x = (", round(x[1], 4), ", ", round(x[2], 4),
                            ") in direction v = (",
                            round(v[1], 4), ", ", round(v[2], 4), ")."))
  }
  else{
    p <- p + ggtitle("Weighted Kernel and Log-Concave Density Estimation")
  }
  return(p)
}

#' Plot reweighted Density
#'
#' Plots the univariate log-concave density reweighted with the Gaussian kernel.
#'
#' @param z Projected data [n]
#' @param w Weights [n]
#' @param h Bandwidth
#'
#' @export
PlotReweightedDensity <- function(z, w, h){
  tmp <- FindReweightedMode(z, w, h)
  logcon <- tmp$logcon
  res2 <- data.frame(x = logcon$x, y = exp(logcon$phi + logcon$x^2 / h^2), w = w)
  p <- ggplot() +
    geom_line(aes(x, y), data = res2) +
    geom_vline(xintercept = tmp$m, color = "red")
  # Add smoothed reweigthed density
  smooth_f <- evaluateLogConDens(logcon$x, logcon)[ , 4]
  p <- p + geom_line(aes(x, y),
                     data.frame(x = logcon$x, y = smooth_f),
                     color = 'blue')
  return(p)
}

#' Plot smoothed density
#'
#' @param z Projected data [n]
#' @param w Weights [n]
#' @param gam Smoothing bandwidth
#'
#' @export
PlotSmoothedDensity <- function(z, w){
  res <- activeSetLogCon(z, w = w)
  xs <- seq(min(z), max(z), length.out = 500)
  res2 <- data.frame(evaluateLogConDens(xs, res))
  p <- ggplot(res2, aes(x = xs)) +
    geom_line(aes(y = smooth.density)) +
    ggtitle('Smoothed Density Estimator')
  return(p)
}

#' Arrow Plot
#'
#' Plots an arrow between each starting and endpoint
#'
#' @param grid Matrix of starting points [d, n].
#' @param ridge Matrix of end points [d, n].
#' @export
PlotArrows <- function(grid, ridge) {
  df <- data.frame(x = grid[1, ], y = grid[2, ], xend = ridge[1, ], yend = ridge[2, ])
  p <- PlotData(data) + geom_segment(aes(x, y, xend = xend, yend = yend),
                             color = "red", data = df, arrow = arrow(length = unit(0.2, "cm")))
  return(p)
}

#' Arrow geom
#'
#' Plots an arrow between each starting and endpoint
#'
#' @param tib A tibble containing X1_start, X2_start, X1, X2
#' @export
geom_Arrows <- function(tib, ...) {
  geom_segment(aes(X1_start, X2_start, xend = X1, yend = X2),
               data = tib, arrow = arrow(length = unit(0.2, "cm")),
               ...)
}


#' True Ridge of Circle Data
#'
#' Adds the true Ridge of Circle data to the plot
#'
#' @export
geom_CircleRidge <- function(r, sigma, ...){
  radius_ridge <- circle_radius(r, sigma)
  if(radius_ridge == 0) warning("Ridge is equal (0, 0)")
  geom_path(aes(x, y), data = tibble(
           x = radius_ridge * cos(seq(0, 2 * pi, length.out = 100)),
           y = radius_ridge * sin(seq(0, 2 * pi, length.out = 100))),
           ...)
}

#' Plot post intervals
#'
#' @param interval Tibble with the variables: X1_start,
#' X2_start, X1_end, X2_end
#' @export
geom_PostInt <- function(interval, ...){
  geom_segment(aes(x = X1_start,
                   y = X2_start,
                   xend = X1_end,
                   yend = X2_end), data = interval, ...)
}

#' Plots the log-concave density of a given step
#'
#' @param list_element Contains one list element generated by the "save"
#' option in "LogConRidge".
#'
#' @export
Plot1Step <- function(list_element){
  p <- PlotDensityBoth(list_element$lcd$x, list_element$lcd$w,
                       list_element$x_old, list_element$v)
  return(p)
}

#' Plot Step Points
#'
#' @param steps List of steps generated by the "save" option in "LogConRidge".
#'
#' @export
geom_StepPoints <- function(steps, ...){
  tib <- GetPath(steps)
  geom_point(aes(x = X1_start, y = X2_start), data = tib, ...)
}

#' Plot step Path
#'
#' @param steps List of steps generated by the "save" option in "LogConRidge".
#'
#' @export
geom_StepPath <- function(steps, ...){
  tib <- GetPath(steps)
  geom_segment(aes(x = X1_start, y = X2_start, xend = X1_end, yend = X2_end),
               data = tib, ...)
}

#' Get Path tibble
#'
#' @param steps List of steps generated by the "save" option in "LogConRidge".
#' @return Tibble with columns X1_start, X2_start, X1_end, X2_end and each row
#' corresponds to one step in the LCRS.
#'
#' @export
GetPath <- function(steps){
  tib <- tibble(X1_start = numeric(),
                X2_start = numeric(),
                X1_end   = numeric(),
                X2_end   = numeric())
  for (step in steps){
    tib <- tib %>% add_row(X1_start = step$x_old[1],
                           X2_start = step$x_old[2],
                           X1_end   = step$x_new[1],
                           X2_end   = step$x_new[2])
  }
  return(tib)
}
