#' Gamma distributed on circle
#'
#' @param n Sample size
#' @param seed Change the seed for RNG
#' @export
GammaCircle <- function(n = 200, shape = 5, rate = 1, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  phi <- runif(n, 0, 2 * pi)
  r   <- rgamma(n, shape, rate)
  return(rbind(
    X1 = cos(phi) * r,
    X2 = sin(phi) * r
  ))
}

#' Gamma on a line
#'
#' 2d-data with X1 uniformlly distributed on [-5, 5] and X2 gamma distributed
#' with certain shape and rate. Note that the true ridge is the line with X2 =
#' (shape - 1) / rate
#'
#' @param n Sample size
#' @param seed Change the seed for RNG
#' @export
GammaLine <- function(n, shape = 2, rate = 4, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  return(rbind(
    X1 = runif(n, -5, 5),
    X2 = rgamma(n, shape, rate)
  ))
}

#' Gamma tiple line
#'
#' The true ridgeset of the data is approximatly at X2 = 1/4, 5/4, 4 + 1/4.
#'
#' @param n Sample size
#' @param seed Change the seed for RNG
#' @export
GammaTripleLine <- function(n, shape = 2, rate = 4, distance = c(0, 1, 4), seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  which <- sample(distance, n, replace = TRUE)
  return(rbind(
    X1 = runif(n, -10, 10),
    X2 = rgamma(n, shape, rate) + which
  ))
}

#' Simulate Sinus data
#'
#' Generates data distributed along a sinus line with Gaussian noise having
#' standard deviation \code{sd}.
#'
#' @param n Sample size
#' @param sd Standard deviation
#' @param seed Change the seed for RNG
#' @return Matrix with rows \code{X1, X2}. [2, n]
#' @export
SimSinus <- function(n, int = c(-2 * pi, 2 * pi), sd = 0.2, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  u <- runif(n, int[1], int[2])
  return(rbind(
    X1 = u + rnorm(n, 0, sd),
    X2 = sin(u) + rnorm(n, 0, sd)
  ))
}

#' Create Tibble
#'
#' Creates a Tibble of the estimated ridge points
#'
#' The Tibble has the following columns:
#' X1_start, X2_start, X1, X2, h, n, method, data, seed, date, time
#'
#' @param grid Matrix of the starting points [d, n]
#' @param result Matrix of the limits [d, n]
#' @param h Used Bandwith
#' @param method Used method
#' @param data_name Name of the used data generation
#' @param seed Used seed for the data generation
#' @param alpha Used alpha for weight vector
#' @export
GetTibble <- function(grid, result, h, n, method, data_name, seed, alpha){
    tib <- tibble(
    X1_start = grid[1, ],
    X2_start = grid[2, ],
    X1 = result[1, ],
    X2 = result[2, ],
    h  = h,
    n  = n,
    method = method,
    data = data_name,
    seed = seed,
    date = format(Sys.Date(), '%Y-%m-%d'),
    time = format(Sys.time(), '%H:%M:%S'),
    alpha = alpha
  )
  return(tib)
}

#' Save tibble
#'
#' Save the tibble under given name
#'
#' @param tibble Tibble to save
#' @param name Name of the file
#' @param add If \code{TRUE}, the data is added to the file, oterhwise the file
#'   is overwritten.
#' @export
SaveTibble <- function(tib, name, add = TRUE){
  write_csv(tib, name, append = add, col_names = !add)
}

#' Load Tibble
#'
#' Loads the tibble from file
#'
#'@param file Name of the file
#'@export
LoadTibble <- function(file_name){
  tib <- read_csv(file_name,
                 col_types = cols(
                    X1_start = col_double(),
                    X2_start = col_double(),
                    X1 = col_double(),
                    X2 = col_double(),
                    h = col_double(),
                    n = col_integer(),
                    method = col_factor(levels = c('lc', 'slc', 'rlc', 'srlc', 'scms')),
                    data = col_character(),
                    seed = col_integer(),
                    date = col_date(format = '%Y-%m-%d'),
                    time = col_time(format = '%H:%M:%S'),
                    alpha = col_double()
                    )
                  )
  return(tib)
}
