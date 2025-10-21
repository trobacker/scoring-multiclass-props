## functions for scoring exploration

## define a density function to use dirichlet values for colors
dd_func <- function(a, b, c, alpha)
  brms::ddirichlet(cbind(a, b, c), alpha)

## define a density function for implied multinomial
mn_func <- function(a, b, c, size, theta)
  mc2d::dmultinomial(cbind(a*size, b*size, c*size), size = size, prob = theta)

## define a function to draw one dirichlet sample and store as a vector
draw_one_dirichlet <- function(n, alpha)
  as.vector(brms::rdirichlet(n, alpha = alpha))

## make colors transparent, from https://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
col_r <- function(col, alpha = NULL) {
  MAX_COL <- 255
  col_mtx <- t(col2rgb(col, alpha = FALSE))

  if (is.null(alpha)) {
    rgb(red = col_mtx, maxColorValue = MAX_COL)
  } else {
    col_alf <- MAX_COL * alpha
    rgb(red = col_mtx, alpha = col_alf, maxColorValue = MAX_COL)
  }
}

## define a function to enumerate the sample space of a multinomial
enumerate_multinomial_sample_space <- function(n){
  ## fix number of categories at 3
  m <- 3
  nrows <- choose(n+m-1, m-1)
  dat <- matrix(nrow = nrows, ncol=3)
  dat[,1] <- rep(0:n, times = (n+1):1)
  first_zeroes_idx <- c(1, 1 + cumsum(1 + n:1))
  for(i in 0:n){
    ## indices of each i in column 2
    idx <- first_zeroes_idx + i
    ## remove elements of idx that are greater than nrows
    idx <- idx[1:(n-i+1)]
    dat[idx,2] <- i
  }
  dat[,3] <- n - dat[,1] - dat[,2]
  ## checks
  if(any(dat[,1] + dat[,2] + dat[,3] != n))
    stop("Error in enumeration, at least one row != n")
  if(any(dat[,1] < 0) || any(dat[,2] < 0) || any(dat[,3] < 0))
    stop("Error in enumeration, one row has a value < 0")
  if(any(dat[,1] > n) || any(dat[,2] > n) || any(dat[,3] > n))
    stop("Error in enumeration, one row has a value > n")
  return(dat)
}

## define a function to plot the implied multinomial distribution
plot_implied_multinomial <- function(n, theta, c=NULL, opacity = 0.5){
  ## n the total size of the multinomial
  ## theta vector of probabilities of each category
  ## c optional, observed counts
  if(!is.null(c))
    if(sum(c) != n)
      stop("Error: sum of observed counts != n")

  ## plot base graph
  TernaryPlot(alab = "Variant A count \u2192",
              blab = "\u2190 Variant B count ",
              clab = "Variant C count \u2192",
              region = Ternary:::ternRegionDefault/(100/n),
              grid.lines = 5,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  ## enuermate the multinomial sample space
  dat <- enumerate_multinomial_sample_space(n)

  # Cut our reflectance data into categories
  # densities <- cut(mc2d::dmultinomial(dat, prob = theta), 10)
  # # Define a colour spectrum
  # mySpectrum <- rowCol <- rev(hcl.colors(10, palette = "viridis"))
  # Assign each data point a colour from the spectrum
  # pointCol <- mySpectrum[densities]

  # Define a size range
  maxSize <- 3 # Size of largest point, in plotting units
  sizes <- mc2d::dmultinomial(dat, prob = theta)
  pointSize <- sizes * maxSize / max(sizes) + 0.0001

  ## add sample space points to graph
  AddToTernary(graphics::points,
               dat,
               pch=20,
               cex = pointSize,
               col = col_r("black", alpha = opacity))

  ## add observation to graph
  if(!is.null(c)) {
    AddToTernary(graphics::points,
                 c,
                 pch=4,
                 col="red", lwd=2)
  }
}

## define a function to plot the implied multinomial mixture
plot_implied_multinomial_mixture <- function(n, thetas, c=NULL, opacity = 0.5){
  ## n the total size of the multinomial
  ## theta vector of probabilities of each category
  ## c optional, observed counts
  if(!is.null(c))
    if(sum(c) != n)
      stop("Error: sum of observed counts != n")

  ## plot base graph
  TernaryPlot(alab = "Variant A count \u2192",
              blab = "\u2190 Variant B count ",
              clab = "Variant C count \u2192",
              main = "Mixture of Implied Multinomials",
              region = Ternary:::ternRegionDefault/(100/n),
              grid.lines = 5,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  ## enuermate the multinomial sample space
  dat <- enumerate_multinomial_sample_space(n)

  # Cut our reflectance data into categories
  # densities <- cut(mc2d::dmultinomial(dat, prob = theta), 10)

  # # Define a colour spectrum
  # mySpectrum <- rowCol <- rev(hcl.colors(10, palette = "viridis"))
  # Assign each data point a colour from the spectrum
  # pointCol <- mySpectrum[densities]

  ## determine mixture probabilities
  sizes <- sapply(FUN = function(x) mc2d::dmultinomial(dat, size = n, prob = x),
                  X = thetas) |>
    rowMeans()

  # Define a size range
  maxSize <- 3 # Size of largest point, in plotting units
  pointSize <- sizes * maxSize / max(sizes) + 0.0001 ## add small number because 0 pointsize plots something

  ## add sample space points to graph
  AddToTernary(graphics::points,
               dat,
               pch=20,
               cex = pointSize,
               col = col_r("black", alpha = opacity))

  ## add observation to graph
  if(!is.null(c)) {
    AddToTernary(graphics::points,
                 c,
                 pch=4,
                 col="red", lwd=2)
  }
}


## define a function to plot a dirichlet multinomial distribution
plot_dirichlet_multinomial <- function(n, alpha, c=NULL, opacity = 0.5){
  require(extraDistr)
  ## n the total size of the multinomial
  ## c optional, observed counts
  if(!is.null(c))
    if(sum(c) != n)
      stop("Error: sum of observed counts != n")

  ## plot base graph
  TernaryPlot(alab = "Variant A count \u2192",
              blab = "\u2190 Variant B count ",
              clab = "Variant C count \u2192",
              main = "Dirichlet-Multinomial",
              region = Ternary:::ternRegionDefault/(100/n),
              grid.lines = 5,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  ## enuermate the multinomial sample space
  dat <- enumerate_multinomial_sample_space(n)

  # Cut our reflectance data into categories
  # densities <- cut(mc2d::dmultinomial(dat, prob = theta), 10)
  # # Define a colour spectrum
  # mySpectrum <- rowCol <- rev(hcl.colors(10, palette = "viridis"))
  # Assign each data point a colour from the spectrum
  # pointCol <- mySpectrum[densities]

  # Define a size range
  maxSize <- 3 # Size of largest point, in plotting units
  sizes <- ddirmnom(dat, size = n, alpha = alpha)
  pointSize <- sizes * maxSize / max(sizes) + 0.0001

  ## add sample space points to graph
  AddToTernary(graphics::points,
               dat,
               pch=20,
               cex = pointSize,
               col = col_r("black", alpha = opacity))

  ## add observation to graph
  if(!is.null(c)) {
    AddToTernary(graphics::points,
                 c,
                 pch=4,
                 col="red", lwd=2)
  }
}
