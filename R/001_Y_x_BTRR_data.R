#' Bayesian Tensor Response Regression Data Generation
#'
#' @param tens_dim A scalar (e.g. 2 or 3) for the dimension of the generated response data for a given time and subject
#' @param subjects A scalar with the sample size of the generated data
#' @param n.time A scalar with the number of time steps for all subjects
#' @param margin_size A numeric vector of length *tens_dim* with the margin sizes of the response tensor
#' @param CNR The contrast-to-noise ratio, as defined by Welvaert, M. and Rosseel, Y., 2013
#' @param obs.var A value for the observational variance
#'
#' @return A list of length two with elements, *dat* and *truth*. The 'dat' object is a list with two elements, *Y* and *x*. *Y* is the tensor response in which the last two array margins represent the time step and the subject number, respectively. *x* is a single vector of length *n.time* that applies as the covariate in the tensor response regression. *truth* is the true coefficient tensor used in the data generation.
#' @export
#'
#' @examples
Y_x_BTRR_data <- function(tens_dim,subjects,n.time,margin_size,CNR,obs.var){
  require(neuRosim)

  margin_size[margin_size < 5] <- 5
  betas <- CNR*sqrt(obs.var)*specifyregion(dim = margin_size,
                                           coord = margin_size * runif(length(margin_size)),
                                           radius = min(round(max(margin_size)*0.1),sample(3:5,1)),
                                           form = "sphere") # Create a signal region
  x <- canonicalHRF(1:n.time,param=list(
    a1 = round(n.time*0.12), # Delay of response relative to onset
    a2 = round(n.time*0.5), # Delay of undershoot relative to onset
    b1 = 2, # Dispersion of response
    b2 = 1, # Dispersion of undershoot
    c = 0.5 # Scale of undershoot
  ))

  y <- sapply(seq(subjects), function(each_subject){
    array(unlist(sapply(x,function(each_x){
      bdim <- dim(betas)
      betas*each_x + array(rnorm(prod(bdim),0,sqrt(obs.var)),dim = bdim)
    },simplify = F)),dim = c(dim(betas),n.time))
  },simplify = "array")

  list(dat=list(Y=y,x=x),truth=list(betas=betas))

}
