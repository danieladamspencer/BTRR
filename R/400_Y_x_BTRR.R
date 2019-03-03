#' Generalized-Dimension Bayesian Tensor Response Regression with Gaussian Graphical Model
#'
#' @param input A list with at least two arguments, *Y* and *x*.The last two dimensions in the array *Y* should represent the time and subject, respectively. The first *D* dimensions of *Y* will represent the tensor dimensions at each time for each subject. *x* should either be a *n.time* by *subject* matrix, or a vector of length *n.time*.
#' @param n.iter The number of iterations that the MCMC will complete in total
#' @param n.burn The number of MCMC iterations that will be discarded as a burn-in
#' @param ranks The number of ranks for the model. See Guhaniyogi and Spencer (2018) for the full details.
#' @param hyperparameters A data frame with one row and six variables, with variables a.lambda, b.lambda, a.tau, b.tau, a.sig, and b.sig for the prior specifications. Default values are 3, 3^(1/(2*D)), D - 1, ranks^((1 / D) - 1), 1, and -log(0.95), respectively.
#' @param save_after The number of iterations after which to have the intermediate results saved. This can be helpful for debugging. The default does not save anything befor the end of the MCMC.
#' @param save_llik Should the log-likelihood be saved? This defaults to `TRUE`.
#'
#' @return A list of length 9 containing posterior draws for $\beta_{j,r}$, $W$, $\lambda_{j,r}$, $\Phi$, $\sigma_y^2$, and $\alpha$. The log-likeihood and acceptance ratio for the Metropolis-Hastings step are also included.
#' @export
#'
#' @examples
#' BTRR_data <- Y_x_BTRR_data(3, 20, 100, c(30,30,30), 1, 1)
#' BTRR_results <- BTR_Y_x(input = temp$dat, 10, 0, 1)
#'
BTR_Y_x <- function(input, n.iter, n.burn, ranks,
                             hyperparameters = NULL, save_after = NULL, save_llik = TRUE) {
  require(tidyverse)
  #  > Read in the data appropriately ----
  lowest_dim_margin <- min(head(dim(input$Y),-2))
  n <- tail(dim(input$Y),1)
  p <- head(dim(input$Y),-2)
  D <- length(dim(input$Y)) - 2

  if(lowest_dim_margin < ranks) stop("The rank of your model cannot be larger than your smallest tensor dimension. Try a lower model rank.")
  # x <- input$x # Stored as a vector of length T or as a matrix with dimensions T x n
  if(is.vector(input$x)) input$x <- tcrossprod(input$x,rep(1,n))
  TT <- dim(input$x)[1] # Number of time steps

  # > Load necessary packages ----
  require(GIGrvg)
  require(abind)
  require(dlm)
  require(statmod)
  require(truncnorm)
  require(Rcpp)
  require(RcppArmadillo)

  # > Functions ----
  source("R/900_misc.R")

  # > Standardize the response and center the covariate ----
  # input$Y <- (input$Y - mean(input$Y)) / sd(input$Y)
  # input$x <- input$x - mean(input$x)
  # This standardization standardizes by voxel and subject
  input$Y <- apply(input$Y,c(seq(D),D+2),function(y){
    if(sd(y) > 0) {(y - mean(y)) / sd(y)}else{rep(0,length(y))}
  }) %>%
    aperm(c(seq(D)+1,1,D+2))


  # > Hyperparameters ----
  if(is.null(hyperparameters)) {
    hyperparameters <- data.frame(
      a.tau = D - 1,
      b.tau = ranks^((1 / D) - 1), # These values are from Guhaniyogi et al [2017]
      a.lambda = 3,
      b.lambda = 3^(1/(2*D)), # These values are from Guhaniyogi et al [2017]
      a.sig = 1,
      b.sig = -log(0.95) # These values are from Guhaniyogi et al [2017]
    )
  }

  a.tau <- hyperparameters$a.tau
  b.tau <- hyperparameters$b.tau
  a.lambda <- hyperparameters$a.lambda
  b.lambda <- hyperparameters$b.lambda
  a.sig <- hyperparameters$a.sig
  b.sig <- hyperparameters$b.sig

  # > Set Storage ----
  results <-
    list(
      B = list(),
      W = list(),
      lambda = list(),
      Phi = list(),
      tau = vector(mode = "numeric",length = n.iter),
      llik = vector(mode = "numeric",length = n.iter),
      sig2y = vector(mode = "numeric",length = n.iter),
      alpha = vector(mode = "numeric",length = n.iter)
    )

  # > Set Initials ----
  # B <- # This is the MLE for voxel-by-voxel regression
  #   apply(input$Y, seq(D), function(each_voxel) {
  #     c(RcppArmadillo::fastLmPure(as.matrix(c(input$x)),c(each_voxel))$coefficients)
  #   })

  ## The initial values are set at random
  betas <-
    sapply(seq(D),function(each_dim){
      out <- matrix(rnorm(p[each_dim]*ranks),p[each_dim],ranks)
      apply(out,2,function(each_rank){each_rank / norm(as.matrix(each_rank))})
    },simplify = FALSE)

  tau <- 1 # Assuming unit value
  alpha.g <-
    seq(ranks ^ (-2), ranks ^ (-.1), length.out = 10) # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
  lambda <- matrix(1, ranks, 2)
  W <-
    sapply(betas, function(each_dim) {
      array(1, dim = dim(each_dim))
    }, simplify = FALSE)
  phi <- rep(1 / ranks, ranks)

  sig2y <- 1 # Observational variance estimate
  if(ranks > 1){
    Xi <- rep(.6,ranks - 1)
    # if(!is.matrix(Xi)) Xi <- matrix(Xi,nrow = 1)
    cov_Metro <- 0.01 * diag(ranks - 1)
  }else{
    Xi <- 1
  }
accept <- 0


  # > Run MCMC ----
beginning_of_sampler <- proc.time()[3]
  # pb <- txtProgressBar(min = 1, max = n.iter, style = 3)
  for (s in 1:n.iter) {
    # >> Griddy-Gibbs ----
      # Set up to sample alpha IF RANKS > 1
      if(ranks > 1){
        M = 4 # Number of reference sets per grid value of alpha
        l.weights <- sapply(alpha.g, function(proposed) {
          bw <- mapply(function(b, w) {
            abind(b, w, along = 3)
          },
          b = betas,
          w = W,
          SIMPLIFY = F)
          # p <- sapply(bw, function(each_dim) {
          #   dim(each_dim)[1]
          # })
          chi <- sapply(bw, function(each_dim) {
            apply(each_dim, 2, function(each_position) {
              t(each_position[, 1]) %*% diag(1 / each_position[, 2]) %*% each_position[, 1]
            })
          })
          ### INCLUDE RANK 1 change:
          if(ranks == 1){
            chi <- sum(chi)
          }else{
            chi <- apply(chi, 1, sum)
          }

          ## Draw Phi proposals
          ##### Phi under a stick-breaking prior
          old_Xi_g <- Xi
          phi.l <- sapply(seq(M),function(m){
            new_Xi_g <- c(old_Xi_g + cov_Metro%*%rnorm(ranks - 1))
            while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
              new_Xi_g <- c(old_Xi_g + cov_Metro%*%rnorm(ranks - 1))
            }
            new_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
              stick_break_log_posterior(new_Xi_g,cr, betas,W,tau,proposed)
            }))
            old_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
              stick_break_log_posterior(old_Xi_g,cr, betas,W,tau,proposed)
            }))
            if(exp(new_post_dens - old_post_dens) > runif(1)) old_Xi_g <- new_Xi_g
            stick_values(old_Xi_g)
          })

          ## Draw tau proposals
          ### ANOTHER RANK 1 CHANGE
          if(ranks == 1){
            chi2 <- chi / phi.l
          }else{
            chi2 <- apply(phi.l, 2, function(each_proposal) {
              chi / each_proposal
            })
            chi2 <- colSums(chi2)
          }
          tau.l <- rgig(M, a.tau - ranks * sum(p)/2, chi2, 2 * b.tau)
          refs <- list(phi = phi.l, tau = tau.l)
          ## Evaluate the densities
          lik.mean.tensor <-
            composeParafac(betas) %o% input$x
          l.lik <-
            sum(dnorm(c(input$Y), c(lik.mean.tensor), sqrt(sig2y), log = T)) # Log-likelihood
          if(ranks == 1){
            l.bdens <- apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
              # Log prior density for all betas
              sapply(each_proposal[-1], function(each_rank_phi) {
                sum(unlist(sapply(bw, function(each_dim) {
                  apply(each_dim, 2, function(each_rank_bw) {
                    dnorm(
                      each_rank_bw[, 1],
                      0,
                      each_proposal[1] * each_rank_phi * each_rank_bw[, 2],
                      log = T
                    )
                  })
                })))
              })
            })
          }else{
            l.bdens <-
              colSums(apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
                # Log prior density for all betas
                sapply(each_proposal[-1], function(each_rank_phi) {
                  sum(unlist(sapply(bw, function(each_dim) {
                    apply(each_dim, 2, function(each_rank_bw) {
                      dnorm(
                        each_rank_bw[, 1],
                        0,
                        each_proposal[1] * each_rank_phi * each_rank_bw[, 2],
                        log = T
                      )
                    })
                  })))
                })
              }))
          }

          l.tau <-
            dgamma(refs$tau, a.tau, b.tau, log = T) # Log prior density for tau

          ### RANK 1 CHANGE:
          if(ranks == 1){
            l.phi <- sapply(refs$phi,function(each_proposal){
              lgamma(ranks * proposed) - ranks * lgamma(proposed) + sum((rep(proposed, ranks) - 1) * log(each_proposal))
            })
          }else{
            l.phi <-
              apply(refs$phi, 2, function(each_proposal) {
                lgamma(ranks * proposed) - ranks * lgamma(proposed) + sum((rep(proposed, ranks) - 1) * log(each_proposal))
              })
          }
          # Log prior density for phi
          apply(cbind(l.phi, l.tau, l.bdens), 1, sum) + l.lik
        })
        mean.lweights <- apply(l.weights, 2, mean)
        weights <- exp(mean.lweights - max(mean.lweights))
        alpha <- sample(alpha.g, 1, prob = weights)
      }else{
        alpha <- 0
      }

      # >> Draw phi ----
      bg <- betas
      wg <- W
      ch <- mapply(function(b, w) {
        apply(abind(b, w, along = 3), 2, function(each_rank) {
          crossprod(each_rank[, 1], diag(1 / each_rank[, 2])) %*% each_rank[, 1]
        })
      }, b = bg, w = wg)

      ##### Phi under a stick-breaking prior
      old_Xi_g <- Xi
      if(ranks == 1){
        phi.g <- 1
      }else{
      accept <- results$accept
      new_Xi_g <- old_Xi_g + c(cov_Metro%*%rnorm(ranks - 1))
      while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
        new_Xi_g <- old_Xi_g + c(cov_Metro%*%rnorm(ranks - 1))
      }
      new_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
        stick_break_log_posterior(new_Xi_g,cr, betas,W,tau,alpha)
      }))
      old_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
        stick_break_log_posterior(old_Xi_g,cr, betas,W,tau,alpha)
      }))
      if(exp(new_post_dens - old_post_dens) > runif(1)){
        old_Xi_g <- new_Xi_g
        accept <- accept + 1
      }
      phi.g <- stick_values(old_Xi_g)
      }
      ##### Phi under a stick-breaking prior

      # >> Draw tau ----
      xi <- a.tau - ranks * sum(p)/2
      ### RANK 1 ADJUSTMENT
      if(ranks == 1){
        chi <- sum(ch)
      }else{
        chi <- sum(apply(ch, 1, sum) / phi.g)
      }
      tau.g <- rgig(1, xi, chi, 2 * b.tau)

      # >> Draw lambda ----
      sumabsb <- sapply(bg, function(each_dim) {
        apply(each_dim, 2, function(each_rank) {
          sum(abs(each_rank))
        })
      })
      if(ranks == 1){
        lambda.g <- sapply(sumabsb,function(each_dim){
          rgamma(ranks,
                 a.lambda + p,
                 b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
        })
        lambda.g <- t(lambda.g)
      }else{
        lambda.g <-
          apply(sumabsb, 2, function(each_dim) {
            rgamma(ranks,
                   a.lambda + p,
                   b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
          })
      }

      # >> Draw omega ----
      omega.g <- sapply(seq(D), function(each_dim) {
        sapply(seq(ranks), function(each_rank) {
          ch <-
            sapply(bg[[each_dim]][, each_rank], function(each_value) {
              each_value ^ 2 / (tau.g * phi.g[each_rank])
            })
          rgig(p[each_dim], 0.5, ch, lambda.g[each_rank, each_dim])
        }, simplify = T)
      }, simplify = F)

      # >> Draw beta ----
      # This next part has to be done sequentially, so the for loops are unavoidable
      for (each_dim in seq(D)) {
        for (each_rank in 1:ranks) {
          if(ranks == 1){
            expected <- array(0,dim=c(p,TT,n))
          }else{
            expected <- composeParafac(lapply(bg, function(each_dim_I) {each_dim_I[, -each_rank,drop = FALSE]})) %o% input$x
          }
          y.til <- (input$Y - expected) %>%
            apply((D+1):(D+2),mode_k_matriz,each_dim) %>%
            array(dim=c(p[each_dim],prod(p[-each_dim]),TT,n))
          betas_by_rank <-
            sapply(seq(ranks), function(each_rank) {
              sapply(seq(D), function(each_dim) {
                bg[[each_dim]][, each_rank]
              },simplify = FALSE)
            },simplify = FALSE) # This restructuring allows for the calculation
          vec_outer_other_betas <- c(Reduce(outer,betas_by_rank[[each_rank]][-each_dim]))
          var <-  1 /  (sum(input$x^2 %o% vec_outer_other_betas^2) / sig2y + (1 / (tau.g * phi.g[each_rank]) / diag(omega.g[[each_dim]][,each_rank]) ))
          mean_beta <- var %*% apply(y.til,1,function(dim_margin){
            sum(sapply(seq(TT),function(each_time){
              sapply(seq(n),function(each_subject){
                input$x[each_time,each_subject] * vec_outer_other_betas * dim_margin[,each_time,each_subject]
              })
            })) / sig2y
          })
          bg[[each_dim]][, each_rank] <-
            rnorm(p[each_dim], mean_beta, sqrt(diag(var)))
          if (each_dim > 1 && each_rank > 1) {
            bg[[each_dim]][1, each_rank] <-
              rtruncnorm(
                1,
                b = bg[[each_dim]][1, (each_rank - 1)],
                mean = (mean_beta)[1],
                sd = sqrt(var)[1]
              )
          }
        }
      }


    al = alpha
    phi = phi.g
    tau = tau.g
    omega = omega.g
    betas = bg
    lambda = lambda.g
    Xi = old_Xi_g

    # >> Draw sig2y ----
    sum_sq_diff <- (input$Y - composeParafac(betas) %o% input$x)^2 %>% sum

    sig2y <- 1/rgamma(1,
                      a.sig + n*TT*sum(prod(p))/2,
                      b.sig + (1/2)*sum(sum_sq_diff))

    W <- omega


    # >> Get the log-likelihood ----
    if(save_llik == TRUE){
      llik <- -0.5*log(2*pi*sig2y)*n*TT*sum(prod(p)) - 0.5/sig2y * sum(sum_sq_diff)
      results$llik[s] <- llik
    }


    results$B[[s]] <- betas
    results$W[[s]] <- W
    results$lambda[[s]] <- lambda
    results$Phi[[s]] <- phi
    results$tau[s] <- tau
    results$sig2y[s] <- sig2y
    results$alpha[s] <- al
    results$accept <- accept

    # save(results,s,file="C:/Users/Dan/Desktop/temp.Rdata")

    # if(save_results & s %% 10 == 0) save(results,file = paste0("intermediate_results",format(Sys.Date(),"%Y%m%d"),".RData"))

    if(!is.null(save_after)){
      if(s %% save_after == 0){
        save(results,file = paste0("../Tests/400_",D,"D_rank_",ranks,"_first_",s,"_samples_",format(Sys.Date(),"%Y%m%d"),".RData"))
      }
    }

    # setTxtProgressBar(pb, s)
    if(s %% ceiling(n.iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Rank = ", ranks,
          " Iteration # ",s," of ",n.iter,
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - beginning_of_sampler, "  seconds #####\n",
          "##### Estimated time remaining: ", ((proc.time()[3] - beginning_of_sampler)/s)*(n.iter - s)," seconds #####\n"
        )
      )
    }

  } # End sampler

  if(n.burn >0){
    results$B <- results$B[-(1:n.burn)]
    results$W <- results$W[-(1:n.burn)]
    results$lambda <- results$lambda[-(1:n.burn)]
    results$Phi <- results$Phi[-(1:n.burn)]
    results$tau <- results$tau[-(1:n.burn)]
    results$sig2y <- results$sig2y[-(1:n.burn)]
  }

  results$accept <- results$accept / n.iter

  results
} # End MCMC function
