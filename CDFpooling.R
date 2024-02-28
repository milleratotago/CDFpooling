# Routines to be used for CDF pooling analyses based on ex-Gaussian RT distributions.
# Notes:
# - The ex-Gaussian parameter values assumed here are mu, sigma, and tau.
# - mu & sigma are the parameters of the normal component, and tau is the mean of the exponential.

# These routines use the ex-Gaussian probability density & other
# functions from the gamlss library, so it should be installed
# (if not already present) and sourced so that its functions
# will be recognized.
if (!require("gamlss")) {  install.packages("gamlss") }
library(gamlss.dist)

pexGAUSbnd <- function(x, mu = 5, sigma = 1, nu = 1) {
  # This function mostly just calls gamlss's pexGAUS, but
  # that function occasionally produces values outside the
  # 0--1 range. This function brings those back to 0 or 1.
  p <- pexGAUS(x, mu=mu, sigma=sigma, nu=nu)
  p[p<0] <- 0
  p[p>1] <- 1
  return(p)
}

exG_lnlike <- function(parms, rts, rt_bounds = NULL) {
  # Compute the log-likelihood of a vector of RT values under a specific
  # set of POSITIVE ex-Gaussian parameter values mu, sigma, and tau.
  #
  # The optional parameter rt_bounds is a 2-position vector giving
  # the lower and upper bounds of RTs, which would be appropriate
  # if the RTs had to fall within a certain interval (e.g., if
  # RTs below 150 were discarded and RTs had to be less than 2000
  # because the trial ended then).  This seems to have negligible
  # effect unless there are a lot of RTs near one of the boundaries.
  mu <- abs(parms[1])
  sigma <- abs(parms[2])
  tau <- abs(parms[3])
  
  pdf_values <- dexGAUS(rts, mu = mu, sigma = sigma, nu = tau)
  
  if (!is.null(rt_bounds)) {
    # Adjust likelihood to allow for the bounds on RTs
    lower_bound <- rt_bounds[1]
    upper_bound <- rt_bounds[2]
    if ( (min(rts) < lower_bound) || (max(rts) > upper_bound) ) {
      stop("Likelihood cannot be computed if any RTs are outside the bounds")
    }
    lower_prob <- pexGAUSbnd(lower_bound, mu = mu, sigma = sigma, nu = tau)
    upper_prob <- 1 - pexGAUSbnd(upper_bound, mu = mu, sigma = sigma, nu = tau)
    lost_prob <- lower_prob + upper_prob
    pdf_values <- pdf_values / (1 - lost_prob)
  }
  
  lnlikes <- log(pdf_values)
  return(sum(lnlikes))
}

exG_estimate <- function(rts, start_prop_var_in_tau = c(0.2, 0.4, 0.6, 0.8),
                         rt_bounds = NULL,
                         gr = NULL,
                         method = "BFGS",
                         lower = -Inf, upper = Inf,
                         control = list(),
                         hessian = FALSE) {
  # Estimate the ex-Gaussian parameters mu, sigma, and tau
  # for a given set of rts (by maximum likelihood).
  #
  # start_prop_var_in_tau is a value or a vector of values
  # in the range 0-1.  optim() is called once for each
  # value of start_prop_var_in_tau, and the starting
  # parameter values passed to optim() are determined by
  # these 0-1 value(s) as shown below.
  #
  # For an explanation of rt_bounds, look at the fn exG_lnlike()
  #
  # The remaining parameters are passed to the optim() function,
  # so check that function for their meanings.  Except:
  #   control$fnscale must negative to optimize/maximize likelihood,
  #   so we make sure that it is.
  control$fnscale <- -1
  
  rts_mean <- mean(rts)
  rts_variance <- var(rts)
  
  n_optim_calls <- length(start_prop_var_in_tau)
  
  # This loop calls optim() multiple times with different starting points,
  # and it saves the best outcome as the final result
  for (i_call in 1:n_optim_calls) {
    # Here the starting parameter values are computed
    # based on the guessed proportions of the RT variance
    # in the exponential (tau) component
    start_tau_variance <- rts_variance * start_prop_var_in_tau[i_call]
    start_normal_variance <- rts_variance - start_tau_variance  # The remaining variance is in sigma
    start_tau <- sqrt(start_tau_variance)  # the mean of an exponential with the specified amount of variance
    start_sigma <- sqrt(start_normal_variance)  # the normal sigma
    start_mu <- rts_mean - start_tau            # start mu so that mu+tau == observed mean
    start_parms <- c(start_mu, start_sigma, start_tau)
    
    # Call optim to do the hard work.
    est_results <- optim(start_parms, exG_lnlike,
                         rts=rts, rt_bounds = rt_bounds,
                         gr=gr, method=method, lower = lower, upper = upper,
                         control = control, hessian = hessian)
    
    # Save the results if they are the first or if they are better
    # than the previous best.
    if (i_call==1) {
      final_results <- est_results
      final_results$start_prop_tau_variance <- start_prop_var_in_tau[i_call]
    } else if (est_results$value>final_results$value) {
      final_results <- est_results
      final_results$start_prop_tau_variance <- start_prop_var_in_tau[i_call]
    }
    # print(paste(i_call,final_results$value))  # for debugging
  } # end for i_call
  final_results$par <- abs(final_results$par)
  return(final_results)
}

CDFpooling_main <- function(rt_df,
                            min_trials_required = 50,
                            start_prop_var_in_tau = c(0.1, 0.3, 0.5, 0.7),
                            rt_bounds = NULL,
                            gr = NULL,
                            method = "BFGS",
                            lower = -Inf, upper = Inf,
                            control = list(), hessian = FALSE) {
  
  uniq_subs <- unique(rt_df$sub)   # NOTE: Take subs in order of file
  uniq_conds <- sort(unique(rt_df$cond))  # but sort conditions
  n_subs <- length(uniq_subs)
  n_conds <- length(uniq_conds)
  
  est_results <- vector(mode="list", length = n_subs*n_conds)
  
  n_trials <- nrow(rt_df)
  cdf <- rep(NaN,n_trials)
  
  # Check to see whether rt_df has a variable 'excluded'.
  # If not, create one and set it to false for all trials
  # so that all trials will be included.
  if ( !('excluded' %in% names(rt_df)) ) {
    rt_df$excluded <- rep(FALSE,n_trials)
    print("Analysis includes all trials because no variable 'excluded' is present.")
  }
  results_computed <- 0
  for (sub in uniq_subs) {
    for (cond in uniq_conds) {
      selected_rows <- (rt_df$sub == sub) & (rt_df$cond == cond) & !(rt_df$excluded)
      selected_rts <- rt_df$rt[selected_rows]
      n_trials_selected <- length(selected_rts)
      if (n_trials_selected < min_trials_required) {
        warning(paste("Skipping estimation for subject",sub,"condition",cond,"because there are only",n_trials_selected,"trials."))
        est_result <- list(estimated=FALSE)
      } else {
        est_result <- exG_estimate(selected_rts,
                                   start_prop_var_in_tau = start_prop_var_in_tau,
                                   rt_bounds = rt_bounds,
                                   gr = gr,
                                   method = method,
                                   lower = lower, upper = upper,
                                   control = control, hessian = hessian)
        est_result$estimated <- TRUE
        # Compute the CDFs of these RTs within their estimated ex-Gaussian
        this_mu <- est_result$par[1]
        this_sigma <- est_result$par[2]
        this_tau <- est_result$par[3]
        these_cdfs <- pexGAUSbnd(selected_rts, mu=this_mu, sigma=this_sigma, nu=this_tau)
        cdf[selected_rows] <- these_cdfs
      }
      est_result$sub <- sub
      est_result$cond <- cond
      est_result$n_trials <- n_trials_selected
      results_computed <- results_computed + 1
      est_results[[results_computed]] <- est_result
    } # for cond
  } # for sub
  
  # Everything has now been computed.

  # The trial-by-trial to-be-pooled CDFs are in the vector 'cdf'.

  # The output of optim(), including the best parameter estimates for each
  # subject/condition combination are in est_results.  The following function
  # converts that to a data frame for ease of accessing the info.
  est_results_df <- est_results_to_df(est_results)
  
  # These two structures are returned in a two-item list:
  return( list(cdf=cdf, est_results_df=est_results_df) )
}

# The next 2 functions are used to extract the values in the list est_results
# to a data frame that is more convenient to process.
myExtract <- function(x, field_wanted, idx=NULL) {
  # This is a helper function to extract info from a field in 'x'
  if (field_wanted %in% names(x)) {
    field_val <- x[[field_wanted]]
    if (!is.null(idx)) {
      field_val <- field_val[idx]
    } # if is.null
  } else {
    field_val <- NA
  }
  return(field_val)
}

est_results_to_df <- function(est_results) {
  subs <- sapply(est_results, FUN=myExtract, field_wanted = "sub")
  conds <- sapply(est_results, FUN=myExtract, field_wanted = "cond")
  mus <- sapply(est_results, FUN=myExtract, field_wanted = "par", idx=1)
  sigmas <- sapply(est_results, FUN=myExtract, field_wanted = "par", idx=2)
  taus <- sapply(est_results, FUN=myExtract, field_wanted = "par", idx=3)
  n_trialss <- sapply(est_results, FUN=myExtract, field_wanted = "n_trials")
  estimateds <- sapply(est_results, FUN=myExtract, field_wanted = "estimated")
  values <- sapply(est_results, FUN=myExtract, field_wanted = "value")
  convergences <- sapply(est_results, FUN=myExtract, field_wanted = "convergence")
  messages <- sapply(est_results, FUN=myExtract, field_wanted = "message")
  start_prop_tau_variances <- sapply(est_results, FUN=myExtract, field_wanted = "start_prop_tau_variance")
  count_functions <- sapply(est_results, FUN=myExtract, field_wanted = "counts", idx=1)
  count_gradients <- sapply(est_results, FUN=myExtract, field_wanted = "counts", idx=2)
  messages[sapply(messages,is.null)] <- ""
  messages <- unlist(messages)
  est_results_df <- data.frame(sub=subs,
                               cond=conds,
                               n_trials=n_trialss,
                               estimated=estimateds,
                               mu=mus,
                               sigma=sigmas,
                               tau=taus,
                               value=values,
                               convergence=convergences,
                               start_prop_tau_variance=start_prop_tau_variances,
                               count_function=count_functions,
                               count_gradient=count_gradients,
                               message=messages)
  return(est_results_df)
}
