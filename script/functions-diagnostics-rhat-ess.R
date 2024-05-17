# Various functions taken from:
#   https://github.com/avehtari/rhat_ess

# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Trustees of Columbia University
# Copyright (C) 2018, 2019 Aki Vehtari, Paul Bürkner
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Rank uniformization
#'
#' Compute rank uniformization for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, uniformize
#' ranks to scale [1/(2S), 1-1/(2S)], where S is the the number of values.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of rank uniformized values with the same
#'     size as input.
u_scale <- function(x) {
  r <- rank(x, ties.method = 'average')
  u <- backtransform_ranks(r)
  u[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    # output should have the input dimension
    u <- array(u, dim = dim(x), dimnames = dimnames(x))
  }
  u
}

#' Backtransformation of ranks
#' @param r array of ranks
#' @param c fractional offset; defaults to c = 3/8 as recommend by Bloom (1985)
backtransform_ranks <- function(r, c = 3/8) {
  S <- length(r)
  (r - c) / (S - 2 * c + 1)
}

#' Effective sample size
#'
#' Compute the effective sample size estimate for a sample of several chains
#' for one parameter. For split-ESS, call this with split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size.
#' 
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
ess_rfun <- function(sims) {
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  # chains = M and n_samples = N in the paper
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  if (n_samples < 3L || anyNA(sims)) {
    return(NA)
  }
  if (any(!is.finite(sims))) {
    return(NaN)
  }
  if (rstan:::is_constant(sims)) {
    return(NA)
  }
  # acov[t,m] = s_m^2 \rho_{m,t} in the paper
  acov <- lapply(seq_len(chains), function(i) rstan:::autocovariance(sims[, i]))
  acov <- do.call(cbind, acov) * n_samples / (n_samples - 1)
  chain_mean <- apply(sims, 2, mean)
  # mean_var = W in the paper
  mean_var <- mean(acov[1, ]) 
  # var_plus = \hat{var}^{+} in the paper
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1)
    var_plus <- var_plus + var(chain_mean)
  
  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_samples)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while (t < nrow(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
         (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t
  # this is used in the improved estimate (see below)
  if (rho_hat_even>0)
    rho_hat_t[max_t + 1] <- rho_hat_even
  
  # Geyer's initial monotone sequence
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
        rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
      rho_hat_t[t + 2] = rho_hat_t[t + 1];
    }
  }
  # nominal sample size S=MN
  S <- chains * n_samples
  # Geyer's truncated estimate is
  #   tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
  # We use an improved estimate, which is equivalent to taking average
  # of truncation with lag max_t and with max_t+1 and which reduces
  # variance in antithetic case
  tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
  # Safety check for negative values and with max ess equal to ess*log10(ess)
  tau_hat <- max(tau_hat, 1/log10(S))
  ess <- S / tau_hat
  ess
}

#' Plot change in ESS (plot_change_ess): Generate ESS
#' @param fit, array, posterior samples
#' @param par, character, parameter to plot
#' @param breaks, numeric, iteration breaks for plotting
#' @param yaxis, character, specifies plotting of absolute or relative ESS
#' @export
plot_change_ess <- function(fit, par, breaks = seq(0.1, 1, 0.05), 
                            yaxis = c("absolute", "relative")) {
  
  
  if (length(par) != 1L) {
    stop("'par' should be of length 1.")
  }
  yaxis <- match.arg(yaxis)
  if (inherits(fit, "stanfit")) {
    if (!is.character(par)) {
      par <- names(fit)[par]
    }
    sims <- as.array(fit, pars = par)[, , 1]
  } else {
    if (!is.character(par)) {
      par <- dimnames(fit)[[3]][par]
    }
    sims <- fit[, , par, drop = TRUE]
  }
  
  iter_breaks <- round(breaks * NROW(sims))
  nbreaks <- length(iter_breaks)
  bulk_seff <- tail_seff <- bulk_reff <- 
    tail_reff <- rep(NA, length(nbreaks))
  for (i in seq_along(iter_breaks)) {
    sims_i <- sims[seq_len(iter_breaks[i]), ]
    nsamples <- prod(dim(sims_i))
    bulk_seff[i] <- rstan:::ess_rfun(
      rstan:::z_scale(rstan:::split_chains(sims_i)))
    tail_seff[i] <- ess_tail(sims_i)
    bulk_reff[i] <- bulk_seff[i] / nsamples
    tail_reff[i] <- tail_seff[i] / nsamples
  }
  df <- data.frame(
    breaks = breaks,
    ndraws = iter_breaks * NCOL(sims),
    seff = c(bulk_seff, tail_seff),
    reff = c(bulk_reff, tail_reff), 
    type = rep(c("bulk", "tail"), each = nbreaks)
  )
  #blues <- bayesplot::color_scheme_get(scheme = "blue", i = c(4, 2))
  #blues <- unname(unlist(blues))
  if (yaxis == "absolute") {
    out <- ggplot(df, aes(ndraws, seff, color = type)) +
      ylab("ESS") +
      geom_hline(yintercept = 0, linetype = 1) +
      geom_hline(yintercept = 400, linetype = 2) +
      labs(color = 'ESS Type', 
           title = paste(par, sep = ''))
  } else if (yaxis == "relative") {
    out <- ggplot(df, aes(ndraws, reff, color = type)) +
      ylab("Relative efficiency") +
      geom_hline(yintercept = 0, linetype = 2) 
  }
  out +  
    geom_line() +
    geom_point() +
    xlab("Total number of draws") #+
    #scale_colour_manual(values = blues) 
}

plot_quantile_ess <- function(fit, par, nalpha = 20, rank = TRUE) {
  if (length(par) != 1L) {
    stop("'par' should be of length 1.")
  }
  if (inherits(fit, "stanfit")) {
    if (!is.character(par)) {
      par <- names(fit)[par]
    }
    sims <- as.array(fit, pars = par)[, , 1]
    params <- as.data.frame(fit, pars = par)
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- do.call(rbind, sampler_params)[, 'divergent__']
    max_depth <- attr(fit@sim$samples[[1]], "args")$control$max_treedepth
    treedepths <- do.call(rbind, sampler_params)[, 'treedepth__']
    params$divergent <- divergent
    params$max_depth <- (treedepths == max_depth) * 1
    params$urank <- u_scale(params[, par])
    params$value <- params[, par]
  } else {
    if (!is.character(par)) {
      par <- dimnames(fit)[[3]][par]
    }
    sims <- fit[, , par, drop = TRUE]
    params <- data.frame(value = as.vector(sims))
    params$divergent <- 0
    params$max_depth <- 0
    params$urank <- u_scale(params$value)
  }
  
  # compute quantile Seff
  delta <- 1 / nalpha
  alphas <- seq(delta, 1 - delta, by = delta)
  zsseffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    I <- sims <= quantile(sims, alpha)
    zsseffs[i] <- ess_rfun(rstan:::split_chains(I))
  }
  S <- prod(dim(I))
  
  # create the plot
  df <- data.frame(
    quantile = seq(delta, 1 - delta, by = delta), 
    value = quantile(params$value, seq(delta, 1 - delta, by = delta)),
    zsseff = zsseffs
  )
  ymax <- max(S, round(max(zsseffs, na.rm = TRUE) * 1.15, 1))
  xname <- if (rank) "quantile" else "value"
  xrug <- if (rank) "urank" else "value"
  out <- ggplot(data = df, aes_string(x = xname, y = "zsseff")) +
    geom_point() + 
    geom_hline(yintercept = c(0, 1)) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = round(0.25*S)), 
      limits = c(0, ymax)
    ) +
    geom_rug(
      data = params[params$divergent == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "red"
    ) +
    geom_rug(
      data = params[params$max_depth == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "orange"
    ) +
    ylab("ESS for quantiles") +
    labs(title = paste(par, sep = ''))
  if (rank) {
    out <- out +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      xlab('Quantile')
  } else {
    out <- out + xlab(par)
  }
  out
}
