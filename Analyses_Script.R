library(bruceR)
library(TreeBUGS)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Fitting the MPT models and running the simulations takes time.
# FALSE loads the saved results instead, TRUE runs everything again.
refit <- FALSE

# Data----
# Exp 1 = DA at retrieval, Exp 2 = DA at encoding
df_E1 <- import("Exp1_Retrieval/data.xlsx")
df_E1$DA <- factor(df_E1$DA, levels = c(0,1,2), labels = c("FA","Addition","Division"))
df_E1$Probe <- factor(df_E1$Probe, levels = c(0,1,2), labels = c("Intact","Related","Unrelated"))

df_E2 <- import("Exp2_Encoding/data.xlsx")
df_E2$DA <- factor(df_E2$DA, levels = c("EncodewithFA","EncodewithAddition","EncodewithDivision"),
                   labels = c("FA","Addition","Division"))
df_E2$Probe <- factor(df_E2$Probe, levels = c("Intact","Related","Unrelated"))

memory_E1 <- df_E1 %>%
  group_by(Subject, DA, Probe) %>%
  summarise(ACC_mean = mean(ACC_Memory), RT_mean = mean(RT_Memory), .groups = "drop")

memory_E2 <- df_E2 %>%
  group_by(Subject, DA, Probe) %>%
  summarise(ACC_mean = mean(ACC_Memory), RT_mean = mean(RT_Memory), .groups = "drop")

# for the cross-experiment comparison: DA cost = FA minus each DA condition
memory_all <- bind_rows(mutate(memory_E1, Testcontext = "Retrieval"),
                        mutate(memory_E2, Testcontext = "Encoding"))

da_cost <- function(d, dv) {
  d %>%
    select(Subject, DA, Probe, Testcontext, all_of(dv)) %>%
    pivot_wider(names_from = DA, values_from = all_of(dv)) %>%
    mutate(cost_Addition = FA - Addition, cost_Division = FA - Division) %>%
    select(-FA, -Addition, -Division) %>%
    pivot_longer(cols = c(cost_Addition, cost_Division), names_to = "Load",
                 names_prefix = "cost_", values_to = "cost")
}

# concurrent arithmetic task
sec_E1 <- import("Exp1_Retrieval/data_secondary.xlsx")
sec_E1$DA <- factor(sec_E1$DA, levels = c("Addition","Division"))
sec_E1$Probe <- factor(sec_E1$Probe, levels = c("Baseline","Intact","Related","Unrelated"))
sec_E1$RT_Number <- as.numeric(sec_E1$RT_Number)

sec_E2 <- import("Exp2_Encoding/data_secondary.xlsx")
sec_E2$DA <- factor(sec_E2$DA, levels = c("Addition","Division"))
sec_E2$Probe <- factor(sec_E2$Probe, levels = c("Baseline","Encoding"))
sec_E2$RT_Number <- as.numeric(sec_E2$RT_Number)

number_E1 <- sec_E1 %>%
  group_by(Subject, DA, Probe) %>%
  summarise(ACC_mean = mean(ACC_Number), RT_mean = mean(RT_Number), .groups = "drop")

number_E2 <- sec_E2 %>%
  group_by(Subject, DA, Probe) %>%
  summarise(ACC_mean = mean(ACC_Number), RT_mean = mean(RT_Number), .groups = "drop")


# Effects of DA on Associative Recognition Accuracy----

## Exp 1----
MANOVA(data = memory_E1, dv = "ACC_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2)

## Exp 2----
MANOVA(data = memory_E2, dv = "ACC_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

## Cross-experiment----
MANOVA(data = da_cost(memory_all, "ACC_mean"), dv = "cost",
       within = c("Load", "Probe"), between = "Testcontext",
       subID = "Subject", digits = 2) %>%
  EMMEANS("Testcontext", by = "Load", digits = 2) %>%
  EMMEANS("Load", digits = 2)


# Effects of DA on Response Latency in the Associative Recognition Task----

## Exp 1----
MANOVA(data = memory_E1, dv = "RT_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

## Exp 2----
MANOVA(data = memory_E2, dv = "RT_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

## Cross-experiment----
MANOVA(data = da_cost(memory_all, "RT_mean"), dv = "cost",
       within = c("Load", "Probe"), between = "Testcontext",
       subID = "Subject", digits = 2, sph.correction = "GG") %>%
  EMMEANS("Testcontext", digits = 2)


# MPT analysis----

## Exp 1----
if (refit) {
  data_mpt_E1 <- import("Exp1_Retrieval/data_mpt.xlsx") %>% subset(select = -c(1))
  model_E1 <- traitMPT(eqnfile = "Exp1_Retrieval/DAretrieval_mpt.eqn", data = data_mpt_E1,
                       modelfilename = "Exp1_Retrieval/mpt_DAretrieval.jags",
                       parEstFile = "Exp1_Retrieval/mpt_DAretrieval.txt",
                       posteriorFile = "Exp1_Retrieval/mpt_DAretrieval.RData",
                       n.iter = 65000*5, n.chain = 3, n.adapt = 10000*5,
                       n.burnin = 15000*5, n.thin = 10, ppp = 1000)
} else {
  load("Exp1_Retrieval/mpt_DAretrieval.RData")
  model_E1 <- fittedModel
}

## Exp 2----
if (refit) {
  data_mpt_E2 <- import("Exp2_Encoding/data_mpt.xlsx") %>% subset(select = -c(1))
  model_E2 <- traitMPT(eqnfile = "Exp2_Encoding/DAencoding_mpt.eqn", data = data_mpt_E2,
                       modelfilename = "Exp2_Encoding/mpt_DAencoding.jags",
                       parEstFile = "Exp2_Encoding/mpt_DAencoding.txt",
                       posteriorFile = "Exp2_Encoding/mpt_DAencoding.RData",
                       n.iter = 65000*5, n.chain = 3, n.adapt = 10000*5,
                       n.burnin = 15000*5, n.thin = 10, ppp = 1000)
} else {
  load("Exp2_Encoding/mpt_DAencoding.RData")
  model_E2 <- fittedModel
}

## Differences between attention conditions----
# posterior of cond1 - cond2 for each parameter, group level
# diff_ratio is the proportion of samples with cond1 > cond2
compare_params <- function(model, params, cond1, cond2, file) {
  res <- data.frame()
  for (p in params) {
    p1 <- paste0(p, cond1)
    p2 <- paste0(p, cond2)
    tp <- transformedParameters(model, list(paste0("diff=", p1, "-", p2),
                                            paste0("p=", p1, ">", p2)), level = "group")
    s <- summary(tp)
    res <- rbind(res, data.frame(parameter = p,
                                 diff_mean = s$statistics["diff", "Mean"],
                                 diff_ratio = s$statistics["p", "Mean"],
                                 diff_lower = s$quantiles["diff", "2.5%"],
                                 diff_upper = s$quantiles["diff", "97.5%"]))
  }
  export(res, file)
  res
}

params <- c("Vi", "Vr", "Gi", "Gr", "a", "b")

compare_params(model_E1, params, "_fa", "_addition", "Exp1_Retrieval/mpt_DAretrieval_FA_Addition.xlsx")
compare_params(model_E1, params, "_fa", "_division", "Exp1_Retrieval/mpt_DAretrieval_FA_Division.xlsx")
compare_params(model_E1, params, "_addition", "_division", "Exp1_Retrieval/mpt_DAretrieval_Addition_Division.xlsx")

compare_params(model_E2, params, "_fa", "_addition", "Exp2_Encoding/mpt_DAencoding_FA_Addition.xlsx")
compare_params(model_E2, params, "_fa", "_division", "Exp2_Encoding/mpt_DAencoding_FA_Division.xlsx")
compare_params(model_E2, params, "_addition", "_division", "Exp2_Encoding/mpt_DAencoding_Addition_Division.xlsx")

## Cross-experiment----
# is the DA effect on a parameter (FA minus DA) different at encoding vs. retrieval?
# The two models were fitted separately, so the posterior samples are just paired by row.
mpt_mean_col <- function(M, label) {
  idx <- match(label, M$mptInfo$thetaUnique)
  if (is.na(idx)) stop("Parameter '", label, "' not found in model's thetaUnique.")
  paste0("mean[", idx, "]")
}

interaction_contrast <- function(M_enc, M_ret, param, da_label,
                                 fa_label = "fa", digits = 4) {
  enc_mat <- as.matrix(M_enc$runjags$mcmc)
  ret_mat <- as.matrix(M_ret$runjags$mcmc)
  draw <- function(M, m, cond) m[, mpt_mean_col(M, paste0(param, "_", cond))]
  enc_eff <- draw(M_enc, enc_mat, fa_label) - draw(M_enc, enc_mat, da_label)
  ret_eff <- draw(M_ret, ret_mat, fa_label) - draw(M_ret, ret_mat, da_label)
  diff    <- enc_eff - ret_eff
  draws <- list(enc_eff, ret_eff, diff)
  m   <- vapply(draws, mean, numeric(1))
  s   <- vapply(draws, sd,   numeric(1))
  lo  <- vapply(draws, function(x) unname(quantile(x, 0.025)), numeric(1))
  hi  <- vapply(draws, function(x) unname(quantile(x, 0.975)), numeric(1))
  pgt <- vapply(draws, function(x) mean(x > 0), numeric(1))
  out <- data.frame(
    Contrast = c(paste0("Encoding: ",  param, "_", fa_label, " - ", param, "_", da_label),
                 paste0("Retrieval: ", param, "_", fa_label, " - ", param, "_", da_label),
                 "Encoding minus Retrieval"),
    Mean    = round(m,   digits),
    SD      = round(s,   digits),
    CI_2.5  = round(lo,  digits),
    CI_97.5 = round(hi,  digits),
    P_gt_0  = round(pgt, digits),
    Estimate = sprintf("%.2f [%.2f, %.2f]", m, lo, hi),
    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# same thing for every parameter x DA condition, in one table
interaction_contrast_table <- function(M_enc, M_ret, params,
                                       da_labels = c("addition", "division"),
                                       fa_label = "fa", digits = 4) {
  grid <- expand.grid(param = params, da = da_labels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rows <- Map(function(p, d) {
    res <- interaction_contrast(M_enc, M_ret, p, d, fa_label = fa_label, digits = digits)
    cbind(Parameter = p, Load = d, res, stringsAsFactors = FALSE)
  }, grid$param, grid$da)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

posterior_cross <- interaction_contrast_table(model_E2, model_E1, params)
openxlsx::write.xlsx(posterior_cross, "Cross-Experiment_Analysis/Encoding_vs_Retrieval_interaction.xlsx")

# Savage-Dickey Bayes factor for the interaction contrast (delta = 0)
# Quantifies evidence that the DA effect is the SAME at encoding and retrieval.
# Induced prior: the model sets mean[s] <- phi(mu[s]) with mu[s] ~ N(0,1), so each
# group mean has an exact Uniform(0,1) prior, independent across parameters and
# across the two (independently fit) models. Hence the interaction contrast
#   delta = (enc_fa - enc_da) - (ret_fa - ret_da) = U1 - U2 - U3 + U4
# has an analytic prior density of 2/3 at delta = 0.
# Because tested and nuisance parameters are prior-independent, the naive
# Savage-Dickey ratio is valid here (Heck, 2018 correction term = 1).
#   BF01 = posterior_density(0) / prior_density(0)
interaction_bf <- function(M_enc, M_ret, param, da_label, fa_label = "fa",
                           rope = 0.05, method = c("logspline", "normal", "kde")) {
  method <- match.arg(method)
  enc_mat <- as.matrix(M_enc$runjags$mcmc)
  ret_mat <- as.matrix(M_ret$runjags$mcmc)
  dr <- function(M, m, cond) m[, mpt_mean_col(M, paste0(param, "_", cond))]
  enc_eff <- dr(M_enc, enc_mat, fa_label) - dr(M_enc, enc_mat, da_label)
  ret_eff <- dr(M_ret, ret_mat, fa_label) - dr(M_ret, ret_mat, da_label)
  delta   <- enc_eff - ret_eff
  prior0  <- 2 / 3                    # exact induced-prior ordinate at delta = 0
  post0 <- switch(method,
    normal    = dnorm(0, mean(delta), sd(delta)),
    kde       = { kd <- density(delta, n = 4096); approx(kd$x, kd$y, 0)$y },
    logspline = { fit <- logspline::logspline(delta); logspline::dlogspline(0, fit) })
  bf01 <- post0 / prior0
  data.frame(
    Parameter = param, Load = da_label,
    post_dens0 = round(post0, 4), prior_dens0 = round(prior0, 4),
    BF01 = round(bf01, 3), BF10 = round(1 / bf01, 3),
    P_in_ROPE = round(mean(abs(delta) < rope), 4),
    P_gt_0    = round(mean(delta > 0), 4),
    stringsAsFactors = FALSE)
}

interaction_bf(model_E2, model_E1, "Vi", "division") # BF01 = 5.064


# Simulated Replications of the Modeling Results----

# group-level estimates of the fitted model (mu, sigma, rho on the probit scale)
sim_pars <- function(model) {
  gp <- model$summary$groupParameters
  nm <- model$mptInfo$thetaNames$Parameter
  mu <- gp$mu[, "Mean"]
  sigma <- gp$sigma[, "Mean"]
  names(mu) <- names(sigma) <- nm
  rho <- matrix(gp[["rho.matrix"]], nrow = 18, ncol = 18, byrow = FALSE)
  dimnames(rho) <- list(nm, nm)
  list(mu = mu, sigma = sigma, rho = rho, mu_prob = pnorm(mu), param_names = nm)
}

# one simulated data set, 12 trials per probe type per attention condition
sim_data <- function(n_subj, n_items_per_cond = 12, mu, sigma, rho, eqnfile, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # same order as in JAGS
  nm <- c("a_addition", "a_division", "a_fa",
          "b_addition", "b_division", "b_fa",
          "Gi_addition", "Gi_division", "Gi_fa",
          "Gr_addition", "Gr_division", "Gr_fa",
          "Vi_addition", "Vi_division", "Vi_fa",
          "Vr_addition", "Vr_division", "Vr_fa")
  stopifnot(length(mu) == 18, length(sigma) == 18)
  mu <- mu[nm]
  sigma <- sigma[nm]
  dimnames(rho) <- list(nm, nm)

  # tree names come from the eqn file
  trees <- gsub("^T_", "", unique(readEQN(eqnfile)$Tree))
  numItems <- setNames(rep(n_items_per_cond, length(trees)), trees)

  gen <- genTraitMPT(N = n_subj, numItems = numItems, eqnfile = eqnfile,
                     mu = mu, sigma = sigma, rho = rho)
  gen$data
}

# does the comparison (e.g. "Vi_fa - Vi_division") go in the expected direction,
# and does its CI exclude 0?
check_effect <- function(group_draws, comparison, expected_direction = "positive", ci_prob = 0.95) {
  lower_p <- (1 - ci_prob) / 2
  upper_p <- 1 - lower_p
  draws <- eval(parse(text = comparison), envir = as.data.frame(group_draws))
  mean_diff <- mean(draws)
  ci_lower <- as.numeric(quantile(draws, lower_p))
  ci_upper <- as.numeric(quantile(draws, upper_p))

  if (expected_direction == "positive") {
    direction_met <- mean_diff > 0
    ci_met <- ci_lower > 0
    p_wrong_sign <- mean(draws < 0)
  } else {
    direction_met <- mean_diff < 0
    ci_met <- ci_upper < 0
    p_wrong_sign <- mean(draws > 0)
  }

  list(comparison = comparison, mean_diff = mean_diff, ci_lower = ci_lower, ci_upper = ci_upper,
       p_wrong_sign = p_wrong_sign, direction_met = direction_met, ci_met = ci_met,
       both_met = direction_met & ci_met)
}

# Null simulation: null_from is target = source, e.g. c(Vi_division = "Vi_fa").
# Each target gets the mean and SD of its source. The correlations can't just be
# copied over (the result is usually not a valid correlation matrix), so rho is
# averaged with a version of itself in which targets and sources are swapped.
# The two conditions then have the same joint distribution.
set_null <- function(pars, null_from) {
  targets <- names(null_from)
  sources <- unname(null_from)
  nm <- names(pars$mu)
  stopifnot(!anyDuplicated(c(targets, sources)), all(c(targets, sources) %in% nm))

  pars$mu[targets] <- pars$mu[sources]
  pars$sigma[targets] <- pars$sigma[sources]

  swapped <- setNames(nm, nm)
  swapped[targets] <- sources
  swapped[sources] <- targets
  rho <- pars$rho
  rho_swapped <- rho[swapped, swapped]
  dimnames(rho_swapped) <- dimnames(rho)
  pars$rho <- (rho + rho_swapped) / 2

  pars$mu_prob <- pnorm(pars$mu)
  pars
}

# The whole loop. Data are generated from gen_eqnfile (all three conditions) and
# fitted with eqnfile, after dropping the columns matching drop_cols.
# A fit that does not converge is thrown away and that simulation is redone with
# a new data set, up to max_attempts fits in total.
# Simulations whose files are already in output_dir are loaded instead of refitted,
# so several R sessions can work on different sim_ids at the same time.
run_recovery <- function(model, comparison, eqnfile, gen_eqnfile, output_dir, drop_cols,
                         true_value = NULL, null_from = NULL, seed = NULL,
                         n_sim = 100, sim_ids = seq_len(n_sim),
                         n_subj = 50, n_items_per_cond = 12,
                         expected_direction = "positive", ci_prob = 0.95,
                         rhat_threshold = 1.05, max_attempts = 2 * length(sim_ids),
                         n.iter = 40000, n.chain = 4, n.adapt = 10000,
                         n.burnin = 10000, n.thin = 10) {

  if (!is.null(seed)) set.seed(seed)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  cat(sprintf("%d simulations, n = %d, %s, %.0f%% CI\n", n_sim, n_subj, comparison, ci_prob * 100))

  pars <- sim_pars(model)
  if (!is.null(null_from)) pars <- set_null(pars, null_from)

  sim_results <- data.frame(
    sim_id = integer(n_sim), converged = logical(n_sim), max_rhat = numeric(n_sim),
    mean_diff = numeric(n_sim), ci_lower = numeric(n_sim), ci_upper = numeric(n_sim),
    p_wrong_sign = numeric(n_sim), direction_met = logical(n_sim), ci_met = logical(n_sim),
    both_met = logical(n_sim), n_attempts = integer(n_sim), stringsAsFactors = FALSE)

  total_attempts <- 0
  for (i in sim_ids) {
    cat(sprintf("Simulation %d/%d\n", i, n_sim))

    sim_num <- sprintf("%02d", i)
    modelfile <- file.path(output_dir, paste0("sim", sim_num, ".jags"))
    parEstFile <- file.path(output_dir, paste0("sim", sim_num, ".txt"))
    drawsFile <- file.path(output_dir, paste0("sim", sim_num, "_draws.rds"))

    attempt <- 0
    converged <- FALSE
    max_rhat <- NA

    # already done in an earlier session
    if (file.exists(drawsFile) && file.exists(parEstFile)) {
      par_rows <- grep("^mean_", readLines(parEstFile), value = TRUE)
      rhat_values <- as.numeric(vapply(strsplit(trimws(par_rows), "\\s+"),
                                       `[`, character(1), 9))
      group_draws <- readRDS(drawsFile)
      attempt <- attr(group_draws, "attempt")
      if (is.null(attempt)) attempt <- 1
      total_attempts <- total_attempts + attempt

      max_rhat <- max(rhat_values, na.rm = TRUE)
      # the txt file rounds Rhat to 3 decimals, use the exact value saved with the draws
      if (!is.null(attr(group_draws, "max_rhat"))) max_rhat <- attr(group_draws, "max_rhat")
      converged <- max_rhat < rhat_threshold
      cat(sprintf("  loaded from saved files (attempt %d): max Rhat = %.3f%s\n",
                  attempt, max_rhat, if (converged) "" else " (not converged, discarded)"))
    }

    while (!converged && total_attempts < max_attempts) {
      attempt <- attempt + 1
      total_attempts <- total_attempts + 1

      # every attempt has its own seed, so a retry gets a new data set
      d <- sim_data(n_subj = n_subj, n_items_per_cond = n_items_per_cond,
                    mu = pars$mu, sigma = pars$sigma, rho = pars$rho, eqnfile = gen_eqnfile,
                    seed = if (!is.null(seed)) seed + (attempt - 1) * n_sim + i else NULL)
      d <- d[, !grepl(drop_cols, colnames(d))]

      model_sim <- traitMPT(eqnfile = eqnfile, data = d,
                            modelfilename = modelfile, parEstFile = parEstFile,
                            n.iter = n.iter, n.chain = n.chain, n.adapt = n.adapt,
                            n.burnin = n.burnin, n.thin = n.thin,
                            ppp = 0, monitorIndividual = FALSE)

      rhat_values <- model_sim[["summary"]][["groupParameters"]][["mean"]][, "Rhat"]

      # only the group-level mean[] samples are kept, chains stacked, to save space
      nm_fit <- model_sim$mptInfo$thetaNames$Parameter
      mean_cols <- sprintf("mean[%d]", seq_along(nm_fit))
      group_draws <- do.call(rbind, lapply(model_sim$runjags$mcmc,
                                           function(ch) ch[, mean_cols, drop = FALSE]))
      colnames(group_draws) <- nm_fit
      attr(group_draws, "attempt") <- attempt
      attr(group_draws, "max_rhat") <- max(rhat_values, na.rm = TRUE)
      saveRDS(group_draws, drawsFile)

      rm(d, model_sim)
      gc()

      max_rhat <- max(rhat_values, na.rm = TRUE)
      converged <- max_rhat < rhat_threshold
      cat(sprintf("  attempt %d: max Rhat = %.3f%s\n", attempt, max_rhat,
                  if (converged) "" else " (not converged, discarded)"))
    }

    if (converged) {
      chk <- check_effect(group_draws, comparison, expected_direction, ci_prob)
      sim_results[i, ] <- list(
        sim_id = i, converged = TRUE, max_rhat = max_rhat,
        mean_diff = chk$mean_diff, ci_lower = chk$ci_lower, ci_upper = chk$ci_upper,
        p_wrong_sign = chk$p_wrong_sign, direction_met = chk$direction_met,
        ci_met = chk$ci_met, both_met = chk$both_met, n_attempts = attempt)
      rm(chk, group_draws)
    } else {
      # ran out of attempts
      if (attempt > 0) {
        sim_results[i, ] <- list(
          sim_id = i, converged = FALSE, max_rhat = max_rhat,
          mean_diff = NA, ci_lower = NA, ci_upper = NA, p_wrong_sign = NA,
          direction_met = FALSE, ci_met = FALSE, both_met = FALSE, n_attempts = attempt)
        rm(group_draws)
      }
      warning(sprintf("max_attempts (%d) reached at simulation %d.", max_attempts, i))
      break
    }
  }

  sim_results <- sim_results[sim_results$sim_id %in% sim_ids, ]
  n_requested <- length(sim_ids)

  # summary over the converged fits
  ok <- sim_results[sim_results$converged, ]
  n_converged <- nrow(ok)
  n_discarded <- total_attempts - n_converged

  cat(sprintf("\nConverged: %d/%d (%d fits in total, %d discarded)\n",
              n_converged, n_requested, total_attempts, n_discarded))

  if (n_converged > 0) {
    if (!is.null(true_value)) {
      covered <- ok$ci_lower <= true_value & ok$ci_upper >= true_value
      n_covered <- sum(covered)
      coverage_rate <- mean(covered)
    } else {
      n_covered <- NA
      coverage_rate <- NA
    }

    summary_stats <- list(
      n_simulations = n_requested,
      n_subj = n_subj,
      ci_prob = ci_prob,
      n_converged = n_converged,
      n_failed = n_requested - n_converged,
      n_attempts = total_attempts,
      n_discarded = n_discarded,
      convergence_rate = n_converged / total_attempts,
      direction_met_rate = mean(ok$direction_met),
      ci_met_rate = mean(ok$ci_met),
      both_met_rate = mean(ok$both_met),
      # CI excludes 0 in either direction = false positives in the null simulation
      ci_excludes_zero_rate = mean(ok$ci_lower > 0 | ok$ci_upper < 0),
      p_direction_mean = mean(1 - ok$p_wrong_sign),
      mean_diff_mean = mean(ok$mean_diff),
      mean_diff_sd = sd(ok$mean_diff),
      true_value = if (is.null(true_value)) NA else true_value,
      n_covered = n_covered,
      coverage_rate = coverage_rate)
  } else {
    summary_stats <- list(
      n_simulations = n_requested, n_subj = n_subj, ci_prob = ci_prob,
      n_converged = 0, n_failed = n_requested,
      n_attempts = total_attempts, n_discarded = n_discarded, convergence_rate = 0,
      direction_met_rate = NA, ci_met_rate = NA, both_met_rate = NA,
      ci_excludes_zero_rate = NA, p_direction_mean = NA,
      mean_diff_mean = NA, mean_diff_sd = NA,
      true_value = if (is.null(true_value)) NA else true_value,
      n_covered = NA, coverage_rate = NA)
    warning("No simulations converged.")
  }

  list(simulation_results = sim_results, summary_stats = summary_stats)
}

## Exp 1----
# Vi under FA vs. division; the addition condition is left out of the refits
pa_E1 <- sim_pars(model_E1)

if (refit) {
  true_E1 <- run_recovery(model_E1, comparison = "Vi_fa - Vi_division",
                          eqnfile = "Exp1_Retrieval/DAretrieval_mpt_simulation.eqn",
                          gen_eqnfile = "Exp1_Retrieval/DAretrieval_mpt.eqn",
                          output_dir = "Exp1_Retrieval/ModelSimulation", drop_cols = "_Addition$",
                          true_value = unname(pa_E1$mu_prob["Vi_fa"] - pa_E1$mu_prob["Vi_division"]),
                          seed = 2026)
  saveRDS(true_E1, "Exp1_Retrieval/ModelSimulation/Undertrue_recovery_result_001_100.rds")

  # no true difference between FA and division
  null_E1 <- run_recovery(model_E1, comparison = "Vi_fa - Vi_division",
                          eqnfile = "Exp1_Retrieval/DAretrieval_mpt_simulation.eqn",
                          gen_eqnfile = "Exp1_Retrieval/DAretrieval_mpt.eqn",
                          output_dir = "Exp1_Retrieval/ModelSimulation", drop_cols = "_Addition$",
                          true_value = 0,
                          null_from = c(a_division = "a_fa", b_division = "b_fa",
                                        Gi_division = "Gi_fa", Gr_division = "Gr_fa",
                                        Vi_division = "Vi_fa", Vr_division = "Vr_fa"),
                          seed = 2026)
  saveRDS(null_E1, "Exp1_Retrieval/ModelSimulation/Undernull_recovery_result_001_100.rds")
} else {
  true_E1 <- readRDS("Exp1_Retrieval/ModelSimulation/Undertrue_recovery_result_001_100.rds")
  null_E1 <- readRDS("Exp1_Retrieval/ModelSimulation/Undernull_recovery_result_001_100.rds")
}
print(true_E1$summary_stats)
print(null_E1$summary_stats)

## Exp 2----
# Gi under FA vs. addition; the division condition is left out of the refits
pa_E2 <- sim_pars(model_E2)

if (refit) {
  true_E2 <- run_recovery(model_E2, comparison = "Gi_fa - Gi_addition",
                          eqnfile = "Exp2_Encoding/DAencoding_mpt_simulation.eqn",
                          gen_eqnfile = "Exp2_Encoding/DAencoding_mpt.eqn",
                          output_dir = "Exp2_Encoding/ModelSimulation", drop_cols = "Division$",
                          true_value = unname(pa_E2$mu_prob["Gi_fa"] - pa_E2$mu_prob["Gi_addition"]),
                          seed = 2026)
  saveRDS(true_E2, "Exp2_Encoding/ModelSimulation/Undertrue_recovery_result_001_100.rds")

  null_E2 <- run_recovery(model_E2, comparison = "Gi_fa - Gi_addition",
                          eqnfile = "Exp2_Encoding/DAencoding_mpt_simulation.eqn",
                          gen_eqnfile = "Exp2_Encoding/DAencoding_mpt.eqn",
                          output_dir = "Exp2_Encoding/ModelSimulation", drop_cols = "Division$",
                          true_value = 0,
                          null_from = c(a_addition = "a_fa", b_addition = "b_fa",
                                        Gi_addition = "Gi_fa", Gr_addition = "Gr_fa",
                                        Vi_addition = "Vi_fa", Vr_addition = "Vr_fa"),
                          seed = 2028)
  saveRDS(null_E2, "Exp2_Encoding/ModelSimulation/Undernull_recovery_result_001_100.rds")
} else {
  true_E2 <- readRDS("Exp2_Encoding/ModelSimulation/Undertrue_recovery_result_001_100.rds")
  null_E2 <- readRDS("Exp2_Encoding/ModelSimulation/Undernull_recovery_result_001_100.rds")
}
print(true_E2$summary_stats)
print(null_E2$summary_stats)


# Concurrent Arithmetic Task Performance----

## Exp 1----
MANOVA(data = number_E1, dv = "ACC_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

MANOVA(data = number_E1, dv = "RT_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2, sph.correction = "GG") %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

## Exp 2----
MANOVA(data = number_E2, dv = "ACC_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

MANOVA(data = number_E2, dv = "RT_mean", within = c("DA", "Probe"),
       subID = "Subject", digits = 2, sph.correction = "GG") %>%
  EMMEANS("DA", digits = 2) %>%
  EMMEANS("Probe", digits = 2) %>%
  EMMEANS("DA", by = "Probe", digits = 2) %>%
  EMMEANS("Probe", by = "DA", digits = 2)

## Cross-experiment----
# cost of doing the arithmetic task together with the memory task, relative to baseline
number_cost <- bind_rows(mutate(number_E1, Testcontext = "Retrieval"),
                         mutate(number_E2, Testcontext = "Encoding")) %>%
  mutate(Probe_collapsed = ifelse(Probe == "Baseline", "Baseline", "DA_phase")) %>%
  group_by(Subject, DA, Testcontext, Probe_collapsed) %>%
  summarise(ACC_mean = mean(ACC_mean), RT_mean = mean(RT_mean), .groups = "drop") %>%
  pivot_wider(names_from = Probe_collapsed, values_from = c(ACC_mean, RT_mean)) %>%
  mutate(cost_ACC = ACC_mean_Baseline - ACC_mean_DA_phase,
         cost_RT  = RT_mean_DA_phase - RT_mean_Baseline) %>%
  select(Subject, Load = DA, Testcontext, cost_ACC, cost_RT)

MANOVA(data = number_cost, dv = "cost_ACC", within = "Load", between = "Testcontext",
       subID = "Subject", digits = 2) %>%
  EMMEANS("Testcontext", by = "Load", digits = 2) %>%
  EMMEANS("Testcontext", digits = 2)

MANOVA(data = number_cost, dv = "cost_RT", within = "Load", between = "Testcontext",
       subID = "Subject", digits = 2) %>%
  EMMEANS("Testcontext", by = "Load", digits = 2)
