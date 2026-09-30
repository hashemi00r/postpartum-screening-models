# Inputs from Wang et al. (2025), Table 1; costs in 2024 HKD.
config <- list(seed = 20260514L, wtp = 422191, months = 9L)

base <- c(
  prevalence = 0.197, epds_sens = 0.81, epds_spec = 0.88,
  whooley_sens = 0.95, whooley_spec = 0.65,
  tau_no_improve = 0.650, cbt_no_improve = 0.312,
  u_depressed = 0.630, u_well = 0.860
)
input_sd <- c(
  prevalence = 0.015, epds_sens = 0.015, epds_spec = 0.008,
  whooley_sens = 0.014, whooley_spec = 0.025,
  tau_no_improve = 0.009, cbt_no_improve = 0.008,
  u_depressed = 0.015, u_well = 0.010
)
unit_cost <- c(epds = 121.93, whooley = 58.90, diagnosis = 1260,
               cbt = 2000, psychiatrist = 1260, nurse = 1550, allied = 1730)
shared_treatment <- 3 * unit_cost['psychiatrist'] +
  4 * unit_cost['nurse'] + unit_cost['allied']
strategies <- list(
  TAU = list(screen = unit_cost['epds'], sens = 'epds_sens',
             spec = 'epds_spec', failure = 'tau_no_improve',
             tp = unit_cost['epds'] + unit_cost['diagnosis'] + shared_treatment),
  Whooley = list(screen = unit_cost['whooley'], sens = 'whooley_sens',
                 spec = 'whooley_spec', failure = 'cbt_no_improve',
                 tp = unit_cost['whooley'] + unit_cost['diagnosis'] +
                   8 * unit_cost['cbt'] + shared_treatment)
)

branch_probabilities <- function(prevalence, sensitivity, specificity) {
  cbind(TP = prevalence * sensitivity,
        FN = prevalence * (1 - sensitivity),
        FP = (1 - prevalence) * (1 - specificity),
        TN = (1 - prevalence) * specificity)
}

evaluate_strategy <- function(x, strategy) {
  p <- strategies[[strategy]]
  b <- branch_probabilities(x$prevalence, x[[p$sens]], x[[p$spec]])
  cost <- b[, 'TP'] * p$tp + b[, 'FP'] * (p$screen + unit_cost['diagnosis']) +
    (b[, 'FN'] + b[, 'TN']) * p$screen
  expected_utility <- b[, 'TP'] * (x[[p$failure]] * x$u_depressed +
    (1 - x[[p$failure]]) * x$u_well) + b[, 'FN'] * x$u_depressed +
    (b[, 'FP'] + b[, 'TN']) * x$u_well
  data.frame(cost = unname(cost), qaly = unname(expected_utility * config$months / 12))
}

evaluate_tree <- function(x) {
  tau <- evaluate_strategy(x, 'TAU')
  whooley <- evaluate_strategy(x, 'Whooley')
  data.frame(cost_tau = tau$cost, cost_whooley = whooley$cost,
             qaly_tau = tau$qaly, qaly_whooley = whooley$qaly,
             delta_cost = whooley$cost - tau$cost,
             delta_qaly = whooley$qaly - tau$qaly)
}

beta_draw <- function(n, mu, sigma) {
  concentration <- mu * (1 - mu) / sigma^2 - 1
  stopifnot(is.finite(concentration), concentration > 0)
  rbeta(n, mu * concentration, (1 - mu) * concentration)
}

sample_inputs <- function(n) {
  x <- lapply(names(base), function(k) beta_draw(n, base[[k]], input_sd[[k]]))
  names(x) <- names(base)
  as.data.frame(x)
}

nmb <- function(delta_cost, delta_qaly, wtp = config$wtp) {
  wtp * delta_qaly - delta_cost
}
