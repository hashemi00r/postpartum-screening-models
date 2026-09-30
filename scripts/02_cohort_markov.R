source('R/model.R')
set.seed(config$seed)
dir.create('results/markov', recursive = TRUE, showWarnings = FALSE)

months <- config$months
states <- c('NoDep', 'DepTreated', 'DepUntreated', 'Remission')

run_markov <- function(x, strategy, trace = FALSE) {
  p <- strategies[[strategy]]
  b <- branch_probabilities(x$prevalence, x[[p$sens]], x[[p$spec]])[1, ]
  current <- c(NoDep = b['FP'] + b['TN'], DepTreated = b['TP'],
               DepUntreated = b['FN'], Remission = 0)
  names(current) <- states
  cost <- evaluate_strategy(x, strategy)$cost
  remission <- 1 - x[[p$failure]]^(1 / months)
  utility <- c(x$u_well, x$u_depressed, x$u_depressed, x$u_well)
  history <- matrix(NA_real_, months + 1L, 4L, dimnames = list(NULL, states))
  history[1, ] <- current
  qaly <- 0
  for (m in seq_len(months)) {
    transition <- diag(4)
    # Treatment ends after the last cycle.
    transition[2, 2] <- if (m == months) 0 else 1 - remission
    transition[2, 3] <- if (m == months) 1 - remission else 0
    transition[2, 4] <- remission
    next_state <- drop(current %*% transition)
    qaly <- qaly + sum((current + next_state) / 2 * utility) / 12
    current <- next_state
    history[m + 1L, ] <- current
  }
  stopifnot(all(abs(rowSums(history) - 1) < 1e-10))
  result <- list(cost = cost, qaly = qaly)
  if (trace) result$trace <- data.frame(month = 0:months, history)
  result
}

base_tau <- run_markov(as.list(base), 'TAU', trace = TRUE)
base_whooley <- run_markov(as.list(base), 'Whooley', trace = TRUE)
base_case <- data.frame(strategy = c('TAU', 'Whooley + CBT'),
                        cost = c(base_tau$cost, base_whooley$cost),
                        qaly = c(base_tau$qaly, base_whooley$qaly))
inputs <- sample_inputs(10000)
psa <- t(vapply(seq_len(nrow(inputs)), function(i) {
  x <- as.list(inputs[i, ])
  a <- run_markov(x, 'TAU')
  b <- run_markov(x, 'Whooley')
  c(delta_cost = b$cost - a$cost, delta_qaly = b$qaly - a$qaly)
}, numeric(2)))
ceac <- data.frame(wtp = seq(0, 800000, by = 10000))
ceac$p_whooley <- vapply(ceac$wtp, function(v)
  mean(nmb(psa[, 'delta_cost'], psa[, 'delta_qaly'], v) > 0), numeric(1))
summary <- data.frame(delta_cost = base_whooley$cost - base_tau$cost,
                      delta_qaly = base_whooley$qaly - base_tau$qaly,
                      icer = (base_whooley$cost - base_tau$cost) /
                        (base_whooley$qaly - base_tau$qaly),
                      psa_p_cost_effective = mean(nmb(psa[, 'delta_cost'],
                                                      psa[, 'delta_qaly']) > 0))
write.csv(base_case, 'results/markov/base_case.csv', row.names = FALSE)
write.csv(summary, 'results/markov/summary.csv', row.names = FALSE)
write.csv(ceac, 'results/markov/ceac.csv', row.names = FALSE)
write.csv(data.frame(strategy = rep(base_case$strategy, each = months + 1L),
                      rbind(base_tau$trace, base_whooley$trace)),
          'results/markov/state_trace.csv', row.names = FALSE)
png('results/markov/remission.png', width = 900, height = 600)
plot(0:months, base_tau$trace$Remission, type = 'l', lwd = 2,
     ylim = range(c(base_tau$trace$Remission, base_whooley$trace$Remission)),
     xlab = 'Month', ylab = 'Proportion in remission')
lines(0:months, base_whooley$trace$Remission, lty = 2, lwd = 2)
legend('topleft', c('TAU', 'Whooley + CBT'), lty = 1:2, lwd = 2, bty = 'n')
dev.off()
print(summary, row.names = FALSE)
