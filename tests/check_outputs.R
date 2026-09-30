source('R/model.R')

x <- as.list(base)
for (strategy in names(strategies)) {
  p <- strategies[[strategy]]
  b <- branch_probabilities(x$prevalence, x[[p$sens]], x[[p$spec]])
  stopifnot(all(b >= 0), isTRUE(all.equal(as.numeric(rowSums(b)), 1)))
}

set.seed(config$seed)
draws <- sample_inputs(1000)
stopifnot(all(as.matrix(draws) > 0), all(as.matrix(draws) < 1))
sample_branches <- branch_probabilities(draws$prevalence,
                                        draws$epds_sens, draws$epds_spec)
stopifnot(all(sample_branches >= 0),
          all(abs(rowSums(sample_branches) - 1) < 1e-12))

same_utility <- x
same_utility$u_depressed <- same_utility$u_well <- 0.7
outcomes <- evaluate_tree(same_utility)
stopifnot(isTRUE(all.equal(outcomes$qaly_tau, 0.7 * config$months / 12)),
          isTRUE(all.equal(outcomes$qaly_whooley, 0.7 * config$months / 12)))

required <- c('results/tree/base_case.csv', 'results/tree/summary.csv',
              'results/markov/summary.csv', 'results/markov/state_trace.csv',
              'results/microsimulation/paired_summary.csv')
stopifnot(all(file.exists(required)))

analytic <- evaluate_tree(x)
reported <- read.csv('results/tree/base_case.csv')
stopifnot(isTRUE(all.equal(reported$delta_cost, analytic$delta_cost)),
          isTRUE(all.equal(reported$delta_qaly, analytic$delta_qaly)))
trace <- read.csv('results/markov/state_trace.csv')
stopifnot(all(abs(rowSums(trace[, c('NoDep', 'DepTreated',
                                      'DepUntreated', 'Remission')]) - 1) < 1e-10))
for (s in unique(trace$strategy)) {
  pathway <- trace[trace$strategy == s, ]
  remission <- pathway$Remission
  stopifnot(all(diff(remission) >= -1e-10))
  p <- strategies[[if (s == 'TAU') 'TAU' else 'Whooley']]
  expected_remission <- base['prevalence'] * base[[p$sens]] *
    (1 - base[[p$failure]])
  stopifnot(abs(tail(remission, 1) - expected_remission) < 1e-10)
  stopifnot(abs(tail(pathway$DepTreated, 1)) < 1e-10)
}
paired <- read.csv('results/microsimulation/paired_summary.csv')
stopifnot(abs(paired$estimate[1] - analytic$delta_cost) < 3 * paired$mc_se[1],
          abs(paired$estimate[2] - analytic$delta_qaly) < 3 * paired$mc_se[2])
cat('Checks passed.\n')
