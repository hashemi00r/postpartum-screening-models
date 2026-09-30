source('R/model.R')
set.seed(config$seed)
dir.create('results/microsimulation', recursive = TRUE, showWarnings = FALSE)

n <- 100000L
has_depression <- runif(n) < base['prevalence']
screen_draw <- runif(n)
response_draw <- runif(n)

simulate <- function(strategy) {
  p <- strategies[[strategy]]
  positive <- ifelse(has_depression, screen_draw < base[[p$sens]],
                     screen_draw < 1 - base[[p$spec]])
  branch <- ifelse(has_depression, ifelse(positive, 'TP', 'FN'),
                   ifelse(positive, 'FP', 'TN'))
  still_depressed <- branch == 'FN' |
    (branch == 'TP' & response_draw < base[[p$failure]])
  cost <- ifelse(branch == 'TP', p$tp,
                 ifelse(branch == 'FP', p$screen + unit_cost['diagnosis'], p$screen))
  qaly <- ifelse(still_depressed, base['u_depressed'], base['u_well']) *
    config$months / 12
  data.frame(branch, cost = unname(cost), qaly = unname(qaly))
}

tau <- simulate('TAU')
whooley <- simulate('Whooley')
delta_cost <- whooley$cost - tau$cost
delta_qaly <- whooley$qaly - tau$qaly
paired <- data.frame(measure = c('delta_cost', 'delta_qaly', 'incremental_nmb'),
                     estimate = c(mean(delta_cost), mean(delta_qaly),
                                  mean(nmb(delta_cost, delta_qaly))),
                     mc_se = c(sd(delta_cost), sd(delta_qaly),
                               sd(nmb(delta_cost, delta_qaly))) / sqrt(n))
paired$lower_95 <- paired$estimate - 1.96 * paired$mc_se
paired$upper_95 <- paired$estimate + 1.96 * paired$mc_se
benchmark <- evaluate_tree(as.list(base))
validation <- data.frame(
  strategy = c('TAU', 'Whooley + CBT'),
  analytic_cost = c(benchmark$cost_tau, benchmark$cost_whooley),
  simulated_cost = c(mean(tau$cost), mean(whooley$cost)),
  analytic_qaly = c(benchmark$qaly_tau, benchmark$qaly_whooley),
  simulated_qaly = c(mean(tau$qaly), mean(whooley$qaly))
)
branches <- data.frame(strategy = rep(c('TAU', 'Whooley + CBT'), each = 4L),
                       branch = rep(c('TP', 'FN', 'FP', 'TN'), 2L),
                       proportion = c(prop.table(table(factor(tau$branch,
                         levels = c('TP', 'FN', 'FP', 'TN')))),
                         prop.table(table(factor(whooley$branch,
                         levels = c('TP', 'FN', 'FP', 'TN'))))))
write.csv(paired, 'results/microsimulation/paired_summary.csv', row.names = FALSE)
write.csv(validation, 'results/microsimulation/validation.csv', row.names = FALSE)
write.csv(branches, 'results/microsimulation/branches.csv', row.names = FALSE)
print(paired, row.names = FALSE)
print(validation, row.names = FALSE)
