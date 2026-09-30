source('R/model.R')
set.seed(config$seed)

dir.create('results/tree', recursive = TRUE, showWarnings = FALSE)
base_case <- evaluate_tree(as.list(base))
psa <- sample_inputs(10000)
psa <- cbind(psa, evaluate_tree(psa))
psa$incremental_nmb <- nmb(psa$delta_cost, psa$delta_qaly)

ceac <- data.frame(wtp = seq(0, 800000, by = 10000))
ceac$p_whooley <- vapply(ceac$wtp, function(v)
  mean(nmb(psa$delta_cost, psa$delta_qaly, v) > 0), numeric(1))
summary <- data.frame(
  measure = c('delta_cost', 'delta_qaly', 'icer', 'p_cost_effective'),
  base_case = c(base_case$delta_cost, base_case$delta_qaly,
                base_case$delta_cost / base_case$delta_qaly,
                as.numeric(nmb(base_case$delta_cost, base_case$delta_qaly) > 0)),
  psa_mean = c(mean(psa$delta_cost), mean(psa$delta_qaly),
               mean(psa$delta_cost) / mean(psa$delta_qaly),
               mean(psa$incremental_nmb > 0))
)
write.csv(base_case, 'results/tree/base_case.csv', row.names = FALSE)
write.csv(summary, 'results/tree/summary.csv', row.names = FALSE)
write.csv(ceac, 'results/tree/ceac.csv', row.names = FALSE)

png('results/tree/ceac.png', width = 900, height = 600)
plot(ceac$wtp / 1000, ceac$p_whooley, type = 'l', lwd = 2,
     col = 'steelblue4', ylim = c(0, 1),
     xlab = 'Willingness to pay (thousand HKD per QALY)',
     ylab = 'Probability Whooley + CBT is cost effective')
abline(v = config$wtp / 1000, lty = 2, col = 'grey50')
dev.off()
png('results/tree/cost_effectiveness_plane.png', width = 900, height = 600)
plot(psa$delta_qaly, psa$delta_cost, pch = 16, cex = 0.35,
     col = rgb(0.1, 0.3, 0.5, 0.25),
     xlab = 'Incremental QALYs', ylab = 'Incremental cost (HKD)')
abline(a = 0, b = config$wtp, lty = 2)
dev.off()
print(summary, row.names = FALSE)
