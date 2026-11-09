# ---- 0. Setup ----
library(CASdatasets)
library(MASS)              # fitdistr (fits distributions by maximum likelihood)
options(scipen = 999)      # print numbers without scientific notation

data(freMTPL2sev)          # one row per claim: IDpol, ClaimAmount
sev <- freMTPL2sev

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Look at the claim amounts ----
dim(sev)                                   # expect 26444 x 2
summary(sev$ClaimAmount)                   # original units
sum(sev$ClaimAmount <= 0)                  # any zero or negative amounts?
quantile(sev$ClaimAmount, c(0.5, 0.9, 0.95, 0.99, 0.999, 1))
head(sort(sev$ClaimAmount, decreasing = TRUE), 10)          # the 10 largest claims
head(sort(table(sev$ClaimAmount), decreasing = TRUE), 10)   # most repeated amounts

# Histogram on a log scale (claim amounts are very skewed)
hist(log(sev$ClaimAmount[sev$ClaimAmount > 0]), breaks = 60,
     main = "Claim amounts (log scale)", xlab = "log(claim amount)")
dev.copy(png, "figures/severity_hist_log.png", width = 800, height = 500)
dev.off()

# Save a summary of section 1
qs <- unname(quantile(sev$ClaimAmount, c(0.9, 0.99, 0.999)))
sec1 <- data.frame(
  item  = c("claims", "min", "median", "mean", "90th percentile",
            "99th percentile", "99.9th percentile", "max", "amounts <= 0"),
  value = c(nrow(sev), min(sev$ClaimAmount), median(sev$ClaimAmount),
            mean(sev$ClaimAmount), qs[1], qs[2], qs[3],
            max(sev$ClaimAmount), sum(sev$ClaimAmount <= 0))
)
write.csv(sec1, "results/severity_section1_summary.csv", row.names = FALSE)

quartiles_check <- quantile(sev$ClaimAmount, c(0.25, 0.5, 0.75))
print(quartiles_check)
write.csv(data.frame(item = names(quartiles_check), value = as.numeric(quartiles_check)),
          "results/severity_quartile_check.csv", row.names = FALSE)

# ---- 2. Fit gamma and log-normal ----
# We work in THOUSANDS (amount / 1000) because the fitting routine converges
# more reliably on smaller numbers. Remember this when reading the results.
x <- sev$ClaimAmount / 1000
x <- x[x > 0]              # both distributions need positive values

fit_gamma <- fitdistr(x, "gamma")
fit_lnorm <- fitdistr(x, "lognormal")
fit_gamma
fit_lnorm
AIC(fit_gamma, fit_lnorm)  # lower = better

# Pull out the parameters and save them
gs    <- unname(fit_gamma$estimate["shape"])
gr    <- unname(fit_gamma$estimate["rate"])
ln_mu <- unname(fit_lnorm$estimate["meanlog"])
ln_sd <- unname(fit_lnorm$estimate["sdlog"])

sev_models <- data.frame(
  model     = c("Gamma", "LogNormal"),
  param_1   = c(gs, ln_mu),            # gamma shape / log-normal meanlog
  param_2   = c(gr, ln_sd),            # gamma rate  / log-normal sdlog
  fitted_mean_thousands   = c(gs / gr, exp(ln_mu + ln_sd^2 / 2)),
  observed_mean_thousands = mean(x),
  AIC       = c(AIC(fit_gamma), AIC(fit_lnorm))
)
print(sev_models)
write.csv(sev_models, "results/severity_models.csv", row.names = FALSE)

mean_gap <- data.frame(
  model = c("Gamma", "LogNormal"),
  fitted_mean = c(gs / gr, exp(ln_mu + ln_sd^2 / 2)) * 1000,
  observed_mean = mean(x) * 1000,
  pct_off = (c(gs / gr, exp(ln_mu + ln_sd^2 / 2)) - mean(x)) / mean(x) * 100
)
print(mean_gap)
write.csv(mean_gap, "results/severity_mean_gap.csv", row.names = FALSE)

# ---- 3. Compare quantiles: what each model predicts vs what we observed ----
probs <- c(0.5, 0.9, 0.95, 0.99, 0.999)
qtab <- data.frame(
  prob      = probs,
  observed  = as.numeric(quantile(x, probs)),
  gamma     = qgamma(probs, shape = gs, rate = gr),
  lognormal = qlnorm(probs, meanlog = ln_mu, sdlog = ln_sd)
)
print(round(qtab, 2))
write.csv(qtab, "results/severity_quantile_comparison.csv", row.names = FALSE)

# ---- 4. QQ plots (points on the line = good fit; both axes on log scale) ----
p   <- ppoints(length(x))
emp <- sort(x)
q_g <- qgamma(p, shape = gs, rate = gr)
q_l <- qlnorm(p, meanlog = ln_mu, sdlog = ln_sd)

par(mfrow = c(1, 2))
plot(q_g, emp, log = "xy", main = "Gamma",
     xlab = "model quantile (thousands)", ylab = "observed (thousands)")
abline(0, 1, col = "red")
plot(q_l, emp, log = "xy", main = "Log-normal",
     xlab = "model quantile (thousands)", ylab = "observed (thousands)")
abline(0, 1, col = "red")
dev.copy(png, "figures/severity_qq.png", width = 1000, height = 500)
dev.off()
par(mfrow = c(1, 1))

# ---- 5. Mean excess plot: does the tail need its own model? ----
# For each threshold u, the average amount by which claims above u exceed u.
# A curve that keeps rising suggests a heavy tail.
u  <- as.numeric(quantile(x, seq(0.05, 0.99, by = 0.01)))
me <- sapply(u, function(t) mean(x[x > t] - t))

plot(u, me, type = "l", xlab = "threshold u (thousands)",
     ylab = "mean excess", main = "Mean excess plot")
dev.copy(png, "figures/severity_mean_excess.png", width = 800, height = 500)
dev.off()

write.csv(data.frame(threshold_thousands = u, mean_excess_thousands = me),
          "results/severity_mean_excess.csv", row.names = FALSE)
