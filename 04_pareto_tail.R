library(CASdatasets)
library(MASS)
options(scipen = 999)

data(freMTPL2sev)
sev <- freMTPL2sev
x   <- sev$ClaimAmount / 1000     # work in thousands

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 0. Threshold ----
threshold_th <- 10   # = $10,000, in thousands 

is_large <- x > threshold_th
p_large  <- mean(is_large)         # share of ALL claims that are "large"
n_large  <- sum(is_large)
n_small  <- sum(!is_large)

cat("Threshold: $10,000\n")
cat("Claims above threshold:", n_large, "of", length(x),
    sprintf("(%.2f%%)\n", 100 * p_large))

write.csv(data.frame(item = c("threshold", "n_total", "n_large", "n_small", "p_large"),
                     value = c(10000, length(x), n_large, n_small, p_large)),
          "results/tail_threshold_summary.csv", row.names = FALSE)

# ---- 1. Refit log-normal to the BODY (small claims only) ----
body_fit <- fitdistr(x[!is_large], "lognormal")
body_mu  <- unname(body_fit$estimate["meanlog"])
body_sd  <- unname(body_fit$estimate["sdlog"])
body_fit

# ---- 2. Fit the GPD to the TAIL (exceedances over the threshold) ----
excess <- x[is_large] - threshold_th      # how far large claims exceed $10,000

# GPD density (standard "peaks over threshold" model):
#   xi > 0  -> heavy (Pareto-type) tail, the shape we expect from the mean
#              excess plot 
#   xi = 0  -> exponential tail
#   xi < 0  -> short, bounded tail 
gpd_negloglik <- function(par, y) {
  sigma <- par[1]; xi <- par[2]
  if (sigma <= 0) return(1e10)
  n <- length(y)
  if (abs(xi) < 1e-8) {
    return(-sum(-log(sigma) - y / sigma))
  }
  z <- 1 + xi * y / sigma
  if (any(z <= 0)) return(1e10)
  -(-n * log(sigma) - (1 / xi + 1) * sum(log(z)))
}

gpd_fit   <- optim(c(mean(excess), 0.1), gpd_negloglik, y = excess, method = "Nelder-Mead")
gpd_sigma <- gpd_fit$par[1]
gpd_xi    <- gpd_fit$par[2]
cat("GPD scale (sigma):", gpd_sigma, "\n")
cat("GPD shape (xi):   ", gpd_xi, "\n")   # xi > 0 confirms a heavy tail

write.csv(data.frame(model = c("LogNormal (body)", "GPD (tail)"),
                     param_1 = c(body_mu, gpd_sigma),
                     param_2 = c(body_sd, gpd_xi)),
          "results/tail_model_params.csv", row.names = FALSE)

# ---- 3. Check the GPD fit against the observed exceedances ----
qgpd <- function(p, sigma, xi) {
  if (abs(xi) < 1e-8) -sigma * log(1 - p) else sigma / xi * ((1 - p)^(-xi) - 1)
}
probs <- c(0.5, 0.9, 0.95, 0.99)
gpd_check <- data.frame(
  prob     = probs,
  observed = as.numeric(quantile(excess, probs)),
  gpd      = qgpd(probs, gpd_sigma, gpd_xi)
)
print(gpd_check)      # "observed" and "gpd" columns should be close
write.csv(gpd_check, "results/tail_gpd_quantile_check.csv", row.names = FALSE)

# QQ plot for the tail fit
p   <- ppoints(length(excess))
plot(qgpd(p, gpd_sigma, gpd_xi), sort(excess), log = "xy",
     main = "GPD fit to claims above $10,000",
     xlab = "model quantile (thousands over threshold)",
     ylab = "observed (thousands over threshold)")
abline(0, 1, col = "red")
dev.copy(png, "figures/tail_gpd_qq.png", width = 700, height = 500)
dev.off()

# ---- 4. Simulate individual claim sizes from the spliced model ----
rgpd <- function(n, sigma, xi) qgpd(runif(n), sigma, xi)

simulate_claim_sizes <- function(n) {
  large <- runif(n) < p_large
  out   <- numeric(n)
  if (sum(!large) > 0) out[!large] <- rlnorm(sum(!large), meanlog = body_mu, sdlog = body_sd)
  if (sum(large)  > 0) out[large]  <- threshold_th + rgpd(sum(large), gpd_sigma, gpd_xi)
  out * 1000   # back to dollars
}

# Sanity check: simulate many individual claims and compare to the real data
sim_claims <- simulate_claim_sizes(50000)
claim_check <- data.frame(
  stat     = c("mean", "p90", "p95", "p99", "p999"),
  observed = c(mean(sev$ClaimAmount), quantile(sev$ClaimAmount, c(.9, .95, .99, .999))),
  spliced  = c(mean(sim_claims), quantile(sim_claims, c(.9, .95, .99, .999)))
)
print(claim_check)
write.csv(claim_check, "results/tail_claim_size_check.csv", row.names = FALSE)
# The spliced model's mean and upper quantiles should now sit much closer to
# the observed values than the plain log-normal model.

# ---- 5. Rerun the aggregate loss simulation with the spliced severity ----
# Frequency parameters (negative binomial)
nb_theta <- 0.5703
nb_rate  <- 0.0746

simulate_one_year_tail <- function(n_policies, theta, rate) {
  N <- rnbinom(1, size = n_policies * theta, mu = n_policies * rate)
  if (N == 0) return(0)
  sum(simulate_claim_sizes(N))
}

set.seed(1234)
n_policies <- 677991
n_sims     <- 1000

S_tail <- replicate(n_sims, simulate_one_year_tail(n_policies, nb_theta, nb_rate))

sim_summary_tail <- data.frame(
  stat  = c("mean", "sd", "variance", "min", "median",
            "p90", "p95", "p99", "p995", "max"),
  value = c(mean(S_tail), sd(S_tail), var(S_tail), min(S_tail), median(S_tail),
            quantile(S_tail, 0.90), quantile(S_tail, 0.95), quantile(S_tail, 0.99),
            quantile(S_tail, 0.995), max(S_tail))
)
print(sim_summary_tail)
write.csv(sim_summary_tail, "results/sim_portfolio_tail_summary.csv", row.names = FALSE)

hist(S_tail, breaks = 40, main = "Simulated total loss, spliced (LogNormal + GPD) model",
     xlab = "Total loss (S)")
abline(v = mean(S_tail), col = "red", lwd = 2)
dev.copy(png, "figures/sim_portfolio_tail_hist.png", width = 800, height = 500)
dev.off()

# ---- 6. Compare: plain log-normal vs. spliced (LogNormal + GPD) ----
# Requires results/sim_portfolio_summary.csv 
plain <- tryCatch(read.csv("results/sim_portfolio_summary.csv"), error = function(e) NULL)
if (!is.null(plain)) {
  comparison <- data.frame(
    stat        = sim_summary_tail$stat,
    lognormal   = plain$value,
    spliced     = sim_summary_tail$value
  )
  comparison$pct_higher <- round((comparison$spliced - comparison$lognormal) /
                                    comparison$lognormal * 100, 2)
  print(comparison)
  write.csv(comparison, "results/sim_model_comparison.csv", row.names = FALSE)
} else {
  cat("Could not find results/sim_portfolio_summary.csv \n")
}
