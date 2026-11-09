library(MASS)
options(scipen = 999)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

set.seed(1234)   # makes the random simulation reproducible: same numbers every run

# ---- 0. Plug in the fitted parameters  ----

# Frequency: negative binomial, fit on ClaimNb ~ offset(log(Exposure))
nb_theta  <- 0.5703          # from results/frequency_models.csv
nb_rate   <- 0.0746          # annual claims per policy-year (NegBin fitted rate)

# Severity: log-normal, fit on ClaimAmount / 1000 (thousands)
ln_mu <- -0.06107862          # from results/severity_models.csv
ln_sd <- 1.13000160            # from results/severity_models.csv

# ---- 1. Simulate one year for a SMALL block of policies ----
# A block of policies, each with 1 year of exposure, so total exposure = n_policies.
n_policies <- 677991
n_sims     <- 1000           # number of times we repeat the whole year

# glm.nb parameterizes the negative binomial by mu (mean) and theta (size).
# For n_policies with 1 year of exposure each, total expected claims = n_policies * nb_rate.
mu_block <- n_policies * nb_rate

simulate_one_year <- function(mu, theta, rate, ln_mu, ln_sd) {
  N <- rnbinom(1, size = n_policies * theta, mu = n_policies * rate)      # total claims for the block this year
  if (N == 0) return(0)
  X <- rlnorm(N, meanlog = ln_mu, sdlog = ln_sd) * 1000   # claim sizes, back in original units
  sum(X)
}

S <- replicate(n_sims, simulate_one_year(n_policies, nb_theta, nb_rate, ln_mu, ln_sd))

# ---- 2. Summarize the simulated aggregate loss ----
sim_summary <- data.frame(
  stat  = c("mean", "sd", "variance", "min", "median",
            "p90", "p95", "p99", "p995", "max"),
  value = c(mean(S), sd(S), var(S), min(S), median(S),
            quantile(S, 0.90), quantile(S, 0.95), quantile(S, 0.99),
            quantile(S, 0.995), max(S))
)
print(sim_summary)
write.csv(sim_summary, "results/sim_portfolio_summary.csv", row.names = FALSE)

# ---- 3. Check the simulation against the theoretical formulas ----
# E(S)   = E(N) * E(X)
# Var(S) = E(N) * Var(X) + Var(N) * E(X)^2
EX   <- exp(ln_mu + ln_sd^2 / 2) * 1000                      # mean claim size
VarX <- (exp(ln_sd^2) - 1) * exp(2 * ln_mu + ln_sd^2) * 1000^2  # variance of claim size

EN   <- n_policies * nb_rate                         # mean number of claims
VarN <- EN + EN^2 / (n_policies * nb_theta)       # negative binomial variance

E_S_theory   <- EN * EX
Var_S_theory <- EN * VarX + VarN * EX^2

check <- data.frame(
  quantity = c("E(S)", "Var(S)"),
  simulated   = c(mean(S), var(S)),
  theoretical = c(E_S_theory, Var_S_theory)
)
check$pct_diff <- round((check$simulated - check$theoretical) / check$theoretical * 100, 2)
print(check)
write.csv(check, "results/sim_portfolio_check.csv", row.names = FALSE)
# These should be close (within a few percent). A big gap usually means a
# units mistake (thousands vs. dollars) or a typo in the parameters above.

# ---- 4. Look at the simulated distribution ----
hist(S, breaks = 60, main = "Simulated total loss, full portfolio (Log-Normal)",
     xlab = "Total loss (S)")
abline(v = mean(S), col = "red", lwd = 2)
dev.copy(png, "figures/sim_portfolio_hist.png", width = 800, height = 500)
dev.off()

# ---- 5. Scale up to the whole portfolio ----
# Once the block-level simulation looks right, re-run with:
#   n_policies <- 677991      # or nrow(freq) if freq is loaded
#   n_sims     <- 1000        # fewer sims needed since each one is more stable
# A portfolio-wide run can take longer, since each simulated year draws many
# more claims. Time a single run first with system.time(simulate_one_year(...))
# before committing to 1000+ repeats.
