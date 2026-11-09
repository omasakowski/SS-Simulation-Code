# ---- 0. Setup (run the install lines ONCE, then comment them out) ----
# install.packages(c("xts", "sp", "zoo", "MASS"))
# install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type = "source")
# (If that repo fails, the package's GitHub page lists alternate repos.)

library(CASdatasets)
library(MASS)          # for glm.nb (negative binomial)

data(freMTPL2freq)     # one row per policy: ClaimNb, Exposure, driver/car info
data(freMTPL2sev)      # one row per claim: IDpol, ClaimAmount

freq <- freMTPL2freq
sev  <- freMTPL2sev

# In some versions of the package ClaimNb is stored as class "table";
# converting to a plain number avoids odd errors later.
freq$ClaimNb <- as.numeric(freq$ClaimNb)

# Folders for saving results (created if they don't exist)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Check the data ----
dim(freq)                             # roughly 678,000 policies
dim(sev)                              # roughly 26,000 claim records
summary(freq$Exposure)                # look at the min and max: any exposure over 1 year?
table(freq$ClaimNb)                   # look at the largest claim counts: any outliers?

sec1 <- data.frame(
  item  = c("policies in freq", "records in sev", "sum of ClaimNb",
            "policies with 0 claims", "policies with 3+ claims", "max ClaimNb",
            "min exposure", "median exposure", "mean exposure", "max exposure",
            "policies with exposure > 1"),
  value = c(nrow(freq), nrow(sev), sum(freq$ClaimNb),
            sum(freq$ClaimNb == 0), sum(freq$ClaimNb >= 3), max(freq$ClaimNb),
            min(freq$Exposure), median(freq$Exposure), mean(freq$Exposure),
            max(freq$Exposure), sum(freq$Exposure > 1))
)
write.csv(sec1, "results/section1_data_summary.csv", row.names = FALSE)

claim_counts <- as.data.frame(table(ClaimNb = freq$ClaimNb))
write.csv(claim_counts, "results/claim_count_table.csv", row.names = FALSE)

# just look at these (no need to save): the policies with the most claims
freq[freq$ClaimNb >= 5, c("IDpol", "ClaimNb", "Exposure")]

# ---- 2. Do the two files agree? (known issue with freMTPL2) ----
sum(freq$ClaimNb)                     # claims counted in the frequency file
nrow(sev)                             # claims that have an amount in the severity file
# (these will NOT match, and that is expected: write down both numbers)

idx <- match(sev$IDpol, freq$IDpol)   # which policy each severity record belongs to
freq$SevNb <- tabulate(idx, nbins = nrow(freq))   # claims per policy per the severity file
sum(is.na(idx))                       # severity records with no matching policy
sum(freq$ClaimNb != freq$SevNb)       # policies where the two files disagree
sum(freq$ClaimNb > 0 & freq$SevNb == 0)  # counted claims with no amount on file
# Do NOT delete these policies for now: doing so would remove many policies that
# had claims and bias the claim rate downward.

checks <- data.frame(
  item  = c("policies in freq", "records in sev", "sum of ClaimNb",
            "sev records with no matching policy",
            "policies where files disagree",
            "counted claims with no amount"),
  value = c(nrow(freq), nrow(sev), sum(freq$ClaimNb),
            sum(is.na(idx)),
            sum(freq$ClaimNb != freq$SevNb),
            sum(freq$ClaimNb > 0 & freq$SevNb == 0))
)
write.csv(checks, "results/section2_file_agreement.csv", row.names = FALSE)

# ---- 3. Annual claim rate and dispersion ----
Y <- freq$ClaimNb                     # claims per policy
E <- freq$Exposure                    # years observed per policy

m   <- sum(Y) / sum(E)                # claims per policy-year
v   <- sum((Y - m * E)^2) / sum(E)    # exposure-weighted variance
phi <- v / m                          # phi = 1 means exactly Poisson-like
cat("Annual claim rate =", m, "\n")
cat("Dispersion phi    =", phi, "\n")

write.csv(data.frame(item  = c("annual claim rate", "dispersion phi"),
                     value = c(m, phi)),
          "results/rate_and_dispersion.csv", row.names = FALSE)

# ---- 4. Fit Poisson and negative binomial (no covariates yet) ----
# offset(log(Exposure)) accounts for policies observed for less than a year
fit_pois <- glm(ClaimNb ~ 1 + offset(log(Exposure)),
                family = poisson, data = freq)
fit_nb   <- glm.nb(ClaimNb ~ 1 + offset(log(Exposure)),
                   data = freq)

summary(fit_pois)
summary(fit_nb)
AIC(fit_pois, fit_nb)                 # lower = better fit (after penalizing complexity)

b <- unname(c(coef(fit_pois), coef(fit_nb)))
frequency_models <- data.frame(
  model       = c("Poisson", "NegBin"),
  intercept   = b,
  annual_rate = exp(b),
  theta       = c(NA, fit_nb$theta),
  AIC         = c(AIC(fit_pois), AIC(fit_nb))
)
write.csv(frequency_models, "results/frequency_models.csv", row.names = FALSE)

# ---- 5. Observed vs expected number of policies with 0, 1, 2, 3+ claims ----
mu_p  <- fitted(fit_pois)             # expected claims per policy under Poisson
mu_n  <- fitted(fit_nb)               # expected claims per policy under neg. binomial
theta <- fit_nb$theta                 # neg. binomial dispersion parameter

ks  <- 0:2
obs <- c(sapply(ks, function(k) sum(Y == k)), sum(Y >= 3))

exp_p <- sapply(ks, function(k) sum(dpois(k, mu_p)))
exp_p <- c(exp_p, length(Y) - sum(exp_p))              # 3+ = everything left over

exp_n <- sapply(ks, function(k) sum(dnbinom(k, size = theta, mu = mu_n)))
exp_n <- c(exp_n, length(Y) - sum(exp_n))

freq_table <- data.frame(claims   = c("0", "1", "2", "3+"),
                         observed = obs,
                         poisson  = round(exp_p, 1),
                         negbin   = round(exp_n, 1))
print(freq_table)

# ---- 6. Save results  ----
write.csv(freq_table, "results/frequency_observed_vs_expected.csv", row.names = FALSE)
saveRDS(list(rate = m, phi = phi, theta = theta,
             aic_pois = AIC(fit_pois), aic_nb = AIC(fit_nb),
             n_claims_freq = sum(Y), n_claims_sev = nrow(sev)),
        "results/frequency_summary.rds")
