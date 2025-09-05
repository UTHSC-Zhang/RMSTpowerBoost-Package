# data-raw/make_datasets.R
# Generate small demo datasets for inclusion. List-only v2 API.

# ensure package functions are available
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
}

dir.create("data-raw", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

covs_base <- list(
  list(name="age", type="continuous", dist="normal",
       params=list(mean=62, sd=10),
       transform=c("center(60)","scale(10)")),
  list(name="sex", type="categorical", dist="bernoulli",
       params=list(p=0.45))
)

# AFT–Weibull (n=200)
rec_weib <- recipe_quick_aft(
  n = 200, tau = 24, model = "aft_weibull",
  baseline = list(shape = 1.3, scale = 12),
  treat_effect = -0.25, covariates = covs_base,
  target_censoring = 0.25, allocation = "1:1", seed = 2025
)
rmst_weibull_200 <- simulate_from_recipe(rec_weib, seed = 2025)

# PH–Weibull with gamma frailty by 3‑level cluster (n=300)
covs_f <- c(covs_base, list(
  list(name="cluster", type="categorical", dist="categorical",
       params=list(prob=c(0.3,0.4,0.3), labels=c("c1","c2","c3")))
))
rec_phw <- list(
  n = 300, covariates = list(defs = covs_f),
  treatment = list(assignment="randomization", allocation="1:1"),
  event_time = list(model="ph_weibull",
                    baseline=list(shape=1.2, scale=10),
                    effects=list(intercept=0, treatment=-0.3, covariates=c(age=0.01)),
                    frailty=list(type="gamma", var=0.5, group="cluster"),
                    tau=24),
  censoring = list(mode="target_overall", target=0.25, admin_time=36)
)
rec_phw <- validate_recipe(rec_phw)
ph_weibull_frailty_300 <- simulate_from_recipe(rec_phw, seed = 31415)

# AFT–loglogistic small (n=100)
rec_ll <- list(
  n = 100, covariates = list(defs = covs_base),
  treatment = list(assignment="randomization", allocation="1:1"),
  event_time = list(model="aft_loglogistic",
                    baseline=list(shape=1.6, scale=10),
                    effects=list(intercept=0, treatment=-0.2, covariates=c(age=0.01)),
                    tau=12),
  censoring = list(mode="target_overall", target=0.20, admin_time=24)
)
rec_ll <- validate_recipe(rec_ll)
aft_loglogistic_100 <- simulate_from_recipe(rec_ll, seed = 7)

# Save .rda objects (for package)
if (!requireNamespace("usethis", quietly = TRUE)) {
  stop("Please install 'usethis' to save package data: install.packages('usethis')")
}
usethis::use_data(rmst_weibull_200, ph_weibull_frailty_300, aft_loglogistic_100,
                  overwrite = TRUE, compress = "xz")

# Also write CSV/TXT mirrors to inst/extdata
write_one <- function(df, stem) {
  utils::write.csv(df, file.path("inst/extdata", paste0(stem, ".csv")), row.names = FALSE)
  utils::write.table(df, file.path("inst/extdata", paste0(stem, ".txt")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}
write_one(rmst_weibull_200,     "rmst_weibull_200")
write_one(ph_weibull_frailty_300,  "ph_weibull_frailty_300")
write_one(aft_loglogistic_100,     "aft_loglogistic_100")

message("Done. Saved data/*.rda and inst/extdata/*.csv/*.txt (v2).")
