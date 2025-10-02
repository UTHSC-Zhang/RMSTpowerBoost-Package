

if (!requireNamespace("RMSTpowerBoost", quietly = TRUE)) {
   stop("Please install/load the package that provides validate_recipe()/simulate_from_recipe().")
}
library(RMSTpowerBoost)

`%||%` <- function(x, y) if (is.null(x)) y else x

dir.create("checks", showWarnings = FALSE, recursive = TRUE)
dir.create("data",   showWarnings = FALSE, recursive = TRUE)

# --- Common covariate definitions ------------------------------------------------
covs_base <- list(
   list(name = "age", type = "continuous", dist = "normal",
        params = list(mean = 62, sd = 10),
        transform = c("center(60)", "scale(10)")),
   list(name = "sex", type = "categorical", dist = "bernoulli",
        params = list(p = 0.45))
)

# --- Scenario recipes (key â†’ list(seed, recipe)) --------------------------------
scenarios <- list(

   # AFT â€” lognormal
   aft_lognormal_L12_n150 = list(
      seed = 1001L,
      recipe = list(
         n = 150,
         covariates = list(defs = covs_base),
         treatment  = list(assignment = "randomization", allocation = "1:1"),
         event_time = list(
            model    = "aft_lognormal",
            baseline = list(mu = 2.8, sigma = 0.6),
            effects  = list(intercept = 0, treatment = -0.25, covariates = list(age = 0.01)),
            L        = 12
         ),
         censoring = list(mode = "target_overall", target = 0.25, admin_time = 30)
      )
   ),

   # AFT â€” Weibull
   aft_weibull_L24_n200 = list(
      seed = 1002L,
      recipe = list(
         n = 200,
         covariates = list(defs = covs_base),
         treatment  = list(assignment = "randomization", allocation = "1:1"),
         event_time = list(
            model    = "aft_weibull",
            baseline = list(shape = 1.3, scale = 12),
            effects  = list(intercept = 0, treatment = -0.25, covariates = list(age = 0.01)),
            L        = 24
         ),
         censoring = list(mode = "target_overall", target = 0.25, admin_time = 36)
      )
   ),

   # PH â€” Weibull
   ph_weibull_L24_n300 = list(
      seed = 1003L,
      recipe = list(
         n = 300,
         covariates = list(defs = covs_base),
         treatment  = list(assignment = "randomization", allocation = "1:1"),
         event_time = list(
            model    = "ph_weibull",
            baseline = list(shape = 1.2, scale = 10),
            effects  = list(intercept = 0, treatment = -0.30, covariates = list(age = 0.012)),
            L        = 24
         ),
         censoring = list(mode = "target_overall", target = 0.25, admin_time = 36)
      )
   ),

   # PH â€” piecewise exponential
   ph_pwexp_L18_n250 = list(
      seed = 1004L,
      recipe = list(
         n = 250,
         covariates = list(defs = covs_base),
         treatment  = list(assignment = "randomization", allocation = "1:1"),
         event_time = list(
            model    = "ph_pwexp",
            baseline = list(rates = c(0.08, 0.05, 0.03), cuts = c(6, 14)),
            effects  = list(intercept = 0, treatment = -0.35, covariates = list(age = 0.015)),
            L        = 18
         ),
         censoring = list(mode = "target_overall", target = 0.25, admin_time = 30)
      )
   )
)

# --- Helper to collect rich metadata from data + recipe -------------------------
.collect_meta <- function(data, recipe, files, seed) {
   has_arm <- "arm" %in% names(data)
   Lval <- attr(data, "L") %||% attr(data, "tau") %||% (recipe$event_time$L %||% NA_real_)
   ac   <- attr(data, "achieved_censoring")
   if (is.null(ac) || is.na(ac)) ac <- mean(data$status == 0, na.rm = TRUE)
   list(
      n                   = nrow(data),
      L                   = Lval,
      event_rate          = mean(data$status == 1, na.rm = TRUE),
      achieved_censoring  = ac,
      n_treat             = if (has_arm) sum(data$arm == 1L) else NA_integer_,
      n_control           = if (has_arm) sum(data$arm == 0L) else NA_integer_,
      model               = recipe$event_time$model,
      baseline            = recipe$event_time$baseline,
      effects             = recipe$event_time$effects,
      treatment           = recipe$treatment %||% NULL,
      censoring           = recipe$censoring %||% NULL,
      covariates          = if (!is.null(recipe$covariates$defs))
         lapply(recipe$covariates$defs, function(d) d[c("name","type","dist","params","transform")])
      else NULL,
      seed_used           = seed,
      files               = files,
      created_at          = as.character(Sys.time())
   )
}

# --- Generate data, write CSV/RDA, build manifest --------------------------------
rows  <- list()
metas <- list()

for (key in names(scenarios)) {
   sc  <- scenarios[[key]]
   rec <- validate_recipe(sc$recipe)
   dat <- simulate_from_recipe(rec, seed = sc$seed)

   # (1) CSV â†’ checks/
   csv_path <- file.path("checks", paste0(key, ".csv"))
   utils::write.csv(dat, csv_path, row.names = FALSE)

   # (2) RDA â†’ data/  (use base save; ensure object exists by that name)
   rda_path <- file.path("data", paste0(key, ".rda"))
   assign(key, dat, envir = environment())
   save(list = key, file = rda_path, compress = "xz", envir = environment())

   # (3) Metadata
   meta <- .collect_meta(data = dat, recipe = rec, files = list(csv = csv_path, rda = rda_path), seed = sc$seed)
   metas[[key]] <- meta

   rows[[length(rows) + 1L]] <- data.frame(
      key = key,
      n   = meta$n,
      L   = meta$L,
      event_rate = meta$event_rate,
      achieved_censoring = meta$achieved_censoring,
      file_csv = csv_path,
      file_rda = rda_path,
      stringsAsFactors = FALSE
   )
}

manifest <- if (length(rows)) do.call(rbind, rows) else data.frame()
manifest$meta <- unname(metas)

saveRDS(manifest, file.path("checks", "manifest.rds"))
utils::write.csv(
   transform(manifest, meta = NULL),  # human-readable
   file.path("checks", "manifest.csv"),
   row.names = FALSE
)

cat(
   "Done.\n",
   "CSV files â†’ ", normalizePath("checks"), "\n",
   "RDA files â†’ ", normalizePath("data"), "\n",
   "Manifest   â†’ ", normalizePath(file.path("checks", "manifest.rds")), "\n",
   sep = ""
)
