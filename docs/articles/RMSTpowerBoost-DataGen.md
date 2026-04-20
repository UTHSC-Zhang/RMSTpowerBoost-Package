# Data Generation: Mechanisms and Examples

## Quick start with `rmst.sim()`

For common settings — AFT or PH baseline, random censoring, and standard
covariate distributions —
[`rmst.sim()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md)
generates data from one function call. Use
[`covar_continuous()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md),
[`covar_binary()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_binary.md),
and
[`covar_categorical()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_categorical.md)
to define covariates directly.

``` r
df <- rmst.sim(
  n              = 200,
  model          = "aft_lognormal",
  baseline       = list(mu = 2.5, sigma = 0.6),
  treat_effect   = -0.25,
  covariates     = list(
    covar_continuous("age",  mean = 62, sd = 10),
    covar_binary("sex",      p = 0.45),
    covar_continuous("x",    dist = "lognormal", meanlog = 0, sdlog = 0.5)
  ),
  target_censoring = 0.25,
  L    = 36,
  seed = 42
)

print(df)
── Simulated RMST Dataset ────────────────────────────────────
  Model          : AFT Log-Normal 
  N (total)      : 200  |  Arm 0: 104  |  Arm 1: 96
  Events         : 147 (73.5%)
  Censored       : 53 (26.5%)
  Truncation time: 36  (stored as attribute)

First 6 rows:
       time status arm      age sex         x
1  8.173606      1   1 75.70958   0 0.9976923
2 12.521599      0   0 56.35302   0 1.4624617
3 17.160187      1   1 65.63128   1 1.0196867
4 15.663621      1   1 68.32863   0 1.4441719
5  6.677792      0   0 66.04268   1 0.9293812
6 16.085942      0   0 60.93875   0 0.9714712
```

The result is a data frame of class `rmst_simdata`. Pass it directly to
[`rmst.power()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md)
or
[`rmst.ss()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md):

``` r

rmst.power(Surv(time, status) ~ age + sex + x,
           data = df, arm = "arm",
           sample_sizes = c(100, 200, 300), L = 36)
```

For designs requiring explicit covariate-dependent censoring or
stratified assignment, use the full recipe system described below.

------------------------------------------------------------------------

## What this covers

We show how to construct a **recipe** (covariates (X), treatment (A),
event time (Y), and censoring (C)), and then simulate datasets with
[`simulate_from_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).
We also show batch generation via
[`generate_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/generate_recipe_sets.md)
and how to read compact metadata back with
[`load_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/load_recipe_sets.md).

The simulator takes a plain **list** as the recipe. No YAML is required.

### Recipe skeleton

``` r
list(
  n = 300,
  covariates = list(defs = list(/* see Covariates */)),
  treatment  = list(/* see Treatment */),
  event_time = list(/* see Event-time engines */),
  censoring  = list(/* see Censoring */),
  seed = 42
)
```

------------------------------------------------------------------------

## Covariates

Each covariate has a `name`, `type` (`"continuous"` or `"categorical"`),
a `dist` with `params`, and optional `transform` steps applied after
generation (e.g., `"center(a)"`, `"scale(b)"`).

**Available distributions**

| Type        | `dist` (string) | Parameters                                |
|-------------|-----------------|-------------------------------------------|
| continuous  | `normal`        | `mean`, `sd`                              |
| continuous  | `lognormal`     | `meanlog`, `sdlog`                        |
| continuous  | `gamma`         | `shape`, `scale`                          |
| continuous  | `weibull`       | `shape`, `scale`                          |
| continuous  | `uniform`       | `min`, `max`                              |
| continuous  | `beta`          | `shape1`, `shape2`                        |
| continuous  | `t`             | `df`                                      |
| categorical | `bernoulli`     | `p` (probability of 1)                    |
| categorical | `categorical`   | `prob = c(...)`, `labels = c(...)` (opt.) |
| categorical | `ordinal`       | `prob = c(...)`, `labels = c(...)` (opt.) |

**Example: define covariates**

``` r

covs <- list(
  list(name="age",   type="continuous",  dist="normal",     params=list(mean=62, sd=10),
       transform=c("center(60)","scale(10)")),
  list(name="sex",   type="categorical", dist="bernoulli",  params=list(p=0.45)),
  list(name="stage", type="categorical", dist="ordinal",
       params=list(prob=c(0.3,0.5,0.2), labels=c("I","II","III"))),
  list(name="x",     type="continuous",  dist="lognormal",  params=list(meanlog=0, sdlog=0.6))
)
```

------------------------------------------------------------------------

## Treatment

Choose one `assignment`:

| Assignment | Key fields | Meaning |
|----|----|----|
| `"randomization"` | `allocation = "a:b"` | Bernoulli with probability (p_1 = a/(a+b)). |
| `"stratified"` | `allocation`, `stratify_by = c(...)` | Same allocation **within** each stratum defined by listed **categorical** covariates. |
| `"logistic_ps"` | `ps_model = list(formula, beta)` | Treatment probability is `plogis(eta)` from the user-specified model. |

**Examples**

Randomization:

``` r

tr_rand <- list(assignment="randomization", allocation="1:1")
```

Stratified by `"stage"`:

``` r

tr_strat <- list(assignment="stratified", allocation="2:1", stratify_by=c("stage"))
```

Logistic propensity:

``` r

tr_ps <- list(
  assignment = "logistic_ps",
  ps_model  = list(
    formula = "~ 1 + x + sex",
    beta    = c(-0.3, 1.2, -0.6)  # (Intercept), x, sex
  )
)
```

------------------------------------------------------------------------

## Event-time engines

Let `Z` be treatment (0/1), `X` be covariates, and `eta` be the **linear
predictor** (defined in `effects`, below). Supported engines and
baseline parameterizations:

| Model (user-facing) | `model` value | Baseline parameters | Notes |
|----|----|----|----|
| AFT Lognormal | `"aft_lognormal"` | `mu`, `sigma` | `log(T) = mu + eta + sigma * epsilon`, with `epsilon ~ N(0, 1)`. |
| AFT Weibull | `"aft_weibull"` | `shape`, `scale` | Baseline survival `S0(t) = exp(-(t/scale)^shape)` with an AFT shift via `eta`. |
| AFT Log-Logistic | `"aft_loglogistic"` | `shape`, `scale` | Log-logistic AFT model with `T = scale * exp(eta) * (U/(1-U))^(1/shape)`. |
| PH Exponential | `"ph_pwexp"` | `rates = c(rate)`, `cuts = numeric(0)` | Piecewise exponential with a **single** segment is exponential. |
| PH Weibull | `"ph_weibull"` | `shape`, `scale` | Proportional hazards with a Weibull baseline. |
| PH Gompertz | `"ph_gompertz"` | `rate`, `gamma` | Hazard `h(t) = rate * exp(gamma * t)`. |
| PH Piecewise Exp | `"ph_pwexp"` | `rates = c(r1,r2,...)`, `cuts = c(...)` | Segment-specific hazard is `r_s * exp(eta)`. |

**Effects and linear predictor**

Specify effects on the appropriate scale (AFT: log-time; PH:
log-hazard):

``` r

effects = list(
  intercept  = 0,                      # default is 0
  treatment  = -0.25,
  covariates = list(age = 0.01, sex = -0.2)
  # or: formula="~ age + sex", beta=c(0.01, -0.2)
)
```

> `effects$covariates` must be a **named list** (e.g.,
> `list(age=0.01)`), not a named vector from
> [`c()`](https://rdrr.io/r/base/c.html).

------------------------------------------------------------------------

## Censoring

We assume a **single censoring mechanism** that may depend on observed
covariates (but not on the unobserved event time after conditioning on
those covariates). Two user-facing modes are supported:

| Mode | Fields | Semantics |
|----|----|----|
| `"target_overall"` | `target`, `admin_time` | Finds an exponential random-censoring rate so the overall censoring fraction is approximately `target`, with an administrative floor at `admin_time` (if provided). |
| `"explicit"` | Any of: `administrative = list(time=...)`; `random = list(dist="exponential", params=list(rate=...))`; `dependent = list(formula="~ ...", base=..., beta=c(...))` | Compose administrative, random, and **covariate-dependent** censoring directly. The `dependent` block defines a covariate-dependent **log-rate** model; it does **not** denote competing risks. |

> **Note.** No cause-specific or competing-risk censoring is generated
> here; the simulator’s `dependent` option is strictly
> **covariate-dependent** censoring under a single mechanism.

### Examples

#### Target overall censoring

``` r

cz_target <- list(
  mode        = "target_overall",
  target      = 0.25,   # desired overall censoring fraction
  admin_time  = 36      # optional administrative cutoff (same units as event time)
)
```

#### Explicit mix (administrative + random)

``` r

cz_explicit <- list(
  mode = "explicit",
  administrative = list(time = 36),                                      # admin cutoff
  random         = list(dist = "exponential", params = list(rate = 0.02))# random censoring
)
```

#### Explicit **covariate-dependent** censoring

The dependent block uses a **log-rate** model:
``` math
\log \lambda_c(X) = \log(b) + X^\top \beta,
```
where `b > 0` is the *baseline rate*. Implementation uses `base` as a
**rate**, and `beta` on the **log-rate** scale. When `base` is supplied,
the **formula must not include an intercept**.

``` r

cz_dep <- list(
  mode = "explicit",
  dependent = list(
    # sex is coded numeric 0/1 in these recipes; if yours is a factor, use I(as.numeric(sex)-1)
    formula = "~ 0 + age + sex",   # NO intercept because we use 'base' as the rate
    base    = 0.03,                # positive baseline rate
    beta    = c(0.012, 0.35)       # (age, sex) on log-rate scale
  ),
  administrative = list(time = 48) # optional admin cutoff
)
```

> If you prefer a log-intercept style, omit `base` and include an
> intercept in the formula and in `beta`,
> e.g. `formula="~ 1 + age + sex"; beta=c(-3.5, 0.012, 0.35)`.

------------------------------------------------------------------------

## Using a censoring recipe inside a full simulation

``` r

# Covariates (sex will be numeric 0/1 here)
covs <- list(
  list(name="age", type="continuous",  dist="normal",    params=list(mean=62, sd=10)),
  list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.45))
)

# Censoring: single covariate-dependent mechanism using a positive base rate
cz_dep <- list(
  mode = "explicit",
  dependent = list(
    formula = "~ 0 + age + sex",  # NO intercept; sex is numeric 0/1
    base    = 0.03,               # baseline censoring rate (>0)
    beta    = c(0.012, 0.35)      # (age, sex) on log-rate scale
  ),
  administrative = list(time = 48)
)

# Quick sanity check (optional)
mm   <- model.matrix(~ 0 + age + sex, data.frame(age=c(50,70), sex=c(0,1)))
rate <- cz_dep$dependent$base * exp(mm %*% matrix(cz_dep$dependent$beta, ncol=1))
stopifnot(all(rate > 0 & is.finite(rate)))

# Recipe
rec <- list(
  n = 300,
  covariates = list(defs = covs),
  treatment  = list(assignment="randomization", allocation="1:1"),
  event_time = list(
    model    = "ph_pwexp",
    baseline = list(rates = c(0.05), cuts = numeric(0)), # exponential baseline
    effects  = list(intercept=0, treatment=-0.3, covariates=list(age=0.01, sex=-0.2))
  ),
  censoring  = cz_dep,
  seed = 2025
)

# Simulate
dat <- simulate_from_recipe(validate_recipe(rec))

# Output
knitr::kable(
  data.frame(achieved_censoring = attr(dat, "achieved_censoring")),
  caption = "Achieved censoring rate"
)
```

| achieved_censoring |
|-------------------:|
|               0.51 |

Achieved censoring rate {.table}

``` r

knitr::kable(utils::head(dat, 8), caption = "First 8 rows of simulated data")
```

|       time | status | arm |      age | sex |
|-----------:|-------:|----:|---------:|----:|
|  0.0749714 |      0 |   0 | 68.20757 |   0 |
|  1.1019457 |      0 |   0 | 62.35641 |   1 |
| 14.3347821 |      0 |   1 | 69.73154 |   0 |
|  1.0472994 |      1 |   1 | 74.72489 |   0 |
|  1.7066127 |      1 |   0 | 65.70975 |   0 |
|  0.8646883 |      1 |   0 | 60.37146 |   0 |
|  0.1501047 |      1 |   1 | 65.97112 |   0 |
|  5.0730373 |      1 |   1 | 61.20011 |   0 |

First 8 rows of simulated data {.table}

------------------------------------------------------------------------

## Worked examples

We now build full recipes and call
[`simulate_from_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).
We report realized censoring via `attr(dat, "achieved_censoring")`.

### Example 1 — AFT Lognormal

``` r

# Covariates
covs1 <- list(
  list(name="age",   type="continuous",  dist="normal",     params=list(mean=62, sd=10),
       transform=c("center(60)","scale(10)")),
  list(name="sex",   type="categorical", dist="bernoulli",  params=list(p=0.45)),
  list(name="stage", type="categorical", dist="ordinal",
       params=list(prob=c(0.3,0.5,0.2), labels=c("I","II","III"))),
  list(name="x",     type="continuous",  dist="lognormal",  params=list(meanlog=0, sdlog=0.6))
)

# Recipe
rec1 <- list(
  n = 300,
  covariates = list(defs = covs1),
  treatment  = list(assignment="randomization", allocation="1:1"),
  event_time = list(model="aft_lognormal",
                    baseline=list(mu=3.0, sigma=0.6),
                    effects=list(intercept=0, treatment=-0.25,
                                 covariates=list(age=0.01, sex=-0.2, x=0.05))),
  censoring  = list(mode="target_overall", target=0.25, admin_time=36),
  seed = 11
)

# Simulate
dat1 <- simulate_from_recipe(validate_recipe(rec1))

# Display: data (first rows)
knitr::kable(utils::head(dat1, 8), caption = "Example 1 — AFT Lognormal: first 8 rows")
```

|      time | status | arm |        age | sex | stage |         x |
|----------:|-------:|----:|-----------:|----:|:------|----------:|
|  8.929004 |      0 |   0 | -0.3910311 |   0 | II    | 1.2192376 |
| 19.531620 |      1 |   1 |  0.2265944 |   1 | II    | 0.6122083 |
| 11.902048 |      1 |   0 | -1.3165531 |   0 | III   | 1.2334828 |
| 18.770511 |      1 |   0 | -1.1626533 |   0 | II    | 0.5784534 |
| 16.584294 |      1 |   1 |  1.3784892 |   0 | II    | 0.2619596 |
|  9.377759 |      1 |   1 | -0.7341513 |   0 | II    | 1.1455235 |
| 12.079110 |      1 |   0 |  1.5236056 |   1 | III   | 1.7865600 |
| 17.642809 |      0 |   1 |  0.8249178 |   1 | II    | 1.4493842 |

Example 1 — AFT Lognormal: first 8 rows {.table}

``` r


# Display: key attributes (compact)
attr_tbl1 <- data.frame(
  attribute = c("achieved_censoring", "n_rows", "n_cols", "class"),
  value = c(
    as.character(attr(dat1, "achieved_censoring")),
    nrow(dat1),
    ncol(dat1),
    paste(class(dat1), collapse = ", ")
  )
)
knitr::kable(attr_tbl1, caption = "Example 1 — Attributes")
```

| attribute          | value             |
|:-------------------|:------------------|
| achieved_censoring | 0.226666666666667 |
| n_rows             | 300               |
| n_cols             | 7                 |
| class              | data.frame        |

Example 1 — Attributes {.table}

### Example 2 — AFT Weibull

``` r

# Start from Example 1 and tweak the event-time engine
rec2 <- rec1
rec2$event_time <- list(
  model    = "aft_weibull",
  baseline = list(shape=1.3, scale=12),
  effects  = list(intercept=0, treatment=-0.20, covariates=list(age=0.008, x=0.04))
)

dat2 <- simulate_from_recipe(validate_recipe(rec2), seed=12)

# Display: data (first rows)
knitr::kable(utils::head(dat2, 8), caption = "Example 2 — AFT Weibull: first 8 rows")
```

|       time | status | arm |        age | sex | stage |         x |
|-----------:|-------:|----:|-----------:|----:|:------|----------:|
|  8.0059365 |      1 |   0 | -1.2805676 |   0 | I     | 0.4280147 |
| 22.1334913 |      1 |   0 |  1.7771695 |   0 | II    | 2.1502917 |
|  3.1817150 |      1 |   0 | -0.7567445 |   0 | II    | 0.9171049 |
|  0.8247902 |      1 |   1 | -0.7200052 |   1 | II    | 1.5528655 |
|  3.9411737 |      1 |   1 | -1.7976421 |   0 | II    | 0.6726879 |
| 12.4980918 |      1 |   1 | -0.0722960 |   0 | II    | 2.2757744 |
|  5.7430478 |      0 |   1 | -0.1153487 |   0 | I     | 0.3663151 |
|  1.0425923 |      1 |   1 | -0.4282552 |   0 | III   | 0.9247253 |

Example 2 — AFT Weibull: first 8 rows {.table}

``` r


# Display: key attributes (compact)
attr_tbl2 <- data.frame(
  attribute = c("achieved_censoring", "n_rows", "n_cols", "class"),
  value = c(
    as.character(attr(dat2, "achieved_censoring")),
    nrow(dat2),
    ncol(dat2),
    paste(class(dat2), collapse = ", ")
  )
)
knitr::kable(attr_tbl2, caption = "Example 2 — Attributes")
```

| attribute          | value             |
|:-------------------|:------------------|
| achieved_censoring | 0.233333333333333 |
| n_rows             | 300               |
| n_cols             | 7                 |
| class              | data.frame        |

Example 2 — Attributes {.table}

### Example 3 — PH Exponential (single segment)

``` r

rec3 <- list(
  n = 400,
  covariates = list(defs = covs1),
  treatment  = list(assignment="randomization", allocation="1:1"),
  event_time = list(
    model    = "ph_pwexp",
    baseline = list(rates = c(0.05), cuts = numeric(0)),
    effects  = list(intercept=0, treatment=-0.3, covariates=list(age=0.01, x=0.03))
  ),
  censoring  = list(mode="target_overall", target=0.20, admin_time=30),
  seed = 13
)

dat3 <- simulate_from_recipe(validate_recipe(rec3))

# Display: data (first rows)
knitr::kable(utils::head(dat3, 8), caption = "Example 3 — PH Exponential: first 8 rows")
```

|     time | status | arm |        age | sex | stage |         x |
|---------:|-------:|----:|-----------:|----:|:------|----------:|
| 30.00000 |      0 |   0 |  0.7543269 |   1 | III   | 2.4109639 |
| 12.97451 |      1 |   1 | -0.0802719 |   1 | II    | 1.0548136 |
| 30.00000 |      0 |   0 |  1.9751634 |   1 | I     | 0.3014883 |
| 10.72582 |      1 |   1 |  0.3873201 |   0 | I     | 2.3418948 |
| 16.35756 |      1 |   1 |  1.3425261 |   1 | I     | 1.8945147 |
| 26.61325 |      1 |   0 |  0.6155261 |   1 | III   | 1.0956426 |
| 30.00000 |      0 |   1 |  1.4295066 |   0 | II    | 1.0172397 |
| 30.00000 |      0 |   1 |  0.4366797 |   0 | II    | 1.9923261 |

Example 3 — PH Exponential: first 8 rows {.table}

``` r


# Display: a small summary + attributes
summ3 <- data.frame(
  Metric = c("n", "Events", "Censoring rate", "Mean time", "Median time"),
  Value  = c(nrow(dat3),
             sum(dat3$status == 1, na.rm = TRUE),
             round(mean(dat3$status == 0, na.rm = TRUE), 3),
             round(mean(dat3$time,   na.rm = TRUE), 2),
             round(stats::median(dat3$time, na.rm = TRUE), 2))
)
knitr::kable(summ3, caption = "Example 3 — Summary")
```

| Metric         |   Value |
|:---------------|--------:|
| n              | 400.000 |
| Events         | 297.000 |
| Censoring rate |   0.258 |
| Mean time      |  16.440 |
| Median time    |  15.130 |

Example 3 — Summary {.table}

``` r


attr_tbl3 <- data.frame(
  attribute = c("achieved_censoring", "n_rows", "n_cols", "class"),
  value = c(
    as.character(attr(dat3, "achieved_censoring")),
    nrow(dat3),
    ncol(dat3),
    paste(class(dat3), collapse = ", ")
  )
)
knitr::kable(attr_tbl3, caption = "Example 3 — Attributes")
```

| attribute          | value      |
|:-------------------|:-----------|
| achieved_censoring | 0.2575     |
| n_rows             | 400        |
| n_cols             | 7          |
| class              | data.frame |

Example 3 — Attributes {.table}

### Example 4 — PH Piecewise Exponential (multi-segment)

``` r

rec4 <- list(
  n = 500,
  covariates = list(defs = list(
    list(name="age", type="continuous",  dist="normal",    params=list(mean=60, sd=8)),
    list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.5)),
    list(name="x",   type="continuous",  dist="lognormal", params=list(meanlog=0, sdlog=0.5))
  )),
  treatment  = list(assignment="randomization", allocation="1:1"),
  event_time = list(
    model    = "ph_pwexp",
    baseline = list(rates=c(0.10, 0.06, 0.03), cuts=c(6, 18)),
    effects  = list(intercept=0, treatment=-0.4, covariates=list(age=0.01, x=0.03))
  ),
  censoring  = list(mode="target_overall", target=0.25, admin_time=30),
  seed = 123
)

dat4 <- simulate_from_recipe(validate_recipe(rec4))

# Display: data (first rows)
knitr::kable(utils::head(dat4, 8), caption = "Example 4 — PH Piecewise Exponential (multi-segment): first 8 rows")
```

|       time | status | arm |      age | sex |         x |
|-----------:|-------:|----:|---------:|----:|----------:|
| 17.1446483 |      1 |   1 | 55.51619 |   0 | 2.1580717 |
|  0.4796944 |      1 |   1 | 58.15858 |   1 | 0.9466222 |
|  7.2330213 |      1 |   1 | 72.46967 |   0 | 1.2914109 |
|  2.4702411 |      1 |   0 | 60.56407 |   1 | 1.1129109 |
|  9.3914480 |      0 |   1 | 61.03430 |   1 | 0.9111385 |
|  1.9382470 |      1 |   0 | 73.72052 |   0 | 0.9415791 |
|  1.0707884 |      0 |   1 | 63.68733 |   1 | 1.6593354 |
|  6.3131555 |      1 |   1 | 49.87951 |   0 | 0.9041780 |

Example 4 — PH Piecewise Exponential (multi-segment): first 8 rows
{.table}

``` r


# Display: summary + attributes
summ4 <- data.frame(
  Metric = c("n", "Events", "Censoring rate", "Mean time", "Median time"),
  Value  = c(nrow(dat4),
             sum(dat4$status == 1, na.rm = TRUE),
             round(mean(dat4$status == 0, na.rm = TRUE), 3),
             round(mean(dat4$time,   na.rm = TRUE), 2),
             round(stats::median(dat4$time, na.rm = TRUE), 2))
)
knitr::kable(summ4, caption = "Example 4 — Summary")
```

| Metric         |   Value |
|:---------------|--------:|
| n              | 500.000 |
| Events         | 372.000 |
| Censoring rate |   0.256 |
| Mean time      |   6.120 |
| Median time    |   3.650 |

Example 4 — Summary {.table}

``` r


attr_tbl4 <- data.frame(
  attribute = c("achieved_censoring", "n_rows", "n_cols", "class"),
  value = c(
    as.character(attr(dat4, "achieved_censoring")),
    nrow(dat4),
    ncol(dat4),
    paste(class(dat4), collapse = ", ")
  )
)
knitr::kable(attr_tbl4, caption = "Example 4 — Attributes")
```

| attribute          | value      |
|:-------------------|:-----------|
| achieved_censoring | 0.256      |
| n_rows             | 500        |
| n_cols             | 6          |
| class              | data.frame |

Example 4 — Attributes {.table}

------------------------------------------------------------------------

## Generate data based on formula (event & censoring)

This example uses `~` formulas for both event-time covariate effects and
covariate-dependent censoring. We keep the treatment effect as a
separate scalar to avoid design ambiguity.

``` r

# Covariates (sex numeric 0/1)
covs_f <- list(
  list(name="age", type="continuous",  dist="normal",    params=list(mean=62, sd=10)),
  list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.45))
)

trt_f <- list(assignment="randomization", allocation="1:1")

# Event-time engine: AFT lognormal with formula-based covariate effects
evt_f <- list(
  model    = "aft_lognormal",
  baseline = list(mu = 3.1, sigma = 0.55),
  effects  = list(
    treatment  = -0.25,                 # effect of treatment (log-time scale)
    formula    = "~ 0 + age + sex",     # covariate effects via formula (no intercept)
    beta       = c(0.010, -0.20)        # (age, sex) -- match order in formula
  )
)

# Censoring: covariate-dependent via log-rate model (single mechanism)
cz_f <- list(
  mode = "explicit",
  dependent = list(
    formula = "~ 0 + age + sex",  # NO intercept; sex numeric 0/1
    base    = 0.03,               # positive rate
    beta    = c(0.012, 0.35)      # (age, sex) on log-rate scale
  ),
  administrative = list(time = 48)  # optional admin cutoff
)

# Full recipe
rec_f <- list(
  n = 350,
  covariates = list(defs = covs_f),
  treatment  = trt_f,
  event_time = evt_f,
  censoring  = cz_f,
  seed = 777
)

# Simulate
dat_f <- simulate_from_recipe(validate_recipe(rec_f))

# Display
knitr::kable(
  data.frame(achieved_censoring = attr(dat_f, "achieved_censoring")),
  caption = "Achieved censoring rate (formula-based example)"
)
```

| achieved_censoring |
|-------------------:|
|          0.9028571 |

Achieved censoring rate (formula-based example) {.table}

``` r

knitr::kable(utils::head(dat_f, 8), caption = "First 8 rows — formula-based example")
```

|       time | status | arm |      age | sex |
|-----------:|-------:|----:|---------:|----:|
|  0.5200119 |      0 |   1 | 66.89786 |   1 |
| 12.4499453 |      0 |   0 | 58.01459 |   0 |
| 18.2041119 |      0 |   0 | 67.10836 |   0 |
| 10.3316633 |      0 |   0 | 58.01188 |   1 |
|  0.7918592 |      0 |   1 | 78.38686 |   0 |
|  9.0102499 |      0 |   1 | 68.21274 |   1 |
|  1.0401903 |      0 |   1 | 64.02704 |   1 |
|  0.7415469 |      0 |   0 | 73.08938 |   0 |

First 8 rows — formula-based example {.table}

------------------------------------------------------------------------

## Batch generation with metadata

For simulation studies, write multiple scenarios and formats together.
The writer creates a `manifest.rds` with a **list-column** `meta`
describing each dataset. The loader reattaches attributes when reading
back.

``` r
# Build a base recipe here (self-contained)
base <- validate_recipe(list(
  n = 200,
  covariates = list(defs = list(
    list(name="age", type="continuous",  dist="normal", params=list(mean=62, sd=10)),
    list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.45)),
    list(name="x",   type="continuous",  dist="lognormal", params=list(meanlog=0, sdlog=0.6))
  )),
  treatment  = list(assignment="randomization", allocation="1:1"),
  event_time = list(
    model    = "aft_weibull",
    baseline = list(shape=1.3, scale=12),
    effects  = list(intercept=0, treatment=-0.20, covariates=list(age=0.008, x=0.04))
  ),
  censoring  = list(mode="target_overall", target=0.25, admin_time=36),
  seed = 3026
))

out_dir <- file.path(tempdir(), "rmstss-manifest-demo")
unlink(out_dir, recursive = TRUE, force = TRUE)

man <- generate_recipe_sets(
  base_recipe = base,
  vary = list(n = c(200, 400),
              "event_time.effects.treatment" = c(-0.15, -0.25)),
  out_dir  = out_dir,
  formats  = c("rds","csv"),
  n_reps   = 1,
  seed_base = 2025
)

# Inspect the first row's compact metadata (fields only; no file paths)
m <- readRDS(file.path(out_dir, "manifest.rds"))
names(m)
 [1] "scenario_id"                     "rep"                            
 [3] "seed"                            "achieved_censoring"             
 [5] "n"                               "file_txt"                       
 [7] "file_csv"                        "file_rds"                       
 [9] "file_rdata"                      "p__n"                           
[11] "p__event_time.effects.treatment" "meta"                           
if ("meta" %in% names(m) && length(m$meta[[1]]) > 0) {
  list(model = m$meta[[1]]$model,
       baseline = m$meta[[1]]$baseline,
       effects = m$meta[[1]]$effects,
       achieved_censoring = m$meta[[1]]$achieved_censoring,
       n = m$meta[[1]]$n)
} else {
  "Manifest is minimal (older run); use rebuild_manifest() to enrich."
}
$model
[1] "aft_weibull"

$baseline
$baseline$shape
[1] 1.3

$baseline$scale
[1] 12


$effects
$effects$intercept
[1] 0

$effects$treatment
[1] -0.15

$effects$covariates
$effects$covariates$age
[1] 0.008

$effects$covariates$x
[1] 0.04



$achieved_censoring
[1] 0.26

$n
[1] 200

# Load datasets back
sets <- load_recipe_sets(file.path(out_dir, "manifest.rds"))
attr(sets[[1]]$data, "achieved_censoring")
[1] 0.26
str(sets[[1]]$meta)
List of 19
 $ dataset_id        : chr "sc001_r01"
 $ scenario_id       : int 1
 $ rep               : int 1
 $ seed_used         : int 3026
 $ n                 : int 200
 $ n_treat           : int 95
 $ n_control         : int 105
 $ event_rate        : num 0.74
 $ achieved_censoring: num 0.26
 $ model             : chr "aft_weibull"
 $ baseline          :List of 2
  ..$ shape: num 1.3
  ..$ scale: num 12
 $ effects           :List of 3
  ..$ intercept : num 0
  ..$ treatment : num -0.15
  ..$ covariates:List of 2
  .. ..$ age: num 0.008
  .. ..$ x  : num 0.04
 $ treatment         :List of 2
  ..$ assignment: chr "randomization"
  ..$ allocation: chr "1:1"
 $ censoring         :List of 3
  ..$ mode      : chr "target_overall"
  ..$ target    : num 0.25
  ..$ admin_time: num 36
 $ covariates        :List of 3
  ..$ :List of 4
  .. ..$ name  : chr "age"
  .. ..$ type  : chr "continuous"
  .. ..$ dist  : chr "normal"
  .. ..$ params:List of 2
  .. .. ..$ mean: num 62
  .. .. ..$ sd  : num 10
  ..$ :List of 4
  .. ..$ name  : chr "sex"
  .. ..$ type  : chr "categorical"
  .. ..$ dist  : chr "bernoulli"
  .. ..$ params:List of 1
  .. .. ..$ p: num 0.45
  ..$ :List of 4
  .. ..$ name  : chr "x"
  .. ..$ type  : chr "continuous"
  .. ..$ dist  : chr "lognormal"
  .. ..$ params:List of 2
  .. .. ..$ meanlog: num 0
  .. .. ..$ sdlog  : num 0.6
 $ allocation        : chr "1:1"
 $ params            :List of 2
  ..$ n                           : num 200
  ..$ event_time.effects.treatment: num -0.15
 $ files             :List of 4
  ..$ txt  : chr NA
  ..$ csv  : chr "C:\\Users\\arnab\\AppData\\Local\\Temp\\RtmpA7al8y/rmstss-manifest-demo/sc1_r1.csv"
  ..$ rds  : chr "C:\\Users\\arnab\\AppData\\Local\\Temp\\RtmpA7al8y/rmstss-manifest-demo/sc1_r1.rds"
  ..$ rdata: chr NA
 $ created_at        : chr "2026-04-19 22:31:18.833457"
```

------------------------------------------------------------------------

## Reproducibility tips

Set `seed` in the recipe or pass `seed=` to
[`simulate_from_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).

> **Assumptions.** The censoring hazard depends on observed covariates
> via a user-specified single mechanism. There are no competing risks.
> Datasets contain only `time`, `status`, treatment, and covariates - no
> dependent-censoring indicators.
