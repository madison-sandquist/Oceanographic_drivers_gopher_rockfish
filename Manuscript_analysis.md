Oceanographic drivers of gopher rockfish recruitment — manuscript
analyses
================
Madison Sandquist
2026-10-14

``` r
# Libraries needed
library(tidyverse) # data cleaning
library(mgcv) # GAM pacage 
library(r4ss) # to get stock assessment deviations
library(DHARMa) # model diagonistics
library(patchwork)   # multi-panel assembly
library(corrplot) # correlation plot
# font size for all figures 
BASE_SIZE <- 16
theme_set(theme_bw(base_size = BASE_SIZE))
```

------------------------------------------------------------------------

# 1. Recruitment deviations

``` r
model_path <- "D:/GPHR_BYEL_POST_STAR_BASE_MODEL" # file from the PFMC website
goph2019 <- r4ss::SS_output(model_path, printstats = FALSE, verbose = FALSE)

# Recruitment deviations from the recruit table
goph_dev <- goph2019$recruit %>%
  filter(Yr >= 1987, Yr <= 2018) %>%
  transmute(recruit_year = Yr, dev)

# Standard error of each recruitment deviation.
rec_dev_se <- goph2019$parameters %>%
  filter(grepl("RecrDev", Label)) %>%
  transmute(recruit_year = as.integer(sub(".*_(\\d{4})$", "\\1", Label)),
            dev_se        = Parm_StDev) %>%
  filter(!is.na(recruit_year))

# joining
goph_recruit <- goph_dev %>%
  left_join(rec_dev_se, by = "recruit_year")
head(goph_recruit)
```

    ##   recruit_year        dev   dev_se
    ## 1         1987 -0.2370630 0.473074
    ## 2         1988 -0.0614859 0.516267
    ## 3         1989  0.1687480 0.606723
    ## 4         1990  0.2624540 0.770210
    ## 5         1991  1.4722400 0.436820
    ## 6         1992  0.1976540 0.702007

## Figure 1 — Spawning output, recruitment, and the stock–recruit relationship

``` r
# recruitment table 
rec <- goph2019$recruit
sb_col  <- intersect(c("SpawnBio", "spawn_bio", "SpawnBioAll"), names(rec))[1]
rec_col <- intersect(c("pred_recr", "predicted_recruits"),      names(rec))[1]
rec <- rec %>%
  mutate(SpawnOut = .data[[sb_col]],
         Age0     = .data[[rec_col]]) %>% 
  filter(Yr > 1987) %>% 
  filter(Yr < 2019)

rec_ts <- rec %>%
  filter(!is.na(Yr),
         !era %in% c("VIRG", "Virg", "Virgin", "INIT", "Init")) %>% 
  filter(Yr > 1987) %>% 
  filter(Yr < 2019)

# years for the plot
FIG_YRS <- c(1987, 2018)

# making beverton holt sr curve
pr  <- goph2019$parameters
R0  <- exp(pr$Value[grep("SR_LN\\(R0\\)", pr$Label)][1])
h   <-      pr$Value[grep("SR_BH_steep",  pr$Label)][1]
dq  <- goph2019$derived_quants
SB0 <- dq$Value[grep("^SSB_(Virgin|[Uu]nfished)$", dq$Label)][1]
if (is.na(SB0)) SB0 <- rec$SpawnOut[rec$era %in% c("VIRG", "Virg", "Virgin")][1]

sr_curve <- tibble(SSB = seq(0, max(rec$SpawnOut, na.rm = TRUE) * 1.02,
                             length.out = 300)) %>%
  mutate(Recruitment = (4 * h * R0 * SSB) / (SB0 * (1 - h) + SSB * (5 * h - 1)))

# spawning outout and recruits 
df_ts <- filter(rec_ts, between(Yr, FIG_YRS[1], FIG_YRS[2]))
coef <- max(df_ts$SpawnOut, na.rm = TRUE) / max(df_ts$Age0, na.rm = TRUE)

# attach the recruitment-deviation SE so the deviations can be drawn with +/- SE
df_ts <- df_ts %>%
  left_join(dplyr::select(goph_recruit, recruit_year, dev_se),
            by = c("Yr" = "recruit_year"))

# map recruitment deviations onto the primary (spawning-biomass) axis 
sb_rng  <- range(df_ts$SpawnOut, na.rm = TRUE)
dev_rng <- range(c(df_ts$dev - df_ts$dev_se, df_ts$dev + df_ts$dev_se), na.rm = TRUE)
a_dev   <- diff(sb_rng) / diff(dev_rng)
b_dev   <- sb_rng[1] - a_dev * dev_rng[1]
to_sb   <- function(d) a_dev * d + b_dev      # deviation -> primary-axis units

sb_colour  <- "black"      
age_colour <- "steelblue"   

pAB <- ggplot(df_ts, aes(x = Yr)) +
  geom_hline(yintercept = to_sb(0), linetype = "dashed", colour = "grey70") +
  geom_line(aes(y = SpawnOut, colour = "Spawning output")) +
  geom_errorbar(aes(ymin = to_sb(dev - dev_se), ymax = to_sb(dev + dev_se),
                    colour = "Recruitment deviation"),
                width = 0.4, linewidth = 0.4, na.rm = TRUE) +
  geom_line(aes(y = to_sb(dev),  colour = "Recruitment deviation"), na.rm = TRUE) +
  geom_point(aes(y = to_sb(dev), colour = "Recruitment deviation"),
             size = 1.4, na.rm = TRUE) +
  geom_vline(xintercept = 1995, colour = "red", linetype = "dashed") +
  scale_y_continuous(
    name     = "Spawning Biomass",
    sec.axis = sec_axis(~ (. - b_dev) / a_dev, name = "Recruitment deviation")
  ) +
  scale_colour_manual(
    values = c("Spawning Biomass (tons)" = sb_colour, "Recruitment deviation" = age_colour),
    name   = NULL
  ) +
  labs(x = "Year") +
  theme_classic() +
  theme(
    legend.position    = "top",
    axis.title.y.left  = element_text(colour = sb_colour),
    axis.title.y.right = element_text(colour = age_colour),
    axis.text.y.right  = element_text(colour = age_colour),
    axis.line.y.right  = element_line(colour = age_colour),
    plot.title         = element_text(hjust = 0.5)
  )

# SR curve
pC <- ggplot(sr_curve, aes(SSB, Recruitment)) +
  geom_line() +
  geom_point(data = filter(rec, Yr > 1987),
             aes(SpawnOut, Age0, colour = Yr)) +
  scale_color_gradientn(colours = c("blue", "purple", "red"),
                        name = "Year") +
  labs(x = "Spawning Biomass", y = "Age-0") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

# assmble into one figure 
Figure1 <- (pAB | pC) + plot_annotation(tag_levels = "A")
Figure1
```

![](Manuscript_analysis_files/figure-gfm/figure1-stock-recruit-1.png)<!-- -->

``` r
ggsave("Figures/Figure1_stock_recruit.jpeg", Figure1,
       width = 9, height = 4, dpi = 300)
```

# 2. Oceanographic data

``` r
# ROMS
roms <- read_csv("Data/monthly_averages_allyear_1995_2020_ALL_41to34p4.csv") %>%
  transmute(Year, Month = match(Month, month.abb),
            Temperature = Temperature_avg, DO = DO_avg, pH = pH_avg,
            Chla = Chla_avg)

# CUTI
cuti <- read.csv("Data/CUTI_daily.csv") %>%
  select(year, month, X35N, X36N, X37N, X38N, X39N, X40N, X41N) %>%
  pivot_longer(starts_with("X"), names_to = "lat", values_to = "cuti") %>%
  group_by(Year = year, Month = month) %>%
  summarise(CUTI = mean(cuti, na.rm = TRUE), .groups = "drop")

# ONI
oni <- read.csv("Data/ONI_data.csv")
names(oni)[1] <- "Year"                     
oni <- oni %>%
  pivot_longer(all_of(month.abb), names_to = "mon", values_to = "ONI") %>%
  transmute(Year, Month = match(mon, month.abb), ONI)

# PDO
pdo <- read.csv("Data/pdo.timeseries.sstens.csv")
names(pdo) <- c("Date", "PDO")
pdo <- pdo %>%
  mutate(Date = as.Date(Date)) %>%
  transmute(Year = year(Date), Month = month(Date), PDO)

# Combine
monthly <- roms %>%
  full_join(cuti, by = c("Year", "Month")) %>%
  full_join(oni,  by = c("Year", "Month")) %>%
  full_join(pdo,  by = c("Year", "Month")) %>%
  mutate(recruit_year = Year)

# STI
sti_raw <- read.csv("Data/sti.csv")
names(sti_raw)[1] <- "Year"
sti <- sti_raw %>%
  pivot_longer(-Year, names_to = "lat", values_to = "v") %>%
  mutate(lat = as.integer(gsub("^X", "", lat))) %>%   
  filter(between(lat, 35, 41)) %>%
  group_by(recruit_year = as.integer(round(Year))) %>%
  summarise(STI = mean(v, na.rm = TRUE), .groups = "drop")

# TUMI
tumi_raw <- read.csv("Data/tumi.csv")
names(tumi_raw)[1] <- "Year"
tumi <- tumi_raw %>%
  pivot_longer(-Year, names_to = "lat", values_to = "v") %>%
  mutate(lat = as.integer(gsub("^X", "", lat))) %>%  
  filter(between(lat, 35, 41)) %>%
  group_by(recruit_year = as.integer(round(Year))) %>%
  summarise(TUMI = mean(v, na.rm = TRUE), .groups = "drop")
```

## 2a. Figure 2B Driver overview time series

``` r
panel_levels <- c(
  "Dissolved Oxygen (mg/L)", "pH", "Temperature (°C)", "CUTI (m²/s)",
  "ONI (°C)", "PDO", "STI (day of year)", "TUMI", "Chlorophyll (mg/m³)"
)

# Six monthly drivers 
monthly_overview <- monthly %>%
  filter(Year >= 1995, Year <= 2018) %>%
  mutate(Date = as.Date(sprintf("%d-%02d-01", Year, Month))) %>%
  select(Date, DO, pH, Temperature, CUTI, ONI, PDO, Chla) %>%
  pivot_longer(-Date, names_to = "driver", values_to = "value") %>%
  mutate(panel = dplyr::recode(driver,
    DO          = "Dissolved Oxygen (mg/L)",
    pH          = "pH",
    Temperature = "Temperature (°C)",
    CUTI        = "CUTI (m²/s)",
    ONI         = "ONI (°C)",
    PDO         = "PDO", 
    Chla = "Chlorophyll (mg/m³)"
  )) %>%
  mutate(panel = factor(panel, levels = panel_levels)) %>%
  filter(!is.na(value))

# Two annual indices -> one point per year 
annual_overview <- bind_rows(
  sti  %>% transmute(year = recruit_year, panel = "STI (day of year)", value = STI),
  tumi %>% transmute(year = recruit_year, panel = "TUMI",              value = TUMI)
) %>%
  filter(year >= 1995, year <= 2018, !is.na(value)) %>%
  mutate(
    Date  = as.Date(sprintf("%d-07-01", year)),
    panel = factor(panel, levels = panel_levels)
  )

# shading band for gopher reproductive season
shade_overview <- tibble(year = 1995:2018) %>%
  mutate(
    xmin = as.Date(sprintf("%d-01-01", year)),
    xmax = as.Date(sprintf("%d-08-31", year))
  )

fig_driver_overview <- ggplot() +
  geom_rect(
    data = shade_overview,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "lightblue", alpha = 0.4
  ) +
  geom_line(
    data = monthly_overview, aes(Date, value),
    color = "black", linewidth = 0.5
  ) +
  geom_line(
    data = annual_overview, aes(Date, value),
    color = "grey40", linewidth = 0.4
  ) +
  geom_point(
    data = annual_overview, aes(Date, value),
    color = "black", size = 1.5
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_x_date(
    date_breaks = "3 years", date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

fig_driver_overview
```

![](Manuscript_analysis_files/figure-gfm/driver-overview-1.png)<!-- -->

``` r
ggsave("Figures/driver_overview_timeseries.jpeg", fig_driver_overview,
       width = 11, height = 9, dpi = 300)
```

------------------------------------------------------------------------

# 3. Stage-windowed, standardized predictors

``` r
# function to calculate window means for each hypothesis below
win <- function(var, months, name) {
  monthly %>%
    filter(Month %in% months) %>%
    group_by(recruit_year) %>%
    summarise(across(all_of(var), ~ mean(.x, na.rm = TRUE), .names = "{name}"),
              .groups = "drop")
}

predictors <- goph_recruit %>%
  # Fertilization / gestation (Jan-Mar, Yr 1)
  left_join(win("Temperature", 1:3,  "Tegg"),   by = "recruit_year") %>%
  left_join(win("DO",          1:3,  "DOegg"),  by = "recruit_year") %>%
  left_join(win("pH",          1:3,  "pHegg"),  by = "recruit_year") %>%
  # Larval release (Mar-Apr, Yr 1)
  left_join(win("Temperature", 3:4,  "Tlarv"),  by = "recruit_year") %>%
  left_join(win("ONI",         3:4,  "ONIlarv"),by = "recruit_year") %>%
  left_join(win("PDO",         3:4,  "PDOlarv"),by = "recruit_year") %>%
  # Pelagic phase (Apr-Jul, Yr 1)
  left_join(win("Temperature", 4:6,  "Tpld"),   by = "recruit_year") %>%
  left_join(win("PDO",         4:6,  "PDOpld"), by = "recruit_year") %>%
  left_join(win("CUTI",        4:6,  "CUTIpld"),by = "recruit_year") %>%
  left_join(win("ONI",         4:6,  "ONIpld"), by = "recruit_year") %>%
  left_join(win("Chla",        4:6,  "Chlapld"),by = "recruit_year") %>%
  # Benthic settlement (Jul-Sep, Yr 1)
  left_join(win("PDO",         7:9,  "PDOben"), by = "recruit_year") %>%
  left_join(win("CUTI",        7:9,  "CUTIben"),by = "recruit_year") %>%
  left_join(win("ONI",         7:9,  "ONIben"), by = "recruit_year") %>%
  left_join(win("Chla",        7:9,  "Chlaben"),by = "recruit_year") %>%
  # Annual upwelling-phenology indices 
  left_join(sti,  by = "recruit_year") %>%
  left_join(tumi, by = "recruit_year") %>%
  filter(recruit_year > 1994)

pred_names <- setdiff(names(predictors), c("recruit_year", "dev", "dev_se"))

# Standardize every predictor to mean 0 
pred_std <- predictors %>%
  mutate(across(all_of(pred_names), ~ as.numeric(scale(.))))
```

# 4. Null benchmark - from SR curve

``` r
# Fit a null GAM model with only an intercept
null_gam_model <- gam(dev ~ 1, data = pred_std)

# Calculate null predictions (intercept only)
pred_std <- pred_std %>%
  mutate(sr_null = 0, mean_null = predict(null_gam_model))

# Caculate Root Mean Squared Error (RMSE) for the null models (full dataset)
rmse_sr_full <- sqrt(mean((pred_std$sr_null - pred_std$dev)^2, na.rm = TRUE))
rmse_mean_full <- sqrt(mean((pred_std$mean_null - pred_std$dev)^2, na.rm = TRUE))

# Generate null model predictions for the LOOCV timeframe
null_lfo_data <- pred_std %>%
  mutate(sr_null = 0, mean_null = predict(null_gam_model)) %>%
  filter(recruit_year >= 2013 & recruit_year < 2018)

# Calculate RMSE for the null models (LOOCV timeframe)
rmse_sr_lfo <- sqrt(mean((null_lfo_data$sr_null - null_lfo_data$dev)^2, na.rm = TRUE))
rmse_mean_lfo <- sqrt(mean((null_lfo_data$mean_null - null_lfo_data$dev)^2, na.rm = TRUE))
```

# 5. Univariate capped GAMs with LOO-CV (Table 2, Figure 3)

``` r
max_covariates <- 1
available_covariates <- pred_names   # univariate screen: every candidate driver, one at a time

# Generate all possible combinations of covariates, with max var set at one it should just list the covariates
all_combinations <- unlist(lapply(1:max_covariates, function(i) {
  combn(available_covariates, i, simplify = FALSE)
}), recursive = FALSE)

length(all_combinations)   # number of admissible candidate models
```

    ## [1] 17

``` r
#### --- LOOCV and model fitting --- ####
model_list       <- list()
model_results    <- data.frame()
predicted_values <- data.frame()

# Loop through each admissible combination, fit and evaluate a GAM
for (i in seq_along(all_combinations)) {

  vars <- all_combinations[[i]]

  d <- pred_std %>%
    select(dev, recruit_year, all_of(vars)) %>%
    filter(complete.cases(.))

  # GAM formula
  smooth_terms   <- paste0("s(", vars, ", k = 3)", collapse = " + ")
  formula_string <- paste("dev ~", smooth_terms)

  # Leave-one-out predictions, one held-out year at a time
  unique_years    <- unique(d$recruit_year)
  loo_predictions <- rep(NA_real_, nrow(d))

  for (j in seq_along(unique_years)) {
    training_data <- d %>% filter(recruit_year != unique_years[j])
    testing_data  <- d %>% filter(recruit_year == unique_years[j])

    gam_loo_model <- gam(as.formula(formula_string), data = training_data, method = "REML")
    loo_predictions[which(d$recruit_year == unique_years[j])] <-
      predict(gam_loo_model, newdata = testing_data)
  }

  # Refit on all years for deviance explained / R-squared
  full_gam_model         <- gam(as.formula(formula_string), data = d, method = "REML")
  full_model_predictions <- predict(full_gam_model)

  # RMSEs on the deviation scale
  rmse_loo        <- sqrt(mean((d$dev - loo_predictions)^2,        na.rm = TRUE))
  rmse_full_model <- sqrt(mean((d$dev - full_model_predictions)^2, na.rm = TRUE))

  r_squared          <- summary(full_gam_model)$r.sq
  deviance_explained <- summary(full_gam_model)$dev.expl

  # variable names so every row has the same columns
  padded_variables <- c(vars, rep(NA, max_covariates - length(vars)))
  var_cols <- setNames(as.list(padded_variables), paste0("var", seq_len(max_covariates)))

  # Store results 
  model_list[[i]] <- full_gam_model
  model_results <- bind_rows(model_results, data.frame(
    ModelID    = i,
    AIC        = AIC(full_gam_model),
    RMSE_loo   = round(rmse_loo, 3),
    RMSE       = round(rmse_full_model, 3),
    rsq_full   = round(r_squared, 2),
    dev.ex     = round(deviance_explained, 4),
    rmse_imp   = (rmse_sr_full - rmse_loo) / rmse_sr_full,
    rmse_ratio = rmse_loo / rmse_sr_full,
    var_cols
  ))

  # Store leave-one-out predictions
  predicted_values <- bind_rows(predicted_values, data.frame(
    ModelID = i,
    pred    = loo_predictions,
    year    = d$recruit_year,
    var_cols
  ))

  print(paste("Model", i, "of", length(all_combinations), "completed."))
}
```

    ## [1] "Model 1 of 17 completed."
    ## [1] "Model 2 of 17 completed."
    ## [1] "Model 3 of 17 completed."
    ## [1] "Model 4 of 17 completed."
    ## [1] "Model 5 of 17 completed."
    ## [1] "Model 6 of 17 completed."
    ## [1] "Model 7 of 17 completed."
    ## [1] "Model 8 of 17 completed."
    ## [1] "Model 9 of 17 completed."
    ## [1] "Model 10 of 17 completed."
    ## [1] "Model 11 of 17 completed."
    ## [1] "Model 12 of 17 completed."
    ## [1] "Model 13 of 17 completed."
    ## [1] "Model 14 of 17 completed."
    ## [1] "Model 15 of 17 completed."
    ## [1] "Model 16 of 17 completed."
    ## [1] "Model 17 of 17 completed."

``` r
n_models <- nrow(model_results)   # total admissible candidate model fits
n_models
```

    ## [1] 17

``` r
loo_gam <- function(form, d) {
  n <- nrow(d); err <- numeric(n)
  for (i in seq_len(n)) {
    mi <- gam(form, data = d[-i, ], method = "REML")
    err[i] <- d$dev[i] - as.numeric(predict(mi, d[i, ]))
  }
  list(rmse = sqrt(mean(err^2)),
       r2   = 1 - sum(err^2) / sum((d$dev - mean(d$dev))^2))
}

uni <- map_dfr(pred_names, function(p) {
  d <- pred_std %>% select(dev, x = all_of(p)) %>% filter(complete.cases(.))
  m <- gam(dev ~ s(x, k = 3), data = d, method = "REML")
  s <- summary(m)
  lo <- loo_gam(dev ~ s(x, k = 3), d)
  null_sr <- sqrt(mean(d$dev^2, na.rm = TRUE))   # S-R null: deviations vs 0
  tibble(
    predictor = p, n = nrow(d),
    edf      = round(s$s.table[1, "edf"], 2),
    F        = round(s$s.table[1, "F"], 2),
    p_value  = signif(s$s.table[1, "p-value"], 3),
    dev_expl = round(s$dev.expl, 3),
    AIC      = round(AIC(m), 1),
    LOO_RMSE = round(lo$rmse, 3),
    LOO_R2   = round(lo$r2, 3),
    null_sr   = round(null_sr, 3),
    RMSE_impr = round((null_sr - lo$rmse) / null_sr, 3)
  )
}) %>% arrange(desc(LOO_R2))

knitr::kable(uni, caption = "Univariate GAMs")
```

| predictor |   n |  edf |     F | p_value | dev_expl |  AIC | LOO_RMSE | LOO_R2 | null_sr | RMSE_impr |
|:----------|----:|-----:|------:|--------:|---------:|-----:|---------:|-------:|--------:|----------:|
| PDOlarv   |  24 | 1.00 | 14.11 | 0.00109 |    0.391 | 16.0 |    0.325 |  0.274 |   0.397 |     0.180 |
| DOegg     |  24 | 1.00 |  8.28 | 0.00875 |    0.273 | 20.2 |    0.348 |  0.169 |   0.397 |     0.123 |
| TUMI      |  24 | 1.41 |  5.61 | 0.00967 |    0.353 | 18.8 |    0.351 |  0.156 |   0.397 |     0.117 |
| PDOpld    |  24 | 1.00 |  8.39 | 0.00837 |    0.276 | 20.1 |    0.352 |  0.150 |   0.397 |     0.114 |
| Tlarv     |  24 | 1.00 |  5.54 | 0.02790 |    0.201 | 22.5 |    0.371 |  0.056 |   0.397 |     0.066 |
| pHegg     |  24 | 1.61 |  5.57 | 0.02860 |    0.306 | 20.8 |    0.372 |  0.050 |   0.397 |     0.062 |
| CUTIben   |  24 | 1.58 |  6.82 | 0.01470 |    0.346 | 19.3 |    0.373 |  0.048 |   0.397 |     0.062 |
| STI       |  24 | 1.00 |  4.55 | 0.04440 |    0.171 | 23.4 |    0.374 |  0.042 |   0.397 |     0.059 |
| ONIpld    |  24 | 1.00 |  4.21 | 0.05230 |    0.161 | 23.7 |    0.379 |  0.016 |   0.397 |     0.046 |
| Tegg      |  24 | 1.00 |  5.54 | 0.02790 |    0.201 | 22.5 |    0.382 |  0.000 |   0.397 |     0.038 |
| Chlapld   |  24 | 1.00 |  3.84 | 0.06290 |    0.149 | 24.0 |    0.384 | -0.009 |   0.397 |     0.034 |
| Chlaben   |  24 | 1.00 |  2.62 | 0.12000 |    0.106 | 25.2 |    0.388 | -0.031 |   0.397 |     0.023 |
| ONIlarv   |  24 | 1.00 |  2.79 | 0.10900 |    0.112 | 25.0 |    0.391 | -0.048 |   0.397 |     0.015 |
| CUTIpld   |  24 | 1.00 |  0.17 | 0.68400 |    0.008 | 27.7 |    0.404 | -0.119 |   0.397 |    -0.017 |
| PDOben    |  24 | 1.12 |  5.08 | 0.02350 |    0.242 | 21.7 |    0.405 | -0.125 |   0.397 |    -0.020 |
| Tpld      |  24 | 1.00 |  1.34 | 0.25900 |    0.057 | 26.5 |    0.408 | -0.143 |   0.397 |    -0.028 |
| ONIben    |  24 | 1.00 |  0.89 | 0.35500 |    0.039 | 26.9 |    0.413 | -0.170 |   0.397 |    -0.040 |

Univariate GAMs

``` r
# How many hypotheses are significant at p < 0.05 ? In manuscript
uni %>%
  summarise(
    total_models       = n(),
    significant_models = sum(p_value < 0.05, na.rm = TRUE),
    percent_significant = 100 * mean(p_value < 0.05, na.rm = TRUE)
  )
```

    ## # A tibble: 1 × 3
    ##   total_models significant_models percent_significant
    ##          <int>              <int>               <dbl>
    ## 1           17                 10                58.8

## 5a. Observed deviation dataset

``` r
obs_dev <- pred_std %>% select(recruit_year, dev) %>% filter(!is.na(dev))
```

## 5b. Figure 3B — deviations with the top univariate GAM fits

``` r
N_TOP <- 6
top_preds <- head(uni$predictor, N_TOP)

fitted_long <- map_dfr(top_preds, function(p) {
  d <- pred_std %>% select(recruit_year, dev, x = all_of(p)) %>% filter(complete.cases(.))
  m <- gam(dev ~ s(x, k = 3), data = d, method = "REML")
  tibble(recruit_year = d$recruit_year, predictor = p, fitted = as.numeric(fitted(m)))
}) %>% mutate(predictor = factor(predictor, levels = top_preds))

fig_uni_ts <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(data = obs_dev, aes(recruit_year, dev), size = 2.2, color = "grey15") +
  geom_line(data = fitted_long, aes(recruit_year, fitted, color = predictor),
            linewidth = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Recruitment year", y = "Recruitment deviation",
       color = "Univariate GAM fit") +
  theme_classic(base_size = BASE_SIZE)
fig_uni_ts
```

![](Manuscript_analysis_files/figure-gfm/uni-timeseries-1.png)<!-- -->

``` r
ggsave("Figures/univariate_gam_timeseries.jpeg", fig_uni_ts, width = 9, height = 5, dpi = 300)
```

## 5c. Figure 3A — driver–deviation relationships

``` r
pts_long <- map_dfr(top_preds, function(p) {
  pred_std %>% select(dev, value = all_of(p)) %>% filter(complete.cases(.)) %>%
    mutate(predictor = p)
}) %>% mutate(predictor = factor(predictor, levels = top_preds))

line_long <- map_dfr(top_preds, function(p) {
  d <- pred_std %>% select(dev, value = all_of(p)) %>% filter(complete.cases(.))
  m <- gam(dev ~ s(value, k = 3), data = d, method = "REML")
  grid <- tibble(value = seq(min(d$value), max(d$value), length.out = 100))
  tibble(predictor = p, value = grid$value, fit = as.numeric(predict(m, grid)))
}) %>% mutate(predictor = factor(predictor, levels = top_preds))

fig_uni_rel <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(data = pts_long, aes(value, dev, color = predictor),
             size = 1.7, alpha = 0.8) +
  geom_line(data = line_long, aes(value, fit, color = predictor), linewidth = 1) +
  facet_wrap(~ predictor, ncol = 3, scales = "free_x") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Standardized driver value", y = "Recruitment deviation") +
  theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "none")
fig_uni_rel
```

![](Manuscript_analysis_files/figure-gfm/uni-relationships-1.png)<!-- -->

``` r
ggsave("Figures/univariate_gam_relationships.jpeg", fig_uni_rel,
       width = 9, height = 5.5, dpi = 300)
```

------------------------------------------------------------------------

# 6. Collinearity screening (Figure S2 & S3)

``` r
cormat <- cor(pred_std[pred_names], use = "pairwise.complete.obs")
round(cormat, 2)
```

    ##          Tegg DOegg pHegg Tlarv ONIlarv PDOlarv  Tpld PDOpld CUTIpld ONIpld
    ## Tegg     1.00  0.73  0.83  0.92    0.75    0.68  0.75   0.46   -0.36   0.53
    ## DOegg    0.73  1.00  0.91  0.70    0.45    0.52  0.37   0.38   -0.34   0.28
    ## pHegg    0.83  0.91  1.00  0.76    0.59    0.63  0.53   0.47   -0.44   0.35
    ## Tlarv    0.92  0.70  0.76  1.00    0.70    0.65  0.83   0.43   -0.25   0.47
    ## ONIlarv  0.75  0.45  0.59  0.70    1.00    0.75  0.55   0.63   -0.35   0.80
    ## PDOlarv  0.68  0.52  0.63  0.65    0.75    1.00  0.46   0.90   -0.42   0.63
    ## Tpld     0.75  0.37  0.53  0.83    0.55    0.46  1.00   0.28   -0.34   0.46
    ## PDOpld   0.46  0.38  0.47  0.43    0.63    0.90  0.28   1.00   -0.49   0.66
    ## CUTIpld -0.36 -0.34 -0.44 -0.25   -0.35   -0.42 -0.34  -0.49    1.00  -0.37
    ## ONIpld   0.53  0.28  0.35  0.47    0.80    0.63  0.46   0.66   -0.37   1.00
    ## Chlapld -0.61 -0.43 -0.47 -0.71   -0.60   -0.73 -0.62  -0.71    0.48  -0.60
    ## PDOben   0.28  0.12  0.23  0.18    0.44    0.60  0.11   0.75   -0.28   0.68
    ## CUTIben -0.13 -0.17 -0.17  0.03   -0.15   -0.35  0.08  -0.41    0.34  -0.35
    ## ONIben  -0.10 -0.23 -0.27 -0.16    0.06    0.12  0.00   0.26   -0.08   0.61
    ## Chlaben  0.30  0.11  0.15  0.26    0.22    0.45  0.41   0.41   -0.36   0.27
    ## STI      0.60  0.77  0.74  0.56    0.44    0.39  0.29   0.33   -0.46   0.25
    ## TUMI    -0.41 -0.59 -0.52 -0.38   -0.34   -0.53 -0.22  -0.57    0.61  -0.40
    ##         Chlapld PDOben CUTIben ONIben Chlaben   STI  TUMI
    ## Tegg      -0.61   0.28   -0.13  -0.10    0.30  0.60 -0.41
    ## DOegg     -0.43   0.12   -0.17  -0.23    0.11  0.77 -0.59
    ## pHegg     -0.47   0.23   -0.17  -0.27    0.15  0.74 -0.52
    ## Tlarv     -0.71   0.18    0.03  -0.16    0.26  0.56 -0.38
    ## ONIlarv   -0.60   0.44   -0.15   0.06    0.22  0.44 -0.34
    ## PDOlarv   -0.73   0.60   -0.35   0.12    0.45  0.39 -0.53
    ## Tpld      -0.62   0.11    0.08   0.00    0.41  0.29 -0.22
    ## PDOpld    -0.71   0.75   -0.41   0.26    0.41  0.33 -0.57
    ## CUTIpld    0.48  -0.28    0.34  -0.08   -0.36 -0.46  0.61
    ## ONIpld    -0.60   0.68   -0.35   0.61    0.27  0.25 -0.40
    ## Chlapld    1.00  -0.48    0.05  -0.22   -0.46 -0.36  0.49
    ## PDOben    -0.48   1.00   -0.50   0.58    0.35  0.18 -0.32
    ## CUTIben    0.05  -0.50    1.00  -0.31   -0.28 -0.26  0.63
    ## ONIben    -0.22   0.58   -0.31   1.00    0.18 -0.27 -0.11
    ## Chlaben   -0.46   0.35   -0.28   0.18    1.00  0.08 -0.42
    ## STI       -0.36   0.18   -0.26  -0.27    0.08  1.00 -0.63
    ## TUMI       0.49  -0.32    0.63  -0.11   -0.42 -0.63  1.00

``` r
lev <- rownames(cormat)[order(tolower(rownames(cormat)))] # alphabetical axis order
cm <- as.data.frame(as.table(cormat)) %>%
  setNames(c("v1", "v2", "r")) %>%
  mutate(v1 = factor(v1, levels = lev),
         v2 = factor(v2, levels = rev(lev)),
         admissible = abs(r) < 0.3 & as.character(v1) != as.character(v2))

# number of correlations admissible into the multivariate models
cm_true <- cm %>% filter(admissible)

fig_corr <- ggplot(cm, aes(v1, v2, fill = r)) +
  geom_tile(color = "grey92") +
  geom_tile(data = filter(cm, admissible), color = "black", linewidth = 1, fill = NA) +
  geom_text(aes(label = sub("^(-?)0", "\\1", formatC(r, format = "f", digits = 2))),
            size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  coord_fixed()
fig_corr
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("Figures/driver_correlation.jpeg", fig_corr, width = 8, height = 7, dpi = 300)
```

``` r
ok  <- abs(cormat) < 0.3
drv <- rownames(cormat)
n_by_size <- sapply(1:3, \(k) sum(combn(drv, k, \(s) all(ok[s, s][upper.tri(ok[s, s])]))))
c(singles = n_by_size[1], pairs = n_by_size[2], triples = n_by_size[3], total = sum(n_by_size))
```

    ## singles   pairs triples   total 
    ##      17      42      21      80

``` r
## correlations of just top 6 univariate predictors
pred_std_top_6 <- pred_std %>% 
  select(PDOlarv, DOegg, TUMI, PDOpld, Tlarv, pHegg)

cormat <- cor(pred_std_top_6, use = "pairwise.complete.obs")

round(cormat, 2)
```

    ##         PDOlarv DOegg  TUMI PDOpld Tlarv pHegg
    ## PDOlarv    1.00  0.52 -0.53   0.90  0.65  0.63
    ## DOegg      0.52  1.00 -0.59   0.38  0.70  0.91
    ## TUMI      -0.53 -0.59  1.00  -0.57 -0.38 -0.52
    ## PDOpld     0.90  0.38 -0.57   1.00  0.43  0.47
    ## Tlarv      0.65  0.70 -0.38   0.43  1.00  0.76
    ## pHegg      0.63  0.91 -0.52   0.47  0.76  1.00

``` r
corrplot(cormat, 
         method = "circle", 
         type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         addCoef.col = "white") 
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# 7. Multivariate model selection (Table 3)

``` r
#### --- Evaluating Correlations Between Covariates --- ####
MAX_VARS <- 3   # up to this many predictors per model
cor_thresh <- 0.3

# Correlation matrix over the candidate predictors (complete cases only)
predictors_for_cor <- pred_std %>%
  select(all_of(pred_names)) %>%
  filter(complete.cases(.))
correlation_matrix <- cor(predictors_for_cor)

# Convert the correlation matrix to long format 
long_correlation_data <- correlation_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("xvar") %>%
  pivot_longer(!xvar, names_to = "yvar", values_to = "correlation")

# Identify correlated variable pairs
correlated_pairs <- long_correlation_data %>%
  filter(abs(correlation) >= cor_thresh, xvar != yvar) %>%
  select(xvar, yvar)

# Turn each correlated pair into a character vector c(xvar, yvar) for filtering
combinations_to_omit <- lapply(seq_len(nrow(correlated_pairs)), function(i) {
  as.character(correlated_pairs[i, ])
})

#### --- 4. Model specification and filtering combinations --- ####
# Candidate covariates and model size
available_covariates <- pred_names
max_covariates       <- MAX_VARS      # up to this many predictors per model

# Generate all possible combinations of covariates
all_combinations <- unlist(lapply(1:max_covariates, function(i) {
  combn(available_covariates, i, simplify = FALSE)
}), recursive = FALSE)

# Ensure unique combinations and sort for consistency
all_combinations <- unique(lapply(all_combinations, sort))

# Check whether a combination contains a highly correlated pair
is_partial_match <- function(combination, omit_list) {
  any(sapply(omit_list, function(omit_pair) all(omit_pair %in% combination)))
}

# Remove combinations containing highly correlated variables
filtered_combinations <- all_combinations[
  !sapply(all_combinations, is_partial_match, omit_list = combinations_to_omit)]

length(filtered_combinations)   # number of admissible candidate models
```

    ## [1] 80

``` r
#### --- 5. LOOCV and model fitting --- ####
model_list       <- list()
model_results    <- data.frame()
predicted_values <- data.frame()

# Loop through each admissible combination, fit and evaluate a GAM
for (i in seq_along(filtered_combinations)) {

  vars <- filtered_combinations[[i]]

  # Complete-case data for just this model's columns (so missing covariates
  # in a year don't break the leave-one-out fold assignment)
  d <- pred_std %>%
    select(dev, recruit_year, all_of(vars)) %>%
    filter(complete.cases(.))

  # Construct the GAM formula
  smooth_terms   <- paste0("s(", vars, ", k = 3)", collapse = " + ")
  formula_string <- paste("dev ~", smooth_terms)

  # Leave-one-out predictions, one held-out year at a time
  unique_years    <- unique(d$recruit_year)
  loo_predictions <- rep(NA_real_, nrow(d))

  for (j in seq_along(unique_years)) {
    training_data <- d %>% filter(recruit_year != unique_years[j])
    testing_data  <- d %>% filter(recruit_year == unique_years[j])

    gam_loo_model <- gam(as.formula(formula_string), data = training_data, method = "REML")
    loo_predictions[which(d$recruit_year == unique_years[j])] <-
      predict(gam_loo_model, newdata = testing_data)
  }

  # Refit on all years for deviance explained / R-squared
  full_gam_model         <- gam(as.formula(formula_string), data = d, method = "REML")
  full_model_predictions <- predict(full_gam_model)

  # RMSEs on the deviation scale
  rmse_loo        <- sqrt(mean((d$dev - loo_predictions)^2,        na.rm = TRUE))
  rmse_full_model <- sqrt(mean((d$dev - full_model_predictions)^2, na.rm = TRUE))

  r_squared          <- summary(full_gam_model)$r.sq
  deviance_explained <- summary(full_gam_model)$dev.expl

  # Pad variable names so every row has the same columns
  padded_variables <- c(vars, rep(NA, max_covariates - length(vars)))
  var_cols <- setNames(as.list(padded_variables), paste0("var", seq_len(max_covariates)))

  # Store results (improvement measured out-of-sample vs the S-R null)
  model_list[[i]] <- full_gam_model
  model_results <- bind_rows(model_results, data.frame(
    ModelID    = i,
    AIC        = AIC(full_gam_model),
    RMSE_loo   = round(rmse_loo, 3),
    RMSE       = round(rmse_full_model, 3),
    rsq_full   = round(r_squared, 2),
    dev.ex     = round(deviance_explained, 4),
    rmse_imp   = (rmse_sr_full - rmse_loo) / rmse_sr_full,
    rmse_ratio = rmse_loo / rmse_sr_full,
    var_cols
  ))

  # Store leave-one-out predictions
  predicted_values <- bind_rows(predicted_values, data.frame(
    ModelID = i,
    pred    = loo_predictions,
    year    = d$recruit_year,
    var_cols
  ))

  print(paste("Model", i, "of", length(filtered_combinations), "completed."))
}
```

    ## [1] "Model 1 of 80 completed."
    ## [1] "Model 2 of 80 completed."
    ## [1] "Model 3 of 80 completed."
    ## [1] "Model 4 of 80 completed."
    ## [1] "Model 5 of 80 completed."
    ## [1] "Model 6 of 80 completed."
    ## [1] "Model 7 of 80 completed."
    ## [1] "Model 8 of 80 completed."
    ## [1] "Model 9 of 80 completed."
    ## [1] "Model 10 of 80 completed."
    ## [1] "Model 11 of 80 completed."
    ## [1] "Model 12 of 80 completed."
    ## [1] "Model 13 of 80 completed."
    ## [1] "Model 14 of 80 completed."
    ## [1] "Model 15 of 80 completed."
    ## [1] "Model 16 of 80 completed."
    ## [1] "Model 17 of 80 completed."
    ## [1] "Model 18 of 80 completed."
    ## [1] "Model 19 of 80 completed."
    ## [1] "Model 20 of 80 completed."
    ## [1] "Model 21 of 80 completed."
    ## [1] "Model 22 of 80 completed."
    ## [1] "Model 23 of 80 completed."
    ## [1] "Model 24 of 80 completed."
    ## [1] "Model 25 of 80 completed."
    ## [1] "Model 26 of 80 completed."
    ## [1] "Model 27 of 80 completed."
    ## [1] "Model 28 of 80 completed."
    ## [1] "Model 29 of 80 completed."
    ## [1] "Model 30 of 80 completed."
    ## [1] "Model 31 of 80 completed."
    ## [1] "Model 32 of 80 completed."
    ## [1] "Model 33 of 80 completed."
    ## [1] "Model 34 of 80 completed."
    ## [1] "Model 35 of 80 completed."
    ## [1] "Model 36 of 80 completed."
    ## [1] "Model 37 of 80 completed."
    ## [1] "Model 38 of 80 completed."
    ## [1] "Model 39 of 80 completed."
    ## [1] "Model 40 of 80 completed."
    ## [1] "Model 41 of 80 completed."
    ## [1] "Model 42 of 80 completed."
    ## [1] "Model 43 of 80 completed."
    ## [1] "Model 44 of 80 completed."
    ## [1] "Model 45 of 80 completed."
    ## [1] "Model 46 of 80 completed."
    ## [1] "Model 47 of 80 completed."
    ## [1] "Model 48 of 80 completed."
    ## [1] "Model 49 of 80 completed."
    ## [1] "Model 50 of 80 completed."
    ## [1] "Model 51 of 80 completed."
    ## [1] "Model 52 of 80 completed."
    ## [1] "Model 53 of 80 completed."
    ## [1] "Model 54 of 80 completed."
    ## [1] "Model 55 of 80 completed."
    ## [1] "Model 56 of 80 completed."
    ## [1] "Model 57 of 80 completed."
    ## [1] "Model 58 of 80 completed."
    ## [1] "Model 59 of 80 completed."
    ## [1] "Model 60 of 80 completed."
    ## [1] "Model 61 of 80 completed."
    ## [1] "Model 62 of 80 completed."
    ## [1] "Model 63 of 80 completed."
    ## [1] "Model 64 of 80 completed."
    ## [1] "Model 65 of 80 completed."
    ## [1] "Model 66 of 80 completed."
    ## [1] "Model 67 of 80 completed."
    ## [1] "Model 68 of 80 completed."
    ## [1] "Model 69 of 80 completed."
    ## [1] "Model 70 of 80 completed."
    ## [1] "Model 71 of 80 completed."
    ## [1] "Model 72 of 80 completed."
    ## [1] "Model 73 of 80 completed."
    ## [1] "Model 74 of 80 completed."
    ## [1] "Model 75 of 80 completed."
    ## [1] "Model 76 of 80 completed."
    ## [1] "Model 77 of 80 completed."
    ## [1] "Model 78 of 80 completed."
    ## [1] "Model 79 of 80 completed."
    ## [1] "Model 80 of 80 completed."

``` r
n_models <- nrow(model_results)   # total admissible candidate model fits
n_models
```

    ## [1] 80

``` r
# Arrange results by AIC and add rank and delta AIC
model_results <- model_results %>%
  arrange(AIC) %>%
  mutate(rankmod = seq_along(ModelID)) %>%
  mutate(deltaAIC = AIC - min(AIC))

# Arrange results by RMSE (LOOCV) and AIC
results_arr_rmse <- arrange(model_results, RMSE_loo)
results_arr_aic <- arrange(model_results, AIC)

# Filter for top models based on Delta AIC
summary_top_models <- model_results %>% filter(deltaAIC < 4)
summary_top_models
```

    ##   ModelID       AIC RMSE_loo  RMSE rsq_full dev.ex  rmse_imp rmse_ratio    var1
    ## 1      69  9.067528    0.296 0.216     0.59 0.6787 0.2554248  0.7445752 Chlaben
    ## 2      65 11.479503    0.334 0.230     0.54 0.6360 0.1581497  0.8418503 Chlaben
    ## 3      24 12.146324    0.311 0.254     0.50 0.5585 0.2164670  0.7835330 CUTIben
    ## 4      33 12.475904    0.325 0.249     0.51 0.5732 0.1817910  0.8182090 CUTIben
    ## 5      60 12.693219    0.317 0.231     0.53 0.6333 0.2002058  0.7997942 Chlaben
    ## 6      63 12.840278    0.330 0.238     0.52 0.6100 0.1683725  0.8316275 Chlaben
    ##      var2  var3 rankmod deltaAIC
    ## 1 CUTIben Tlarv       1 0.000000
    ## 2 CUTIben pHegg       2 2.411975
    ## 3   DOegg  <NA>       3 3.078796
    ## 4   Tlarv  <NA>       4 3.408375
    ## 5 CUTIben  Tegg       5 3.625690
    ## 6 CUTIben DOegg       6 3.772750

``` r
# Write top models to a CSV file
write_csv(summary_top_models, "Figures/TopModels.csv")
```

# 8. Permutation test of selection skill (Figure S3)

``` r
set.seed(42)
N_PERM <- 400   

# Precompute each admissible combo's complete-case data + formula ONCE. Only the
# response is shuffled, so predictors/formulas are fixed. Reuses the SAME `combos`
combo_spec <- lapply(filtered_combinations, function(vars) {
  rows <- which(complete.cases(pred_std[, vars, drop = FALSE]))
  list(rows = rows,
       data = pred_std[rows, vars, drop = FALSE],
       form = as.formula(paste("dev ~", paste0("s(", vars, ", k = 3)", collapse = " + "))))
})

# Deviance explained of the AIC-selected best model for an arbitrary response.
select_best_devexpl <- function(dev_vec) {
  best_aic <- Inf; best_de <- NA_real_
  for (cs in combo_spec) {
    d <- cs$data; d$dev <- dev_vec[cs$rows]
    m <- try(gam(cs$form, data = d, method = "REML"), silent = TRUE)
    if (inherits(m, "try-error")) next
    if (isTRUE(m$aic < best_aic)) {
      best_aic <- m$aic
      best_de  <- (m$null.deviance - m$deviance) / m$null.deviance
    }
  }
  best_de
}

# Observed skill: deviance explained of the AIC-best model from the real selection.
obs_best_devexpl <- model_results %>% slice_min(AIC, n = 1, with_ties = FALSE) %>% pull(dev.ex)

# Null distribution: shuffle the deviations, re-select, record best-model dev. explained.
null_devexpl <- numeric(N_PERM)
pb <- txtProgressBar(min = 0, max = N_PERM, style = 3)
```

    ##   |                                                                              |                                                                      |   0%

``` r
for (b in seq_len(N_PERM)) {
  null_devexpl[b] <- select_best_devexpl(sample(pred_std$dev))
  setTxtProgressBar(pb, b)
}
```

    ##   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

``` r
close(pb)
```

``` r
null_devexpl <- null_devexpl[is.finite(null_devexpl)]

perm_p    <- (1 + sum(null_devexpl >= obs_best_devexpl)) / (1 + length(null_devexpl))
fpr_at_55 <- mean(null_devexpl >= 0.55)

tibble(
  observed_dev_expl     = round(obs_best_devexpl, 3),
  null_median_dev_expl  = round(median(null_devexpl), 3),
  null_95th_dev_expl    = round(quantile(null_devexpl, 0.95), 3),
  empirical_p           = signif(perm_p, 3),
  prop_null_ge_0.55     = round(fpr_at_55, 3),
  n_perm                = length(null_devexpl)
) %>%
  knitr::kable()
```

| observed_dev_expl | null_median_dev_expl | null_95th_dev_expl | empirical_p | prop_null_ge_0.55 | n_perm |
|---:|---:|---:|---:|---:|---:|
| 0.679 | 0.37 | 0.588 | 0.0175 | 0.085 | 400 |

``` r
proj_models <- list(
  "Tlarv + CUTIben + Chlaben" = c("Tlarv", "CUTIben", "Chlaben"),
  "DOegg + CUTIben + Chlaben" = c("DOegg", "CUTIben", "Chlaben"),
  "pHegg + CUTIben + Chlaben" = c("pHegg", "CUTIben", "Chlaben"))

proj_de <- imap_dfr(proj_models, function(vars, nm) {
  d <- pred_std %>% select(dev, all_of(vars)) %>% filter(complete.cases(.))
  f <- as.formula(paste("dev ~", paste0("s(", vars, ", k = 3)", collapse = " + ")))
  de <- summary(gam(f, data = d, method = "REML"))$dev.expl
  tibble(model = nm, dev_expl = de,
         emp_p = (1 + sum(null_devexpl >= de)) / (1 + length(null_devexpl)))
})

model_cols <- c("Tlarv + CUTIben + Chlaben" = "#8DA0CB",
                "DOegg + CUTIben + Chlaben" = "#66C2A5",
                "pHegg + CUTIben + Chlaben" = "#FC8D62")  
proj_de   <- proj_de %>% mutate(model = factor(model, levels = names(model_cols)))
lab_vec   <- setNames(
  sprintf("%s\n(dev. expl. = %.2f, p = %.3f)", proj_de$model, proj_de$dev_expl, proj_de$emp_p),
  as.character(proj_de$model))

fig_perm <- ggplot(tibble(dev_expl = null_devexpl), aes(dev_expl)) +
  geom_histogram(bins = 30, fill = "grey80", color = "grey55") +
  geom_vline(data = proj_de, aes(xintercept = dev_expl, color = model), linewidth = 1) +
  scale_color_manual(values = model_cols, labels = lab_vec, name = "Projected model") +
  labs(x = "Deviance explained (null = best model under selection)",
       y = "Permutations", title = "") +
  theme_classic(base_size = BASE_SIZE) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.key.height = unit(1.1, "lines"))
fig_perm
```

![](Manuscript_analysis_files/figure-gfm/permutation-figure-1.png)<!-- -->

``` r
ggsave("Figures/permutation_null.jpeg", fig_perm,
       width = 7.5, height = 4.5, dpi = 300)
```

------------------------------------------------------------------------

# 9. Sensitivity to early, low-precision recruitment years (1995–1997)

``` r
# Refit and re-select the univariate and multivariate GAMs after dropping low SE years

# --- univariate screen on an arbitrary data subset ---
uni_screen <- function(d_full) {
  map_dfr(pred_names, function(p) {
    d <- d_full %>% select(dev, x = all_of(p)) %>% filter(complete.cases(.))
    m <- gam(dev ~ s(x, k = 3), data = d, method = "REML")
    tibble(predictor = p,
           dev_expl  = summary(m)$dev.expl,
           AIC       = AIC(m),
           LOO_R2    = loo_gam(dev ~ s(x, k = 3), d)$r2)
  }) %>% arrange(desc(dev_expl))
}

# --- multivariate selection on an arbitrary data subset (same candidate set) ---
mv_screen <- function(d_full) {
  map_dfr(seq_along(filtered_combinations), function(i) {
    vars <- filtered_combinations[[i]]
    d  <- d_full %>% select(dev, recruit_year, all_of(vars)) %>% filter(complete.cases(.))
    fs <- as.formula(paste("dev ~", paste0("s(", vars, ", k = 3)", collapse = " + ")))
    m  <- gam(fs, data = d, method = "REML")
    tibble(model = paste(vars, collapse = " + "), n_cov = length(vars),
           AIC = AIC(m), dev_expl = summary(m)$dev.expl)
  }) %>% arrange(AIC) %>% mutate(rank = row_number(), deltaAIC = AIC - min(AIC))
}

# --- leave-year-out scenarios ---
drop_sets <- list(
  "Full data"      = integer(0),
  "Drop 1995"      = 1995,
  "Drop 1996"      = 1996,
  "Drop 1997"      = 1997,
  "Drop 1995-1997" = c(1995, 1996, 1997)
)

# --- univariate: best driver per scenario ---
sens_uni <- imap_dfr(drop_sets, function(yrs, nm) {
  d_full <- pred_std %>% filter(!recruit_year %in% yrs)
  best   <- uni_screen(d_full) %>% slice_max(dev_expl, n = 1, with_ties = FALSE)
  tibble(scenario = nm, n_years = nrow(d_full),
         top_predictor = best$predictor,
         dev_expl = round(best$dev_expl, 3), LOO_R2 = round(best$LOO_R2, 3))
})

# --- multivariate: best model + top-3 per scenario ---
sens_mv_all <- imap_dfr(drop_sets, function(yrs, nm) {
  d_full <- pred_std %>% filter(!recruit_year %in% yrs)
  mv_screen(d_full) %>% mutate(scenario = nm)
})
sens_mv_best <- sens_mv_all %>% filter(rank == 1) %>%
  transmute(scenario, top_model = model, dev_expl = round(dev_expl, 3), AIC = round(AIC, 1))
sens_mv_top3 <- sens_mv_all %>% filter(rank <= 3) %>%
  transmute(scenario, rank, model, dev_expl = round(dev_expl, 3), deltaAIC = round(deltaAIC, 2))

# --- the three manuscript models: rank + dev. explained in every scenario ---
paper_models <- list("Tlarv + CUTIben + Chlaben" = c("Tlarv","CUTIben","Chlaben"),
                     "DOegg + CUTIben + Chlaben" = c("DOegg","CUTIben","Chlaben"),
                     "pHegg + CUTIben + Chlaben" = c("pHegg","CUTIben","Chlaben"))
paper_keys <- sapply(paper_models, function(v) paste(sort(v), collapse = " + "))
sens_paper <- imap_dfr(drop_sets, function(yrs, nm) {
  d_full <- pred_std %>% filter(!recruit_year %in% yrs)
  r <- mv_screen(d_full)
  imap_dfr(paper_keys, function(key, pname) {
    row <- r %>% filter(model == key)
    tibble(scenario = nm, paper_model = pname, rank = row$rank,
           dev_expl = round(row$dev_expl, 3), deltaAIC = round(row$deltaAIC, 2))
  })
})
sens_paper_devexpl <- sens_paper %>%
  select(paper_model, scenario, dev_expl) %>%
  pivot_wider(names_from = scenario, values_from = dev_expl)
sens_paper_rank <- sens_paper %>%
  select(paper_model, scenario, rank) %>%
  pivot_wider(names_from = scenario, values_from = rank)
```

## 9a. Combined leave-year-out sensitivity table (CSV export)

``` r
# Self-contained LOO-CV RMSE helper (same logic as loo_gam above; redefined here
# so this chunk runs on its own as long as `pred_std` exists).
.loo_rmse <- function(form, d) {
  n <- nrow(d); err <- numeric(n)
  for (i in seq_len(n)) {
    mi     <- gam(form, data = d[-i, ], method = "REML")
    err[i] <- d$dev[i] - as.numeric(predict(mi, d[i, ]))
  }
  sqrt(mean(err^2))
}

# The three manuscript models
sens_models <- list(
  "Tlarv + CUTIben + Chlaben" = c("Tlarv", "CUTIben", "Chlaben"),
  "DOegg + CUTIben + Chlaben" = c("DOegg", "CUTIben", "Chlaben"),
  "pHegg + CUTIben + Chlaben" = c("pHegg", "CUTIben", "Chlaben")
)

# Leave-year-out scenarios ("Full data" is the baseline the changes are measured against)
sens_drop_sets <- list(
  "Full data"      = integer(0),
  "Drop 1995"      = 1995,
  "Drop 1996"      = 1996,
  "Drop 1997"      = 1997,
  "Drop 1995-1997" = c(1995, 1996, 1997)
)

# Fit one model on one scenario's data; return the three metrics
sens_fit_one <- function(vars, drop_years) {
  d <- pred_std %>%
    dplyr::select(dev, recruit_year, dplyr::all_of(vars)) %>%
    dplyr::filter(complete.cases(.), !recruit_year %in% drop_years)
  f <- as.formula(paste("dev ~", paste0("s(", vars, ", k = 3)", collapse = " + ")))
  m <- gam(f, data = d, method = "REML")
  tibble(dev_expl = summary(m)$dev.expl,
         AIC      = AIC(m),
         RMSE_loo = .loo_rmse(f, d))
}

# Every model x scenario combination
sens_raw <- tidyr::expand_grid(Model    = names(sens_models),
                               Scenario = names(sens_drop_sets)) %>%
  mutate(res = purrr::map2(sens_models[Model], sens_drop_sets[Scenario],
                           ~ sens_fit_one(.x, .y))) %>%
  tidyr::unnest(res)

# Full-data baseline per model, then change = scenario - baseline
sens_baseline <- sens_raw %>%
  filter(Scenario == "Full data") %>%
  select(Model, base_dev = dev_expl, base_AIC = AIC, base_RMSE = RMSE_loo)

sensitivity_change <- sens_raw %>%
  left_join(sens_baseline, by = "Model") %>%
  filter(Scenario != "Full data") %>%
  transmute(Model, Scenario,
            delta_dev_expl = round(dev_expl - base_dev, 3),
            delta_AIC      = round(AIC      - base_AIC, 2),
            delta_RMSE_loo = round(RMSE_loo - base_RMSE, 3)) %>%
  mutate(Model    = factor(Model,    levels = names(sens_models)),
         Scenario = factor(Scenario, levels = names(sens_drop_sets)[-1])) %>%
  arrange(Model, Scenario)

# Export the single table
readr::write_csv(sensitivity_change, "Figures/Sensitivity_leave_year_out.csv")

knitr::kable(
  sensitivity_change,
  col.names = c("Model", "Scenario",
                "Change in deviance explained", "Change in AIC", "Change in RMSE (LOO-CV)"))
```

| Model | Scenario | Change in deviance explained | Change in AIC | Change in RMSE (LOO-CV) |
|:---|:---|---:|---:|---:|
| Tlarv + CUTIben + Chlaben | Drop 1995 | -0.015 | -1.73 | -0.014 |
| Tlarv + CUTIben + Chlaben | Drop 1996 | -0.019 | 0.69 | 0.004 |
| Tlarv + CUTIben + Chlaben | Drop 1997 | 0.009 | 0.87 | 0.012 |
| Tlarv + CUTIben + Chlaben | Drop 1995-1997 | -0.034 | 0.07 | 0.003 |
| DOegg + CUTIben + Chlaben | Drop 1995 | -0.047 | 0.02 | 0.007 |
| DOegg + CUTIben + Chlaben | Drop 1996 | -0.020 | 0.76 | 0.014 |
| DOegg + CUTIben + Chlaben | Drop 1997 | 0.015 | -0.81 | -0.009 |
| DOegg + CUTIben + Chlaben | Drop 1995-1997 | -0.050 | -1.05 | -0.013 |
| pHegg + CUTIben + Chlaben | Drop 1995 | -0.027 | 0.57 | 0.024 |
| pHegg + CUTIben + Chlaben | Drop 1996 | -0.021 | 1.10 | 0.011 |
| pHegg + CUTIben + Chlaben | Drop 1997 | 0.015 | -0.05 | 0.004 |
| pHegg + CUTIben + Chlaben | Drop 1995-1997 | -0.041 | 0.23 | 0.003 |

------------------------------------------------------------------------

# 10. Leave-future-out validation (Table S2, Figure 3)

``` r
# ---- Model 1: Tlarv + CUTIben + Chlaben ----
d_m1 <- pred_std %>%
  select(recruit_year, dev, Tlarv, CUTIben, Chlaben) %>%
  filter(complete.cases(.)) %>%
  arrange(recruit_year)

# Training = all but the last 5 years; prediction = the final 5 years
horizon  <- 5
n_m1     <- nrow(d_m1)
train_m1 <- d_m1[1:(n_m1 - horizon), ]
test_m1  <- d_m1[(n_m1 - horizon + 1):n_m1, ]

# Fit the GAM on the training years
form_m1 <- dev ~ s(Tlarv, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
gam_m1  <- gam(form_m1, data = train_m1, method = "REML")

# Forecast the held-out final 5 years (obs, model prediction, and null)
pred_m1 <- data.frame(
  model        = "Tlarv + CUTIben + Chlaben",
  recruit_year = test_m1$recruit_year,
  obs          = test_m1$dev,
  pred         = as.numeric(predict.gam(gam_m1, newdata = test_m1)),
  null         = mean(train_m1$dev)
)

# Forecast skill for this model
rmse_m1  <- sqrt(mean((pred_m1$obs - pred_m1$pred)^2))
nrmse_m1 <- sqrt(mean((pred_m1$obs - pred_m1$null)^2))
skill_m1 <- data.frame(
  model      = "Tlarv + CUTIben + Chlaben",
  n_fore     = nrow(pred_m1),
  first      = min(pred_m1$recruit_year),
  last       = max(pred_m1$recruit_year),
  LFO_RMSE   = round(rmse_m1, 3),
  null_RMSE  = round(nrmse_m1, 3),
  RMSE_skill = round((nrmse_m1 - rmse_m1) / nrmse_m1, 3),
  R2         = round(1 - sum((pred_m1$obs - pred_m1$pred)^2) /
                         sum((pred_m1$obs - mean(pred_m1$obs))^2), 3)
)
skill_m1
```

    ##                       model n_fore first last LFO_RMSE null_RMSE RMSE_skill
    ## 1 Tlarv + CUTIben + Chlaben      5  2014 2018    0.378     0.561      0.325
    ##      R2
    ## 1 0.345

``` r
d_m2 <- pred_std %>%
  select(recruit_year, dev, pHegg, CUTIben, Chlaben) %>%
  filter(complete.cases(.)) %>%
  arrange(recruit_year)

horizon  <- 5
n_m2     <- nrow(d_m2)
train_m2 <- d_m2[1:(n_m2 - horizon), ]
test_m2  <- d_m2[(n_m2 - horizon + 1):n_m2, ]

form_m2 <- dev ~ s(pHegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
gam_m2  <- gam(form_m2, data = train_m2, method = "REML")

pred_m2 <- data.frame(
  model        = "pHegg + CUTIben + Chlaben",
  recruit_year = test_m2$recruit_year,
  obs          = test_m2$dev,
  pred         = as.numeric(predict.gam(gam_m2, newdata = test_m2)),
  null         = mean(train_m2$dev)
)

rmse_m2  <- sqrt(mean((pred_m2$obs - pred_m2$pred)^2))
nrmse_m2 <- sqrt(mean((pred_m2$obs - pred_m2$null)^2))
skill_m2 <- data.frame(
  model      = "pHegg + CUTIben + Chlaben",
  n_fore     = nrow(pred_m2),
  first      = min(pred_m2$recruit_year),
  last       = max(pred_m2$recruit_year),
  LFO_RMSE   = round(rmse_m2, 3),
  null_RMSE  = round(nrmse_m2, 3),
  RMSE_skill = round((nrmse_m2 - rmse_m2) / nrmse_m2, 3),
  R2         = round(1 - sum((pred_m2$obs - pred_m2$pred)^2) /
                         sum((pred_m2$obs - mean(pred_m2$obs))^2), 3)
)
skill_m2
```

    ##                       model n_fore first last LFO_RMSE null_RMSE RMSE_skill
    ## 1 pHegg + CUTIben + Chlaben      5  2014 2018     0.41     0.561      0.268
    ##      R2
    ## 1 0.231

``` r
d_m3 <- pred_std %>%
  select(recruit_year, dev, DOegg, CUTIben, Chlaben) %>%
  filter(complete.cases(.)) %>%
  arrange(recruit_year)

horizon  <- 5
n_m3     <- nrow(d_m3)
train_m3 <- d_m3[1:(n_m3 - horizon), ]
test_m3  <- d_m3[(n_m3 - horizon + 1):n_m3, ]

form_m3 <- dev ~ s(DOegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
gam_m3  <- gam(form_m3, data = train_m3, method = "REML")

pred_m3 <- data.frame(
  model        = "DOegg + CUTIben + Chlaben",
  recruit_year = test_m3$recruit_year,
  obs          = test_m3$dev,
  pred         = as.numeric(predict.gam(gam_m3, newdata = test_m3)),
  null         = mean(train_m3$dev)
)

rmse_m3  <- sqrt(mean((pred_m3$obs - pred_m3$pred)^2))
nrmse_m3 <- sqrt(mean((pred_m3$obs - pred_m3$null)^2))
skill_m3 <- data.frame(
  model      = "DOegg + CUTIben + Chlaben",
  n_fore     = nrow(pred_m3),
  first      = min(pred_m3$recruit_year),
  last       = max(pred_m3$recruit_year),
  LFO_RMSE   = round(rmse_m3, 3),
  null_RMSE  = round(nrmse_m3, 3),
  RMSE_skill = round((nrmse_m3 - rmse_m3) / nrmse_m3, 3),
  R2         = round(1 - sum((pred_m3$obs - pred_m3$pred)^2) /
                         sum((pred_m3$obs - mean(pred_m3$obs))^2), 3)
)
skill_m3
```

    ##                       model n_fore first last LFO_RMSE null_RMSE RMSE_skill
    ## 1 DOegg + CUTIben + Chlaben      5  2014 2018    0.412     0.561      0.265
    ##      R2
    ## 1 0.225

``` r
lfo_predictions <- bind_rows(pred_m1, pred_m2, pred_m3)
lfo_skill       <- bind_rows(skill_m1, skill_m2, skill_m3)

knitr::kable(lfo_skill, caption = "Leave-future-out forecast skill (final 5 years)")
```

| model                     | n_fore | first | last | LFO_RMSE | null_RMSE | RMSE_skill |    R2 |
|:--------------------------|-------:|------:|-----:|---------:|----------:|-----------:|------:|
| Tlarv + CUTIben + Chlaben |      5 |  2014 | 2018 |    0.378 |     0.561 |      0.325 | 0.345 |
| pHegg + CUTIben + Chlaben |      5 |  2014 | 2018 |    0.410 |     0.561 |      0.268 | 0.231 |
| DOegg + CUTIben + Chlaben |      5 |  2014 | 2018 |    0.412 |     0.561 |      0.265 | 0.225 |

Leave-future-out forecast skill (final 5 years)

## 10a. Figure 3 — in-sample fits with one-year-ahead forecasts

``` r
## full-data fits (gam_m1/m2/m3 were training-only), one per model.
fit_full_m1 <- gam(form_m1, data = d_m1, method = "REML")
p1 <- predict(fit_full_m1, se.fit = TRUE)
insample_m1 <- tibble(model = "Tlarv + CUTIben + Chlaben",
                      recruit_year = d_m1$recruit_year,
                      fit = as.numeric(p1$fit),
                      lo  = as.numeric(p1$fit - 2 * p1$se.fit),
                      hi  = as.numeric(p1$fit + 2 * p1$se.fit))

fit_full_m2 <- gam(form_m2, data = d_m2, method = "REML")
p2 <- predict(fit_full_m2, se.fit = TRUE)
insample_m2 <- tibble(model = "pHegg + CUTIben + Chlaben",
                      recruit_year = d_m2$recruit_year,
                      fit = as.numeric(p2$fit),
                      lo  = as.numeric(p2$fit - 2 * p2$se.fit),
                      hi  = as.numeric(p2$fit + 2 * p2$se.fit))

fit_full_m3 <- gam(form_m3, data = d_m3, method = "REML")
p3 <- predict(fit_full_m3, se.fit = TRUE)
insample_m3 <- tibble(model = "DOegg + CUTIben + Chlaben",
                      recruit_year = d_m3$recruit_year,
                      fit = as.numeric(p3$fit),
                      lo  = as.numeric(p3$fit - 2 * p3$se.fit),
                      hi  = as.numeric(p3$fit + 2 * p3$se.fit))

insample_fit <- bind_rows(insample_m1, insample_m2, insample_m3)

## Dots = static-block LOO forecasts for 2014-2018 (from the lfo-combine chunk)
lfo_fore <- lfo_predictions


ord <- c("DOegg + CUTIben + Chlaben",
         "pHegg + CUTIben + Chlaben",
         "Tlarv + CUTIben + Chlaben")
model_cols <- c("DOegg + CUTIben + Chlaben" = "#66C2A5",   # green
                "pHegg + CUTIben + Chlaben" = "#FC8D62",   # orange
                "Tlarv + CUTIben + Chlaben" = "#8DA0CB")   # purple

insample_fit <- mutate(insample_fit, model = factor(model, levels = ord))
lfo_fore     <- mutate(lfo_fore,     model = factor(model, levels = ord))

fig4 <- ggplot() +
  geom_ribbon(data = insample_fit,
              aes(recruit_year, ymin = lo, ymax = hi, fill = model), alpha = 0.20) +
  geom_line(data = insample_fit,
            aes(recruit_year, fit, color = model), linewidth = 0.9) +
  geom_point(data = insample_fit,
             aes(recruit_year, fit, color = model), alpha = .5) +
  geom_line(data = obs_dev, aes(recruit_year, dev),
            linetype = "dashed", color = "grey15") +
  geom_point(data = lfo_fore, aes(recruit_year, pred),
             shape = 21, fill = "red", color = "red",
             size = 2, stroke = 0.6) +
  scale_color_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  facet_wrap(~ model, nrow = 1) +
  labs(x = "Recruitment year", y = "Recruitment deviation") +
  theme_classic(base_size = BASE_SIZE) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none")
fig4
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggsave("Figures/fig4.jpeg", fig4, width = 10, height = 3.8, dpi = 300)
```

------------------------------------------------------------------------

# 11. Three top models: fits and partial residuals

``` r
md_tlarv <- pred_std %>% select(recruit_year, dev, Tlarv, CUTIben, Chlaben) %>% filter(complete.cases(.))
md_doegg <- pred_std %>% select(recruit_year, dev, DOegg, CUTIben, Chlaben) %>% filter(complete.cases(.))
md_pheg  <- pred_std %>% select(recruit_year, dev, pHegg, CUTIben, Chlaben) %>% filter(complete.cases(.))

m_tlarv <- gam(dev ~ s(Tlarv, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_tlarv, method = "REML")
m_doegg <- gam(dev ~ s(DOegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_doegg, method = "REML")
m_pheg  <- gam(dev ~ s(pHegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_pheg,  method = "REML")
summary(m_tlarv)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## dev ~ s(Tlarv, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.10863    0.05072  -2.142    0.046 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##              edf Ref.df      F  p-value    
    ## s(Tlarv)   1.156  1.284 11.542 0.003463 ** 
    ## s(CUTIben) 1.793  1.954 10.151 0.000986 ***
    ## s(Chlaben) 1.839  1.972  2.832 0.104982    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.594   Deviance explained = 67.9%
    ## -REML = 8.6208  Scale est. = 0.061735  n = 24

``` r
summary(m_doegg)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## dev ~ s(DOegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  -0.1086     0.0552  -1.968   0.0641 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##              edf Ref.df     F p-value   
    ## s(DOegg)   1.000  1.000 9.231 0.00676 **
    ## s(CUTIben) 1.734  1.927 5.837 0.01275 * 
    ## s(Chlaben) 1.609  1.845 0.735 0.44532   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.519   Deviance explained =   61%
    ## -REML =  9.666  Scale est. = 0.07314   n = 24

``` r
summary(m_pheg)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## dev ~ s(pHegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.10863    0.05371  -2.022   0.0579 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##              edf Ref.df      F p-value   
    ## s(pHegg)   1.000  1.000 10.050 0.00530 **
    ## s(CUTIben) 1.807  1.961  6.973 0.00574 **
    ## s(Chlaben) 1.803  1.960  2.168 0.15921   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.545   Deviance explained = 63.6%
    ## -REML = 9.6268  Scale est. = 0.069246  n = 24

``` r
# Partial-residual tables for each smooth term in a fitted GAM.
partial_tables <- function(m, d, vars, n_grid = 100) {

  term_fit  <- predict(m, type = "terms")      # each smooth's contribution
  work_res  <- residuals(m, type = "working")  # = response residuals for Gaussian/identity

  # Partial residuals = smooth contribution + model residual
  points <- map_dfr(vars, function(v) {
    tibble(term = v,
           x    = d[[v]],
           pr   = as.numeric(term_fit[, paste0("s(", v, ")")] + work_res),
           year = d$recruit_year)
  })

  # Smooth curve and band for each term, with the other drivers held at their mean
  lines <- map_dfr(vars, function(v) {
    grid <- as.data.frame(lapply(d[vars], function(z) rep(mean(z, na.rm = TRUE), n_grid)))
    grid[[v]] <- seq(min(d[[v]], na.rm = TRUE), max(d[[v]], na.rm = TRUE), length.out = n_grid)

    p   <- predict(m, newdata = grid, type = "terms", se.fit = TRUE)
    col <- paste0("s(", v, ")")
    tibble(term = v,
           x    = grid[[v]],
           fit  = as.numeric(p$fit[, col]),
           se   = as.numeric(p$se.fit[, col])) %>%
      mutate(lo = fit - 2 * se, hi = fit + 2 * se) %>%
      select(-se)
  })

  list(points = points, lines = lines)
}
```

``` r
partial_plot <- function(m, d, vars, title = "") {
  pa <- partial_tables(m, d, vars)
  pa$points <- pa$points %>% mutate(term = factor(term, levels = vars))
  pa$lines  <- pa$lines  %>% mutate(term = factor(term, levels = vars))

  subt <- sprintf("deviance explained = %.0f%%,  adj. R² = %.2f",
                  100 * summary(m)$dev.expl, summary(m)$r.sq)

  lab <- if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = pa$points, aes(x, pr, label = year),
                             size = 2.4, max.overlaps = 30, seed = 1, color = "grey35")
  else geom_text(data = pa$points, aes(x, pr, label = year),
                 size = 2.4, vjust = -0.6, color = "grey35")

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(data = pa$lines, aes(x, ymin = lo, ymax = hi), fill = "grey85") +
    geom_line(data = pa$lines, aes(x, fit), color = "steelblue", linewidth = 1) +
    geom_point(data = pa$points, aes(x, pr), size = 1.5, color = "grey25") +
    lab +
    facet_wrap(~ term, scales = "free_x") +
    labs(x = "Standardized driver value", y = "Partial residual",
         title = "", subtitle = subt) +
    theme_bw(base_size = BASE_SIZE)
}
```

``` r
p_tlarv <- partial_plot(m_tlarv, md_tlarv, c("Tlarv",  "CUTIben", "Chlaben"))
p_doegg <- partial_plot(m_doegg, md_doegg, c("DOegg",  "CUTIben", "Chlaben"))
p_pheg  <- partial_plot(m_pheg,  md_pheg,  c("pHegg",  "CUTIben", "Chlaben"))

no_x <- theme(axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())

full_fig <- (p_tlarv + no_x) / (p_doegg + no_x) / p_pheg
full_fig
```

![](Manuscript_analysis_files/figure-gfm/partial-resid-full-1.png)<!-- -->

``` r
ggsave("Figures/partial_resid_full_3x3.jpeg", full_fig,
       width = 9, height = 11, dpi = 300)
```

------------------------------------------------------------------------

# 12. Sign-prediction skill (Figure S8)

``` r
# list of top three models 
model_list <- list("DOegg + CUTIben + CHLAben" = m_doegg,
                   "Tlarv + CUTIben + CHLAben" = m_tlarv,
                   "pHegg + CUTIben + CHLAben" = m_pheg)
data_list  <- list("DOegg + CUTIben + CHLAben" = md_doegg,
                   "Tlarv + CUTIben + CHLAben" = md_tlarv,
                   "pHegg + CUTIben + CHLAben" = md_pheg)

fit3 <- imap_dfr(model_list, function(m, nm)
  tibble(recruit_year = data_list[[nm]]$recruit_year, model = nm,
         fitted = as.numeric(fitted(m))))

# Sign-skill uses the deviations from 1996 on (matches the fitted model years).
obs_dev_sign <- pred_std %>% select(recruit_year, dev) %>% filter(!is.na(dev))
```

``` r
# Univariate GAM fitted values for each driver used in the top models
uni_params <- c("DOegg", "Tlarv", "pHegg", "CUTIben", "Chlaben")
uni_fit <- map_dfr(uni_params, function(p) {
  d <- pred_std %>% select(recruit_year, dev, x = all_of(p)) %>% filter(complete.cases(.))
  m <- gam(dev ~ s(x, k = 3), data = d, method = "REML")
  tibble(recruit_year = d$recruit_year, model = p, type = "Univariate",
         fitted = as.numeric(fitted(m)))
})

fit_all <- bind_rows(fit3 %>% mutate(type = "Multi-driver"), uni_fit)

# Classify observed and predicted deviations: "increased" (>= 0) or "poor" (< 0)
sign_skill <- fit_all %>%
  left_join(obs_dev_sign, by = "recruit_year") %>%
  mutate(obs_sign  = if_else(dev    > 0, "increased", "poor"),
         pred_sign = if_else(fitted > 0, "increased", "poor"),
         correct   = obs_sign == pred_sign)

sign_summary <- bind_rows(
  sign_skill %>% group_by(type, model, observed = obs_sign) %>%
    summarise(n = n(), n_correct = sum(correct),
              accuracy = round(mean(correct), 2), .groups = "drop"),
  sign_skill %>% group_by(type, model) %>%
    summarise(observed = "overall", n = n(), n_correct = sum(correct),
              accuracy = round(mean(correct), 2), .groups = "drop")
) %>%
  mutate(type = factor(type, levels = c("Univariate", "Multi-driver"))) %>%
  arrange(type, model, observed)

knitr::kable(sign_summary)
```

| type         | model                     | observed  |   n | n_correct | accuracy |
|:-------------|:--------------------------|:----------|----:|----------:|---------:|
| Univariate   | CUTIben                   | increased |   8 |         5 |     0.62 |
| Univariate   | CUTIben                   | overall   |  24 |        21 |     0.88 |
| Univariate   | CUTIben                   | poor      |  16 |        16 |     1.00 |
| Univariate   | Chlaben                   | increased |   8 |         3 |     0.38 |
| Univariate   | Chlaben                   | overall   |  24 |        15 |     0.62 |
| Univariate   | Chlaben                   | poor      |  16 |        12 |     0.75 |
| Univariate   | DOegg                     | increased |   8 |         4 |     0.50 |
| Univariate   | DOegg                     | overall   |  24 |        16 |     0.67 |
| Univariate   | DOegg                     | poor      |  16 |        12 |     0.75 |
| Univariate   | Tlarv                     | increased |   8 |         4 |     0.50 |
| Univariate   | Tlarv                     | overall   |  24 |        15 |     0.62 |
| Univariate   | Tlarv                     | poor      |  16 |        11 |     0.69 |
| Univariate   | pHegg                     | increased |   8 |         4 |     0.50 |
| Univariate   | pHegg                     | overall   |  24 |        15 |     0.62 |
| Univariate   | pHegg                     | poor      |  16 |        11 |     0.69 |
| Multi-driver | DOegg + CUTIben + CHLAben | increased |   8 |         6 |     0.75 |
| Multi-driver | DOegg + CUTIben + CHLAben | overall   |  24 |        21 |     0.88 |
| Multi-driver | DOegg + CUTIben + CHLAben | poor      |  16 |        15 |     0.94 |
| Multi-driver | Tlarv + CUTIben + CHLAben | increased |   8 |         5 |     0.62 |
| Multi-driver | Tlarv + CUTIben + CHLAben | overall   |  24 |        19 |     0.79 |
| Multi-driver | Tlarv + CUTIben + CHLAben | poor      |  16 |        14 |     0.88 |
| Multi-driver | pHegg + CUTIben + CHLAben | increased |   8 |         6 |     0.75 |
| Multi-driver | pHegg + CUTIben + CHLAben | overall   |  24 |        21 |     0.88 |
| Multi-driver | pHegg + CUTIben + CHLAben | poor      |  16 |        15 |     0.94 |

``` r
sign_plot_df <- filter(sign_summary, observed != "overall")

fig_sign <- ggplot(sign_plot_df, aes(model, accuracy, fill = observed)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(increased = "#4C78A8", poor = "#E45756")) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(~ type, scales = "free_x", space = "free") +
  labs(x = NULL, y = "Sign-prediction accuracy", fill = "Observed year", title = "") +
  theme_bw(base_size = BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(fig_sign)
```

![](Manuscript_analysis_files/figure-gfm/sign-skill-figure-1.png)<!-- -->

``` r
ggsave("Figures/sign_prediction_accuracy.jpeg", fig_sign, width = 10, height = 4, dpi = 300)
```

------------------------------------------------------------------------

# 13. Future projections to 2100 (Figures 4 & 5, Table S3)

``` r
# The three top models refit for projection (identical to fits ealier but renamed for clarity
md_tlarv3 <- pred_std %>% select(recruit_year, dev, Tlarv, CUTIben, Chlaben) %>% filter(complete.cases(.))
md_doegg3 <- pred_std %>% select(recruit_year, dev, DOegg, CUTIben, Chlaben) %>% filter(complete.cases(.))
md_pheg3  <- pred_std %>% select(recruit_year, dev, pHegg, CUTIben, Chlaben) %>% filter(complete.cases(.))

m_tlarv3 <- gam(dev ~ s(Tlarv, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_tlarv3, method = "REML")
m_doegg3 <- gam(dev ~ s(DOegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_doegg3, method = "REML")
m_pheg3  <- gam(dev ~ s(pHegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_pheg3,  method = "REML")
```

``` r
proj_vars <- c("Tlarv", "DOegg", "pHegg", "CUTIben", "Chlaben")
hist_stats <- map_dfr(proj_vars, function(v) {
  x <- predictors[[v]]; mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  z <- (x - mu) / sdv
  tibble(driver = v, mu = mu, sd = sdv,
         zmin = min(z, na.rm = TRUE), zmax = max(z, na.rm = TRUE))
})
zt    <- function(x, v) { s <- hist_stats[hist_stats$driver == v, ]; (x - s$mu) / s$sd }
inenv <- function(z, v) { s <- hist_stats[hist_stats$driver == v, ]; z >= s$zmin & z <= s$zmax }

# Shared ESM colour palette (used in every ESM figure)
ESM_colors <- c(GFDL = "#798E87", HADL = "#C27D38", IPSL = "gold3")
```

``` r
future_ocean <- bind_rows(
  read.csv("Data/GFDL_monthly_averages_allyear_2000_2100_ALL_41to34p4.csv") %>% mutate(ESM = "GFDL"),
  read.csv("Data/HADL_monthly_averages_allyear_2000_2100_ALL_41to34p4.csv") %>% mutate(ESM = "HADL"),
  read.csv("Data/IPSL_monthly_averages_allyear_2000_2100_ALL_41to34p4.csv") %>% mutate(ESM = "IPSL")
) %>% mutate(Mnum = match(Month, month.abb)) 

future_cuti <- read.csv("Data/future_cuti_all_models.csv") %>%
  filter(between(latitude, 35, 41)) %>%
  group_by(ESM = toupper(model), Year = as.integer(year), Mnum = as.integer(month)) %>%
  summarise(CUTI = mean(upwelling_index, na.rm = TRUE), .groups = "drop")

f_Tlarv <- future_ocean %>% filter(Mnum %in% 3:4) %>% group_by(ESM, Year) %>%
  summarise(Tlarv = mean(Temperature_avg, na.rm = TRUE), .groups = "drop")
f_DOegg <- future_ocean %>% filter(Mnum %in% 1:3) %>% group_by(ESM, Year) %>%
  summarise(DOegg = mean(DO_avg, na.rm = TRUE), .groups = "drop")
f_pHegg <- future_ocean %>% filter(Mnum %in% 1:3) %>% group_by(ESM, Year) %>%
  summarise(pHegg = mean(pH_avg, na.rm = TRUE), .groups = "drop")
f_CUTIben <- future_cuti %>% filter(Mnum %in% 7:9) %>% group_by(ESM, Year) %>%
  summarise(CUTIben = mean(CUTI, na.rm = TRUE), .groups = "drop")
f_Chlaben <- future_ocean %>% filter(Mnum %in% 7:9) %>% group_by(ESM, Year) %>%
  summarise(Chlaben = mean(Chla_avg, na.rm = TRUE), .groups = "drop")

# Standardize future drivers with HISTORICAL mean/SD; keep names the models expect
future_std <- f_Tlarv %>%
  full_join(f_DOegg,   by = c("ESM", "Year")) %>%
  full_join(f_pHegg,   by = c("ESM", "Year")) %>%
  full_join(f_CUTIben, by = c("ESM", "Year")) %>%
  full_join(f_Chlaben, by = c("ESM", "Year")) %>%
  filter(Year >= 2000, complete.cases(.)) %>%
  transmute(ESM, Year,
            Tlarv   = zt(Tlarv,   "Tlarv"),
            DOegg   = zt(DOegg,   "DOegg"),
            pHegg   = zt(pHegg,   "pHegg"),
            CUTIben = zt(CUTIben, "CUTIben"),
            Chlaben = zt(Chlaben, "Chlaben"))
```

## 13a. Figure 4 — future ocean conditions by ESM

``` r
ocean_future_long <- future_ocean %>%
  transmute(ESM, Year,
            `DO (mg/L)` = DO_avg, pH = pH_avg, `Temperature (°C)` = Temperature_avg,
            `Chl (mg/m²)` = Chla_avg) %>%
  pivot_longer(-c(ESM, Year), names_to = "variable", values_to = "value") %>%
  bind_rows(future_cuti %>% transmute(ESM, Year, variable = "CUTI (m²/s)", value = CUTI)) %>%
  filter(Year >= 2000, Year <= 2100) %>%
  mutate(period = case_when(
           Year <= 2020 ~ "2000-2020",
           Year <= 2060 ~ "2020-2060",
           TRUE         ~ "2060-2100"),
         period = factor(period, levels = c("2000-2020", "2020-2060", "2060-2100")),
         ESM = factor(ESM, levels = c("GFDL", "IPSL", "HADL")),
         variable = factor(variable,
                           levels = c("DO (mg/L)", "pH", "Temperature (°C)", "CUTI (m²/s)", "Chl (mg/m²)")))

fig_future_box <- ggplot(ocean_future_long, aes(period, value, fill = ESM)) +
  geom_point(aes(color = ESM),
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8),
             size = 0.3, alpha = 0.18) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
               alpha = 0.75, outlier.shape = NA, linewidth = 0.35) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = ESM_colors) +
  scale_color_manual(values = ESM_colors) +
  labs(color = "", fill = "")+
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = BASE_SIZE) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = rel(0.9)),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom") +
  stat_summary(aes(group = ESM, color = ESM), fun = min, geom = "point",
               position = position_dodge(width = 0.8),
               shape = 23, size = 2.2, stroke = 0.3)+
  stat_summary(aes(group = ESM, color = ESM), fun = max, geom = "point",
               position = position_dodge(width = 0.8),
               shape = 23, size = 2.2, stroke = 0.3)
fig_future_box
```

![](Manuscript_analysis_files/figure-gfm/future-ocean-boxplot-1.png)<!-- -->

``` r
ggsave("Figures/future_ocean_boxplot.jpeg", fig_future_box, width = 9, height = 3.6, dpi = 300)
```

``` r
median_table <- ocean_future_long %>%
  group_by(variable, period, ESM) %>%
  summarise(median = median(value, na.rm = TRUE),
            n = sum(!is.na(value)), .groups = "drop") %>%
  arrange(variable, period, ESM)

median_wide <- median_table %>%
  select(-n) %>%
  pivot_wider(names_from = ESM, values_from = median) %>%
  arrange(variable, period)

print(median_wide, n = Inf)
```

    ## # A tibble: 15 × 5
    ##    variable         period      GFDL   IPSL   HADL
    ##    <fct>            <fct>      <dbl>  <dbl>  <dbl>
    ##  1 DO (mg/L)        2000-2020  5.88   5.79   5.84 
    ##  2 DO (mg/L)        2020-2060  5.92   5.74   6.06 
    ##  3 DO (mg/L)        2060-2100  5.42   5.45   6.43 
    ##  4 pH               2000-2020  7.79   7.80   7.80 
    ##  5 pH               2020-2060  7.71   7.70   7.77 
    ##  6 pH               2060-2100  7.49   7.54   7.64 
    ##  7 Temperature (°C) 2000-2020 12.2   12.6   12.6  
    ##  8 Temperature (°C) 2020-2060 13.0   13.3   14.1  
    ##  9 Temperature (°C) 2060-2100 13.4   14.8   15.7  
    ## 10 CUTI (m²/s)      2000-2020  0.800  0.8    0.801
    ## 11 CUTI (m²/s)      2020-2060  0.809  0.797  0.791
    ## 12 CUTI (m²/s)      2060-2100  0.906  0.815  0.756
    ## 13 Chl (mg/m²)      2000-2020  0.325  0.314  0.310
    ## 14 Chl (mg/m²)      2020-2060  0.316  0.311  0.318
    ## 15 Chl (mg/m²)      2060-2100  0.316  0.318  0.321

## 13b. Projected recruitment and extrapolation (Table S3, Figure 5)

``` r
future_pred <- future_std %>%
  mutate(
    `Tlarv + CUTIben + Chlaben` = as.numeric(predict(m_tlarv3, newdata = .)),
    `DOegg + CUTIben + Chlaben` = as.numeric(predict(m_doegg3, newdata = .)),
    `pHegg + CUTIben + Chlaben` = as.numeric(predict(m_pheg3,  newdata = .)),
    env_T = inenv(Tlarv, "Tlarv") & inenv(CUTIben, "CUTIben") & inenv(Chlaben, "Chlaben"),
    env_D = inenv(DOegg, "DOegg") & inenv(CUTIben, "CUTIben") & inenv(Chlaben, "Chlaben"),
    env_P = inenv(pHegg, "pHegg") & inenv(CUTIben, "CUTIben") & inenv(Chlaben, "Chlaben")
  )

proj_long <- future_pred %>%
  pivot_longer(c(`Tlarv + CUTIben + Chlaben`, `DOegg + CUTIben + Chlaben`, `pHegg + CUTIben + Chlaben`),
               names_to = "model", values_to = "pred") %>%
  mutate(in_env = case_when(
    model == "Tlarv + CUTIben + Chlaben" ~ env_T,
    model == "DOegg + CUTIben + Chlaben" ~ env_D,
    model == "pHegg + CUTIben + Chlaben" ~ env_P))

proj_long %>% group_by(model, ESM) %>%
  summarise(pct_extrapolated = round(mean(!in_env) * 100), .groups = "drop") %>%
  pivot_wider(names_from = ESM, values_from = pct_extrapolated) %>%
  knitr::kable()
```

| model                     | GFDL | HADL | IPSL |
|:--------------------------|-----:|-----:|-----:|
| DOegg + CUTIben + Chlaben |   21 |    6 |    4 |
| Tlarv + CUTIben + Chlaben |   47 |   87 |   62 |
| pHegg + CUTIben + Chlaben |   65 |   44 |   55 |

``` r
fig_proj <- ggplot(proj_long, aes(Year, pred, color = ESM)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(shape = in_env), size = 1.5) +
  geom_line(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4, linewidth = 1) +
  scale_color_manual(values = ESM_colors) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "within envelope", `FALSE` = "extrapolated"),
                     name = NULL) +
  facet_wrap(~ model, ncol = 1) +
  labs(color = "", fill = "")+
  labs(x = "Year", y = "Projected recruitment deviation", title = "") +
  theme_bw(base_size = BASE_SIZE)
fig_proj
```

![](Manuscript_analysis_files/figure-gfm/proj-figure-1.png)<!-- -->

``` r
ggsave("Figures/projected_recruitment.jpeg", fig_proj, width = 9, height = 8, dpi = 300)
```

------------------------------------------------------------------------

# 14. Diagnostics for the three projected models

## 14.1 Model 1 — `Tlarv + CUTIben + Chlaben`

``` r
# call model
m_tlarv3 <- gam(dev ~ s(Tlarv, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_tlarv3, method = "REML")
md_tlarv3
```

    ##    recruit_year        dev      Tlarv     CUTIben     Chlaben
    ## 1          1995  0.5680260  0.4145250 -0.88253163 -0.78756255
    ## 2          1996  0.3087000 -0.1208118 -1.02634325  0.40280498
    ## 3          1997  0.0133416 -1.5733468 -1.98895740  1.76812698
    ## 4          1998 -0.1074160  2.1620936 -0.31990106  1.41227483
    ## 5          1999 -0.1049780 -0.7959688 -0.45140085 -0.73879768
    ## 6          2000 -0.2997240 -0.4689689  0.16854843  0.31605284
    ## 7          2001  0.1698870 -0.2491260  1.32398111 -0.59138354
    ## 8          2002 -0.1225090 -0.5697488  0.89472764 -0.82503487
    ## 9          2003  0.3585030  0.7061463  0.29158704  0.87442735
    ## 10         2004 -0.4306280  0.2790103  0.01211821 -1.87044615
    ## 11         2005 -0.4336440  0.6690441  1.14465948 -1.52035414
    ## 12         2006 -0.1256660  0.8713119  1.86007138  0.90010076
    ## 13         2007 -0.5975890 -0.9825600  0.97870337 -0.02148589
    ## 14         2008 -0.3832460 -1.9608454 -0.05417421 -0.13132948
    ## 15         2009 -0.5731490 -1.0635331 -0.49347208 -0.48099160
    ## 16         2010 -0.1621180  0.8697893  0.37146505 -0.05861031
    ## 17         2011 -0.3476180 -0.5651663  0.60604636 -0.51990015
    ## 18         2012 -0.4497780 -0.8020198  0.25999769 -0.16798380
    ## 19         2013 -0.5673630 -0.9114658 -0.07811063 -1.74535361
    ## 20         2014  0.4953730  0.7055824 -1.92779806  1.13825654
    ## 21         2015  0.6440870  0.9991512 -1.59609690  0.65860302
    ## 22         2016 -0.1030920  0.8197790  0.28820990  1.29622122
    ## 23         2017  0.2970000  1.3969598 -0.44792734 -0.28198550
    ## 24         2018 -0.6535000  0.1701686  1.06659775  0.97435073

``` r
# cooks distance
md_tlarv3$cooks.distance <- cooks.distance(m_tlarv3)
md_tlarv3$leverage <-influence.gam(m_tlarv3)
md_tlarv3$residuals.standard <- scale(residuals.gam(m_tlarv3))
n <- length(md_tlarv3$residuals.standard )  # Number of data points
influence  <- md_tlarv3[(md_tlarv3$cooks.distance> 4 / n |md_tlarv3$cooks.distance> 1|md_tlarv3$cooks.distance> 0.5),]

simulationOutput <- simulateResiduals(fittedModel = m_tlarv3)
plot(simulationOutput)
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
dharma_plot <- function() {
  plot(simulationOutput)
}

png("Figures/model_1_res.jpeg", width = 7, height = 5, units = "in", res = 300)
plot(simulationOutput)
dev.off()
```

    ## png 
    ##   2

``` r
cooks_by_year <- ggplot(md_tlarv3, aes(x = recruit_year, y = cooks.distance)) +
  geom_segment(aes(xend = recruit_year, yend = 0), color = "grey40") +   # the stems
  geom_point(size = 2) +
  geom_hline(yintercept = 4 / n, linetype = "dashed", color = "red") +
  annotate("text", x = max(md_tlarv3$recruit_year), y = 4 / n,
           label = "4/n", color = "red", hjust = 1, vjust = -0.4, size = 3) +
  geom_text(data = influence, aes(label = recruit_year),
            color = "red", vjust = -0.9, size = 3) +               # label only influential years
  labs(title = "",
       x = "Recruitment year", y = "Cook's distance (leave-one-out)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cooks_by_year
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
ggsave("Figures/model_1_cooks.jpeg", cooks_by_year, width = 4, height = 4, dpi = 300)

influence_plot <- ggplot(md_tlarv3, aes(x = leverage, y = residuals.standard)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 2 * mean(md_tlarv3$leverage),
             linetype = "dotted", color = "grey60") +              # high-leverage guide (2 * mean hat)
  geom_point(aes(size = cooks.distance), alpha = 0.6) +
  geom_text(aes(label = recruit_year), size = 2.5, vjust = -0.9, color = "grey35") +
  scale_size_continuous(name = "Cook's D") +
  labs(title = "",
       x = "Leverage (hat value)", y = "Studentized residual") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
influence_plot
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
ggsave("Figures/model_1_influence.jpeg", influence_plot, width = 4, height = 4, dpi = 300)
```

## 14.2 Model 6 — `DOegg + CUTIben + Chlaben`

``` r
# call model
md_doegg3 <- pred_std %>% select(recruit_year, dev, DOegg, CUTIben, Chlaben) %>% filter(complete.cases(.))
m_doegg3 <- gam(dev ~ s(DOegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_doegg3, method = "REML")
 

# cooks distance
md_doegg3$cooks.distance <- cooks.distance(m_doegg3)
md_doegg3$leverage <-influence.gam(m_doegg3)
md_doegg3$residuals.standard <- scale(residuals.gam(m_doegg3))
n <- length(md_doegg3$residuals.standard )  # Number of data points
influence  <- md_doegg3[(md_doegg3$cooks.distance> 4 / n |md_doegg3$cooks.distance> 1|md_doegg3$cooks.distance> 0.5),]


simulationOutput <- simulateResiduals(fittedModel = m_doegg3)
plot(simulationOutput)
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
dharma_plot <- function() {
  plot(simulationOutput)
}

png("Figures/model_6_res.jpeg", width = 7, height = 5, units = "in", res = 300)
plot(simulationOutput)
dev.off()
```

    ## png 
    ##   2

``` r
cooks_by_year <- ggplot(md_doegg3, aes(x = recruit_year, y = cooks.distance)) +
  geom_segment(aes(xend = recruit_year, yend = 0), color = "grey40") +   # the stems
  geom_point(size = 2) +
  geom_hline(yintercept = 4 / n, linetype = "dashed", color = "red") +
  annotate("text", x = max(md_doegg3$recruit_year), y = 4 / n,
           label = "4/n", color = "red", hjust = 1, vjust = -0.4, size = 3) +
  geom_text(data = influence, aes(label = recruit_year),
            color = "red", vjust = -0.9, size = 3) +               # label only influential years
  labs(title = "",
       x = "Recruitment year", y = "Cook's distance (leave-one-out)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cooks_by_year
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
ggsave("Figures/model_6_cooks.jpeg", cooks_by_year, width = 4, height = 4, dpi = 300)


influence_plot <- ggplot(md_doegg3, aes(x = leverage, y = residuals.standard)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 2 * mean(md_doegg3$leverage),
             linetype = "dotted", color = "grey60") +              # high-leverage guide (2 * mean hat)
  geom_point(aes(size = cooks.distance), alpha = 0.6) +
  geom_text(aes(label = recruit_year), size = 2.5, vjust = -0.9, color = "grey35") +
  scale_size_continuous(name = "Cook's D") +
  labs(title = "",
       x = "Leverage (hat value)", y = "Studentized residual") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
influence_plot
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
ggsave("Figures/model_6_influence.jpeg", influence_plot, width = 4, height = 4, dpi = 300)
```

## 14.3 Model 2 — `pHegg + CUTIben + Chlaben`

``` r
# call model
md_pheg3 <- pred_std %>% select(recruit_year, dev, pHegg, CUTIben, Chlaben) %>% filter(complete.cases(.))
m_pheg3 <- gam(dev ~ s(pHegg, k = 3) + s(CUTIben, k = 3) + s(Chlaben, k = 3), data = md_pheg3, method = "REML")
 

# cooks distance
md_pheg3$cooks.distance <- cooks.distance(m_pheg3)
md_pheg3$leverage <-influence.gam(m_pheg3)
md_pheg3$residuals.standard <- scale(residuals.gam(m_pheg3))
n <- length(md_pheg3$residuals.standard )  # Number of data points
influence  <- md_pheg3[(md_pheg3$cooks.distance> 4 / n |md_pheg3$cooks.distance> 1|md_pheg3$cooks.distance> 0.5),]


simulationOutput <- simulateResiduals(fittedModel = m_pheg3)
plot(simulationOutput)
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
dharma_plot <- function() {
  plot(simulationOutput)
}

png("Figures/model_2_res.jpeg", width = 7, height = 5, units = "in", res = 300)
plot(simulationOutput)
dev.off()
```

    ## png 
    ##   2

``` r
cooks_by_year <- ggplot(md_pheg3, aes(x = recruit_year, y = cooks.distance)) +
  geom_segment(aes(xend = recruit_year, yend = 0), color = "grey40") +   # the stems
  geom_point(size = 2) +
  geom_hline(yintercept = 4 / n, linetype = "dashed", color = "red") +
  annotate("text", x = max(md_pheg3$recruit_year), y = 4 / n,
           label = "4/n", color = "red", hjust = 1, vjust = -0.4, size = 3) +
  geom_text(data = influence, aes(label = recruit_year),
            color = "red", vjust = -0.9, size = 3) +               # label only influential years
  labs(title = "",
       x = "Recruitment year", y = "Cook's distance (leave-one-out)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cooks_by_year
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
ggsave("Figures/model_2_cooks.jpeg", cooks_by_year, width = 4, height = 4, dpi = 300)

influence_plot <- ggplot(md_pheg3, aes(x = leverage, y = residuals.standard)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 2 * mean(md_pheg3$leverage),
             linetype = "dotted", color = "grey60") +              # high-leverage guide (2 * mean hat)
  geom_point(aes(size = cooks.distance), alpha = 0.6) +
  geom_text(aes(label = recruit_year), size = 2.5, vjust = -0.9, color = "grey35") +
  scale_size_continuous(name = "Cook's D") +
  labs(title = "",
       x = "Leverage (hat value)", y = "Studentized residual") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
influence_plot
```

![](Manuscript_analysis_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

``` r
ggsave("Figures/model_2_influence.jpeg", influence_plot, width = 4, height = 4, dpi = 300)
```
