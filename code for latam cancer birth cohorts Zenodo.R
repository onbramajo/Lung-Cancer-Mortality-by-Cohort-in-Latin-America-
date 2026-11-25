# =============================================================================
# Birth-cohort patterns of lung cancer mortality in seven Latin American
# countries, 1997-2023
#
# 
#
# Reproduces all results, figures and tables of the paper.
# Data sources and folder structure are described in README.md.
#
#   0. Packages and settings
#   1. Deaths (WHO Mortality Database)
#   2. Population (UN WPP 2024) and death rates
#   3. Data quality: ill-defined causes (Table S1)
#   4. Age-standardised rates (Figure 1, Table S2)
#   5. Smoothed rates and ROMI (Figure 2, Figure S3)
#   6. Tests for non-linear cohort effects (Table S4)
#   7. Two-factor age-period and age-cohort models (Table S3)
#   8. Figures
#   9. Tables
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES AND SETTINGS
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(ungroup)
library(metR)
library(splines)
library(broom)

# Folders, relative to the project folder (where this script is)
data_dir   <- "data"
output_dir <- "output"
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# Countries, ordered from south to north by the latitude of the capital.
# Brazil starts in 2006, after the large fall in deaths from ill-defined causes.
countries <- tribble(
  ~country,     ~who_code, ~first_year,
  "Uruguay",    2460,      1997,
  "Argentina",  2020,      1997,
  "Chile",      2120,      1997,
  "Brazil",     2070,      2006,
  "Colombia",   2130,      1997,
  "Mexico",     2310,      1998,
  "Cuba",       2150,      2001
)
country_list <- countries$country
sexes        <- c("Females", "Males")

lung     <- "Lung (C33-C34)"
lung_ill <- "Lung + C76-C80"
groups   <- c(lung, lung_ill)

ill_codes <- c("C76", "C77", "C78", "C79", "C80")   # ill-defined and unspecified sites

# WHO world standard population 2000-2025 (Ahmad et al. 2001), per 100,000,
# ages 40-44 ... 85-89 and 90+
std_world <- tibble(age = seq(40, 90, 5),
                    w   = c(6590, 6040, 5370, 4550, 3720, 2960, 2210, 1520, 910, 440, 195))

# Smoothing with ungroup (pclm2D)
ungroup_scale  <- 1000            # deaths are multiplied by this before fitting
ungroup_kr     <- 2
ungroup_lambda <- c(1e7, 1e7)
rho_cap        <- 0.25            # ROMI capped at +/- 25% per year

# Models
last_year_tests  <- 2019          # tests A and B use data up to 2019
reference_cohort <- "1940-1944"
min_cells        <- 5             # cohorts observed in fewer cells are dropped

# Colours and labels
romi_breaks <- c(-Inf, -10, -5, -2, 0, 2, 5, 10, 15, Inf)
romi_labels <- c("< -10", "-10 to -5", "-5 to -2", "-2 to 0",
                 "0 to 2", "2 to 5", "5 to 10", "10 to 15", "> 15")
romi_colours <- c("< -10" = "black", "-10 to -5" = "grey35",
                  "-5 to -2" = "grey60", "-2 to 0" = "grey85",
                  "0 to 2" = "#2166AC", "2 to 5" = "#1A9850",
                  "5 to 10" = "#FFD92F", "10 to 15" = "#F46D43",
                  "> 15" = "#B2182B")
sex_colours <- c(Females = "#B2182B", Males = "#2166AC")

# Order countries and sexes in every table and figure
set_order <- function(d) {
  d %>% mutate(country = factor(country, levels = country_list),
               sex     = factor(sex, levels = sexes))
}


# -----------------------------------------------------------------------------
# 1. DEATHS (WHO MORTALITY DATABASE)
# -----------------------------------------------------------------------------
# Raw ICD-10 files (Morticd10_part1 ... part6) in data/.
# Age columns: Deaths13 = 35-39, Deaths14 = 40-44, ..., Deaths23 = 85-89,
# Deaths24 and Deaths25 = 90+.

who_files <- list.files(data_dir, pattern = "^morticd10_part", full.names = TRUE,
                        ignore.case = TRUE)
who_files <- who_files[!grepl("\\.(zip|pdf|docx?|xlsx?)$", who_files, ignore.case = TRUE)]

who_all <- map_dfr(who_files, function(f) {
  read_csv(f, show_col_types = FALSE,
           col_types = cols(Cause = col_character(), List = col_character(),
                            Admin1 = col_character(), SubDiv = col_character(),
                            Frmat = col_character(), Sex = col_character(),
                            .default = col_guess())) %>%
    filter(Country %in% countries$who_code, Sex %in% c("1", "2"),
           is.na(Admin1), is.na(SubDiv),
           List %in% c("103", "104", "10M"), as.numeric(Frmat) %in% c(0, 1, 2)) %>%
    select(Country, Year, Cause, Sex, Deaths13:Deaths25)
}) %>%
  left_join(countries, by = c("Country" = "who_code")) %>%
  mutate(code = substr(Cause, 1, 3))

deaths_long <- who_all %>%
  filter(Year >= first_year) %>%
  pivot_longer(Deaths13:Deaths25, names_to = "col", values_to = "deaths") %>%
  mutate(age    = pmin(5 * (as.numeric(sub("Deaths", "", col)) - 6), 90),
         deaths = replace_na(deaths, 0),
         sex    = ifelse(Sex == "1", "Males", "Females"))

sum_codes <- function(codes, label) {
  deaths_long %>%
    filter(code %in% codes) %>%
    group_by(country, year = Year, sex, age) %>%
    summarise(deaths = sum(deaths), .groups = "drop") %>%
    mutate(group = label)
}

deaths <- bind_rows(
  sum_codes(c("C33", "C34"), lung),
  sum_codes(c("C33", "C34", ill_codes), lung_ill)
) %>%
  complete(nesting(country, year, sex), age, group, fill = list(deaths = 0))


# -----------------------------------------------------------------------------
# 2. POPULATION (UN WPP 2024) AND DEATH RATES
# -----------------------------------------------------------------------------
# Single-age population files in data/:
#   WPP2024_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx
#   WPP2024_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx
# Sheet "Estimates", population in thousands.

read_wpp <- function(file, sex_label) {
  raw <- read_excel(file, sheet = "Estimates", col_names = FALSE, col_types = "text",
                    .name_repair = "minimal")
  header_row <- which(raw[[1]] == "Index")[1]
  names(raw) <- as.character(unlist(raw[header_row, ]))
  raw[-(1:header_row), ] %>%
    rename(country = `Region, subregion, country or area *`) %>%
    filter(country %in% country_list) %>%
    select(country, Year, `0`:`100+`) %>%
    pivot_longer(`0`:`100+`, names_to = "age_single", values_to = "pop") %>%
    mutate(age_single = as.numeric(sub("\\+", "", age_single)),
           pop  = 1000 * as.numeric(pop),
           year = as.numeric(Year),
           sex  = sex_label)
}

pop <- bind_rows(
  read_wpp(file.path(data_dir, "WPP2024_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx"), "Males"),
  read_wpp(file.path(data_dir, "WPP2024_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx"), "Females")
) %>%
  filter(age_single >= 35) %>%
  mutate(age = pmin(5 * (age_single %/% 5), 90)) %>%
  group_by(country, year, sex, age) %>%
  summarise(pop = sum(pop), .groups = "drop")

dat <- deaths %>%
  inner_join(pop, by = c("country", "year", "sex", "age")) %>%
  mutate(rate = deaths / pop)

# Years available by country (Uruguay has no data for 2011)
print(dat %>% group_by(country) %>%
        summarise(first = min(year), last = max(year),
                  missing = paste(setdiff(min(year):max(year), unique(year)), collapse = " ")))


# -----------------------------------------------------------------------------
# 3. DATA QUALITY: ILL-DEFINED CAUSES (TABLE S1)
# -----------------------------------------------------------------------------
# Ages 40+, both sexes, all available years.
# pct_ill_defined = deaths coded R00-R99 / all deaths
# pct_c76_c80     = cancer deaths coded C76-C80 / all cancer deaths (C00-C97)

data_quality <- who_all %>%
  filter(grepl("^[A-Z][0-9][0-9]$", code)) %>%                 # single causes, no totals
  pivot_longer(Deaths14:Deaths25, values_to = "deaths") %>%
  group_by(country, year = Year, first_year) %>%
  summarise(all_deaths    = sum(deaths, na.rm = TRUE),
            ill_defined   = sum(deaths[substr(code, 1, 1) == "R"], na.rm = TRUE),
            cancer_deaths = sum(deaths[code %in% sprintf("C%02d", 0:97)], na.rm = TRUE),
            c76_c80       = sum(deaths[code %in% ill_codes], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pct_ill_defined = round(100 * ill_defined / all_deaths, 1),
         pct_c76_c80     = round(100 * c76_c80 / cancer_deaths, 1))

write_csv(data_quality, file.path(output_dir, "data_quality_all_years.csv"))


# -----------------------------------------------------------------------------
# 4. AGE-STANDARDISED RATES (FIGURE 1, TABLE S2)
# -----------------------------------------------------------------------------
# Main standard: WHO world 2000-2025.
# Sensitivity: pooled population of the seven countries in 2010, both sexes.

std_latam <- pop %>%
  filter(year == 2010, age >= 40) %>%
  group_by(age) %>%
  summarise(w = sum(pop), .groups = "drop")

asmr_with <- function(std, label) {
  dat %>%
    filter(age >= 40) %>%
    left_join(std, by = "age") %>%
    group_by(country, sex, year, group) %>%
    summarise(asmr = 1e5 * sum(w * rate) / sum(w), deaths = sum(deaths), .groups = "drop") %>%
    mutate(standard = label)
}

asmr <- bind_rows(asmr_with(std_world, "WHO world"),
                  asmr_with(std_latam, "Latin America 2010"))
write_csv(asmr, file.path(output_dir, "asmr.csv"))

# Sensitivity: % change between the first and last year under each standard
pct_change <- asmr %>%
  filter(group == lung) %>%
  group_by(country, sex, standard) %>%
  summarise(change = 100 * (asmr[year == max(year)] / asmr[year == min(year)] - 1),
            .groups = "drop") %>%
  pivot_wider(names_from = standard, values_from = change)
cat("Largest difference in % change between standards (percentage points):",
    round(max(abs(pct_change$`WHO world` - pct_change$`Latin America 2010`)), 1), "\n")


# -----------------------------------------------------------------------------
# 5. SMOOTHED RATES AND ROMI (FIGURE 2, FIGURE S3)
# -----------------------------------------------------------------------------
# 5-year age groups 35-39 ... 85-89 are split into single ages with the
# two-dimensional penalised composite link model, and single ages 40-84 are
# kept. A missing year (Uruguay 2011) is filled with the mean of the two
# neighbouring years for the fit only; ROMI values that use it are set to NA.
#
# ROMI at age x in year t = -(log m(x, t) - log m(x, t - 1)), in %.
# Positive values mean that mortality falls.

smoothed <- list()
for (cn in country_list) {
  for (s in sexes) {
    for (g in groups) {
      
      d <- dat %>% filter(country == cn, sex == s, group == g, age <= 85)
      
      D <- d %>% select(age, year, deaths) %>%
        pivot_wider(names_from = year, values_from = deaths) %>% arrange(age)
      N <- d %>% select(age, year, pop) %>%
        pivot_wider(names_from = year, values_from = pop) %>% arrange(age)
      
      all_years <- min(d$year):max(d$year)
      gaps <- setdiff(all_years, d$year)
      for (yy in gaps) {
        D[[as.character(yy)]] <- (D[[as.character(yy - 1)]] + D[[as.character(yy + 1)]]) / 2
        N[[as.character(yy)]] <- (N[[as.character(yy - 1)]] + N[[as.character(yy + 1)]]) / 2
      }
      D <- as.matrix(D[, as.character(all_years)])
      N <- as.matrix(N[, as.character(all_years)])
      D[D == 0] <- 0.00001
      
      fit <- pclm2D(x = seq(35, 85, by = 5), y = D * ungroup_scale, nlast = 5, offset = N,
                    control = list(kr = ungroup_kr, lambda = ungroup_lambda))
      
      mx <- fit$fitted / ungroup_scale             # single ages 35-89
      colnames(mx) <- all_years
      
      smoothed[[paste(cn, s, g)]] <- as.data.frame(mx) %>%
        mutate(age = 35:89) %>%
        pivot_longer(-age, names_to = "year", values_to = "mx") %>%
        mutate(year = as.integer(year), country = cn, sex = s, group = g,
               gap_year = year %in% gaps) %>%
        filter(age >= 40, age <= 84)
      cat("Smoothed:", cn, s, g, "\n")
    }
  }
}
smoothed <- bind_rows(smoothed)

romi <- smoothed %>%
  arrange(country, sex, group, age, year) %>%
  group_by(country, sex, group, age) %>%
  mutate(rho = -(log(mx) - log(lag(mx))),
         touches_gap = gap_year | lag(gap_year)) %>%
  filter(year > min(year) + 1, year < max(year)) %>%       # drop edge years
  ungroup() %>%
  filter(!is.na(rho)) %>%
  mutate(romi = 100 * pmax(pmin(rho, rho_cap), -rho_cap),
         romi = ifelse(touches_gap, NA, romi),
         romi_cat = cut(romi, breaks = romi_breaks, labels = romi_labels))

write_csv(romi, file.path(output_dir, "romi.csv"))


# -----------------------------------------------------------------------------
# 6. TESTS FOR NON-LINEAR COHORT EFFECTS (TABLE S4)
# -----------------------------------------------------------------------------
# Unsmoothed data: 5-year age groups 40-44 ... 80-84, single years, to 2019.
# Test A: age + year (factor), with and without a natural spline of cohort.
# Test B: the same comparison, with an age-specific linear trend in year
#         added to both models.
# Quasi-Poisson F-tests. Deviance reduction = share of the deviance of the
# model without cohort that is removed by adding the cohort spline.

test_data <- dat %>%
  filter(age >= 40, age <= 80, year <= last_year_tests) %>%
  mutate(agef = factor(age), yearf = factor(year),
         yr = year - 2008, coh = (year + 0.5) - (age + 2.5))

run_tests <- function(d) {
  qp <- quasipoisson(link = "log")
  m_ap    <- glm(deaths ~ agef + yearf,                             offset = log(pop), family = qp, data = d)
  m_apc   <- glm(deaths ~ agef + yearf + ns(coh, df = 5),           offset = log(pop), family = qp, data = d)
  m_apx   <- glm(deaths ~ agef + yearf + agef:yr,                   offset = log(pop), family = qp, data = d)
  m_apxc  <- glm(deaths ~ agef + yearf + agef:yr + ns(coh, df = 5), offset = log(pop), family = qp, data = d)
  p_value <- function(a, b) anova(a, b, test = "F")$`Pr(>F)`[2]
  reduction <- function(a, b) (deviance(a) - deviance(b)) / deviance(a)
  tibble(dispersion = summary(m_apxc)$dispersion,
         p_A = p_value(m_ap, m_apc),   reduction_A = reduction(m_ap, m_apc),
         p_B = p_value(m_apx, m_apxc), reduction_B = reduction(m_apx, m_apxc))
}

cohort_tests <- list()
for (cn in country_list) {
  for (s in sexes) {
    for (g in groups) {
      d <- test_data %>% filter(country == cn, sex == s, group == g) %>% droplevels()
      cohort_tests[[paste(cn, s, g)]] <- run_tests(d) %>%
        mutate(country = cn, sex = s, group = g, .before = 1)
    }
  }
}
cohort_tests <- bind_rows(cohort_tests)
write_csv(cohort_tests, file.path(output_dir, "cohort_tests.csv"))


# -----------------------------------------------------------------------------
# 7. TWO-FACTOR AGE-PERIOD AND AGE-COHORT MODELS (TABLE S3)
# -----------------------------------------------------------------------------
# Unsmoothed data: 5-year age groups 40-44 ... 80-84, single years,
# 5-year birth cohorts (cohort = year - age).
#  - Age-period model:  age + year (factor)
#  - Age-cohort model:  age + cohort (factor), reference cohort 1940-1944
#  - Drift: age + linear cohort trend (cohort centred on its mean)
#  - Curvature: cohort coefficients minus their inverse-variance weighted
#    linear fit
#  - QAIC: deviance / c_hat + 2 x number of parameters, with c_hat the larger
#    of the two dispersion estimates (at least 1), used for both models
#  - Period effects: relative to the mean of the first three years
# Run twice: all years, and without 2020-2023.

prep_models <- function(d) {
  starts <- seq(1890, 1985, by = 5)
  d <- d %>%
    mutate(deaths      = round(deaths),
           age_group   = factor(paste0(age, "-", age + 4)),
           year_factor = factor(year),
           cohort_year = year - age)
  d$cohort_group <- cut(d$cohort_year, breaks = seq(1890, 1990, by = 5),
                        right = FALSE, labels = paste0(starts, "-", starts + 4))
  cells <- d %>% count(cohort_group)
  d <- d %>%
    filter(cohort_group %in% cells$cohort_group[cells$n >= min_cells]) %>%
    mutate(cohort_group = droplevels(cohort_group),
           cohort_c     = cohort_year - mean(cohort_year))
  d$age_group    <- relevel(d$age_group, ref = "40-44")
  d$cohort_group <- relevel(d$cohort_group, ref = reference_cohort)
  d$year_factor  <- relevel(d$year_factor, ref = as.character(min(d$year)))
  d
}

for (run in c("all_years", "without_2020_2023")) {
  
  comparison <- list(); period_effects <- list()
  curvature <- list(); drift <- list()
  
  for (cn in country_list) {
    for (s in sexes) {
      for (g in groups) {
        
        d <- dat %>% filter(country == cn, sex == s, group == g, age >= 40, age <= 80)
        if (run == "without_2020_2023") d <- d %>% filter(year <= 2019)
        d <- prep_models(d)
        id <- paste(cn, s, g)
        
        m_ap    <- glm(deaths ~ age_group + year_factor + offset(log(pop)),  family = quasipoisson, data = d)
        m_ac    <- glm(deaths ~ age_group + cohort_group + offset(log(pop)), family = quasipoisson, data = d)
        m_drift <- glm(deaths ~ age_group + cohort_c + offset(log(pop)),     family = quasipoisson, data = d)
        
        # Model comparison
        c_hat <- max(summary(m_ap)$dispersion, summary(m_ac)$dispersion, 1)
        qaic  <- function(m) deviance(m) / c_hat + 2 * length(coef(m))
        comparison[[id]] <- tibble(country = cn, sex = s, group = g, c_hat = c_hat,
                                   delta_QAIC = qaic(m_ap) - qaic(m_ac))   # > 0: age-cohort better
        
        # Period effects
        pe <- tidy(m_ap) %>% filter(grepl("year_factor", term)) %>%
          mutate(year = as.numeric(gsub("year_factor", "", term))) %>%
          select(year, estimate, std.error)
        period_effects[[id]] <- bind_rows(tibble(year = min(d$year), estimate = 0, std.error = 0), pe) %>%
          arrange(year) %>%
          mutate(estimate = estimate - mean(estimate[1:3]),
                 ci_lower = estimate - 1.96 * std.error,
                 ci_upper = estimate + 1.96 * std.error,
                 country = cn, sex = s, group = g)
        
        # Cohort effects and curvature
        ce <- tidy(m_ac) %>% filter(grepl("cohort_group", term)) %>%
          mutate(cohort = gsub("cohort_group", "", term)) %>%
          select(cohort, estimate, std.error)
        ce <- bind_rows(tibble(cohort = reference_cohort, estimate = 0, std.error = 0), ce) %>%
          mutate(cohort_mid = as.numeric(substr(cohort, 1, 4)) + 2) %>%
          arrange(cohort_mid)
        weights <- 1 / ce$std.error^2
        weights[is.infinite(weights)] <- max(weights[is.finite(weights)])   # reference cohort
        line <- lm(estimate ~ cohort_mid, data = ce, weights = weights)
        curvature[[id]] <- ce %>%
          mutate(curvature = estimate - fitted(line),
                 ci_lower  = curvature - 1.96 * std.error,
                 ci_upper  = curvature + 1.96 * std.error,
                 country = cn, sex = s, group = g)
        
        # Drift, % per year of birth
        drift[[id]] <- tidy(m_drift) %>% filter(term == "cohort_c") %>%
          transmute(country = cn, sex = s, group = g,
                    drift    = 100 * (exp(estimate) - 1),
                    drift_lo = 100 * (exp(estimate - 1.96 * std.error) - 1),
                    drift_hi = 100 * (exp(estimate + 1.96 * std.error) - 1))
      }
    }
  }
  
  write_csv(bind_rows(comparison),     file.path(output_dir, paste0("model_comparison_", run, ".csv")))
  write_csv(bind_rows(period_effects), file.path(output_dir, paste0("period_effects_", run, ".csv")))
  write_csv(bind_rows(curvature),      file.path(output_dir, paste0("cohort_curvature_", run, ".csv")))
  write_csv(bind_rows(drift),          file.path(output_dir, paste0("drift_", run, ".csv")))
}


# -----------------------------------------------------------------------------
# 8. FIGURES
# -----------------------------------------------------------------------------
# Saved as PNG (600 dpi) and PDF.

theme_paper <- theme_bw(base_size = 24) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 22),
        legend.title     = element_text(size = 22),
        legend.key.size  = unit(1.1, "cm"),
        strip.text       = element_text(size = 22),
        strip.background = element_rect(fill = "grey92"),
        axis.text        = element_text(size = 19, colour = "black"),
        axis.title       = element_text(size = 24),
        panel.grid.minor = element_blank())

save_figure <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}

plot_asmr <- function(d) {
  ggplot(d, aes(x = year, y = asmr, colour = country)) +
    geom_line(linewidth = 1.4) +
    facet_wrap(~ sex, scales = "free_y") +
    scale_x_continuous(breaks = seq(2000, 2020, 5)) +
    labs(x = "Year", y = "Deaths per 100,000, ages 40+\n(age-standardised)", colour = NULL) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 3))) +
    theme_paper
}

plot_romi <- function(d) {
  ggplot(d, aes(x = year, y = age)) +
    geom_tile(aes(fill = romi_cat)) +
    scale_fill_manual(values = romi_colours, drop = FALSE, breaks = romi_labels,
                      na.value = "white", name = "ROMI (%), positive = mortality falls") +
    geom_contour(aes(z = year - age), breaks = seq(1910, 1985, 5),
                 colour = "black", linewidth = 0.2) +
    geom_text_contour(aes(z = year - age), breaks = seq(1920, 1980, 10),
                      skip = 0, stroke = 0.2, size = 6) +
    facet_grid(sex ~ country) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(breaks = seq(2000, 2020, 10)) +
    labs(x = "Year", y = "Age") +
    guides(fill = guide_legend(nrow = 2, title.position = "top")) +
    theme_paper +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

plot_cohort <- function(d, y_title = "Cohort curvature\n(log rate ratio)") {
  ggplot(d, aes(x = cohort_mid, y = curvature, colour = sex)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 1.5, linewidth = 0.8,
                  position = position_dodge(width = 2)) +
    geom_line(linewidth = 1.2, position = position_dodge(width = 2)) +
    geom_point(size = 3.5, position = position_dodge(width = 2)) +
    facet_wrap(~ country, nrow = 2) +
    scale_x_continuous(breaks = seq(1920, 1980, 20)) +
    scale_colour_manual(values = sex_colours, name = NULL) +
    labs(x = "Birth cohort (midpoint of 5-year group)", y = y_title) +
    theme_paper +
    theme(legend.position = c(0.88, 0.22))
}

asmr_f      <- read_csv(file.path(output_dir, "asmr.csv"), show_col_types = FALSE) %>% set_order()
romi_f      <- read_csv(file.path(output_dir, "romi.csv"), show_col_types = FALSE) %>%
  mutate(romi_cat = factor(romi_cat, levels = romi_labels)) %>% set_order()
curv_all    <- read_csv(file.path(output_dir, "cohort_curvature_all_years.csv"), show_col_types = FALSE) %>% set_order()
curv_no_cov <- read_csv(file.path(output_dir, "cohort_curvature_without_2020_2023.csv"), show_col_types = FALSE) %>% set_order()
period_f    <- read_csv(file.path(output_dir, "period_effects_all_years.csv"), show_col_types = FALSE) %>% set_order()

# Figure 1: age-standardised death rates, lung cancer
save_figure(plot_asmr(asmr_f %>% filter(group == lung, standard == "WHO world")),
            "Figure_1", 13, 8)

# Figure 2: ROMI, lung cancer
save_figure(plot_romi(romi_f %>% filter(group == lung)), "Figure_2", 16, 11)

# Figure 3: cohort curvature, lung cancer
save_figure(plot_cohort(curv_all %>% filter(group == lung)), "Figure_3", 15, 10)

# Figure S1: cohort effects before removing the linear trend, lung cancer
save_figure(plot_cohort(curv_all %>% filter(group == lung) %>%
                          mutate(curvature = estimate,
                                 ci_lower  = estimate - 1.96 * std.error,
                                 ci_upper  = estimate + 1.96 * std.error),
                        y_title = "Cohort effect\n(log rate ratio vs 1940-1944)"),
            "Figure_S1", 15, 10)

# Figure S2: cohort curvature, lung cancer + C76-C80
save_figure(plot_cohort(curv_all %>% filter(group == lung_ill)), "Figure_S2", 15, 10)

# Figure S3: ROMI, lung cancer + C76-C80
save_figure(plot_romi(romi_f %>% filter(group == lung_ill)), "Figure_S3", 16, 11)

# Figure S4: cohort curvature, lung cancer, without 2020-2023
save_figure(plot_cohort(curv_no_cov %>% filter(group == lung)), "Figure_S4", 15, 10)

# Figure S5: period effects, lung cancer
p <- period_f %>% filter(group == lung) %>%
  ggplot(aes(x = year, y = estimate, colour = sex, fill = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ country, nrow = 2) +
  scale_x_continuous(breaks = seq(2000, 2020, 10)) +
  scale_colour_manual(values = sex_colours, name = NULL) +
  scale_fill_manual(values = sex_colours, name = NULL) +
  labs(x = "Year", y = "Period effect\n(log rate ratio vs first 3 years)") +
  theme_paper +
  theme(legend.position = c(0.88, 0.22))
save_figure(p, "Figure_S5", 15, 10)


# -----------------------------------------------------------------------------
# 9. TABLES
# -----------------------------------------------------------------------------

# Table S1: all deaths, % ill-defined (R00-R99), cancer deaths, % C76-C80,
# ages 40+, first available year and every 5 years to 2023
table_s1 <- data_quality %>%
  filter(year >= first_year) %>%
  group_by(country) %>%
  filter(year %in% c(min(year), 2000, 2005, 2010, 2015, 2020, 2023)) %>%
  ungroup() %>%
  set_order() %>%
  arrange(country, year) %>%
  select(country, year, all_deaths, pct_ill_defined, cancer_deaths, pct_c76_c80)
write_csv(table_s1, file.path(tab_dir, "Table_S1.csv"))

# Table S2: lung cancer deaths (ages 40+) and age-standardised rates (WHO world
# standard) in the first year, 2010 and 2023
table_s2 <- asmr_f %>%
  filter(group == lung, standard == "WHO world") %>%
  group_by(country, sex) %>%
  mutate(first_year = min(year)) %>%
  filter(year %in% c(first_year, 2010, 2023)) %>%
  mutate(when = case_when(year == first_year ~ "first_year_value",
                          year == 2010 ~ "2010", TRUE ~ "2023"),
         asmr = round(asmr, 1)) %>%
  ungroup() %>%
  select(country, sex, first_year, when, deaths, asmr) %>%
  pivot_wider(names_from = when, values_from = c(deaths, asmr)) %>%
  arrange(sex, country)
write_csv(table_s2, file.path(tab_dir, "Table_S2.csv"))

# Table S3: model comparison (all years and without 2020-2023) and drift
table_s3 <- read_csv(file.path(output_dir, "model_comparison_all_years.csv"), show_col_types = FALSE) %>%
  left_join(read_csv(file.path(output_dir, "model_comparison_without_2020_2023.csv"), show_col_types = FALSE) %>%
              select(country, sex, group, delta_QAIC_without_2020_2023 = delta_QAIC),
            by = c("country", "sex", "group")) %>%
  left_join(read_csv(file.path(output_dir, "drift_all_years.csv"), show_col_types = FALSE),
            by = c("country", "sex", "group")) %>%
  filter(group == lung) %>%
  set_order() %>%
  arrange(sex, country) %>%
  mutate(across(c(c_hat, delta_QAIC, delta_QAIC_without_2020_2023), ~ round(.x, 1)),
         across(starts_with("drift"), ~ round(.x, 2))) %>%
  select(country, sex, c_hat, delta_QAIC, delta_QAIC_without_2020_2023, drift, drift_lo, drift_hi)
write_csv(table_s3, file.path(tab_dir, "Table_S3.csv"))

# Table S4: tests A and B, lung cancer
table_s4 <- cohort_tests %>%
  filter(group == lung) %>%
  set_order() %>%
  arrange(sex, country) %>%
  mutate(across(starts_with("p_"), ~ ifelse(.x < 0.001, "<0.001", sprintf("%.3f", .x))),
         across(starts_with("reduction"), ~ round(100 * .x, 1))) %>%
  select(country, sex, p_A, reduction_A, p_B, reduction_B)
write_csv(table_s4, file.path(tab_dir, "Table_S4.csv"))

# Drift for lung cancer + C76-C80, reported in the text (Uruguay females)
print(read_csv(file.path(output_dir, "drift_all_years.csv"), show_col_types = FALSE) %>%
        filter(country == "Uruguay", sex == "Females"))