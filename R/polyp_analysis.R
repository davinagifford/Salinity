# polyp_analysis.R
#
# analysis of polyp numbers across the experiment.
#
# Created: 20 November 2026
# Created by: Davina Gifford

# load libaries
library(tidyverse)
library(glmmTMB)
library(DHARMa)


# load data
data <- readRDS("var/polyp_tidy.rds")
treatment_data <- readRDS("var/treatment_tidy.rds")

# rename Treatment to salinity level for consistency with the polyp data
treatment_data <- treatment_data %>%
  rename(salinity_level = Treatment)

# # combine into a single df for analysis

data <- data %>%
  left_join(
    treatment_data,
    by = c(
      "salinity_level" = "salinity_level",
      "Day" = "Day"
    )
  )

# data wrangling ------

polyp_data <- data %>%
  group_by(salinity_level, replicate) %>%
  #arrange(Day, .by_group = TRUE) %>%
  mutate(
    dead = replace_na(dead, 0), # replace empty cells with zeroes
    
    # Initial population
    initial_polyps = polyp_number[Day == 1][1], # create a new column with the inital polyp numbers
    
    # Relative abundance since Day 1
    relative_polyp_abundance = if_else(
      initial_polyps > 0,
      polyp_number / initial_polyps,
      NA_real_
    ),
    
    # Previous observation
    previous_polyps = lag(polyp_number),
    days_elapsed = Day - lag(Day),
    
    # Net population change
    population_delta = polyp_number - previous_polyps,
    
    # Daily net growth rate
    growth_rate = if_else(
      previous_polyps > 0 &
        polyp_number > 0 &
        days_elapsed > 0,
      log(polyp_number / previous_polyps) / days_elapsed,
      NA_real_
    ),
    
    # Observed mortality (minimum recorded deaths)
    observed_deaths = dead
  ) %>%
  # budding ratio
  mutate(
    budding_ratio = if_else(
      polyp_number > 0,
      replace_na(buds, 0) / polyp_number,
      NA_real_
    )
  ) %>% 
  ungroup()

# Update polyp data
polyp_data <- polyp_data %>% 
  mutate(across(where(is.character), as.factor),
         plot_day = Day - 2, # create a column for plotting the data that reflects the exposure day as Day Zero
         replicate_id = interaction(salinity_level, replicate) # create unique ids for each replicate
         )


# remove the baseline data for analysis
polyp_data <- polyp_data %>% 
  filter(Day > 1) %>% 
  droplevels()

# visualise the data

ggplot(
  polyp_data,
  aes(
    x = Day,
    y = polyp_number,
    group = replicate_id
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ salinity_level) +
  geom_vline(
    xintercept = c(5, 10),
    linetype = "dashed",
    colour = "grey50"
  ) +
  labs(
    x = "Experimental day",
    y = "Number of polyps",
    title = "Polyp population trajectories"
  ) +
  theme_classic()

# summary statistics
polyp_data %>%
  group_by(salinity_level) %>%
  summarise(
    Mean = mean(polyp_number, na.rm = TRUE),
    Variance = var(polyp_number, na.rm = TRUE),
    Zeros = sum(polyp_number == 0, na.rm = TRUE), # understand number of zeroes in case need to consider zero inflation factors during modelling. 
    Observations = n(),
    .groups = "drop"
  )

# Models -------

# GLMM using poisson as response variable polyp_number is a count variable
# fixed effect of salinity level and phase of experiment, with random effect of replicate id

m1 <- glmmTMB(polyp_number ~ salinity_level * phase + (1|replicate_id), family = poisson, data = polyp_data)
summary(m1)

simulationOutput <- simulateResiduals(fittedModel = m1, plot = FALSE)
plotQQunif(simulationOutput, testUniformity = TRUE, testOutliers = TRUE, testDispersion = TRUE)

# Check convergence and Hessian
m1$fit$convergence
m1$sdr$pdHess

# account for difference in initial polyp numbers. Did not use an offset as expected abundance may not be proportional to initial abundance. 
m2 <- glmmTMB(polyp_number ~ salinity_level * phase + log(initial_polyps) + (1|replicate_id), family = poisson, data = polyp_data)
summary(m2)

simulationOutput2 <- simulateResiduals(fittedModel = m2, plot = FALSE)
plotQQunif(simulationOutput2, testUniformity = TRUE, testOutliers = TRUE, testDispersion = TRUE)

# Check convergence and Hessian
m2$fit$convergence
m2$sdr$pdHess

anova(m1, m2) # compare the two
# m2 performs better, but AIC is still very high. 

# consider separate time points for analysis rather than just using the phase variable 
# becomes a peicewise method.

polyp_data <- polyp_data %>%
  mutate(
    t_exp = pmin(pmax(Day - 1, 0), 5), # time for exposure phase
    t_rec = pmin(pmax(Day - 6, 0), 4), # time for recovery phase
    t_post = pmax(Day - 10, 0) # time for post-recovery phase
  )

m3 <- glmmTMB(polyp_number ~ salinity_level * t_exp + t_rec + t_post + log(initial_polyps) + (1|replicate_id), family = poisson(link = "log"), data = polyp_data)
summary(m3)

simulationOutput3 <- simulateResiduals(fittedModel = m3, plot = FALSE)
plotQQunif(simulationOutput3, testUniformity = TRUE, testOutliers = TRUE, testDispersion = TRUE)

# Check convergence and Hessian
m3$fit$convergence
m3$sdr$pdHess

anova(m2, m3)
# m3 slightly better

# lots of deviation in the diagnostics, so could still be improved. 

# include interactions between all the time points and the salinity levels.

m4 <- glmmTMB(polyp_number ~ salinity_level * (t_exp + t_rec + t_post) + log(initial_polyps) + (1|replicate_id), family = poisson(link = "log"), data = polyp_data) 
summary(m4)

simulationOutput4 <- simulateResiduals(fittedModel = m4, plot = FALSE)
plotQQunif(simulationOutput4, testUniformity = TRUE, testOutliers = TRUE, testDispersion = TRUE)

# Check convergence and Hessian
m4$fit$convergence
m4$sdr$pdHess

anova(m4, m3)
# m4 better

# check assumptions

# Simulate residuals
res <- simulateResiduals(
  fittedModel = m4,
  n = 1000
)

# Diagnostic plots
plot(res)

# Formal diagnostic tests
testDispersion(res)
testZeroInflation(res)
testUniformity(res)

# significant deviation in the KS test. 
# more zeroes than expected

# Example: inspect one replicate
one_id <- levels(polyp_data$replicate_id)[1]

one_rep <- polyp_data %>%
  filter(replicate_id == one_id) %>%
  arrange(Day)

# Plot raw time series
plot(
  one_rep$Day,
  one_rep$polyp_number,
  type = "b",
  xlab = "Day",
  ylab = "Polyp count"
)


# Compare slopes among treatments

# Exposure slopes
exp_slopes <- emtrends(
  m4,
  ~ salinity_level,
  var = "t_exp"
)

summary(exp_slopes)

# Recovery slopes
rec_slopes <- emtrends(
  m4,
  ~ salinity_level,
  var = "t_rec"
)

summary(rec_slopes)

# Post-recovery slopes
post_slopes <- emtrends(
  m4,
  ~ salinity_level,
  var = "t_post"
)

summary(post_slopes)


# Example: exposure slopes
exp_results <- as.data.frame(exp_slopes)

exp_results <- exp_results %>%
  mutate(
    Percent_Change_Per_Day =
      100 * (exp(t_exp.trend) - 1)
  )

print(exp_results)


# Pairwise comparisons during exposure
pairs(exp_slopes, adjust = "tukey")

# Pairwise comparisons during recovery
pairs(rec_slopes, adjust = "tukey")

# Pairwise comparisons during post-recovery
pairs(post_slopes, adjust = "tukey")


# remove 5ppt from the analysis
polyp_data <- polyp_data %>% 
  mutate(
    post_10 = t_post * (salinity_level == "10 ppt"),
    post_15 = t_post * (salinity_level == "15 ppt"),
    post_20 = t_post * (salinity_level == "20 ppt"),
    post_control =
      t_post * (salinity_level == "Control"),
    post_procedural =
      t_post *
      (salinity_level == "Procedural control")
  )

m5 <- glmmTMB(polyp_number ~ salinity_level * (t_exp + t_rec) +
                post_10 + post_15 + post_20 +
                post_control + post_procedural + 
                log(initial_polyps) + (1|replicate_id), 
              family = poisson(link = "log"), data = polyp_data) 
summary(m5)


m5$fit$convergence
m5$sdr$pdHess

anova(m5, m4)

# diagnostics
# Simulated residuals
res_5 <- simulateResiduals(
  fittedModel = m5,
  n = 1000
)

plot(res_5)

testDispersion(res_5)
testZeroInflation(res_5)
testUniformity(res_5)

# Residual patterns across time and treatments
plotResiduals(
  res_5,
  form = polyp_data$Day
)

plotResiduals(
  res_5,
  form = polyp_data$salinity_level
)


testCategorical(
  res_5,
  catPred = polyp_data$salinity_level
)



# test treatment specific dispersion

m6 <- glmmTMB(polyp_number ~ salinity_level * (t_exp + t_rec) +
                post_10 + post_15 + post_20 +
                post_control + post_procedural + 
                log(initial_polyps) + (1|replicate_id), 
              dispformula = ~ salinity_level,
              family = poisson(link = "log"), data = polyp_data) 
summary(m6)

m6$fit$convergence
m6$sdr$pdHess

anova(m5, m6)
AIC(m5, m6)

# no difference, keep m5

# Temporal Autocorrelation

# Test temporal autocorrelation separately
# within each replicate

replicate_ids <- unique(polyp_data$replicate_id)

autocorrelation_results <- lapply(
  replicate_ids,
  function(id) {
    
    idx <- which(
      polyp_data$replicate_id == id
    )
    
    # Extract residuals for this replicate
    res_subset <- recalculateResiduals(
      res_corrected,
      sel = idx
    )
    
    # Test temporal autocorrelation
    testTemporalAutocorrelation(
      res_subset,
      time = polyp_data$Day[idx],
      plot = FALSE
    )
  }
)

# Extract p-values
acf_pvalues <- sapply(
  autocorrelation_results,
  function(x) x$p.value
)

summary(acf_pvalues)

# there are some replicates with temporal autocorrelation
# lets find out which ones

# Number of significant tests
sum(acf_pvalues < 0.05, na.rm = TRUE)

# Proportion of significant tests
mean(acf_pvalues < 0.05, na.rm = TRUE)

# Identify affected replicates
data.frame(
  replicate_id = replicate_ids,
  p_value = acf_pvalues
) %>%
  filter(p_value < 0.05)


# Adjust p-values for multiple testing
acf_adjusted <- p.adjust(
  acf_pvalues,
  method = "BH"
)

# Create summary table
acf_results <- data.frame(
  replicate_id = replicate_ids,
  p_value = acf_pvalues,
  adjusted_p = acf_adjusted
)

# Display smallest p-values
acf_results %>%
  arrange(p_value) %>%
  print()

# Number remaining significant
sum(acf_adjusted < 0.05, na.rm = TRUE)


# compare salinity treatments during exposure

library(emmeans)

# Estimate exposure slopes
exposure_slopes <- emtrends(
  m5,
  ~ salinity_level,
  var = "t_exp"
)

# Display slopes and confidence intervals
summary(
  exposure_slopes,
  infer = c(TRUE, TRUE)
)


exposure_results <- as.data.frame(
  confint(exposure_slopes)
)

exposure_results <- exposure_results %>%
  mutate(
    Percent_Change = 100 * (
      exp(t_exp.trend) - 1
    ),
    Lower_Percent = 100 * (
      exp(asymp.LCL) - 1
    ),
    Upper_Percent = 100 * (
      exp(asymp.UCL) - 1
    )
  )

print(exposure_results)


# Compare each treatment against Control
exposure_contrasts <- contrast(
  exposure_slopes,
  method = "trt.vs.ctrl",
  ref = which(
    levels(model_data2$salinity_level) == "Control"
  ),
  adjust = "dunnettx"
)

summary(exposure_contrasts)



library(dplyr)
library(DHARMa)
library(purrr)

# Identify observation positions for each culture
replicate_indices <- split(
  seq_len(nrow(polyp_data)),
  polyp_data$replicate_id
)

# Test temporal autocorrelation separately for each replicate
autocorrelation_results <- map_dfr(
  names(replicate_indices),
  function(rep_id) {
    
    idx <- replicate_indices[[rep_id]]
    
    idx <- idx[order(polyp_data$Day[idx])]
    
    residual_subset <- recalculateResiduals(
      res_poisson,
      sel = idx
    )
    
    test_result <- testTemporalAutocorrelation(
      residual_subset,
      time = polyp_data$Day[idx],
      plot = FALSE
    )
    
    tibble(
      replicate_id = rep_id,
      DW_statistic = unname(test_result$statistic),
      p_value = test_result$p.value
    )
  }
)

# Adjust for testing 30 cultures
autocorrelation_results <- autocorrelation_results %>%
  mutate(
    adjusted_p = p.adjust(p_value, method = "BH")
  )

print(autocorrelation_results, n = 30)

autocorrelation_results %>%
  summarise(
    Total_Replicates = n(),
    Significant_Unadjusted =
      sum(p_value < 0.05, na.rm = TRUE),
    Significant_Adjusted =
      sum(adjusted_p < 0.05, na.rm = TRUE),
    Minimum_p = min(p_value, na.rm = TRUE),
    Minimum_Adjusted_p =
      min(adjusted_p, na.rm = TRUE)
  )

# sensitivity check to see if PC2 is disproportionally affecting the DHARMa uniformity result

m7 <- update(
  m5,
  data = polyp_data %>%
    filter(
      as.character(replicate_id) !=
        "Procedural control.2"
    ) %>%
    droplevels()
)

# Check convergence
m7$fit$convergence
m7$sdr$pdHess

# Repeat DHARMa diagnostics
set.seed(123)

res_7 <- simulateResiduals(
  m7,
  n = 1000
)

testUniformity(res_7)
testDispersion(res_7)
testZeroInflation(res_7)


# new model with random slopes to fix issues with uniformity


# Centre and scale experimental day
polyp_data <- polyp_data %>%
  mutate(
    Day_scaled = (Day - mean(Day)) / sd(Day)
  )

# Fit Poisson GLMM with replicate-specific time slopes
m8 <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 + Day_scaled | replicate_id),
  
  family = poisson(link = "log"),
  data = polyp_data
)

summary(m8)

# Check convergence
m8$fit$convergence
m8$sdr$pdHess

# Compare AIC values
AIC(m5, m8)

# that produced a correlation of 1 for the intercept/slope correlation

# next try 


m10 <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 | replicate_id) +
    (0 + Day_scaled | replicate_id),
  
  family = poisson(link = "log"),
  data = polyp_data
)

summary(m10)

# Check numerical stability
m10$fit$convergence
m10$sdr$pdHess

# Compare all three Poisson models
AIC(
  m5,
  m8,
  m10
)


# m8 better than m5 but still has issues. 

# do some sensitivity checks



exposure_intercept <- emtrends(
  m5,
  ~ salinity_level,
  var = "t_exp"
)

exposure_slope <- emtrends(
  m8,
  ~ salinity_level,
  var = "t_exp"
)

# Dunnett-adjusted comparisons against Control
contrast(
  exposure_intercept,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)
contrast(
  exposure_slope,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)


recovery_intercept <- emtrends(
  m5,
  ~ salinity_level,
  var = "t_rec"
)

recovery_slope <- emtrends(
  m8,
  ~ salinity_level,
  var = "t_rec"
)

contrast(
  recovery_intercept,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)

contrast(
  recovery_slope,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)



post_contrasts <- function(model) {
  
  b <- fixef(model)$cond
  V <- vcov(model)$cond
  
  treatments <- c(
    "10 ppt",
    "15 ppt",
    "20 ppt",
    "Procedural control"
  )
  
  terms <- c(
    "post_10",
    "post_15",
    "post_20",
    "post_procedural"
  )
  
  results <- lapply(seq_along(terms), function(i) {
    
    L <- setNames(rep(0, length(b)), names(b))
    
    L[terms[i]] <- 1
    L["post_control"] <- -1
    
    estimate <- sum(L * b)
    se <- sqrt(
      as.numeric(t(L) %*% V %*% L)
    )
    
    z <- estimate / se
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    
    data.frame(
      Treatment = treatments[i],
      Estimate = estimate,
      SE = se,
      z = z,
      p = p
    )
  })
  
  results <- dplyr::bind_rows(results)
  
  results$Holm_p <- p.adjust(
    results$p,
    method = "holm"
  )
  
  results
}

# Random-intercept model
post_contrasts(m5)

# Random-slope model
post_contrasts(m8)


# leave one replicate out senstitivity analysis


# refit the model excluding one replicate at a time

replicates <- unique(
  as.character(polyp_data$replicate_id)
)

loo_recovery <- map_dfr(replicates, function(rep_id) {
  
  # Remove all observations from one culture
  reduced_data <- polyp_data %>%
    filter(
      as.character(replicate_id) != rep_id
    ) %>%
    droplevels()
  
  # Refit the original random-intercept Poisson GLMM
  fit <- tryCatch(
    update(
      m5,
      data = reduced_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Omitted_Replicate = rep_id,
      Estimate = NA_real_,
      SE = NA_real_,
      Adjusted_p = NA_real_,
      Converged = FALSE
    ))
  }
  
  converged <- (
    fit$fit$convergence == 0 &&
      isTRUE(fit$sdr$pdHess)
  )
  
  if (!converged) {
    return(tibble(
      Omitted_Replicate = rep_id,
      Estimate = NA_real_,
      SE = NA_real_,
      Adjusted_p = NA_real_,
      Converged = FALSE
    ))
  }
  
  # Estimate recovery-phase slopes
  trends <- emtrends(
    fit,
    ~ salinity_level,
    var = "t_rec"
  )
  
  # Dunnett-adjusted comparisons with Control
  comparisons <- as.data.frame(
    summary(
      contrast(
        trends,
        method = "trt.vs.ctrl",
        ref = "Control",
        adjust = "dunnett"
      )
    )
  )
  
  # Extract 10 ppt vs Control
  result <- comparisons %>%
    filter(contrast == "10 ppt - Control")
  
  tibble(
    Omitted_Replicate = rep_id,
    Estimate = result$estimate,
    SE = result$SE,
    Adjusted_p = result$p.value,
    Converged = TRUE
  )
})

print(loo_recovery, n = 30)

# summaries the sensitivy analysis

loo_recovery_summary <- loo_recovery %>%
  summarise(
    Total_Fits = n(),
    
    Converged_Fits = sum(Converged),
    
    Minimum_Estimate = min(
      Estimate,
      na.rm = TRUE
    ),
    
    Maximum_Estimate = max(
      Estimate,
      na.rm = TRUE
    ),
    
    Median_Estimate = median(
      Estimate,
      na.rm = TRUE
    ),
    
    Significant_Fits = sum(
      Adjusted_p < 0.05,
      na.rm = TRUE
    ),
    
    Minimum_Adjusted_p = min(
      Adjusted_p,
      na.rm = TRUE
    ),
    
    Maximum_Adjusted_p = max(
      Adjusted_p,
      na.rm = TRUE
    )
  )

print(loo_recovery_summary)

# which replicates are most influential

loo_recovery %>%
  mutate(
    Change_From_Original = Estimate - (-0.406932)
  ) %>%
  arrange(
    desc(abs(Change_From_Original))
  ) %>%
  print(n = 30)
