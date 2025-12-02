# Data project scripts for limnology.
# Hierarchical temperature models.

library(mgcv)


# Import and prep lake pulse data.

lakepulse_rene <-
    read_csv(here("Data", "lakepulse_data_curated_mediancorrected_2021-03-23.csv")) |>
    arrange(lake_ID)

zooplankton <-
    read_delim(here("Data", "Zooplanktonbiomass.csv"), delim = ";") |>
    rename(lake_ID = Lake_ID) |>
    arrange(lake_ID)

plankton <- prep_plankton(lakepulse_rene, zooplankton)
plankton

max(plankton$watercolumn_temp) - min(plankton$watercolumn_temp)
summary(plankton$watercolumn_temp)
mean(plankton$watercolumn_temp)
sd(plankton$watercolumn_temp)


# Plankton densities ----

lm1 <- lm(log_zoo_density ~ log_chla, data = plankton)
lm2 <- lm(log_zoo_density ~ log_chla + I(log_chla^2), data = plankton)
gam1 <- gam(log_zoo_density ~ s(log_chla, k = 20), data = plankton)

summary(gam1)
gam.check(gam1)

plankton_plot <- ggplot(
    filter(plankton, !is.na(log_zoo_density)),
    aes(x = log_chla, y = log_zoo_density)
    ) +
    geom_point(mapping = aes(fill = watercolumn_temp), alpha = 0.6, size = 3, shape = 21) +
    geom_line(mapping = aes(y = fitted(lm1)),
    #geom_smooth(method = lm, alpha = 0.2, #se = FALSE,
        color = "black", linewidth = 1.2) +
    geom_line(mapping = aes(y = fitted(gam1)), color = "red", linewidth = 1.2) +
    scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 17) +
    labs(
        x = "Chlorophyll a (log µg/L)",
        y = "Zooplankton biomass (log µg/L)",
        fill = "Temperature (°C)"
    )

plankton_plot

# Test significance for filtered & in vivo measures

gam1 <- gam(log_zoo_density ~ log_chla, data = plankton)
gam2 <- gam(log_zoo_density ~ s(log_chla), data = plankton)

summary(gam1)
summary(gam2) # Significantly linear, though low R2

gam1 <- gam(log_zoo_density ~ log_invivo_chla, data = plankton)
gam2 <- gam(log_zoo_density ~ s(log_invivo_chla), data = plankton)

summary(gam1)
summary(gam2) # Even lower R2 than filtered measure

gam1 <- gam(log_chla ~ log_zoo_density, data = plankton)
gam2 <- gam(log_chla ~ s(log_zoo_density), data = plankton)

summary(gam1)
summary(gam2)


# Take temperature ntiles and nest regressions.

# For 7 (or 10, 11, 12) ntiles:
ntiles <- ncut_temperature(12, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

plot_nest(cut_midtemp, ntiles_estimates)
ntiles_estimates |> arrange(term) |> print(n = 20)


# Then fit simple regression and GAM on nested linear estimates ----

lm7 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(lm7)

gam7 <- gam(estimate ~ s(cut_midtemp, k = 3), data = filter(ntiles_estimates, term == "Effect"))
summary(gam7)

plot(gam7)
gam.check(gam7)

# For 14 ntiles:
ntiles <- ncut_temperature(14, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

lm14 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(lm14)

gam14 <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "Effect"))
gam14 <- gam(estimate ~ s(cut_midtemp, k = 5), method = "REML", data = filter(ntiles_estimates, term == "Effect"))
summary(gam14)

lm14_intercepts <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Intercept"))
gam14_intercepts <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "Intercept"))

summary(lm14_intercepts)
summary(gam14_intercepts)

plot_nest(cut_midtemp, ntiles_estimates, plot_fits = TRUE)
plot(gam14_intercepts, xlab = "Temperature (°C)", ylab = "Smoothed intercept estimate")
plot(gam14, xlab = "Temperature (°C)", ylab = "Smoothed effect estimate")

gam.check(gam14)


# Could fit polynomial regression instead of GAM, which
# is perhaps more sound with so few observations... ----

lm_lin <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))

lm_cubic <- lm(estimate ~ cut_midtemp + I(cut_midtemp^2) + I(cut_midtemp^3), data = filter(ntiles_estimates, term == "Effect"))
summary(lm_cubic) # Final model for effect estimates (on n = 7, n = 11, n = 14).

lm_intercepts <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Intercept"))

lm_quadintercepts <- lm(estimate ~ cut_midtemp + I(cut_midtemp^2), data = filter(ntiles_estimates, term == "Intercept"))
summary(lm_quadintercepts) # Final model for intercept estimates.

# Fitted polynomial regression estimates.
polyplot_nest1(ntiles_estimates)

# Smoothed polynomial errors.
polyplot_nest2(ntiles_estimates)


# Results table

lm1 <- lm(log_zoo_density ~ log_chla, data = plankton)
lm_results <-
    bind_rows(glance(lm1), glance(lm_quadintercepts), glance(lm_cubic)) |>
    mutate(
        name = c("All (linear)", "Intercepts (quadratic)",
        "Effects (cubic)"), .before = r.squared) |>
    select(name, df, p.value, adj.r.squared, nobs)

lm_results
write_csv(lm_results, "12n_nested_regression_outputs.csv")

tidy(lm1)
tidy(lm_quadintercepts)
tidy(lm_cubic)

tidy(gam14)
tidy(gam14_intercepts)


# WIP ----

# Linear test

gam2 <- gam(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(gam2)

AIC(gam14, gam2)
gam.check(gam14)

# Take temperature cuts and nest regressions

pools <- cut_temperature(9, plankton, pool_cuts = TRUE)
pools_estimates <- nest_temperature(cut_midtemp, pools)
plot_nest(cut_midtemp, pools_estimates)
# High uncertainties at temperature extremes

drops <- cut_temperature(7, plankton, pool_cuts = FALSE)
drops_estimates <- nest_temperature(cut_midtemp, drops)
plot_nest(cut_midtemp, drops_estimates)
