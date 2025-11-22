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


# Plot plankton densities.

ggplot(plankton, aes(x = log_chla, y = log_zoo_density)) +
    geom_point(mapping = aes(fill = watercolumn_temp), alpha = 0.8, size = 3, shape = 21) +
    geom_smooth(
        color = "black", linetype = 1,
        method = "lm", se = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 17) +
    labs(
        x = "Chlorophyll a (log µg/L)",
        y = "Zooplankton biomass (log µg/L)",
        fill = "Temperature (°C)"
    )

# Test significance for filtered & in vivo measures

gam1 <- gam(log_zoo_density ~ log_chla, data = plankton)
gam2 <- gam(log_zoo_density ~ s(log_chla), data = plankton)

summary(gam1)
summary(gam2)

gam1 <- gam(log_zoo_density ~ log_invivo_chla, data = plankton)
gam2 <- gam(log_zoo_density ~ s(log_invivo_chla), data = plankton)

summary(gam1)
summary(gam2) # lower R2 than filtered measure


# Take temperature ntiles and nest regressions.
# Then fit simple regression and GAM on nested linear estimates ----

# For 7 ntiles:
ntiles <- ncut_temperature(7, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

plot_nest(cut_midtemp, ntiles_estimates)
ntiles_estimates |> arrange(desc(term)) |> print(n = 20)

lm7 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(lm7)

# For 14 ntiles:
ntiles <- ncut_temperature(14, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

lm14 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(lm14)

gam14 <- gam(estimate ~ s(cut_midtemp), method = "REML", data = filter(ntiles_estimates, term == "Effect"))
gam14 <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "Effect"))
summary(gam14)

lm14_intercepts <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Intercept"))
gam14_intercepts <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "Intercept"))

summary(lm14_intercepts)
summary(gam14_intercepts)

plot_nest(cut_midtemp, ntiles_estimates, plot_fits = TRUE)
plot(gam14, xlab = "Temperature (°C)", ylab = "Smoothed effect estimate")
plot(gam14_intercepts, xlab = "Temperature (°C)", ylab = "Smoothed intercept estimate")

# Linear test

gam2 <- gam(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "Effect"))
summary(gam2)

AIC(gam14, gam2)
gam.check(gam14)


# Take temperature cuts and nest regressions ----

pools <- cut_temperature(9, plankton, pool_cuts = TRUE)
pools_estimates <- nest_temperature(cut_midtemp, pools)
plot_nest(cut_midtemp, pools_estimates)
# High uncertainties at temperature extremes

drops <- cut_temperature(7, plankton, pool_cuts = FALSE)
drops_estimates <- nest_temperature(cut_midtemp, drops)
plot_nest(cut_midtemp, drops_estimates)
