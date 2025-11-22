# Data project scripts for limnology.
# Hierarchical temperature models.


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


# Take temperature ntiles and nest regressions.
# Then fit simple regression and GAM on nested linear estimates.

# For 7 ntiles:
ntiles <- ncut_temperature(7, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

plot_nest(cut_midtemp, ntiles_estimates)
ntiles_estimates |> arrange(desc(term)) |> print(n = 20)

lm7 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "adjusted_log_chla"))
summary(lm7)

# For 14 ntiles:
ntiles <- ncut_temperature(14, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

lm14 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "adjusted_log_chla"))
summary(lm14)

gam14 <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "adjusted_log_chla"))
summary(gam1)
plot_nest(cut_midtemp, ntiles_estimates, plot_fits = TRUE)
plot(gam1, xlab = "Temperature (°C)", ylab = "Smoothed effect estimate")

gam2 <- gam(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "adjusted_log_chla"))
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
