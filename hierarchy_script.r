#
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


# Take temperature ntiles and nest regressions.

ntiles <- ncut_temperature(12, plankton)
ntiles_estimates <- nest_temperature(cut_midtemp, ntiles)

plot_nest(cut_midtemp, ntiles_estimates)
ntiles_estimates |> arrange(desc(term)) |> print(n = 20)

# Fit models on nested estimates.
lm1 <- lm(estimate ~ cut_midtemp, data = filter(ntiles_estimates, term == "adjusted_log_chla"))
summary(lm1)

gam1 <- gam(estimate ~ s(cut_midtemp), data = filter(ntiles_estimates, term == "adjusted_log_chla"))
summary(gam1)

plot_nest(cut_midtemp, ntiles_estimates, plot_fits = TRUE)


# Take temperature cuts and nest regressions ----

pools <- cut_temperature(7, plankton, pool_cuts = TRUE)
drops <- cut_temperature(7, plankton, pool_cuts = FALSE)

pools_estimates <- nest_temperature(cut_midtemp, pools)
drops_estimates <- nest_temperature(cut_midtemp, drops)
