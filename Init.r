# Data project scripts for Limnology.
# This script imports Lake Pulse environment & plankton data.

library(patchwork)
library(tidyverse)
library(here)
library(httpgd)

hgd()
hgd_browse()

theme_set(theme_bw())


# Import Lake Pulse ----

lakepulse_rene <- 
    read_csv(here("lakepulse_data_curated_mediancorrected_2021-03-23.csv"))

zooplankton <- read_delim("Zooplanktonbiomass.csv", delim = ";")
# Zooplankton csv delimits with ";" not ","

lakepulse_rene
str(lakepulse_rene)

zooplankton
str(zooplankton)
ncol(zooplankton) # 102 species found (& site col)

zooplankton <- zooplankton |> arrange(Lake_ID)

# Subset columns of interest.

plankton <- lakepulse_rene |> select(
    lake_ID, lake_name:longitude, ecozone:size_class,
    watershed_area, max_depth,
    rbr_chla_mean_tube, rbr_chla_mean_bottom, chla,
    rbr_temperature_mean_tube, rbr_temperature_mean_watercolumn,
    rbr_temperature_mean_epilimnion, rbr_temperature_mean_hypolimnion,
    tp_tube:tntp_ratio
    )

str(plankton)

plankton |> count(is.na(rbr_temperature_mean_epilimnion))
plankton |> count(is.na(rbr_temperature_mean_hypolimnion))
# All temperatures are value-complete.


# Zooplankton processing ----

zooplankton_pa <-
    mutate(zooplankton, across(!Lake_ID, \(x) replace(x, x != 0, 1))) |>
    rownames_to_column("order")

zooplankton_richness <-
    rowSums(select(zooplankton_pa, !c(order, Lake_ID))) |>
    as_tibble() |>
    rename(zoo_richness = value) |>
    rownames_to_column("order")

richness <- zooplankton_pa |>
    select(order, Lake_ID) |>
    left_join(zooplankton_richness) |>
    rename(lake_ID = Lake_ID) |>
    select(!order)

plankton <- left_join(plankton, richness)

plankton
plankton |> count(is.na(zoo_richness)) # Only 40 NAs in richness

lake_id <- zooplankton |>
    rownames_to_column("order") |>
    select(order, Lake_ID)

lake_biomass <- rowSums(zooplankton[-1]) |>
    as_tibble() |>
    rownames_to_column("order")

zoo_density <- left_join(lake_id, lake_biomass)[-1] |>
    rename(zoo_density = value, lake_ID = Lake_ID)

zoo_density

plankton <-
    left_join(plankton, zoo_density) |>
    rename(
        invivo_chla_surface = rbr_chla_mean_tube,
        invivo_chla_bottom = rbr_chla_mean_bottom,
        filtered_chla_surface = chla,
        surface_temp = rbr_temperature_mean_tube,
        watercolumn_temp = rbr_temperature_mean_watercolumn,
        epilimnion_temp = rbr_temperature_mean_epilimnion,
        hypolimnion_temp = rbr_temperature_mean_hypolimnion
    )

plankton
rm(lake_biomass, lake_id, richness, zoo_density, zooplankton_richness)

#write_csv(plankton, "plankton.csv")


# Helper functions ----

scatterplot_chla <- function(xvar = lake_ID, data = plankton) {

    p1 <- ggplot(data, aes(y = log_filtered_chla_surface, x = {{xvar}})) +
        geom_point(alpha = 0.6) +
        geom_smooth() +
        theme(axis.title.x = element_blank())

    p2 <- ggplot(data, aes(y = log_invivo_chla_surface, x = {{xvar}})) +
        geom_point(alpha = 0.6) +
        geom_smooth()
    
    p3 <- ggplot(data, aes(y = log_invivo_chla_bottom, x = {{xvar}})) +
        geom_point(alpha = 0.6) +
        geom_smooth() +
        theme(axis.title.x = element_blank())

    p1 + p2 + p3

}

scatterplot <- function(xvar = surface_temp, yvar = log_zoo_density, data = plankton) {

    ggplot(data, aes(x = {{xvar}}, y = {{yvar}})) +
        geom_point(alpha = 0.6) #+
        #geom_smooth()

}

# Temperature nester
temperature_nest <- function(nesting_var = surface_temp_bin, input_data = plankton) {

    input_data |>
        nest_by({{nesting_var}}) |>
        mutate(
            mean_chla = mean(data$log_filtered_chla_surface),
            data_c = list(mutate(data, adjusted_log_chla = log_filtered_chla_surface - mean_chla)),
            models = list(lm(log_zoo_density~adjusted_log_chla, data = data_c)),
            tidymodels = list(broom::tidy(models))
        ) |>
        select(!c(data, models)) |>
        unnest(tidymodels)

}

nest_plot <- function(nesting_var = surface_temp_bin, data_nest = surface_nest) {

    effect_plot <- ggplot(filter(data_nest, term == "adjusted_log_chla"), aes(y = {{nesting_var}}, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0)
    
    print(effect_plot)

    med_ntile_temp <- median(data_nest$ntile_midtemp)
    print(med_ntile_temp)

    intercept_effect_plot1 <- ggplot(data_nest, aes(x = {{nesting_var}}, y = estimate)) +
        geom_point(aes(fill = {{nesting_var}}), size = 3.5, shape = 21) +
        geom_smooth(se = FALSE, span = 0.8) +
        geom_smooth(method = lm, color = "black", linewidth = 1, alpha = 0.2) +
        scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = med_ntile_temp) +
        facet_wrap(~term, scales = "free")

    print(intercept_effect_plot1)

    intercept_effect_plot2 <- data_nest |> select({{nesting_var}}, term, estimate) |>
        pivot_wider(names_from = term, values_from = estimate) |>
        ggplot(mapping = aes(x = `(Intercept)`, y = adjusted_log_chla)) +
        geom_point() +
        geom_smooth(method = lm)
    
    print(intercept_effect_plot2)

    data_nest <- filter(data_nest, term == "adjusted_log_chla")

    print(data_nest)
    print(median(data_nest$p.value))

}

# This pools the values from the highest temperature cuts so n > 5 per cut.
cut_nest_plot1 <- function(n = 10, input_data = plankton) {

    temperature_cuts <- input_data |>
        mutate(temperature_range = cut(watercolumn_temp, n)) |>
        group_by(temperature_range) |>
        summarize(
            min = min(watercolumn_temp), max = max(watercolumn_temp),
            interval = max - min, n = n()
        ) |>
        rownames_to_column(var = "temperature_cut") |>
        mutate(temperature_cut = as.numeric(temperature_cut))

    print(temperature_cuts)

    pooled_cuts <- filter(temperature_cuts, n < 5) |> select(temperature_cut) |> t() |> as.vector()
    max_cut <- max(filter(temperature_cuts, n >= 5)$temperature_cut)
    temperature_cuts <- mutate(temperature_cuts, temperature_cut = if_else(
        temperature_cut %in% pooled_cuts,
        max_cut, temperature_cut)) |>
        select(temperature_cut, temperature_range)

    tempcut_plankton <-
        mutate(input_data, temperature_range = cut(watercolumn_temp, n)) |>
        left_join(temperature_cuts, by = "temperature_range")

    tempcut_plankton |> count(temperature_range, temperature_cut) |> print()
    
    tempcut_nest <- temperature_nest(temperature_cut, input_data = tempcut_plankton)
    tempcut_nest |> arrange(term) |> print(n = 40)

    estimate_plot <- ggplot(filter(tempcut_nest, term == "adjusted_log_chla"), aes(y = temperature_cut, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0)

    intercept_estimate_plot1 <- ggplot(tempcut_nest, aes(x = temperature_cut, y = estimate)) +
        geom_point(size = 3) +
        #geom_smooth(se = FALSE, span = 0.8) +
        geom_smooth(se = FALSE, method = lm, color = "red") +
        facet_wrap(~term, scales = "free")

    print(estimate_plot)
    print(intercept_estimate_plot1)

    cut_estimate_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "adjusted_log_chla"))
    print(summary(cut_estimate_lm))

    cut_intercept_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "(Intercept)"))
    print(summary(cut_intercept_lm))

}

# This drops the highest temperature ranges so n > 5 per cut.
cut_nest_plot2 <- function(n = 10, input_data = plankton) {

    temperature_cuts <- input_data |>
        mutate(temperature_range = cut(watercolumn_temp, n)) |>
        group_by(temperature_range) |>
        summarize(
            min = min(watercolumn_temp), max = max(watercolumn_temp),
            interval = max - min, n = n()
        ) |>
        rownames_to_column(var = "temperature_cut") |>
        mutate(temperature_cut = as.numeric(temperature_cut))

    print(temperature_cuts)

    temperature_cuts <-
        filter(temperature_cuts, n > 6) |>
        select(temperature_cut, temperature_range)

    tempcut_plankton <-
        mutate(input_data, temperature_range = cut(watercolumn_temp, n)) |>
        left_join(temperature_cuts, by = "temperature_range")

    tempcut_plankton |> count(temperature_range, temperature_cut) |> print()

    tempcut_plankton <- filter(tempcut_plankton, !is.na(temperature_cut))
    
    tempcut_nest <- temperature_nest(temperature_cut, input_data = tempcut_plankton)
    tempcut_nest |> arrange(term) |> print()

    estimate_plot <- ggplot(filter(tempcut_nest, term == "adjusted_log_chla"), aes(y = temperature_cut, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0)

    intercept_estimate_plot1 <- ggplot(tempcut_nest, aes(x = temperature_cut, y = estimate)) +
        geom_point(size = 3) +
        #geom_smooth(se = FALSE, span = 0.8) +
        geom_smooth(se = FALSE, method = lm, color = "red") +
        facet_wrap(~term, scales = "free")

    print(estimate_plot)
    print(intercept_estimate_plot1)

    cut_estimate_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "adjusted_log_chla"))
    print(summary(cut_estimate_lm))

    cut_intercept_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "(Intercept)"))
    print(summary(cut_intercept_lm))

}

cut_nest_plotgam <- function(n = 10, input_data = plankton) {

    temperature_cuts <- input_data |>
        mutate(temperature_range = cut(watercolumn_temp, n)) |>
        group_by(temperature_range) |>
        summarize(
            min = min(watercolumn_temp), max = max(watercolumn_temp),
            interval = max - min, n = n()
        ) |>
        rownames_to_column(var = "temperature_cut") |>
        mutate(temperature_cut = as.numeric(temperature_cut))

    print(temperature_cuts)

    pooled_cuts <- filter(temperature_cuts, n < 5) |> select(temperature_cut) |> t() |> as.vector()
    max_cut <- max(filter(temperature_cuts, n >= 5)$temperature_cut)
    temperature_cuts <- mutate(temperature_cuts, temperature_cut = if_else(
        temperature_cut %in% pooled_cuts,
        max_cut, temperature_cut)) |>
        select(temperature_cut, temperature_range)

    tempcut_plankton <-
        mutate(input_data, temperature_range = cut(watercolumn_temp, n)) |>
        left_join(temperature_cuts, by = "temperature_range")

    tempcut_plankton |> count(temperature_range, temperature_cut) |> print()
    
    tempcut_nest <- temperature_nest(temperature_cut, input_data = tempcut_plankton)
    tempcut_nest |> arrange(term) |> print(n = 40)

    estimate_plot <- ggplot(filter(tempcut_nest, term == "adjusted_log_chla"), aes(y = temperature_cut, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0)

    intercept_estimate_plot1 <- ggplot(tempcut_nest, aes(x = temperature_cut, y = estimate)) +
        geom_point(size = 3) +
        geom_smooth(se = FALSE, span = 0.8) +
        geom_smooth(se = FALSE, method = lm, color = "red") +
        facet_wrap(~term, scales = "free")

    cut_estimate_lm <- gam(
        estimate ~ temperature_cut,
        data = filter(tempcut_nest, term == "adjusted_log_chla"))

    cut_estimate_gam <- gam(
        estimate ~ s(temperature_cut),
        data = filter(tempcut_nest, term == "adjusted_log_chla"))

    intercept_estimate_plot2 <- ggplot(
            filter(tempcut_nest, term == "adjusted_log_chla"),
            aes(x = temperature_cut, y = estimate)
        ) +
        geom_point(size = 3) +
        geom_line(aes(y = fitted(cut_estimate_lm))) +
        geom_line(aes(y = fitted(cut_estimate_gam)))

    print(estimate_plot)
    print(intercept_estimate_plot1)
    print(intercept_estimate_plot2)

    print(summary(cut_estimate_gam))

    #cut_estimate_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "adjusted_log_chla"))
    print(summary(cut_estimate_lm))

    cut_intercept_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "(Intercept)"))
    print(summary(cut_intercept_lm))

}