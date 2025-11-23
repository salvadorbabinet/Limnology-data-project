# Data project scripts for limnology.
# This script loads helpers for hierarchical temperature models.

prep_plankton <- function(lakepulse, zooplankton) {

    plankton <- lakepulse |> select(
        lake_ID, lake_name:longitude, ecozone:size_class,
        watershed_area, max_depth,
        rbr_chla_mean_tube, #rbr_chla_mean_bottom,
        rbr_temperature_mean_tube, rbr_temperature_mean_watercolumn,
        chla, fraction_agriculture, tp_tube:tntp_ratio
    )

    lake_id <- zooplankton |>
        rownames_to_column("order") |>
        select(order, lake_ID)

    lake_biomass <-
        rowSums(zooplankton[-1]) |>
        as_tibble() |>
        rownames_to_column("order")

    zoo_density <-
        left_join(lake_id, lake_biomass)[-1] |>
        rename(zoo_density = value)

    left_join(plankton, zoo_density) |>
        rename(
            invivo_chla = rbr_chla_mean_tube,
            #invivo_chla_bottom = rbr_chla_mean_bottom,
            #filtered_chla_surface = chla,
            surface_temp = rbr_temperature_mean_tube,
            watercolumn_temp = rbr_temperature_mean_watercolumn,
        ) |>
        mutate(
            log_chla = log(chla),
            log_invivo_chla = log(invivo_chla),
            log_zoo_density = log(zoo_density),
            log_nitrogen = log(tn_tube),
            log_phosphorus = log(tp_tube),
            .after = chla
        ) |>
        relocate(chla, invivo_chla, zoo_density, .after = fraction_agriculture)

}


# Cut temperature by ntile, so same no. of observations
# but different intervals between temperature cuts.
ncut_temperature <- function(n, input_data) {

    temperature_cuts <- input_data |>
        mutate(temperature_ntile = ntile(watercolumn_temp, n)) |>
        group_by(temperature_ntile) |>
        summarize(
            min = min(watercolumn_temp), max = max(watercolumn_temp),
            cut_midtemp = (max + min) / 2, interval = max - min, n = n()
        )

    print(temperature_cuts)
    temperature_cuts <- select(temperature_cuts, temperature_ntile, cut_midtemp)

    mutate(input_data,
            temperature_ntile = ntile(watercolumn_temp, n),
            .after = watercolumn_temp) |>
        left_join(temperature_cuts) |>
        relocate(cut_midtemp, .after = temperature_ntile)

}


# Cut temperature by interval, so same intervals but
# different no. of observations between temperature cuts.
cut_temperature <- function(n, input_data, pool_cuts = TRUE) {

    temperature_cuts <- input_data |>
        mutate(temperature_interval = cut(watercolumn_temp, n)) |>
        group_by(temperature_interval) |>
        summarize(
            min = min(watercolumn_temp), max = max(watercolumn_temp),
            cut_midtemp = (max + min) / 2, interval = max - min, n = n()
        ) |>
        rownames_to_column(var = "temperature_cut") |>
        mutate(temperature_cut = as.numeric(temperature_cut))

    print(temperature_cuts)

    # For high temperature intervals with few observations, pool the
    # upper range into a wider combined interval.
    pooled_cutlevels <-
        filter(temperature_cuts, n < 5) |>
        select(temperature_cut) |>
        t() |> as.vector()
    max_cutlevel <- max(filter(temperature_cuts, n >= 5)$temperature_cut)
    lower_bound <- max(filter(temperature_cuts, n >= 5)$min)
    upper_bound <- max(temperature_cuts$max)

    pooled_cuts <- mutate(temperature_cuts,
        cut_midtemp = if_else(
            temperature_cut %in% pooled_cutlevels,
            (lower_bound + upper_bound) / 2, cut_midtemp),
        temperature_cut = if_else(
            temperature_cut %in% pooled_cutlevels,
            max_cutlevel, temperature_cut),
        cut_midtemp = if_else(
            temperature_cut == max(temperature_cut),
            max(cut_midtemp), cut_midtemp)
        ) |>
        select(temperature_interval, temperature_cut, cut_midtemp)

    # Or keep intervals constant & drop ranges with few observations.
    # This means losing the entire upper temperature extreme.
    dropped_cuts <-
        filter(temperature_cuts, n > 6) |>
        select(temperature_interval, temperature_cut, cut_midtemp)

    if (pool_cuts == TRUE) {

        print(pooled_cuts)

        mutate(input_data,
            temperature_interval = cut(watercolumn_temp, n)) |>
            left_join(pooled_cuts, by = "temperature_interval") |>
            relocate(temperature_cut, cut_midtemp, .after = watercolumn_temp)

    } else if (pool_cuts == FALSE) {

        print(dropped_cuts)

        mutate(input_data,
            temperature_interval = cut(watercolumn_temp, n)) |>
            left_join(dropped_cuts, by = "temperature_interval") |>
            relocate(temperature_cut, cut_midtemp, .after = watercolumn_temp)

    }

}


# Nest data by temperature cut (ntile or fixed interval) and
# fit simple regressions within each bin.
nest_temperature <- function(nesting_var, input_data) {

    input_data |>
        nest_by({{nesting_var}}) |>
        mutate(
            mean_chla = mean(data$log_chla),
            data_c = list(mutate(data, adjusted_log_chla = log_chla - mean_chla)),
            models = list(lm(log_zoo_density~adjusted_log_chla, data = data_c)),
            tidymodels = list(broom::tidy(models))
        ) |>
        select(!c(data, models)) |>
        unnest(tidymodels) |>
        mutate(term = if_else(term == "(Intercept)", "Intercept", "Effect"))

}


# Plot nested estimates.
plot_nest <- function(nesting_var, data_nest, plot_fits = FALSE) {

    effect_plot <- 
        ggplot(
            filter(data_nest, term == "Effect"),
            aes(y = {{ nesting_var }}, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0) +
        labs(
            x = "Effect estimate", y = "Temperature (°C)"
        )

    print(effect_plot)

    intercept_effect_plot1 <-
        ggplot(data_nest, aes(x = {{ nesting_var }}, y = estimate)) +
        geom_point(aes(fill = {{ nesting_var }}), size = 3.5, shape = 21) +
        geom_smooth(method = lm, color = "black", linewidth = 1, alpha = 0.2) +
        scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = median(data_nest$cut_midtemp)) +
        facet_wrap(~term, scales = "free") +
        labs(
            x = "Temperature (°C)", fill = "Temperature (°C)",
            y = "Estimate"
        )

    print(intercept_effect_plot1)

    if (plot_fits == TRUE) {

    estimate_plot2 <-
        ggplot(
            filter(data_nest, term == "Effect"),
            aes(x = {{ nesting_var }}, y = estimate)
        ) +
        geom_point(aes(fill = {{ nesting_var }}), shape = 21, size = 4) +
        geom_line(aes(y = fitted(lm14)), linetype = 2) +
        geom_line(aes(y = fitted(gam14)), color = "red") +
        scale_fill_gradient2(
            low = "blue", mid = "orange", high = "red",
            midpoint = median(data_nest$cut_midtemp)
        ) +
        labs(
            x = "Temperature (°C)", fill = "Temperature (°C)",
            y = "Effect estimate"
        )

    estimate_plot1 <-
        ggplot(
            filter(data_nest, term == "Intercept"),
            aes(x = {{ nesting_var }}, y = estimate)
        ) +
        geom_point(aes(fill = {{ nesting_var }}), shape = 21, size = 4) +
        geom_line(aes(y = fitted(lm14_intercepts)), linetype = 2) +
        geom_line(aes(y = fitted(gam14_intercepts)), color = "red") +
        scale_fill_gradient2(
            low = "blue", mid = "orange", high = "red",
            midpoint = median(data_nest$cut_midtemp)
        ) +
        labs(
            x = "Temperature (°C)", fill = "Temperature (°C)",
            y = "Intercept estimate"
        )

    print(estimate_plot1)
    print(estimate_plot2)

    estimate_plot1 <-
        estimate_plot1 +
        theme(legend.position = "none")

    print(estimate_plot1 + estimate_plot2 + plot_layout(axes = "collect_x"))

    }

}
