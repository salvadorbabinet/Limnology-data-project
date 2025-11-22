# Data project scripts for Limnology.
# This script runs simple models for hypothesis testing.

library(mgcv)

str(plankton)
plankton <- mutate(plankton, lake_ID = factor(lake_ID))


# Hypothesis 1: Phyto- and zooplankton biomass respond to temperature.


# Compare the available chlorophyll a measures in Lake Pulse ----

p1 <- ggplot(plankton, aes(x = filtered_chla_surface)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = invivo_chla_surface)) + geom_histogram()
p3 <- ggplot(plankton, aes(x = invivo_chla_bottom)) + geom_histogram()

p1 + p2 + p3

plankton <-
    mutate(plankton, across(
        invivo_chla_surface:filtered_chla_surface,
        log, .names = "log_{.col}")
    ) |>
    select(!invivo_chla_surface:filtered_chla_surface)

p1 <- ggplot(plankton, aes(x = log_filtered_chla_surface)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = log_invivo_chla_surface)) + geom_histogram()
p3 <- ggplot(plankton, aes(x = log_invivo_chla_bottom)) + geom_histogram()

p1 + p2 + p3 # Log(chla) distributions

# Chla measures to ambient environment

p1 <- ggplot(plankton, aes(x = max_depth)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = surface_temp)) + geom_histogram()
p3 <- ggplot(plankton, aes(x = watercolumn_temp)) + geom_histogram()

p1 + p2 + p3

scatterplot_chla(max_depth)

# Shallow lakes tend to have high chla at surface & bottom samples
plankton |> select(
    lake_ID, max_depth,
    log_invivo_chla_surface, log_invivo_chla_bottom
    ) |>
    arrange(desc(log_invivo_chla_bottom)) |>
    print(n = 50)

p1 <- scatterplot(max_depth, surface_temp)
p2 <- scatterplot(max_depth, watercolumn_temp)

p1 + p2 # Water column temp tracks depth, surface temp much less so

scatterplot_chla(surface_temp)
scatterplot_chla(watercolumn_temp)

# Nutrient load distributions
p1 <- ggplot(plankton, aes(x = tp_tube)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = tn_tube)) + geom_histogram()
p3 <- ggplot(plankton, aes(x = tp_bottom)) + geom_histogram()
p4 <- ggplot(plankton, aes(x = tn_bottom)) + geom_histogram()
p1 + p2 + p3 + p4 + plot_layout(ncol = 2) # 130 missing values

plankton <-
    mutate(plankton, across(
        tp_tube:tn_bottom,
        log, .names = "log_{.col}")
    ) |>
    select(!tp_tube:tn_bottom)

p1 <- ggplot(plankton, aes(x = log_tp_tube)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = log_tn_tube)) + geom_histogram()
p3 <- ggplot(plankton, aes(x = log_tp_bottom)) + geom_histogram(bins = 20)
p4 <- ggplot(plankton, aes(x = log_tn_bottom)) + geom_histogram(bins = 20)
p1 + p2 + p3 + p4 + plot_layout(ncol = 2) # Log nutrient distributions

# Chla measures to nutrient load
scatterplot_chla(log_tp_tube)
scatterplot_chla(log_tn_tube)


# Plot zooplankton density to environment ----

plankton <- mutate(plankton, log_zoo_density = log(zoo_density))

p1 <- ggplot(plankton, aes(x = zoo_density)) + geom_histogram()
p2 <- ggplot(plankton, aes(x = log_zoo_density)) + geom_histogram()
p1 + p2

# Phyto- & zooplankton densities to surface water & water column temperatures.
p1 <- scatterplot(surface_temp)
p2 <- scatterplot(watercolumn_temp) + theme(axis.title.y = element_blank())
p3 <- scatterplot(surface_temp, log_invivo_chla_surface) + theme(axis.title.x = element_blank())
p4 <- scatterplot(watercolumn_temp, log_invivo_chla_surface) + theme(axis.title = element_blank())

p3 + p4 + p1 + p2 + plot_layout(ncol = 2)

p5 <- scatterplot(surface_temp, log_filtered_chla_surface) + theme(axis.title.x = element_blank())
p6 <- scatterplot(watercolumn_temp, log_filtered_chla_surface) + theme(axis.title = element_blank())

# Plot #1 in analysis doc:
p3 + p4 + p5 + p6 + p1 + p2 + plot_layout(ncol = 2)

# Stratified temperatures?
p1 <- scatterplot(epilimnion_temp, log_invivo_chla_surface) + theme(axis.title.x = element_blank())
p2 <- scatterplot(hypolimnion_temp, log_invivo_chla_surface) + theme(axis.title = element_blank())
p3 <- scatterplot(epilimnion_temp, log_filtered_chla_surface) + theme(axis.title.x = element_blank())
p4 <- scatterplot(hypolimnion_temp, log_filtered_chla_surface) + theme(axis.title = element_blank())
p5 <- scatterplot(epilimnion_temp)
p6 <- scatterplot(hypolimnion_temp)

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

# Phyto- to zooplankton biomass
p1 <- scatterplot(log_filtered_chla_surface)
p2 <- scatterplot(log_invivo_chla_surface) + theme(axis.title.y = element_blank())
p3 <- scatterplot(log_invivo_chla_bottom) + theme(axis.title.y = element_blank())

p1 + p2 + p3


# Linear relationships ----

# Filtered chla (color by water column temp)

ggplot(plankton, aes(x = log_filtered_chla_surface, y = log_zoo_density)) +
    geom_point(mapping = aes(fill = watercolumn_temp), alpha = 0.8, size = 2.5, shape = 21) +
    geom_smooth(
        color = "black", linetype = 1,
        method = "lm", se = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 17)

filtered_lm <- lm(log_zoo_density~log_filtered_chla_surface, data = plankton)

par(mfrow = c(2,2))
plot(filtered_lm)
shapiro.test(residuals(filtered_lm)) # Not normally distributed...

summary(filtered_lm)

# In vivo chla (color by surface temp)

# ...tube (surface) sample measure
ggplot(plankton, aes(x = log_invivo_chla_surface, y = log_zoo_density)) +
    geom_point(mapping = aes(fill = watercolumn_temp), alpha = 0.8, size = 2.5, shape = 21) +
    geom_smooth(
        color = "black", linetype = 1,
        method = "lm", se = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 17)

invivo_lm <- lm(log_zoo_density~log_invivo_chla_surface, data = plankton)

plot(invivo_lm)
shapiro.test(residuals(invivo_lm)) # Not normally distributed...

summary(invivo_lm)

# ...bottom sample measure
ggplot(plankton, aes(x = log_invivo_chla_bottom, y = log_zoo_density)) +
    geom_point(mapping = aes(fill = watercolumn_temp), alpha = 0.8, size = 2.5, shape = 21) +
    geom_smooth(
        color = "black", linetype = 1,
        method = "lm", se = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 17)

invivo_bottom_lm <- lm(log_zoo_density~log_invivo_chla_bottom, data = plankton)

plot(invivo_bottom_lm)
shapiro.test(residuals(invivo_bottom_lm)) # Most, but significantly not normal

summary(invivo_bottom_lm)


# Simple regression suggests significant relationships between phytoplankton
# and zooplankton densities. But the R2 is low, and the residuals are not
# normally distributed (but maybe close enough?)

# A GLM may be better suited, and/or incorporating nutrient load & depth.
# GLM would address the normal distribution issue, but not the possibility
# that effects are non-linear... in which case GAM is best.


# This relationship is generally linear, except for some temperature subsets.

gam_data <- filter(plankton, !is.na(log_zoo_density))

gam1 <- gam(log_zoo_density ~ s(log_filtered_chla_surface), data = gam_data)

summary(gam1)

scatterplot(log_filtered_chla_surface, data = gam_data) +
    geom_line(mapping = aes(y = fitted(gam1)))


# Nested linear relation between densities


# Surface temperature ----

n <- 12

surface_temp_bins <- plankton |>
    mutate(surface_temp_bin = ntile(surface_temp, n)) |>
    group_by(surface_temp_bin) |>
    summarize(
        min = min(surface_temp), max = max(surface_temp),
        interval = max - min, n = n()
    )

surface_temp_bins

ncut_plankton <- mutate(plankton, surface_temp_bin = ntile(surface_temp, n))
ncut_plankton |> count(surface_temp_bin)

surface_nest <- temperature_nest(input_data = ncut_plankton)
nest_plot()

surface_temp_plankton |> select(surface_temp_bin, term, estimate) |>
    pivot_wider(names_from = term, values_from = estimate) |>
    ggplot(mapping = aes(x = `(Intercept)`)) +
    geom_histogram(bins = 8)

surface_temp_plankton |> select(surface_temp_bin, term, estimate) |>
    pivot_wider(names_from = term, values_from = estimate) |>
    ggplot(mapping = aes(x = adjusted_log_chla)) +
    geom_histogram(bins = 8)


# Water column temperature

n <- 15 # No. of bins... p value increases too much past n = 20

# Take ntile ----

temperature_cuts <- plankton |>
    mutate(temp_ntile = ntile(watercolumn_temp, n)) |>
    group_by(temp_ntile) |>
    summarize(
        min = min(watercolumn_temp), max = max(watercolumn_temp),
        ntile_midtemp = (min + max) / 2, interval = max - min, n = n()
    )

temperature_cuts
temperature_cuts <- select(temperature_cuts, temp_ntile, ntile_midtemp)

ncut_plankton <- plankton |>
    select(!epilimnion_temp:tntp_ratio) |>
    mutate(temp_ntile = ntile(watercolumn_temp, n), .after = watercolumn_temp) |>
    left_join(temperature_cuts) |>
    relocate(ntile_midtemp, .after = temp_ntile)

ncut_plankton

tempntile_nest <- temperature_nest(ntile_midtemp, input_data = ncut_plankton)
tempntile_nest |> arrange(term) |> print(n = 22)

lm1 <- lm(estimate ~ ntile_midtemp, data = filter(tempntile_nest, term == "adjusted_log_chla"))
summary(lm1)

lm2 <- lm(estimate ~ ntile_midtemp, data = filter(tempntile_nest, term == "(Intercept)"))
summary(lm2)

nest_plot(ntile_midtemp, tempntile_nest)

# Compare midpoint temp to ntile number

tempntile_nest <- temperature_nest(temp_ntile, input_data = ncut_plankton)

lm3 <- lm(estimate ~ temp_ntile, data = filter(tempntile_nest, term == "adjusted_log_chla"))
summary(lm3)

lm4 <- lm(estimate ~ temp_ntile, data = filter(tempntile_nest, term == "(Intercept)"))
summary(lm4)

nest_plot(watercolumn_temp_bin, watercolumn_nest)

ggplot(filter(watercolumn_nest, term == "adjusted_log_chla"), aes(x = watercolumn_temp_bin, y = estimate)) +
    geom_point(size = 3.5) +
    geom_smooth(se = FALSE, method = lm)

estimate_lm <- lm(estimate ~ watercolumn_temp_bin, data = filter(watercolumn_nest, term == "adjusted_log_chla"))
summary(estimate_lm)

estimate_gam <- gam(
    estimate ~ s(watercolumn_temp_bin),
    data = filter(watercolumn_nest, term == "adjusted_log_chla"))
summary(estimate_gam)

ggplot(filter(watercolumn_nest, term == "adjusted_log_chla"), aes(x = watercolumn_temp_bin, y = estimate)) +
    geom_point(size = 3.5) +
    geom_smooth(method = lm, se = FALSE, color = "red") +
    geom_line(aes(y = fitted(estimate_gam)))

temperature_nest(watercolumn_temp_bin) |>
    filter(watercolumn_temp_bin == 6, term == "adjusted_log_chla") |>
    unnest(data_c) |>
    ggplot(mapping = aes(x = log_filtered_chla_surface, y = log_zoo_density)) +
    geom_point(aes(color = watercolumn_temp)) +
    geom_smooth(method = lm)


# Take cut ----

cut_nest_plot1(7) # 11, 15 cuts are good
cut_nest_plot2(7)

tempcut_nest <- temperature_nest(temperature_cut, input_data = tempcut_plankton)
tempcut_nest |> arrange(term) #|> select(temperature_cut, term, estimate, p.value)
tempcut_nest |> group_by(term) |> summarize(med_p = median(p.value))

ggplot(filter(tempcut_nest, term == "adjusted_log_chla", temperature_cut != 9), aes(y = temperature_cut, x = estimate)) +
        geom_point(size = 3.5) +
        geom_errorbar(mapping = aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
        geom_vline(xintercept = 0)

ggplot(filter(tempcut_nest, temperature_cut != 9), aes(x = temperature_cut, y = estimate)) +
        geom_point(size = 3) +
        geom_smooth(se = FALSE, span = 1) +
        geom_smooth(se = FALSE, method = lm, color = "red") +
        facet_wrap(~term, scales = "free")

cutestimate_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "(Intercept)", temperature_cut != 9))
summary(cutestimate_lm)

cutintercept_lm <- lm(estimate ~ temperature_cut, data = filter(tempcut_nest, term == "adjusted_log_chla", temperature_cut != 9))
summary(cutestimate_lm)

tempcut_nest |>
    filter(temperature_cut == 9, term == "adjusted_log_chla") |>
    unnest(data_c) |>
    ggplot(mapping = aes(x = log_filtered_chla_surface, y = log_zoo_density)) +
    geom_point(size = 3, aes(color = watercolumn_temp)) +
    geom_smooth(method = lm)


# Multiple regression ----

lm1 <- lm(log_filtered_chla_surface ~ watercolumn_temp + log_tn_tube, data = plankton)
summary(lm1)

lm2 <- lm(log_zoo_density ~ log_filtered_chla_surface, data = ncut_plankton)
lm2_resid <- rstandard(lm2)

summary(lm2)
plot(lm2_resid ~ ncut_plankton$temp_ntile)

ggplot(plankton, aes(x = log_tn_tube, y = log_filtered_chla_surface)) +
    geom_point()

ggplot(plankton, aes(x = log_tp_tube, y = log_filtered_chla_surface)) +
    geom_point()

lm3 <- lm(log_zoo_density ~ log_invivo_chla_surface * watercolumn_temp, data = plankton)
summary(lm3)

ggplot(plankton, aes(x = log_zoo_density, y = log_filtered_chla_surface)) +
    geom_point()

lm4 <- lm(log_zoo_density ~ log_filtered_chla_surface * watercolumn_temp + max_depth + log_tn_tube, data = plankton)
summary(lm4)
