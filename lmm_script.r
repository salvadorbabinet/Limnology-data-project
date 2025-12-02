# Data project scripts for limnology.
# Linear mixed models and GAMs.

library(mgcv)
library(lme4)


plankton # From hierarchy script

ggplot(plankton, aes(x = surface_temp)) + geom_histogram()
ggplot(plankton, aes(x = watercolumn_temp)) + geom_histogram()
ggplot(plankton, aes(x = log(watercolumn_temp))) + geom_histogram()

modeling_data <- plankton |>
    mutate(
        centered_temp = scale(watercolumn_temp, scale = FALSE),
        scaled_temp = scale(watercolumn_temp))

ggplot(modeling_data, aes(x = centered_temp)) + geom_histogram()
ggplot(modeling_data, aes(x = scaled_temp)) + geom_histogram()


# GAMs ----

par(mfrow = c(1,1))

# Phytoplankton ~ ...

gam1 <- gam(log_chla ~ s(watercolumn_temp, k = 5), method = "REML", data = plankton)
summary(gam1)

gam.check(gam1)
plot(gam1, main = "Phytoplankton",
    xlab = "Temperature (°C)", ylab = "Smoothed effect estimate")

plankton |>
    select(
        log_chla, watercolumn_temp,
        max_depth, log_phosphorus,
        fraction_agriculture
        ) |>
    as.data.frame() |>
    cor()

# Choose N or P, probably P

gam2 <- gam(log_chla ~ s(log_phosphorus, k = 20), data = plankton)
summary(gam2)
gam.check(gam2)
plot(gam2)

lm2 <- lm(log_chla ~ log_phosphorus, data = plankton)
summary(lm2)

ggplot(plankton, aes(x = log_phosphorus, y = log_chla)) +
    geom_point() +
    geom_line(aes(y = fitted(lm2))) +
    geom_line(aes(y = fitted(gam2)), color = "red")

gam3 <- gam(log_chla ~ s(watercolumn_temp, k = 7) + s(log_phosphorus, k = 9), data = plankton)
summary(gam3)
gam.check(gam3)
plot(gam3)

# I think keeping the phosphorus response linear is OK, though
# smoothed phosphorus effect increases R2 and decreases AIC somewhat
gam4 <- gam(log_chla ~ s(watercolumn_temp, k = 7) + log_phosphorus, method = "REML", data = plankton)
summary(gam4)
gam.check(gam4)
plot(gam4, all.terms = TRUE)

gam.resid <- residuals(gam4)
write_csv(as.data.frame(gam.resid), "PhytoplanktonGAMResiduals.csv")

AIC(gam3, gam4)
AIC(gam1, lm2, gam2, gam3, gam4)

# Zooplankton ~ ...

zgam1 <- gam(log_zoo_density ~ s(watercolumn_temp, k = 7), method = "REML", data = plankton)
summary(zgam1)

gam.check(zgam1)
plot(zgam1, main = "Zooplankton",
    xlab = "Temperature (°C)", ylab = "Smoothed effect estimate")

ggplot(filter(plankton, !is.na(log_zoo_density)), aes(x = watercolumn_temp, y = log_zoo_density)) +
    geom_point() +
    geom_line(aes(y = fitted(zgam1)), color = "red")

plankton |>
    select(
        log_zoo_density, watercolumn_temp,
        max_depth, watershed_area,
        log_phosphorus, log_nitrogen,
        fraction_agriculture
        ) |>
    filter(!is.na(log_zoo_density)) |>
    as.data.frame() |>
    cor()

zgam2 <- gam(log_zoo_density ~ s(log_phosphorus), method = "REML", data = plankton)
summary(zgam2)

zlm2 <- lm(log_zoo_density ~ log_phosphorus, data = plankton)
summary(zlm2)

zgam3 <- gam(
    log_zoo_density ~ s(watercolumn_temp, k = 8) + log_phosphorus,
    method = "REML", data = plankton)
summary(zgam3)

gam.check(zgam3)
plot(zgam3, all.terms = TRUE)

zgam4 <- gam(
    log_zoo_density ~ s(watercolumn_temp, k = 8) + s(log_phosphorus),
    method = "REML", data = plankton)
summary(zgam4)

gam.check(zgam4)
plot(zgam4)

AIC(zgam1, zlm2, zgam2, zgam3, zgam4)


# Results table

summary(gam4)
summary(zgam3)

glance(gam4)
glance(zgam3)

tidy(gam4)
tidy(zgam3)


# LMMs ----

lm1 <- lm(log_zoo_density ~ log_chla, data = plankton)
summary(lm1)

lmm1 <- lmer(log_zoo_density ~ log_chla, data = plankton)
summary(lmm1)
