# Data project scripts for limnology.
# Linear mixed models.

library(lme4)
plankton # From hierarchy script

lm1 <- lm(log_zoo_density ~ log_chla, data = plankton)
summary(lm1)

lmm1 <- lmer(log_zoo_density ~ log_chla + (watercolumn_temp | size_class), data = plankton)
summary(lmm1)
