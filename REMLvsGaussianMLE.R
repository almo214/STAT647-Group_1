library(nlme)


df.bond <- data.frame(ingot = rep(1:7, each = 3),
                      metal = rep(c("n", "i", "c"), 7),
                      pres = c(67, 71.9, 72.2,
                               67.5, 68.8, 66.4,
                               76, 82.6, 74.5,
                               72.7, 78.1, 67.3,
                               73.1, 74.2, 73.2,
                               65.8, 70.8, 68.7,
                               75.6, 84.9, 69))
df.bond$ingot <- factor(df.bond$ingot)

# Linear mixed model by MLE
lmm.bond_mle <- lme(pres ~ metal, random = ~1 | ingot, data = df.bond, method = "ML")

# Linear mixed model by REML
lmm.bond_reml <- lme(pres ~ metal, random = ~1 | ingot, data = df.bond, method = "REML")

# Compare AIC/BIC
AIC_mle <- AIC(lmm.bond_mle)
BIC_mle <- BIC(lmm.bond_mle)
AIC_reml <- AIC(lmm.bond_reml)
BIC_reml <- BIC(lmm.bond_reml)

# Print AIC and BIC values
cat("MLE AIC:", AIC_mle, "\n")
cat("MLE BIC:", BIC_mle, "\n")
cat("REML AIC:", AIC_reml, "\n")
cat("REML BIC:", BIC_reml, "\n")

# Summary of MLE model
summary(lmm.bond_mle)

# Summary of REML model
summary(lmm.bond_reml)
