# R script used to conduct the statistics presented in the manuscript

#load packages
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(broom)

# Clear ggplot2# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

# load data set 
data <- read.csv("Gestation_2018.csv")

#### calculate neonate Body condition index (BCI) ####  

# Fit the linear regression model
model <- lm(Mass ~ SVL, data = data)

# Extract the residuals
residuals <- resid(model)

# Assign residuals to the BCI column
data$BCI <- residuals

# Extract the fitted values and residuals
data$fitted <- fitted(model)
data$residuals <- resid(model)

# plot residuals
ggplot(data, aes(x = SVL, y = Mass, shape = Treatment, color = residuals)) +
  geom_point(size = 3) +  # Scatter plot with residuals as color and Treatment as shape
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Regression line
  geom_segment(aes(xend = SVL, yend = fitted), linetype = "dashed", color = "red") +  # Residuals
  labs(x = "Snout-Vent Length (SVL)",
    y = "Mass",
    color = "Residuals (BCI)",
    shape = "Treatment") +
  theme_minimal()


#### Add the MomRelSize column to the data ####
# the relative size is calculated based on ranked ID order, converted so that larger individuals have higher number
data <- data %>%
  mutate(MomRelSize = 100 - as.numeric(substr(as.character(Mother), 
                                              nchar(as.character(Mother)) - 1, 
                                              nchar(as.character(Mother)))))

#write.csv(x = data, file = "Full_data.csv",row.names = FALSE)

#### Run Linear Mixed-effects Model comparisons ####

# For each LMM 3 models were tested using the Akaike information criterion (AIC).
# The full model includes Treatment and relative maternal body size as fixed effects
# and mother identity (specific litter) as a random effect

#MASS

Mass_model_full <- lmer(Mass ~ Treatment + MomRelSize + (1 | Mother), data = data)
summary(Mass_model_full) 

Mass_model_mid <- lmer(Mass ~ Treatment + (1 | Mother), data = data)
summary(Mass_model_mid)

Mass_model_reduced <- lmer(Mass ~ (1 | Mother), data = data)
summary(Mass_model_reduced)

# Compare the 3 models
AIC_Mass <- anova(Mass_model_reduced, Mass_model_mid, Mass_model_full)


#SVL
SVL_model_full <- lmer(SVL ~ Treatment + MomRelSize + (1 | Mother), data = data)
summary(SVL_model_full)

SVL_model_mid <- lmer(SVL ~ Treatment + (1 | Mother), data = data)
summary(SVL_model_mid)

SVL_model_reduced <- lmer(SVL ~ (1 | Mother), data = data)
summary(SVL_model_reduced)

# Compare the 3 models
AIC_SVL <- anova(SVL_model_reduced, SVL_model_mid, SVL_model_full)

#BCI
BCI_model_full <- lmer(BCI ~ Treatment + MomRelSize + (1 | Mother), data = data)
summary(BCI_model_full)

BCI_model_mid <- lmer(BCI ~ Treatment + (1 | Mother), data = data)
summary(BCI_model_mid)

BCI_model_reduced <- lmer(BCI ~ (1 | Mother), data = data)
summary(BCI_model_reduced)

# Compare the 3 models
AIC_BCI <- anova(BCI_model_reduced, BCI_model_mid, BCI_model_full)

#Gestation
Gestation_model_full <- lmer(Gestation ~ Treatment + MomRelSize + (1 | Mother), data = data)
summary(Gestation_model_full)

Gestation_model_mid <- lmer(Gestation ~ Treatment + (1 | Mother), data = data)
summary(Gestation_model_mid)

Gestation_model_reduced <- lmer(Gestation ~ (1 | Mother), data = data)
summary(Gestation_model_reduced)

# Compare the 3 models
AIC_Gestation <- anova(Gestation_model_reduced, Gestation_model_mid, Gestation_model_full)

# summary of all AIc results
AIC <- bind_rows(AIC_Mass, AIC_SVL, AIC_BCI, AIC_Gestation)

#write.csv(x = AIC,file = "AIC_Results.csv", row.names = TRUE)


####### litter based effects ###########

### calculate total and mean traits by litter ### 

data_by_mother <- data %>%
  group_by(Mother) %>%
  summarize(
    TotalMass = sum(Mass),
    AverageMass = mean(Mass),
    Treatment = first(Treatment),
    Gestation = first(Gestation),
    Litter_size = first(Litter_size),
    AverageSVL = mean(SVL),
    AverageBCI = mean(BCI),
    MomRelSize = first(MomRelSize)
  )

# read in DF with stillborn and deformation data
stillborn <- read.csv("StillBorn.csv")

# add that data to the DF
data_by_mother <- left_join(data_by_mother, stillborn, by = "Mother")

#write.csv(data_by_mother, file = "data_by_litter.csv", row.names = FALSE)

#### Run Linear Model ####

# test the effect of Treatment on total mass
tot_mass_model <- lm(TotalMass ~ Treatment, data = data_by_mother)
summary(tot_mass_model) #p-value: 0.1278

# test the effect of Treatment on mean neonate mass per litter
ave_mass_model <- lm(AverageMass ~ Treatment, data = data_by_mother)
summary(ave_mass_model) #p-value: 0.03142*

# test the effect of Treatment on mean neonate SVL per litter
ave_SVL_model <- lm(AverageSVL ~ Treatment, data = data_by_mother)
summary(ave_SVL_model) #p-value: 0.621

# test the effect of Treatment on stillborns per litter
still_model <- lm(stillborn ~ Treatment, data = data_by_mother)
summary(still_model) #p-value: 0.1678

# test the effect of Treatment on deformations per litter
deform_model <- lm(deformed ~ Treatment, data = data_by_mother)
summary(deform_model) #p-value: 0.914

# test the effect of Treatment on litter size
litter_size_model <- lm(Litter_size ~ Treatment, data = data_by_mother)
summary(litter_size_model) # p-value: 0.566

# test the effect of Treatment on mean neonate BCI per litter
ave_BCI_model <- lm(AverageBCI ~ Treatment, data = data_by_mother)
summary(ave_BCI_model) # p-value: 0.00635 **



####### correlations between neonate traits by treatment ##########


# Relationship between neonate SVL and neonate mass by treatment
svl_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(Mass ~ SVL, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic)

print(svl_mass_R2)

# Relationship between neonate SVL and neonate BCI by treatment
svl_bci_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(BCI ~ SVL, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic)

print(svl_bci_R2)

# Relationship between neonate SVL and neonate BCI by treatment
mass_bci_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(BCI ~ Mass, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic)

print(mass_bci_R2)

# create summary table
trait_R2 <- bind_rows(svl_mass_R2, svl_bci_R2, mass_bci_R2)
print(trait_R2)


#### correlations between litter size and neonate traits ####

# Relationship between litter size and neonate mass by treatment
litsize_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(Mass ~ Litter_size, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic, nobs)

print(litsize_mass_R2)

# Relationship between litter size and neonate SVL by treatment
litsize_svl_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(SVL ~ Litter_size, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic, nobs)

print(litsize_svl_R2)

# Relationship between litter size and gestation length by treatment
litsize_gest_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(Gestation ~ Litter_size, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic, nobs)

print(litsize_gest_R2)

# Relationship between litter size and neonate BCI by treatment
litsize_bci_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(BCI ~ Litter_size, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic, nobs)

print(litsize_bci_R2)


#### correlations between Gestation length and neonate traits ####

# Relationship between gestation length and neonate mass by treatment
gest_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(Mass ~ Gestation, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic)

print(gest_mass_R2)

# Relationship between gestation length and neonate SVL by treatment
gest_svl_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(SVL ~ Gestation, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic)

print(gest_svl_R2)


# Relationship between gestation length and neonate BCI by treatment
gest_bci_R2 <- data %>%
  group_by(Treatment) %>%
  do(glance(lm(BCI ~ Gestation, data = .))) %>%
  select(Treatment, r.squared, p.value, df, statistic) #should redo this order...

print(gest_bci_R2)



##### correlations between ranked maternal size and traits (Kendall's tau) ####

# Relationship between relative maternal body size and neonate mass by treatment
Mass_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$Mass, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}


# Apply the function to each treatment group
tau_Mass <- data %>%
  group_by(Treatment) %>%
  do(Mass_tau(.))

print(tau_Mass)

# Relationship between relative maternal body size and neonate BCI by treatment
bci_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$BCI, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_bci <- data %>%
  group_by(Treatment) %>%
  do(bci_tau(.))

print(tau_bci)

# Relationship between relative maternal body size and neonate SVL by treatment
svl_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$SVL, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_svl <- data %>%
  group_by(Treatment) %>%
  do(svl_tau(.))

print(tau_svl)

# Relationship between relative maternal body size and gestation length by treatment
gest_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$Gestation, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_gest <- data %>%
  group_by(Treatment) %>%
  do(gest_tau(.))

print(tau_gest)

# Relationship between relative maternal body size and litter size by treatment
lit_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$Litter_size, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_lit <- data %>%
  group_by(Treatment) %>%
  do(lit_tau(.))

print(tau_lit)

# Relationship between relative maternal body size and total litter mass by treatment
litmass_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$TotalMass, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_litmass <- data_by_mother %>%
  group_by(Treatment) %>%
  do(litmass_tau(.))

print(tau_litmass)

