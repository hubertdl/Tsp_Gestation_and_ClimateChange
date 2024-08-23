# R script used to conduct the statistics and generate figures presented in the manuscript

# Set working directory to a location that contains the following two files:
# Gestation_2018.csv and StillBorn.csv
setwd("")

#load packages
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(broom)
library(MuMIn)
library(broom.mixed)

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

# load data set 
data <- read.csv("Gestation_2018.csv")

# Statistics #### 

### calculate Scaled-mass index (SMI) for each treatment[Sup Table 1] #####

# Separate data by treatment
data_treatment <- subset(data, Treatment == "treatment")
data_control <- subset(data, Treatment == "control")

# Calculate SMI for Treatment
scaling_model1 <- lmer(log(Mass) ~ log(SVL) + (1 | Mother), data = data_treatment)
b1 <- fixef(scaling_model1)["log(SVL)"]
svl_0_1 <- mean(data_treatment$SVL)
data_treatment$SMI <- data_treatment$Mass * (svl_0_1 / data_treatment$SVL) ^ b1

# Calculate SMI for control
scaling_model2 <- lmer(log(Mass) ~ log(SVL) + (1 | Mother), data = data_control)
b2 <- fixef(scaling_model2)["log(SVL)"]
svl_0_2 <- mean(data_control$SVL)
data_control$SMI <- data_control$Mass * (svl_0_2 / data_control$SVL) ^ b2

#add SMI column to DF
data$SMI <- NA 

# Fill in SMI values for Treatment
data$SMI[data$Treatment == "treatment"] <- data_treatment$SMI
# Fill in SMI values for control
data$SMI[data$Treatment == "control"] <- data_control$SMI

### Add Relative maternal size column to the data [Sup Table 1] ####
# the relative size is calculated based on ranked ID order, converted so that larger individuals have higher number
data <- data %>%
  mutate(MomRelSize = 100 - as.numeric(substr(as.character(Mother), 
                                              nchar(as.character(Mother)) - 1, 
                                              nchar(as.character(Mother)))))

#write.csv(x = data, file = "Full_data.csv",row.names = FALSE)

### Linear Mixed-effects Model comparisons (AIC) [Sup Table 2] ####

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

#SMI
SMI_model_full <- lmer(SMI ~ Treatment + MomRelSize + (1 | Mother), data = data)
summary(SMI_model_full)

SMI_model_mid <- lmer(SMI ~ Treatment + (1 | Mother), data = data)
summary(SMI_model_mid)

SMI_model_reduced <- lmer(SMI ~ (1 | Mother), data = data)
summary(SMI_model_reduced)

# Compare the 3 models
AIC_SMI <- anova(SMI_model_reduced, SMI_model_mid, SMI_model_full)

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
AIC <- bind_rows(AIC_Mass, AIC_SVL, AIC_SMI, AIC_Gestation)
AIC

#write.csv(x = AIC,file = "AIC_Results.csv", row.names = TRUE)


### Calculate Litter based effects [Sup Table 3] ####

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
    AverageSMI = mean(SMI),
    MomRelSize = first(MomRelSize)
  )

# read in DF with stillborn and deformation data
stillborn <- read.csv("StillBorn.csv")

# add that data to the DF
data_by_mother <- left_join(data_by_mother, stillborn, by = "Mother")

#write.csv(data_by_mother, file = "data_by_litter.csv", row.names = FALSE)

# Run Linear Models #

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

# test the effect of Treatment on mean neonate SMI per litter
ave_SMI_model <- lm(AverageSMI ~ Treatment, data = data_by_mother)
summary(ave_SMI_model) # p-value: 0.0004693 ***



### Relationship between neonate traits by treatment (LMM)[Sup Table 4] ####

# Relationship between neonate SVL and neonate mass by treatment
svl_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(Mass ~ SVL + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["SVL"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

svl_mass_R2


# Relationship between neonate SVL and neonate SMI by treatment

svl_SMI_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SMI ~ SVL + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["SVL"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

# View the results
print(svl_SMI_R2)


# Relationship between neonate Mass and neonate SMI by treatment
mass_SMI_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SMI ~ Mass + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Mass"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)


mass_SMI_R2


# create summary table
trait_R2 <- bind_rows(svl_mass_R2, svl_SMI_R2, mass_SMI_R2)
print(trait_R2)

### Relationship between litter size and neonate traits (LMM)[Sup Table 4] ####

# Relationship between litter size and neonate mass by treatment
lit_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(Mass ~ Litter_size + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Litter_size"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

lit_mass_R2

# Relationship between litter size and neonate SVL by treatment
lit_SVL_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SVL ~ Litter_size + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Litter_size"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

lit_SVL_R2

# Relationship between litter size and gestation length by treatment
lit_gest_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(Gestation ~ Litter_size + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Litter_size"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

lit_gest_R2

# Relationship between litter size and neonate SMI by treatment
lit_SMI_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SMI ~ Litter_size + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Litter_size"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

lit_SMI_R2


### Relationship between Gestation length and neonate traits (LMM) [Sup Table 4] ####

# Relationship between gestation length and neonate mass by treatment
gest_mass_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(Mass ~ Gestation + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Gestation"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

gest_mass_R2

# Relationship between gestation length and neonate SVL by treatment
gest_SVL_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SVL ~ Gestation + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Gestation"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

gest_SVL_R2


# Relationship between gestation length and neonate SMI by treatment
gest_SMI_R2 <- data %>%
  group_by(Treatment) %>%
  do({
    model <- lmer(SMI ~ Gestation + (1 | Mother), data = .)
    r.squared_vals <- r.squaredGLMM(model)
    
    coefficients <- fixef(model)  # Extract fixed effects
    slope <- coefficients["Gestation"]  # Slope for SVL
    intercept <- coefficients["(Intercept)"]  # Intercept
    
    tibble(
      r.squared = r.squared_vals[1], # Marginal R-squared
      p.value = coef(summary(model))[, "Pr(>|t|)"][2],
      df = summary(model)$coefficients[2, "df"],
      statistic = coef(summary(model))[, "t value"][2],
      slope = slope,
      intercept = intercept
    )
  }) %>%
  select(Treatment, df, statistic, r.squared, p.value, slope, intercept)

gest_SMI_R2


### Trait interaction models (LMM) [Sup Table 5] ####

# neonate ass
Int_model_mass <- lmer(Mass ~ SVL * Gestation * Litter_size * Treatment + (1 | Mother), data = data)
summary(Int_model_mass)

# neonate SVL
Int_model_svl <- lmer(SVL ~ Mass * Gestation * Litter_size * Treatment + (1 | Mother), data = data)
summary(Int_model_svl)

# neonate SMI
Int_model_smi <- lmer(SMI ~ Mass * SVL * Gestation * Litter_size * Treatment + (1 | Mother), data = data)
summary(Int_model_smi)

### correlation between ranked maternal size and traits (Kendall's tau) [Sup Table 6] ####

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

# Relationship between relative maternal body size and neonate SMI by treatment
SMI_tau <- function(df) {
  cor_test <- cor.test(df$MomRelSize, df$SMI, method = "kendall")
  tibble(
    Treatment = unique(df$Treatment),
    tau = cor_test$estimate,
    p.value = cor_test$p.value
  )
}

# Apply the function to each treatment group
tau_SMI <- data %>%
  group_by(Treatment) %>%
  do(SMI_tau(.))

print(tau_SMI)

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


# Figures #####

# Load additional required packages for figures
library(ggpubr)
library(grid)
library(gridExtra)
library(ggpmisc)
library(emmeans)
library(ggsignif)



### Bar plot for birthrate ####

# create data frame for proportion of females that successfully reached praturation for each condition
mdf3 <- data.frame(Treatment=c("control", "treatment"), gave_birth=c((14/40), (12/40)))

#set color objects
GRAY<-c("grey91","grey45")

# Make barplot using percents
birthplot <-  ggplot(mdf3, aes(x=Treatment, y=gave_birth)) +
  geom_bar(stat = "identity", 
           fill=c(GRAY),
           colour="black")+
  ylim(0, 1) + 
  theme_minimal()+
  theme(axis.title.x = element_blank()) +
  ylab("Birth rate") + 
  annotate("text",label="ns", x=1.5, y=0.7, size=4, colour = "black") +
  labs(x = "Treatment") 

# A 2 proportions z.test to compare differences in whether they gave birth or not
prop.test(x= c(12, 14), n= c(40,40), p = NULL, alternative = "two.sided", correct = TRUE)

#birthplot






# Boxplots #

### Box plots for Linear Mixed Effects Models #####

# Mass #

# Fit the linear mixed model
Mass_model <- lmer(Mass ~ Treatment + (1 | Mother), data = data)
summary(Mass_model)

# Perform pairwise comparisons
emm <- emmeans(Mass_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]

# Create the boxplot with ggplot2
p_mass <- ggplot(data, aes(x = Treatment, y = Mass)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Mass (g)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_mass)


# SVL #

# Fit the linear mixed model
SVL_model <- lmer(SVL ~ Treatment + (1 | Mother), data = data)

# Perform pairwise comparisons
emm <- emmeans(SVL_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_SVL <- ggplot(data, aes(x = Treatment, y = SVL)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "SVL (cm)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_SVL)


# SMI #

#Fit the linear mixed model
SMI_model <- lmer(SMI ~ Treatment + (1 | Mother), data = data)

# Perform pairwise comparisons
emm <- emmeans(SMI_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_SMI <- ggplot(data, aes(x = Treatment, y = SMI)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "SMI") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_SMI)


# Gestation #

#Fit the linear mixed model
Gestation_model <- lmer(Gestation ~ Treatment + (1 | Mother), data = data)

# Perform pairwise comparisons
emm <- emmeans(Gestation_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_Gestation <- ggplot(data, aes(x = Treatment, y = Gestation)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Gestation (days)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_Gestation)

### Box plots for Litter effects ######

# total litter mass LM #

#Fit the linear model for total litter mass by treatment
tot_mass_model <- lm(TotalMass ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(tot_mass_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_totmass <- ggplot(data_by_mother, aes(x = Treatment, y = TotalMass)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Total litter mass (g)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_totmass)


# average neonate mass per litter LM #

#Fit the linear model
ave_mass_model <- lm(AverageMass ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(ave_mass_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_avemass <- ggplot(data_by_mother, aes(x = Treatment, y = AverageMass)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Average neonate mass by litter (g)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_avemass)


# average neonate svl per litter LM #

#Fit the linear model
ave_svl_model <- lm(AverageSVL ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(ave_svl_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_avesvl <- ggplot(data_by_mother, aes(x = Treatment, y = AverageSVL)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Average neonate SVL by litter (cm)") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_avesvl)


# average neonate SMI per litter LM #

#Fit the linear model
ave_SMI_model <- lm(AverageSMI ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(ave_SMI_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_aveSMI <- ggplot(data_by_mother, aes(x = Treatment, y = AverageSMI)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Average neonate SMI by litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_aveSMI)


# Stillborns per litter LM #

#Fit the linear model
still_model <- lm(stillborn ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(still_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_still <- ggplot(data_by_mother, aes(x = Treatment, y = stillborn)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Stillbirths per litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_still)


# Deformities per litter LM #

#Fit the linear model
deform_model <- lm(deformed ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(deform_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_deform <- ggplot(data_by_mother, aes(x = Treatment, y = deformed)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Defomred neonates per litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_deform)


# Average litter size LM #

#Fit the linear model
size_model <- lm(Litter_size ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(size_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_size <- ggplot(data_by_mother, aes(x = Treatment, y = Litter_size)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Litter size") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
#print(p_size)



# FIGURE 1: box plots ####

ggarrange(p_Gestation, p_mass, p_SVL, p_SMI, p_totmass, p_still, p_size, birthplot,
          labels= c("A", "B", "C", "D","E", "F", "G", "H"), 
          ncol = 4, nrow = 2)


# Supplemental Figure 1: box plots #### 
ggarrange(p_deform, p_avesvl, p_avemass, p_aveSMI,
          labels= c("A", "B", "C", "D"), 
          ncol = 4, nrow = 1)



### Neonate SVL and neonate mass ####
 

# Assign LMM_SVLM values
LMM_SVLM <- svl_mass_R2

# Define X_Trait_SVLM for "Gestation", Y_Trait_SVLM for "Mass", and X_Var_SVLM for unit label
X_Trait_SVLM <- "SVL"
Y_Trait_SVLM <- "Mass"
X_Var_SVLM <- "(cm)"
Y_Var_SVLM <- "(g)"


# Extract and round values for the control group
slope_C_SVLM <- LMM_SVLM$slope[LMM_SVLM$Treatment == "control"]
intercept_C_SVLM <- LMM_SVLM$intercept[LMM_SVLM$Treatment == "control"]
predefined_r2_C_SVLM <- ifelse(
  round(LMM_SVLM$r.squared[LMM_SVLM$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_SVLM$r.squared[LMM_SVLM$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_SVLM <- LMM_SVLM$slope[LMM_SVLM$Treatment == "treatment"]
intercept_T_SVLM <- LMM_SVLM$intercept[LMM_SVLM$Treatment == "treatment"]
predefined_r2_T_SVLM <- ifelse(
  round(LMM_SVLM$r.squared[LMM_SVLM$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_SVLM$r.squared[LMM_SVLM$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_SVLM <- range(data[[X_Trait_SVLM]][data$Treatment == "control"])
x_limits_T_SVLM <- range(data[[X_Trait_SVLM]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_SVLM <- slope_C_SVLM * x_limits_C_SVLM[1] + intercept_C_SVLM
y_end_C_SVLM <- slope_C_SVLM * x_limits_C_SVLM[2] + intercept_C_SVLM
# Treatment
y_start_T_SVLM <- slope_T_SVLM * x_limits_T_SVLM[1] + intercept_T_SVLM
y_end_T_SVLM <- slope_T_SVLM * x_limits_T_SVLM[2] + intercept_T_SVLM

# Define the position for the R2 annotations
x_annotation_C_SVLM <- max(data[[X_Trait_SVLM]]) * 0.75  # % of the maximum x value
x_annotation_T_SVLM <- max(data[[X_Trait_SVLM]]) * 0.748  # % of the maximum x value
y_annotation_C_SVLM <- max(data[[Y_Trait_SVLM]]) * 0.95  # % of the maximum y value for the first annotation
y_annotation_T_SVLM <- max(data[[Y_Trait_SVLM]]) * 0.92  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_SVLM <- ggplot(data, aes(x = .data[[X_Trait_SVLM]], y = .data[[Y_Trait_SVLM]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_SVLM[1], y = y_start_C_SVLM, xend = x_limits_C_SVLM[2], yend = y_end_C_SVLM), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_SVLM, y = y_annotation_C_SVLM, label = paste(predefined_r2_C_SVLM), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_SVLM[1], y = y_start_T_SVLM, xend = x_limits_T_SVLM[2], yend = y_end_T_SVLM), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_SVLM, X_Var_SVLM), y = paste(Y_Trait_SVLM, Y_Var_SVLM)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_SVLM, y = y_annotation_T_SVLM, label = paste(predefined_r2_T_SVLM), 
           parse = TRUE, color = "#00AFBB")

#print(plot_SVLM)


### Neonate SVL and neonate SMI ####


# Assign LMM_SVLSMI values
LMM_SVLSMI <- svl_SMI_R2

# Define X_Trait_SVLSMI for "Gestation", Y_Trait_SVLSMI for "Mass", and X_Var_SVLSMI for unit label
X_Trait_SVLSMI <- "SVL"
Y_Trait_SVLSMI <- "SMI"
X_Var_SVLSMI <- "(cm)"
Y_Var_SVLSMI <- ""


# Extract and round values for the control group
slope_C_SVLSMI <- LMM_SVLSMI$slope[LMM_SVLSMI$Treatment == "control"]
intercept_C_SVLSMI <- LMM_SVLSMI$intercept[LMM_SVLSMI$Treatment == "control"]
predefined_r2_C_SVLSMI <- ifelse(
  round(LMM_SVLSMI$r.squared[LMM_SVLSMI$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_SVLSMI$r.squared[LMM_SVLSMI$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_SVLSMI <- LMM_SVLSMI$slope[LMM_SVLSMI$Treatment == "treatment"]
intercept_T_SVLSMI <- LMM_SVLSMI$intercept[LMM_SVLSMI$Treatment == "treatment"]
predefined_r2_T_SVLSMI <- ifelse(
  round(LMM_SVLSMI$r.squared[LMM_SVLSMI$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_SVLSMI$r.squared[LMM_SVLSMI$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_SVLSMI <- range(data[[X_Trait_SVLSMI]][data$Treatment == "control"])
x_limits_T_SVLSMI <- range(data[[X_Trait_SVLSMI]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_SVLSMI <- slope_C_SVLSMI * x_limits_C_SVLSMI[1] + intercept_C_SVLSMI
y_end_C_SVLSMI <- slope_C_SVLSMI * x_limits_C_SVLSMI[2] + intercept_C_SVLSMI
# Treatment
y_start_T_SVLSMI <- slope_T_SVLSMI * x_limits_T_SVLSMI[1] + intercept_T_SVLSMI
y_end_T_SVLSMI <- slope_T_SVLSMI * x_limits_T_SVLSMI[2] + intercept_T_SVLSMI

# Define the position for the R2 annotations
x_annotation_C_SVLSMI <- max(data[[X_Trait_SVLSMI]]) * 0.75  # % of the maximum x value
x_annotation_T_SVLSMI <- max(data[[X_Trait_SVLSMI]]) * 0.75  # % of the maximum x value
y_annotation_C_SVLSMI <- max(data[[Y_Trait_SVLSMI]]) * 0.95  # % of the maximum y value for the first annotation
y_annotation_T_SVLSMI <- max(data[[Y_Trait_SVLSMI]]) * 0.92  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_SVLSMI <- ggplot(data, aes(x = .data[[X_Trait_SVLSMI]], y = .data[[Y_Trait_SVLSMI]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_SVLSMI[1], y = y_start_C_SVLSMI, xend = x_limits_C_SVLSMI[2], yend = y_end_C_SVLSMI), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_SVLSMI, y = y_annotation_C_SVLSMI, label = paste(predefined_r2_C_SVLSMI), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_SVLSMI[1], y = y_start_T_SVLSMI, xend = x_limits_T_SVLSMI[2], yend = y_end_T_SVLSMI), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_SVLSMI, X_Var_SVLSMI), y = paste(Y_Trait_SVLSMI, Y_Var_SVLSMI)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_SVLSMI, y = y_annotation_T_SVLSMI, label = paste(predefined_r2_T_SVLSMI), 
           parse = TRUE, color = "#00AFBB")

#print(plot_SVLSMI)


### Neonate mass and neonate SMI ####


# Assign LMM_MSMI values
LMM_MSMI <- mass_SMI_R2

# Define X_Trait_MSMI for "Gestation", Y_Trait_MSMI for "Mass", and X_Var_MSMI for unit label
X_Trait_MSMI <- "Mass"
Y_Trait_MSMI <- "SMI"
X_Var_MSMI <- "(g)"
Y_Var_MSMI <- ""


# Extract and round values for the control group
slope_C_MSMI <- LMM_MSMI$slope[LMM_MSMI$Treatment == "control"]
intercept_C_MSMI <- LMM_MSMI$intercept[LMM_MSMI$Treatment == "control"]
predefined_r2_C_MSMI <- ifelse(
  round(LMM_MSMI$r.squared[LMM_MSMI$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_MSMI$r.squared[LMM_MSMI$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_MSMI <- LMM_MSMI$slope[LMM_MSMI$Treatment == "treatment"]
intercept_T_MSMI <- LMM_MSMI$intercept[LMM_MSMI$Treatment == "treatment"]
predefined_r2_T_MSMI <- ifelse(
  round(LMM_MSMI$r.squared[LMM_MSMI$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_MSMI$r.squared[LMM_MSMI$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_MSMI <- range(data[[X_Trait_MSMI]][data$Treatment == "control"])
x_limits_T_MSMI <- range(data[[X_Trait_MSMI]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_MSMI <- slope_C_MSMI * x_limits_C_MSMI[1] + intercept_C_MSMI
y_end_C_MSMI <- slope_C_MSMI * x_limits_C_MSMI[2] + intercept_C_MSMI
# Treatment
y_start_T_MSMI <- slope_T_MSMI * x_limits_T_MSMI[1] + intercept_T_MSMI
y_end_T_MSMI <- slope_T_MSMI * x_limits_T_MSMI[2] + intercept_T_MSMI

# Define the position for the R2 annotations
x_annotation_C_MSMI <- max(data[[X_Trait_MSMI]]) * 0.55  # % of the maximum x value
x_annotation_T_MSMI <- max(data[[X_Trait_MSMI]]) * 0.55  # % of the maximum x value
y_annotation_C_MSMI <- max(data[[Y_Trait_MSMI]]) * 0.95  # % of the maximum y value for the first annotation
y_annotation_T_MSMI <- max(data[[Y_Trait_MSMI]]) * 0.92  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_MSMI <- ggplot(data, aes(x = .data[[X_Trait_MSMI]], y = .data[[Y_Trait_MSMI]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_MSMI[1], y = y_start_C_MSMI, xend = x_limits_C_MSMI[2], yend = y_end_C_MSMI), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_MSMI, y = y_annotation_C_MSMI, label = paste(predefined_r2_C_MSMI), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_MSMI[1], y = y_start_T_MSMI, xend = x_limits_T_MSMI[2], yend = y_end_T_MSMI), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_MSMI, X_Var_MSMI), y = paste(Y_Trait_MSMI, Y_Var_MSMI)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_MSMI, y = y_annotation_T_MSMI, label = paste(predefined_r2_T_MSMI), 
           parse = TRUE, color = "#00AFBB")

#print(plot_MSMI)


# Figure 2: SVL, Mass, and SMI ####

# Remove the legend title using theme() for each plot
plot_SVLSMI <- plot_SVLSMI + theme(legend.title = element_blank())
plot_SVLM  <- plot_SVLM + theme(legend.title = element_blank())
plot_MSMI  <- plot_MSMI + theme(legend.title = element_blank())


# plot all basic neonate comparisons FIGURE 2 
ggarrange(plot_SVLSMI, plot_SVLM, plot_MSMI, 
          labels= c("A", "B", "C"), 
          ncol = 3, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")



### Litter size and neonate mass by treatment ####


# Assign LMM_LM values
LMM_LM <- lit_mass_R2

# Define X_Trait_LM for "Gestation", Y_Trait_LM for "Mass", and X_Var_LM for unit label
X_Trait_LM <- "Litter_size"
Y_Trait_LM <- "Mass"
X_Var_LM <- ""
Y_Var_LM <- "(g)"


# Extract and round values for the control group
slope_C_LM <- LMM_LM$slope[LMM_LM$Treatment == "control"]
intercept_C_LM <- LMM_LM$intercept[LMM_LM$Treatment == "control"]
predefined_r2_C_LM <- ifelse(
  round(LMM_LM$r.squared[LMM_LM$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LM$r.squared[LMM_LM$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_LM <- LMM_LM$slope[LMM_LM$Treatment == "treatment"]
intercept_T_LM <- LMM_LM$intercept[LMM_LM$Treatment == "treatment"]
predefined_r2_T_LM <- ifelse(
  round(LMM_LM$r.squared[LMM_LM$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LM$r.squared[LMM_LM$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_LM <- range(data[[X_Trait_LM]][data$Treatment == "control"])
x_limits_T_LM <- range(data[[X_Trait_LM]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_LM <- slope_C_LM * x_limits_C_LM[1] + intercept_C_LM
y_end_C_LM <- slope_C_LM * x_limits_C_LM[2] + intercept_C_LM
# Treatment
y_start_T_LM <- slope_T_LM * x_limits_T_LM[1] + intercept_T_LM
y_end_T_LM <- slope_T_LM * x_limits_T_LM[2] + intercept_T_LM

# Define the position for the R2 annotations
x_annotation_C_LM <- max(data[[X_Trait_LM]]) * 0.892  # % of the maximum x value
x_annotation_T_LM <- max(data[[X_Trait_LM]]) * 0.90  # % of the maximum x value
y_annotation_C_LM <- max(data[[Y_Trait_LM]]) * 0.85  # % of the maximum y value for the first annotation
y_annotation_T_LM <- max(data[[Y_Trait_LM]]) * 0.82  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_LM <- ggplot(data, aes(x = .data[[X_Trait_LM]], y = .data[[Y_Trait_LM]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_LM[1], y = y_start_C_LM, xend = x_limits_C_LM[2], yend = y_end_C_LM), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_LM, y = y_annotation_C_LM, label = paste(predefined_r2_C_LM), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_LM[1], y = y_start_T_LM, xend = x_limits_T_LM[2], yend = y_end_T_LM), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_LM, X_Var_LM), y = paste(Y_Trait_LM, Y_Var_LM)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_LM, y = y_annotation_T_LM, label = paste(predefined_r2_T_LM), 
           parse = TRUE, color = "#00AFBB")

#print(plot_LM)



### Litter size and neonate SVL ####


# Assign LMM_LSVL values
LMM_LSVL <- lit_SVL_R2

# Define X_Trait_LSVL for "Gestation", Y_Trait_LSVL for "Mass", and X_Var_LSVL for unit label
X_Trait_LSVL <- "Litter_size"
Y_Trait_LSVL <- "SVL"
X_Var_LSVL <- ""
Y_Var_LSVL <- "(cm)"


# Extract and round values for the control group
slope_C_LSVL <- LMM_LSVL$slope[LMM_LSVL$Treatment == "control"]
intercept_C_LSVL <- LMM_LSVL$intercept[LMM_LSVL$Treatment == "control"]
predefined_r2_C_LSVL <- ifelse(
  round(LMM_LSVL$r.squared[LMM_LSVL$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LSVL$r.squared[LMM_LSVL$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_LSVL <- LMM_LSVL$slope[LMM_LSVL$Treatment == "treatment"]
intercept_T_LSVL <- LMM_LSVL$intercept[LMM_LSVL$Treatment == "treatment"]
predefined_r2_T_LSVL <- ifelse(
  round(LMM_LSVL$r.squared[LMM_LSVL$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LSVL$r.squared[LMM_LSVL$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_LSVL <- range(data[[X_Trait_LSVL]][data$Treatment == "control"])
x_limits_T_LSVL <- range(data[[X_Trait_LSVL]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_LSVL <- slope_C_LSVL * x_limits_C_LSVL[1] + intercept_C_LSVL
y_end_C_LSVL <- slope_C_LSVL * x_limits_C_LSVL[2] + intercept_C_LSVL
# Treatment
y_start_T_LSVL <- slope_T_LSVL * x_limits_T_LSVL[1] + intercept_T_LSVL
y_end_T_LSVL <- slope_T_LSVL * x_limits_T_LSVL[2] + intercept_T_LSVL

# Define the position for the R2 annotations
x_annotation_C_LSVL <- max(data[[X_Trait_LSVL]]) * 0.20  # % of the maximum x value
x_annotation_T_LSVL <- max(data[[X_Trait_LSVL]]) * 0.20  # % of the maximum x value
y_annotation_C_LSVL <- max(data[[Y_Trait_LSVL]]) * 0.80  # % of the maximum y value for the first annotation
y_annotation_T_LSVL <- max(data[[Y_Trait_LSVL]]) * 0.78  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_LSVL <- ggplot(data, aes(x = .data[[X_Trait_LSVL]], y = .data[[Y_Trait_LSVL]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_LSVL[1], y = y_start_C_LSVL, xend = x_limits_C_LSVL[2], yend = y_end_C_LSVL), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_LSVL, y = y_annotation_C_LSVL, label = paste(predefined_r2_C_LSVL), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_LSVL[1], y = y_start_T_LSVL, xend = x_limits_T_LSVL[2], yend = y_end_T_LSVL), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_LSVL, X_Var_LSVL), y = paste(Y_Trait_LSVL, Y_Var_LSVL)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_LSVL, y = y_annotation_T_LSVL, label = paste(predefined_r2_T_LSVL), 
           parse = TRUE, color = "#00AFBB")

#print(plot_LSVL)



### Litter size and gestation length ####


# Assign LMM_LG values
LMM_LG <- lit_gest_R2

# Define X_Trait_LG for "Gestation", Y_Trait_LG for "Mass", and X_Var_LG for unit label
X_Trait_LG <- "Litter_size"
Y_Trait_LG <- "Gestation"
X_Var_LG <- ""
Y_Var_LG <- "(days)"


# Extract and round values for the control group
slope_C_LG <- LMM_LG$slope[LMM_LG$Treatment == "control"]
intercept_C_LG <- LMM_LG$intercept[LMM_LG$Treatment == "control"]
predefined_r2_C_LG <- ifelse(
  round(LMM_LG$r.squared[LMM_LG$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LG$r.squared[LMM_LG$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_LG <- LMM_LG$slope[LMM_LG$Treatment == "treatment"]
intercept_T_LG <- LMM_LG$intercept[LMM_LG$Treatment == "treatment"]
predefined_r2_T_LG <- ifelse(
  round(LMM_LG$r.squared[LMM_LG$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LG$r.squared[LMM_LG$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_LG <- range(data[[X_Trait_LG]][data$Treatment == "control"])
x_limits_T_LG <- range(data[[X_Trait_LG]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_LG <- slope_C_LG * x_limits_C_LG[1] + intercept_C_LG
y_end_C_LG <- slope_C_LG * x_limits_C_LG[2] + intercept_C_LG
# Treatment
y_start_T_LG <- slope_T_LG * x_limits_T_LG[1] + intercept_T_LG
y_end_T_LG <- slope_T_LG * x_limits_T_LG[2] + intercept_T_LG

# Define the position for the R2 annotations
x_annotation_C_LG <- max(data[[X_Trait_LG]]) * 0.85  # % of the maximum x value
x_annotation_T_LG <- max(data[[X_Trait_LG]]) * 0.85  # % of the maximum x value
y_annotation_C_LG <- max(data[[Y_Trait_LG]]) * 0.95  # % of the maximum y value for the first annotation
y_annotation_T_LG <- max(data[[Y_Trait_LG]]) * 0.92  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_LG <- ggplot(data, aes(x = .data[[X_Trait_LG]], y = .data[[Y_Trait_LG]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_LG[1], y = y_start_C_LG, xend = x_limits_C_LG[2], yend = y_end_C_LG), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_LG, y = y_annotation_C_LG, label = paste(predefined_r2_C_LG), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_LG[1], y = y_start_T_LG, xend = x_limits_T_LG[2], yend = y_end_T_LG), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_LG, X_Var_LG), y = paste(Y_Trait_LG, Y_Var_LG)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_LG, y = y_annotation_T_LG, label = paste(predefined_r2_T_LG), 
           parse = TRUE, color = "#00AFBB")

#print(plot_LG)

### Litter size and neonate SMI #### 


# Assign LMM_LSMI values
LMM_LSMI <- lit_SMI_R2

# Define X_Trait_LSMI for "Gestation", Y_Trait_LSMI for "Mass", and X_Var_LSMI for unit label
X_Trait_LSMI <- "Litter_size"
Y_Trait_LSMI <- "SMI"
X_Var_LSMI <- ""
Y_Var_LSMI <- ""


# Extract and round values for the control group
slope_C_LSMI <- LMM_LSMI$slope[LMM_LSMI$Treatment == "control"]
intercept_C_LSMI <- LMM_LSMI$intercept[LMM_LSMI$Treatment == "control"]
predefined_r2_C_LSMI <- ifelse(
  round(LMM_LSMI$r.squared[LMM_LSMI$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LSMI$r.squared[LMM_LSMI$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_LSMI <- LMM_LSMI$slope[LMM_LSMI$Treatment == "treatment"]
intercept_T_LSMI <- LMM_LSMI$intercept[LMM_LSMI$Treatment == "treatment"]
predefined_r2_T_LSMI <- ifelse(
  round(LMM_LSMI$r.squared[LMM_LSMI$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_LSMI$r.squared[LMM_LSMI$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_LSMI <- range(data[[X_Trait_LSMI]][data$Treatment == "control"])
x_limits_T_LSMI <- range(data[[X_Trait_LSMI]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_LSMI <- slope_C_LSMI * x_limits_C_LSMI[1] + intercept_C_LSMI
y_end_C_LSMI <- slope_C_LSMI * x_limits_C_LSMI[2] + intercept_C_LSMI
# Treatment
y_start_T_LSMI <- slope_T_LSMI * x_limits_T_LSMI[1] + intercept_T_LSMI
y_end_T_LSMI <- slope_T_LSMI * x_limits_T_LSMI[2] + intercept_T_LSMI

# Define the position for the R2 annotations
x_annotation_C_LSMI <- max(data[[X_Trait_LSMI]]) * 0.85  # % of the maximum x value
x_annotation_T_LSMI <- max(data[[X_Trait_LSMI]]) * 0.85  # % of the maximum x value
y_annotation_C_LSMI <- max(data[[Y_Trait_LSMI]]) * 0.95  # % of the maximum y value for the first annotation
y_annotation_T_LSMI <- max(data[[Y_Trait_LSMI]]) * 0.92  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_LSMI <- ggplot(data, aes(x = .data[[X_Trait_LSMI]], y = .data[[Y_Trait_LSMI]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_LSMI[1], y = y_start_C_LSMI, xend = x_limits_C_LSMI[2], yend = y_end_C_LSMI), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_LSMI, y = y_annotation_C_LSMI, label = paste(predefined_r2_C_LSMI), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_LSMI[1], y = y_start_T_LSMI, xend = x_limits_T_LSMI[2], yend = y_end_T_LSMI), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_LSMI, X_Var_LSMI), y = paste(Y_Trait_LSMI, Y_Var_LSMI)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_LSMI, y = y_annotation_T_LSMI, label = paste(predefined_r2_T_LSMI), 
           parse = TRUE, color = "#00AFBB")

#print(plot_LSMI)


# Supplemental Figure 2: Litter scatter plots####
# Remove the legend title using theme() for each plot
plot_LM <- plot_LM  + theme(legend.title = element_blank())
plot_LSVL  <- plot_LSVL + theme(legend.title = element_blank())
plot_LG  <- plot_LG + theme(legend.title = element_blank())
plot_LSMI <- plot_LSMI + theme(legend.title = element_blank())

# plot all litter size figures 
ggarrange(plot_LM, plot_LSVL, plot_LG, plot_LSMI, 
          labels= c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")


### Gestation length and neonate mass ####

# Assign LMM_GM values
LMM_GM <- gest_mass_R2

# Define X_Trait_GM for "Gestation", Y_Trait_GM for "Mass", and X_Var_GM for unit label
X_Trait_GM <- "Gestation"
Y_Trait_GM <- "Mass"
X_Var_GM <- "(days)"
Y_Var_GM <- "(g)"


# Extract and round values for the control group
slope_C_GM <- LMM_GM$slope[LMM_GM$Treatment == "control"]
intercept_C_GM <- LMM_GM$intercept[LMM_GM$Treatment == "control"]
predefined_r2_C_GM <- ifelse(
  round(LMM_GM$r.squared[LMM_GM$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GM$r.squared[LMM_GM$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_GM <- LMM_GM$slope[LMM_GM$Treatment == "treatment"]
intercept_T_GM <- LMM_GM$intercept[LMM_GM$Treatment == "treatment"]
predefined_r2_T_GM <- ifelse(
  round(LMM_GM$r.squared[LMM_GM$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GM$r.squared[LMM_GM$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_GM <- range(data[[X_Trait_GM]][data$Treatment == "control"])
x_limits_T_GM <- range(data[[X_Trait_GM]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_GM <- slope_C_GM * x_limits_C_GM[1] + intercept_C_GM
y_end_C_GM <- slope_C_GM * x_limits_C_GM[2] + intercept_C_GM
# Treatment
y_start_T_GM <- slope_T_GM * x_limits_T_GM[1] + intercept_T_GM
y_end_T_GM <- slope_T_GM * x_limits_T_GM[2] + intercept_T_GM

# Define the position for the R2 annotations
x_annotation_C_GM <- max(data[[X_Trait_GM]]) * 0.95  # % of the maximum x value
x_annotation_T_GM <- max(data[[X_Trait_GM]]) * 0.95  # % of the maximum x value
y_annotation_C_GM <- max(data[[Y_Trait_GM]]) * 0.85  # % of the maximum y value for the first annotation
y_annotation_T_GM <- max(data[[Y_Trait_GM]]) * 0.81  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_GM <- ggplot(data, aes(x = .data[[X_Trait_GM]], y = .data[[Y_Trait_GM]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_GM[1], y = y_start_C_GM, xend = x_limits_C_GM[2], yend = y_end_C_GM), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_GM, y = y_annotation_C_GM, label = paste(predefined_r2_C_GM), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_GM[1], y = y_start_T_GM, xend = x_limits_T_GM[2], yend = y_end_T_GM), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_GM, X_Var_GM), y = paste(Y_Trait_GM, Y_Var_GM)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_GM, y = y_annotation_T_GM, label = paste(predefined_r2_T_GM), 
           parse = TRUE, color = "#00AFBB")

#print(plot_GM)


### Gestation length and neonate SVL ####
# Gestation v SVL

# Assign LMM_GSVL values
LMM_GSVL <- gest_SVL_R2

# Define X_Trait_GSVL for "Gestation", Y_Trait_GSVL for "Mass", and X_Var_GSVL for unit label
X_Trait_GSVL <- "Gestation"
Y_Trait_GSVL <- "SVL"
X_Var_GSVL <- "(days)"
Y_Var_GSVL <- "(cm)"


# Extract and round values for the control group
slope_C_GSVL <- LMM_GSVL$slope[LMM_GSVL$Treatment == "control"]
intercept_C_GSVL <- LMM_GSVL$intercept[LMM_GSVL$Treatment == "control"]
predefined_r2_C_GSVL <- ifelse(
  round(LMM_GSVL$r.squared[LMM_GSVL$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GSVL$r.squared[LMM_GSVL$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_GSVL <- LMM_GSVL$slope[LMM_GSVL$Treatment == "treatment"]
intercept_T_GSVL <- LMM_GSVL$intercept[LMM_GSVL$Treatment == "treatment"]
predefined_r2_T_GSVL <- ifelse(
  round(LMM_GSVL$r.squared[LMM_GSVL$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GSVL$r.squared[LMM_GSVL$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_GSVL <- range(data[[X_Trait_GSVL]][data$Treatment == "control"])
x_limits_T_GSVL <- range(data[[X_Trait_GSVL]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_GSVL <- slope_C_GSVL * x_limits_C_GSVL[1] + intercept_C_GSVL
y_end_C_GSVL <- slope_C_GSVL * x_limits_C_GSVL[2] + intercept_C_GSVL
# Treatment
y_start_T_GSVL <- slope_T_GSVL * x_limits_T_GSVL[1] + intercept_T_GSVL
y_end_T_GSVL <- slope_T_GSVL * x_limits_T_GSVL[2] + intercept_T_GSVL

# Define the position for the R2 annotations
x_annotation_C_GSVL <- max(data[[X_Trait_GSVL]]) * 0.95  # % of the maximum x value
x_annotation_T_GSVL <- max(data[[X_Trait_GSVL]]) * 0.95  # % of the maximum x value
y_annotation_C_GSVL <- max(data[[Y_Trait_GSVL]]) * 0.75  # % of the maximum y value for the first annotation
y_annotation_T_GSVL <- max(data[[Y_Trait_GSVL]]) * 0.73  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_GSVL <- ggplot(data, aes(x = .data[[X_Trait_GSVL]], y = .data[[Y_Trait_GSVL]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_GSVL[1], y = y_start_C_GSVL, xend = x_limits_C_GSVL[2], yend = y_end_C_GSVL), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_GSVL, y = y_annotation_C_GSVL, label = paste(predefined_r2_C_GSVL), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_GSVL[1], y = y_start_T_GSVL, xend = x_limits_T_GSVL[2], yend = y_end_T_GSVL), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_GSVL, X_Var_GSVL), y = paste(Y_Trait_GSVL, Y_Var_GSVL)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_GSVL, y = y_annotation_T_GSVL, label = paste(predefined_r2_T_GSVL), 
           parse = TRUE, color = "#00AFBB")

#print(plot_GSVL)

### Gestation length and neonate SMI ####
# Gestation v SMI

# Assign LMM_GSMI values
LMM_GSMI <- gest_SMI_R2

# Define X_Trait_GSMI for "Gestation", Y_Trait_GSMI for "Mass", and X_Var_GSMI for unit label
X_Trait_GSMI <- "Gestation"
Y_Trait_GSMI <- "SMI"
X_Var_GSMI <- "(days)"
Y_Var_GSMI <- ""


# Extract and round values for the control group
slope_C_GSMI <- LMM_GSMI$slope[LMM_GSMI$Treatment == "control"]
intercept_C_GSMI <- LMM_GSMI$intercept[LMM_GSMI$Treatment == "control"]
predefined_r2_C_GSMI <- ifelse(
  round(LMM_GSMI$r.squared[LMM_GSMI$Treatment == "control"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GSMI$r.squared[LMM_GSMI$Treatment == "control"], 3))
)

# Extract and round values for the treatment group
slope_T_GSMI <- LMM_GSMI$slope[LMM_GSMI$Treatment == "treatment"]
intercept_T_GSMI <- LMM_GSMI$intercept[LMM_GSMI$Treatment == "treatment"]
predefined_r2_T_GSMI <- ifelse(
  round(LMM_GSMI$r.squared[LMM_GSMI$Treatment == "treatment"], 3) < 0.01,
  "R^2 == \"< 0.01\"",
  paste0("R^2 == ", round(LMM_GSMI$r.squared[LMM_GSMI$Treatment == "treatment"], 3))
)

# Set range limits for regression line
x_limits_C_GSMI <- range(data[[X_Trait_GSMI]][data$Treatment == "control"])
x_limits_T_GSMI <- range(data[[X_Trait_GSMI]][data$Treatment == "treatment"])

# Calculate the y-values for the endpoints of the lines
# Control
y_start_C_GSMI <- slope_C_GSMI * x_limits_C_GSMI[1] + intercept_C_GSMI
y_end_C_GSMI <- slope_C_GSMI * x_limits_C_GSMI[2] + intercept_C_GSMI
# Treatment
y_start_T_GSMI <- slope_T_GSMI * x_limits_T_GSMI[1] + intercept_T_GSMI
y_end_T_GSMI <- slope_T_GSMI * x_limits_T_GSMI[2] + intercept_T_GSMI

# Define the position for the R2 annotations
x_annotation_C_GSMI <- max(data[[X_Trait_GSMI]]) * 0.95  # % of the maximum x value
x_annotation_T_GSMI <- max(data[[X_Trait_GSMI]]) * 0.95  # % of the maximum x value
y_annotation_C_GSMI <- max(data[[Y_Trait_GSMI]]) * 0.85  # % of the maximum y value for the first annotation
y_annotation_T_GSMI <- max(data[[Y_Trait_GSMI]]) * 0.82  # % of the maximum y value for the second annotation

# Plot scatterplot by treatment with LMM regression line
plot_GSMI <- ggplot(data, aes(x = .data[[X_Trait_GSMI]], y = .data[[Y_Trait_GSMI]], color = Treatment)) +
  geom_point() +
  geom_segment(aes(x = x_limits_C_GSMI[1], y = y_start_C_GSMI, xend = x_limits_C_GSMI[2], yend = y_end_C_GSMI), 
               color = "#F8766D", linetype = "solid", linewidth = 1) +
  annotate("text", x = x_annotation_C_GSMI, y = y_annotation_C_GSMI, label = paste(predefined_r2_C_GSMI), 
           parse = TRUE, color = "#F8766D") +
  geom_segment(aes(x = x_limits_T_GSMI[1], y = y_start_T_GSMI, xend = x_limits_T_GSMI[2], yend = y_end_T_GSMI), 
               color = "#00AFBB", linetype = "solid", linewidth = 1) +
  labs(x = paste(X_Trait_GSMI, X_Var_GSMI), y = paste(Y_Trait_GSMI, Y_Var_GSMI)) +
  theme_minimal() +
  annotate("text", x = x_annotation_T_GSMI, y = y_annotation_T_GSMI, label = paste(predefined_r2_T_GSMI), 
           parse = TRUE, color = "#00AFBB")

#print(plot_GSMI)


# Supplemental Figure 3: Gestation scatter plots ####

# Remove the legend title using theme() for each plot
plot_GM <- plot_GM + theme(legend.title = element_blank())
plot_GSVL  <- plot_GSVL + theme(legend.title = element_blank())
plot_GSMI  <- plot_GSMI + theme(legend.title = element_blank())

# plot all gestation length figures
ggarrange(plot_GM, plot_GSVL, plot_GSMI, 
          labels = c("A", "B", "C"), 
          ncol = 3, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")

### correlation between ranked maternal size and traits (Kendall's tau) #####

# Relative (ranked) size of mother and neonate mass

mom_mass <- 
  ggplot(data, aes(x = MomRelSize, y = Mass, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "Mass (g)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "left",
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_mass)

# Relative (ranked) size of mother and neonate SVL
mom_svl <- 
  ggplot(data, aes(x = MomRelSize, y = SVL, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "SVL (cm)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center", 
           label.y.npc = "bottom",
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_svl)

# Relative (ranked) size of mother and neonate SMI
mom_SMI <-
  ggplot(data, aes(x = MomRelSize, y = SMI, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "SMI") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "left", 
           label.y.npc = "bottom",
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_SMI)

# Relative (ranked) size of mother and litter size
mom_lit <- 
  ggplot(data, aes(x = MomRelSize, y = Litter_size, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "Litter size") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center",
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_lit)

# Relative (ranked) size of mother and gestation length
mom_gest <- 
  ggplot(data, aes(x = MomRelSize, y = Gestation, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "Gestation lenght (days)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01,
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_gest)

# Relative (ranked) size of mother and total litter mass
mom_totmass <- 
  ggplot(data_by_mother, aes(x = MomRelSize, y = TotalMass, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mother relative size",
       y = "Total mass/litter (g)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center",
           show.legend = FALSE) + 
  theme_minimal()

#print(mom_totmass)

# Supplemental Figure 4: Relative maternal size figure ####

# Remove the legend title using theme() for each plot
mom_mass <- mom_mass + theme(legend.title = element_blank())
mom_svl  <- mom_svl + theme(legend.title = element_blank())
mom_SMI  <- mom_SMI + theme(legend.title = element_blank())
mom_lit  <- mom_lit + theme(legend.title = element_blank())
mom_gest  <- mom_gest + theme(legend.title = element_blank())
mom_totmass  <- mom_totmass + theme(legend.title = element_blank())

# plot all mom figures
ggarrange(mom_SMI, mom_mass, mom_svl, mom_gest, mom_lit, mom_totmass,
          labels= c("A", "B", "C", "D", "E", "F"), 
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
