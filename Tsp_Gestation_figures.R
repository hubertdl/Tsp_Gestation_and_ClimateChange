# R script used to generate the figures presented in the manuscript


# Load required packages
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggpmisc)

# Clear ggplot2# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)


# Read in counts of DE receptor transcripts
data <- read.csv("Gestation_2018.csv")
data

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

# Add the MomRelSize column to the data dataframe calculated based on ID order
data <- data %>%
  mutate(MomRelSize = 100 - as.numeric(substr(as.character(Mother), 
                                              nchar(as.character(Mother)) - 1, 
                                              nchar(as.character(Mother)))))


#### create bar plot for birthrate ####

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

birthplot


# Load necessary libraries
library(lme4)
library(dplyr)
library(emmeans)
library(ggsignif)

##### Linear Mixed Effects Models #####

#### Mass ####

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
print(p_mass)


##### SVL ######

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
print(p_SVL)


###### BCI #####

#Fit the linear mixed model
BCI_model <- lmer(BCI ~ Treatment + (1 | Mother), data = data)

# Perform pairwise comparisons
emm <- emmeans(BCI_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_BCI <- ggplot(data, aes(x = Treatment, y = BCI)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "BCI") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
print(p_BCI)


###### Gestation #####

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
print(p_Gestation)

##### Litter effects ######

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
  ) %>%
  mutate(MomRelSize = 100 - as.numeric(substr(as.character(Mother), nchar(as.character(Mother)) - 1, nchar(as.character(Mother)))))

# read in DF with stillborn and deformation data
stillborn <- read.csv("StillBorn.csv")

# add that data to the DF
data_by_mother <- left_join(data_by_mother, stillborn, by = "Mother")

#write.csv(data_by_mother, file = "data_by_litter.csv", row.names = FALSE)


###### total litter mass LM #####

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
print(p_totmass)


###### average neonate mass per litter LM #####

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
print(p_avemass)


###### average neonate svl per litter LM #####

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
print(p_avesvl)


##### average neonate BCI per litter LM #####

#Fit the linear model
ave_bci_model <- lm(AverageBCI ~ Treatment, data = data_by_mother)

# Perform pairwise comparisons
emm <- emmeans(ave_bci_model, ~ Treatment)
pairwise_comparisons <- pairs(emm)
pairwise_comparisons_summary <- summary(pairwise_comparisons)

# Extract p-value for significance annotation
p_value <- pairwise_comparisons_summary$p.value[1]  # Assuming one comparison for simplicity

# Create the boxplot with ggplot2
p_avebci <- ggplot(data_by_mother, aes(x = Treatment, y = AverageBCI)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue", aes(group = Treatment)) +
  labs(title = "",
       x = "Treatment",
       y = "Average neonate BCI by litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
print(p_avebci)


###### average stillborns per litter LM #####

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
       y = "Average stillbirths per litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
print(p_still)


###### average deformities per litter LM #####

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
       y = "Average defomred neonates per litter") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("control", "treatment")), 
              map_signif_level = TRUE, 
              annotations = ifelse(p_value < 0.001, "***", 
                                   ifelse(p_value < 0.01, "**", 
                                          ifelse(p_value < 0.05, "*", "ns"))))
# Print the plot
print(p_deform)


###### average litter size LM #####

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
print(p_size)


# plot all ALL boxplots together
ggarrange(p_Gestation, p_mass, p_SVL, p_BCI, p_totmass, p_still, p_deform, p_size, p_avesvl, p_avemass, p_avebci, birthplot,
          labels= c("A", "B", "C", "D","E", "F", "G", "H", "I", "J", "K", "L"), 
          ncol = 4, nrow = 3)




##### correlations #####

# neonate traits

# Correlation between neonate SVL and neonate mass
svl_mass <- 
  ggplot(data, aes(x = SVL, y = Mass, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -1, vjust = 1) +
  geom_point() +
  labs(x = "SVL (cm)",
    y = "Mass (g)") +
  theme_minimal()

print(svl_mass)

# Correlation between neonate SVL and neonate BCI
svl_bci <- ggplot(data, aes(x = SVL, y = BCI, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -1, vjust = 1) +
  geom_point() +
  labs(x = "SVL (cm)",
    y = "BCI") +
  theme_minimal()

print(svl_bci)

# Correlation between neonate mass and neonate BCI
mass_bci <- ggplot(data, aes(x = Mass, y = BCI, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -1, vjust = 1) +
  geom_point() +
  labs(x = "Mass (g)",
    y = "BCI") +
  theme_minimal()

print(mass_bci)

# plot all basic neonate comparisons
ggarrange(svl_bci, svl_mass, mass_bci, 
          labels= c("A", "B", "C"), 
          ncol = 3, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")




# Litter size and neonate traits

# Correlation between Litter size and neonate mass
lit_mass <- ggplot(data, aes(x = Litter_size, y = Mass, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -5.5, vjust = 1) +
  geom_point() +
  labs(x = "Litter size",
    y = "Mass (g)") +
  theme_minimal()

print(lit_mass)

# Correlation between Litter size and neonate SVL
lit_svl <- ggplot(data, aes(x = Litter_size, y = SVL, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -5.5, vjust = 1) +
  geom_point() +
  labs(x = "Litter size",
    y = "SVL (cm)") +
  theme_minimal()

print(lit_svl)

# Correlation between Litter size and gestation length
lit_ges <- ggplot(data, aes(x = Litter_size, y = Gestation, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -5.5, vjust = 1) +
  geom_point() +
  labs(x = "Litter size",
    y = "Gestaton (days)") +
  theme_minimal()

print(lit_ges)

# Correlation between Litter size and neonate BCI
lit_bci <- ggplot(data, aes(x = Litter_size, y = BCI, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -5.5, vjust = 1) +
  geom_point() +
  labs(x = "Litter size",
    y = "BCI") +
  theme_minimal()

print(lit_bci)

# plot all litter size figures 
ggarrange(lit_mass, lit_svl, lit_ges, lit_bci, 
          labels= c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")


# Gestation length

# Correlation between Gestation length and neonate mass
ges_mass <- ggplot(data, aes(x = Gestation, y = Mass, color = Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -2, vjust = 1) +
  geom_point() +
  labs(x = "Gestation (days)",
    y = "Mass (g)") +
  theme_minimal()

print(ges_mass)

# Correlation between Gestation length and neonate SVL
ges_svl <- ggplot(data, aes(x = Gestation, y = SVL, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -2, vjust = 1) +
  geom_point() +
  labs(x = "Gestation (days)",
    y = "SVL (cm)") +
  theme_minimal()

print(ges_svl)

# Correlation between Gestation length and BCI
ges_bci <- ggplot(data, aes(x = Gestation, y = BCI, color = Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -2, vjust = 1) +
  geom_point() +
  labs(x = "Gestation (days)",
    y = "BCI") +
  theme_minimal()

print(ges_bci)

# Correlation between gestation length and litter size
ges_lit <- ggplot(data, aes(x = Gestation, y = Litter_size, color=Treatment)) +
  stat_poly_line(se = TRUE) +
  stat_poly_eq(aes(label = ..rr.label..), 
               formula = y ~ x, 
               parse = TRUE, 
               hjust = -2, vjust = 1) +
  geom_point() +
  labs(x = "Gestaton (days)",
    y = "Litter size") +
  theme_minimal()

print(ges_lit)

# plot all gestation length figures
ggarrange(ges_mass, ges_svl, ges_bci, ges_lit, 
          labels = c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")

###### correlation between ranked mother size and variable traits using Kendall's Tau (ranked variable with >2 options) #####
# Relative (ranked) size of mother and neonate mass

mom_mass <- 
  ggplot(data, aes(x = MomRelSize, y = Mass, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "Mass (g)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "left",
           show.legend = FALSE) + 
  theme_minimal()

print(mom_mass)

# Relative (ranked) size of mother and neonate SVL
mom_svl <- 
  ggplot(data, aes(x = MomRelSize, y = SVL, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "SVL (cm)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center", 
           label.y.npc = "bottom",
           show.legend = FALSE) + 
  theme_minimal()

print(mom_svl)

# Relative (ranked) size of mother and neonate BCI
mom_bci <-
  ggplot(data, aes(x = MomRelSize, y = BCI, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "BCI") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "left", 
           label.y.npc = "bottom",
           show.legend = FALSE) + 
  theme_minimal()

print(mom_bci)

# Relative (ranked) size of mother and litter size
mom_lit <- 
  ggplot(data, aes(x = MomRelSize, y = Litter_size, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "Litter size") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center",
           show.legend = FALSE) + 
  theme_minimal()

print(mom_lit)

# Relative (ranked) size of mother and gestation length
mom_gest <- 
  ggplot(data, aes(x = MomRelSize, y = Gestation, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "Gestation lenght (days)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01,
           show.legend = FALSE) + 
  theme_minimal()

print(mom_gest)

# Relative (ranked) size of mother and total litter mass
mom_totmass <- 
  ggplot(data_by_mother, aes(x = MomRelSize, y = TotalMass, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Mother relative size",
       y = "Total mass/litter (g)") +
  stat_cor(method = "kendall", 
           cor.coef.name = "tau", 
           p.accuracy = 0.01, 
           label.x.npc = "center",
           show.legend = FALSE) + 
  theme_minimal()

print(mom_totmass)

# plot all mom figures
ggarrange(mom_bci, mom_mass, mom_svl, mom_gest, mom_lit, mom_totmass,
          labels= c("A", "B", "C", "D", "E", "F"), 
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")


