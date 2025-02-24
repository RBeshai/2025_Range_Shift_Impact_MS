# Created by [name redacted for peer review]
# Last edited: 01/09/2025


## Load Packages ---------------------------------------------------------------

library(tidyverse)    #loads tidyr, readr, dplyr, and ggplot
library(viridis)      #enables colorblind-friendly palettes 
library(vegan)        #community data functions 
library(lme4)         #for running linear mixed model statistics
library(lmerTest)     #for determining the statistical significance
library(performance)  #can use to check model assumptions with check_model()
library(ggpubr)       #for visualizing qq plots



## Set and Check Working Directory ---------------------------------------------
# Check the working directory 
getwd()



## Clear the Workspace ---------------------------------------------------------

# Even if there's nothing there
rm(list=ls())



## Code basic functions that I like to use -------------------------------------
# Code the "not in" symbol 
"%!in%" <- Negate("%in%")

# Define a rounding function so that 0.5 rounds up 
round.off <- function (x, digits=0) {
  posneg = sign(x)
  z = trunc(abs(x) * 10 ^ (digits + 1)) / 10
  z = floor(z * posneg + 0.5) / 10 ^ digits
  return(z)
}

# Define a function to calculate a running average
running.avg <- function(x) {cumsum(x)/1:length(x)}



# Read in the data -------------------------------------------------------------

# Load the manipulative experiment dataset 
CageData <- read.csv("ExptCommunityData.csv")



# Edit the data ----------------------------------------------------------------

# Data Processing Steps Occurring Below:
  # (1) Calculating the average number of whelks in the plot & join to main dataset
  # (2) Calculating the running average number of whelks in the plot 
        # E.g. If we observed 10 whelks in Week 2, 8 in Week 4, and 9 in Week 6, 
        # the running average would be 9 whelks in the plot
  # (3) Calculate the average density of whelks in a plot 
  # (4) Create a column for target density (the desired density of whelks in the plot)
  # (5) Calculate diversity metrics
  # (6) Calculate the response variables (change in diversity metrics over time)
 

## (1) Add a column for average number of whelks in plot
# Collect all average whelk values in one table 
Avg_Whelks_All <- CageData %>% 
  # Pull in the data
  group_by(Site, Target_Species, Plot_ID) %>% 
  # Group by species and plot
  summarise(Avg_No_Whelks = round.off(mean(No_Target_Whelks_In_Plot)))
  # Calculate the average number of whelks we observed in the plot

# Join this to the main dataset
CageData <- left_join(CageData, Avg_Whelks_All, by = c("Site", "Target_Species", "Plot_ID"))


## (2) Calculate the running average of whelks in the plots
# Run the calculation
RunAvg_Whelks <- CageData %>%
  # Pull in the data
  select("Site", "Target_Species", "Plot_ID", "No_Target_Whelks_In_Plot", "Days_Deployed") %>% 
  # Subset the data
  group_by(Site, Target_Species, Plot_ID) %>% 
  # Group
  mutate(RunAvg_Whelks = (running.avg(No_Target_Whelks_In_Plot)))
  # Calculate running average using function defined above (line 37)

# Join this to the main dataset
CageData <- left_join(CageData, RunAvg_Whelks, by = c("Site", "Target_Species", "Plot_ID", "No_Target_Whelks_In_Plot", "Days_Deployed"))

# Clean dataframe 
rm(Avg_Whelks_All, RunAvg_Whelks)


## (3) Add a column for average plot density 
CageData <- CageData %>%
  # Pull in the data
  mutate(Avg_Density = Avg_No_Whelks * 10)
  # Create a new column 


## (4) Create a column for "Target Density"
CageData <- CageData %>% 
  # Pull in the data
  ungroup(.) %>%    
  # Ungroup 
  mutate(Target_Density = case_when(CageData$Plot_ID == "No Cage" ~ NA,
                                    CageData$Plot_ID == "Partial" ~ NA,
                                    CageData$Plot_ID == "0 Whelk" ~ 0,
                                    CageData$Plot_ID == "1 Whelk" ~ 10, 
                                    CageData$Plot_ID == "3 Whelk" ~ 30,
                                    CageData$Plot_ID == "6 Whelk" ~ 60,
                                    CageData$Plot_ID == "12 Whelk" ~ 120, 
                                    CageData$Plot_ID == "24 Whelk" ~ 240))
  # These values correlate to the target densities within the cages (whelks/m^2)


## (5) Calculate diversity metrics and desired responses 

# Start by creating a new dataframe to hold diversity metrics
DiversityData <- select(CageData, -starts_with("Dead") & -ends_with("Count")) 
  # This removes all "dead" cover and data collected as counts (rather than cover)
  # allows for all columns to be compared with equivalent units of abundance

# Calculate the diversity metrics in this new dataframe 
DiversityData <- DiversityData %>%
  # Pull in the data
  ungroup() %>% 
  # Ungroup 
  mutate(Shannon_Mex = diversity(x = DiversityData[c(18:55, 57:75)], index = "shannon"), 
         Richness_Mex = specnumber(x = DiversityData[c(18:55, 57:75)]),
         Evenness_Mex = Shannon_Mex/log(Richness_Mex),
         Shannon_As = diversity(x = DiversityData[c(18:54, 56:75)], index = "shannon"),
         Richness_As = specnumber(x = DiversityData[c(18:54, 56:75)]),
         Evenness_As = Shannon_As/log(Richness_As))
  # Create new diversity metrics for each species 
  # Each metric is calculated to exclude the driving species 


## (6) Calculate response variables (the change in diversity metrics between each time point)

# Start by cleaning the diversity dataframe
DiversityData <- DiversityData[, c(3:6, 10:13, 78, 81:86)]
  # This removes all the species information 

# Then calculate the change in the metrics from week 0 to week x
# This will be done in 2 separate dataframes (1 per focal whelk)
DiversityMex <- filter(DiversityData, Target_Species == "M") %>% 
  arrange(Site, Plot_ID, Weeks_Deployed) %>%   
  # Order the dataframe
  group_by(Site, Plot_ID) %>%                  
  # Create the groups
  mutate(Change_Shannon_Mex = (Shannon_Mex - first(Shannon_Mex))) %>%
  mutate(Change_Richness_Mex = (Richness_Mex - first(Richness_Mex))) %>% 
  mutate(Change_Evenness_Mex = (Evenness_Mex - first(Evenness_Mex)))
  # Calculate change in diversity metric 

DiversityAs <- filter(DiversityData, Target_Species == "As") %>% 
  arrange(Site, Plot_ID, Weeks_Deployed) %>%   
  # Order the dataframe
  group_by(Site, Plot_ID) %>%                  
  # Create the groups
  mutate(Change_Shannon_As = (Shannon_As - first(Shannon_As))) %>%
  mutate(Change_Richness_As = (Richness_As - first(Richness_As))) %>% 
  mutate(Change_Evenness_As = (Evenness_As - first(Evenness_As)))
  # Calculate change in diversity metric


################################################################################
## Impacts on Diversity and Richness  ------------------------------------------
################################################################################
# All models will be assessing the full 8-week experiment
# Each focal whelk section will have 4x site level models per response (1 per site),
# across 3 diversity responses (n = 12 site-level models per focal whelk). Each set of 
# site level model will be followed by a global model  

# This script will focus on the linear impacts. To assess non-linear impacts
# replace "RunAvg_Whelks" with "poly(RunAvg_Whelks, 2)" in each lm() function

###########################################
# Acanthinucella 
###########################################

## Acanthinucella --> Shannon Diversity (Site) ---------------------------------
# Define the models 
CM_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))

CMS_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))

Dana_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Check mode assumptions
check_model(CM_AsShannon)
check_model(CMS_AsShannon)
check_model(Dana_AsShannon)
check_model(Scripps_AsShannon)

# Determine significance
CM_AS_Shannon_Summary <- summary(CM_AsShannon)
CMS_AS_Shannon_Summary <-summary(CMS_AsShannon)
Dana_AS_Shannon_Summary <-summary(Dana_AsShannon)
Scripps_AS_Shannon_Summary <-summary(Scripps_AsShannon)

# Clean the workspace 
rm(CM_AsShannon, CMS_AsShannon, Dana_AsShannon, Scripps_AsShannon)

## Acanthinucella --> Shannon Diversity (Global Level)
# Create the model
Global_AsShannon_lm <- lmer(formula = Change_Shannon_As ~ RunAvg_Whelks * Region + (1|Site), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))
  # Dropped Region term from random effects to ensure convergence

# Check model assumptions
check_model(Global_AsShannon_lm)

# Determine significance
summary(Global_AsShannon_lm)
  # No overall effect of whelk or region

# Because no significance was found, try a reduced model looking only at the effect of whelks
# Build the model
Global_AsShannon_lm <- lmer(formula = Change_Shannon_As ~ RunAvg_Whelks + (1|Site), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))
  # Did not include region random effect here to ensure it is a subset of the first LME

# Check model assumptions
check_model(Global_AsShannon_lm)

# Determine significance
As_Shannon_Global_Summary <- summary(Global_AsShannon_lm)

# Clean the workspace
rm(Global_AsShannon_lm)


## Acanthinucella --> Richness (Site Level) ------------------------------------
# Define the models 
CM_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))

CMS_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))

Dana_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Check model assumptions 
check_model(CM_AsRichness)
check_model(CMS_AsRichness)
check_model(Dana_AsRichness)
check_model(Scripps_AsRichness)

# Determine significance
CM_AS_Richness_Summary <- summary(CM_AsRichness)
CMS_AS_Richness_Summary <-summary(CMS_AsRichness)
Dana_AS_Richnessn_Summary <-summary(Dana_AsRichness)
Scripps_Richness_Summary <-summary(Scripps_AsRichness)

# Clear the workspace
rm(CM_AsRichness, CMS_AsRichness, Dana_AsRichness, Scripps_AsRichness)

## Acanthinucella --> Richness (Global Level)
# Define the model
Global_AsRich_lm <- lmer(formula = Change_Richness_As ~ RunAvg_Whelks * Region + (1|Site/Region), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))

# Check model assumptions
check_model(Global_AsRich_lm)

# Determine significance
summary(Global_AsRich_lm)

# Because no significance of whelk effect was found, try a reduced model looking only at the effect of whelks
# Define the reduced model
Global_AsRich_lm <- lmer(formula = Change_Richness_As ~ RunAvg_Whelks + (1|Site/Region), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))
  # Throws a convergence error, but inspection of the output suggests that this is a false error
  # lmerTest has a lower threshhold for reporting convergence than lme4
  # Convergence error ignored because the summary output is nearly identical when using different optimizers

# Check model assumptions
check_model(Global_AsRich_lm)

# Determine significance
As_Richness_Global_Summary <- summary(Global_AsRich_lm)
# Summary: no significance

# Clean the workspace
rm(Global_AsRich_lm)


## Acanthinucella --> Evenness (Site Level) ------------------------------------
# Define the models 
CM_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))

CMS_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))

Dana_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Check model assumptions 
check_model(CM_AsEvenness)
check_model(CMS_AsEvenness)
check_model(Dana_AsEvenness)
check_model(Scripps_AsEvenness)

# Determine significance
CM_AS_Evenness_Summary <- summary(CM_AsEvenness)
CMS_AS_Evennes_Summary <-summary(CMS_AsEvenness)
Dana_AS_Evennes_Summary <-summary(Dana_AsEvenness)
Scripps_Evennes_Summary <-summary(Scripps_AsEvenness)

# Clean workspace 
rm(CM_AsEvenness, CMS_AsEvenness, Dana_AsEvenness, Scripps_AsEvenness)

## Acanthinucella --> Evenness (Global Level)
# Define the model
Global_AsEven_lm <- lmer(formula = Change_Evenness_As ~ RunAvg_Whelks * Region + (1|Site/Region), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))

# Check model assumptions
check_model(Global_AsEven_lm)

# Determine significance
summary(Global_AsEven_lm)

# Because no significance of whelk effect was found, try a reduced model looking only at the effect of whelks
# Define reduced model
Global_AsEven_lm <- lmer(formula = Change_Evenness_As ~ RunAvg_Whelks + (1|Site/Region), data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))

# Check model assumptions
check_model(Global_AsEven_lm)

# Determine significance
As_Evenness_Global_Summary <- summary(Global_AsEven_lm)

# Clean the workspace
rm(Global_AsEven_lm)


###########################################
# Mexacanthina
###########################################

## Mexacanthina --> Shannon Diversity (Site Level) -----------------------------
# Define the models 
Dana_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

PM_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))

CK_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

# Check model assumptions 
check_model(Dana_MexShannon)
check_model(Scripps_MexShannon)
check_model(PM_MexShannon)
check_model(CK_MexShannon)

# Determine significance 
Dana_Mex_Shannon_Summary <- summary(Dana_MexShannon)
Scripps_Mex_Shannon_Summary <- summary(Scripps_MexShannon)
PM_Mex_Shannon_Summary <- summary(PM_MexShannon)
CK_Mex_Shannon_Summary <- summary(CK_MexShannon)

# Clean workspace
rm(CK_MexShannon, PM_MexShannon, Dana_MexShannon, Scripps_MexShannon)

## Mexacanthina --> Shannon Diversity (Global Level)
# Define the model 
Global_MexShannon_lm <- lmer(formula = Change_Shannon_Mex ~ RunAvg_Whelks * Region + (1|Site/Region), data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))

# Check model assumptions
check_model(Global_MexShannon_lm)

# Determine significance 
Mex_Shannon_Global_Summary <- summary(Global_MexShannon_lm)

# Clean the workspace
rm(Global_MexShannon_lm)


## Mexacanthina --> Richness (Site Level) --------------------------------------
# Define the models 
Dana_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

PM_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))

CK_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

# Check model assumptions 
check_model(Dana_MexRichness)
check_model(Scripps_MexRichness)
check_model(PM_MexRichness)
check_model(CK_MexRichness)

# Determine significance
Dana_Mex_Richness_Summary <- summary(Dana_MexRichness)
Scripps_Mex_Richness_Summary <- summary(Scripps_MexRichness)
PM_Mex_Richness_Summary <- summary(PM_MexRichness)
CK_Mex_Richness_Summary <- summary(CK_MexRichness)

# Clean the workspace 
rm(CK_MexRichness, PM_MexRichness, Dana_MexRichness, Scripps_MexRichness)


##  Mexacanthina --> Richness (Global Level)
Global_MexRich_lm <- lmer(formula = Change_Richness_Mex ~ RunAvg_Whelks * Region + (1|Site/Region), data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))
    # This throws a convergence error. However, it can be likely ignored because if we try multiple optimizers
    # we don't see convergence issues in the summary and receive nearly identical outputs. 
    # To test, add "control = lmerControl(optimizer ="nloptwrap"). Can change "bobyqa" to "Nelder_Mead" or
    # "nloptwrap"
    # lmerTest has lower threshold for convergence errors than lme4, meaning that the 
    # error is thrown from from the lmerTest overlay rather than lme4 itself  

# Check model assumptions
check_model(Global_MexRich_lm) 
  # Model assumptions look good 

# Determine significance
Mex_Richness_Global_Summary <- summary(Global_MexRich_lm)

# Clean the workspace
rm(Global_MexRich_lm)


##  Mexacanthina --> Evenness (Site Level) -------------------------------------
# Define the models 
Dana_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))

Scripps_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

PM_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))

CK_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

# Check model assumptions 
check_model(Dana_MexEvenness)
check_model(Scripps_MexEvenness)
check_model(PM_MexEvenness)
check_model(CK_MexEvenness)

# Determine significance 
Dana_Mex_Evenness_Summary <- summary(Dana_MexEvenness)
Scripps_Mex_Evenness_Summary <- summary(Scripps_MexEvenness)
PM_Mex_Evenness_Summary <- summary(PM_MexEvenness)
CK_Mex_Evenness_Summary <- summary(CK_MexEvenness)

# Clean the workspace 
rm(CK_MexEvenness, PM_MexEvenness, Dana_MexEvenness, Scripps_MexEvenness)

##  Mexacanthina --> Evenness (Global Level)
# Define the global model
Global_MexEven_lm <- lmer(formula = Change_Evenness_Mex ~ RunAvg_Whelks * Region + (1|Site/Region), data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8))

# Check model assumptions 
check_model(Global_MexEven_lm)

# Determine significance
Mex_Evenness_Global_Summary <- summary(Global_MexEven_lm)

# Clean workspace
rm(Global_MexEven_lm)
 

###############################################
# MODEL OUTPUTS (Run these to see significance)
###############################################

#################
# Acanthinucella
#################

# Site-level Acanthinucella on Shannon Diversity 
CM_AS_Shannon_Summary 
CMS_AS_Shannon_Summary 
Dana_AS_Shannon_Summary 
Scripps_AS_Shannon_Summary

# Site-level Acanthinucella on Richness 
CM_AS_Richness_Summary 
CMS_AS_Richness_Summary
Dana_AS_Richnessn_Summary 
Scripps_Richness_Summary 

# Site-level Acanthinucella on Evenness 
CM_AS_Evenness_Summary
CMS_AS_Evennes_Summary 
Dana_AS_Evennes_Summary 
Scripps_Evennes_Summary

# Global-level Acanthinucella Impacts
As_Shannon_Global_Summary
As_Richness_Global_Summary
As_Evenness_Global_Summary

#################
# Mexacanthina
#################

# Site-level Mexacanthina on Shannon Diversity 
Dana_Mex_Shannon_Summary 
Scripps_Mex_Shannon_Summary
PM_Mex_Shannon_Summary
CK_Mex_Shannon_Summary

# Site-level Mexacanthina on Richness 
Dana_Mex_Richness_Summary 
Scripps_Mex_Richness_Summary
PM_Mex_Richness_Summary
CK_Mex_Richness_Summary

# Site-level Mexacanthina on Evenness 
Dana_Mex_Evenness_Summary 
Scripps_Mex_Evenness_Summary
PM_Mex_Evenness_Summary
CK_Mex_Evenness_Summary

# Global-level Mexacanthina Impacts
Mex_Shannon_Global_Summary
Mex_Richness_Global_Summary
Mex_Evenness_Global_Summary

###############################################

# Remove Acanthinuella summaries from workspace
rm(CM_AS_Shannon_Summary, CMS_AS_Shannon_Summary, Dana_AS_Shannon_Summary, Scripps_AS_Shannon_Summary)
rm(CM_AS_Richness_Summary, CMS_AS_Richness_Summary, Dana_AS_Richnessn_Summary, Scripps_Richness_Summary)
rm(CM_AS_Evenness_Summary, CMS_AS_Evennes_Summary, Dana_AS_Evennes_Summary, Scripps_Evennes_Summary)
rm(As_Shannon_Global_Summary, As_Richness_Global_Summary, As_Evenness_Global_Summary)

# Remove Mexacanthina summaries from the workspace
rm(Dana_Mex_Shannon_Summary, Scripps_Mex_Shannon_Summary, PM_Mex_Shannon_Summary, CK_Mex_Shannon_Summary)
rm(Dana_Mex_Richness_Summary, Scripps_Mex_Richness_Summary, PM_Mex_Richness_Summary, CK_Mex_Richness_Summary)
rm(Dana_Mex_Evenness_Summary, Scripps_Mex_Evenness_Summary, PM_Mex_Evenness_Summary, CK_Mex_Evenness_Summary)
rm(Mex_Shannon_Global_Summary, Mex_Richness_Global_Summary, Mex_Evenness_Global_Summary)



################################################################################
## Abundance-Impact Figures (Not presented in MS) ------------------------------
################################################################################

###########################################
# Acanthinucella
###########################################
# To create plot - generate plots for Acanthinucella's impact on each diversity metric
# Then combine into groups based on site 
# Finally, combine all site-diversity metric groups into a single overall plot


## Plot Acanthinucella --> Shannon Diversity (SITE LEVEL)
# Mendocino North
CMShannonPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Mendocino North" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_As)) +
  geom_point() +
  geom_smooth(color = "black", linetype = "dashed", method = "lm", se = F) +
  labs(x = " ", y = "Mendocino North") +
  ggtitle(expression(bold(Delta~"Shannon"))) +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Mendocino South
CMSShannonPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Mendocino South" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_As)) +
  geom_point() +
  labs(x = " ", y = "Mendocino South") +
  #ggtitle("Mendocino South") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Dana Point
DanaShannonPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_As)) +
  geom_point() +
  geom_smooth(color = "black", method = "lm", se = F) +
  labs(x = " ", y = "Dana Point") +
  #ggtitle("Dana Point") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Scripps
ScrippsShannonPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_As)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = "Scripps") +
  #ggtitle("Scripps") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size=12, face = "bold", colour = "black"))


## Plot Acanthinucella --> Richness (SITE LEVEL)
# Mendocino North
CMRichPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Mendocino North" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_As)) +
  geom_point() +
  labs(x = " ", y = " ") +
  ggtitle(expression(bold(Delta~"Richness"))) +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Mendocino South
CMSRichPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Mendocino South" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_As)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Dana Point
DanaRichPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_As)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Scripps
ScrippsRichPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), 
                             aes(x = RunAvg_Whelks, y = Change_Richness_As)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))


## Plot Acanthinucella --> Evenness (SITE LEVEL)
# Mendocino North
CMEvenPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Mendocino North" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_As)) +
  geom_point() +
  geom_smooth(color = "black", linetype = "dashed", method = "lm", se = F) +
  labs(x = " ", y = " ") +
  ggtitle(expression(bold(Delta~"Evenness"))) + 
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Mendocino South
CMSEvenPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Mendocino South" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_As)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Dana Point
DanaEvenPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_As)) +
  geom_point() +
  geom_smooth(color = "black", method = "lm", se = F) +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Scripps
ScrippsEvenPlot_As <- ggplot(data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_As)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = " ") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-0.05, 15)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))


## Create a figure that is grouped within site (Shannon, Rich, Evenness in a column), with all 4 sites as a square
# Will need to create 4 arranges, 1 per site, 1 to combine 

CMDivAs <- ggarrange(CMShannonPlot_As, CMRichPlot_As, CMEvenPlot_As,
                     nrow = 1, ncol = 3, align = "h")

CMSDivAs <- ggarrange(CMSShannonPlot_As, CMSRichPlot_As, CMSEvenPlot_As,
                      nrow = 1, ncol = 3, align = "h")

DanaDivAs <- ggarrange(DanaShannonPlot_As, DanaRichPlot_As, DanaEvenPlot_As,
                       nrow = 1, ncol = 3, align = "h")

ScrippsDivAs <- ggarrange(ScrippsShannonPlot_As, ScrippsRichPlot_As, ScrippsEvenPlot_As,
                          nrow = 1, ncol = 3, align = "h")


# Combine all four site-diversity groups into main figure
AllDiversity_As <- ggarrange(CMDivAs, CMSDivAs, DanaDivAs, ScrippsDivAs, 
                             nrow = 4, ncol = 1, align = "v",
                             widths = c(0.8, 1, 1), heights = c(1, 1, 1, 1.2))


# Print the final plot
AllDiversity_As

# Clear the workspace
rm(CMShannonPlot_As, CMRichPlot_As, CMEvenPlot_As)
rm(CMSShannonPlot_As, CMSRichPlot_As, CMSEvenPlot_As)
rm(DanaShannonPlot_As, DanaRichPlot_As, DanaEvenPlot_As)
rm(ScrippsShannonPlot_As, ScrippsRichPlot_As, ScrippsEvenPlot_As)
rm(CMDivAs, CMSDivAs, DanaDivAs, ScrippsDivAs, AllDiversity_As)



###########################################
# Mexacanthina
###########################################
# Plot production process is the same for Mexacanthina

## Plot Mexacanthina --> Shannon Diversity (Site Level)

# Dana Point
DanaShannonPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_Mex)) +
  geom_point() +
  labs(x = " ", y = "Dana Point") +
  ggtitle(expression(bold(Delta~"Shannon"))) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Scripps
ScrippsShannonPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_Mex)) +
  geom_point() +
  labs(x = " ", y = "Scripps") +
  #ggtitle("Scripps") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Punta Morro
PMShannonPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Punta Morro" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_Mex)) +
  geom_point() +
  geom_smooth(color = "black", linetype = "dashed", method = "lm", se = F) +
  labs(x = " ", y = "Punta Morro") +
  #ggtitle("Punta Morro") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.y = element_text(size=12, face = "bold", colour = "black"))

# Campo Kennedy
CKShannonPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Campo Kennedy" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Shannon_Mex)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = "Campo Kennedy") +
  #ggtitle("Campo Kennedy") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.50, 1.50, by = 0.50), limits = c(-0.50, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size=12, face = "bold", colour = "black"))


## Plot Mexacanthina --> Richness (SITE LEVEL)
# Dana Point
DanaRichPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_Mex)) +
  geom_point() +
  labs(x = " ", y = " ") +
  ggtitle(expression(bold(Delta~"Richness"))) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Scripps
ScrippsRichPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_Mex)) +
  geom_point() +
  geom_smooth(color = "black", method = "lm", se = F) +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Punta Morro
PMRichPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Punta Morro" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_Mex)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Campo Kennedy
CKRichPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Campo Kennedy" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Richness_Mex)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 4), limits = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))


## Plot Mexacanthina --> Evenness (SITE LEVEL)
# Dana Point
DanaEvenPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Dana Point" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_Mex)) +
  geom_point() +
  labs(x = " ", y = " ") +
  ggtitle(expression(bold(Delta~"Evenness"))) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Scripps
ScrippsEvenPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Scripps" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_Mex)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Punta Morro
PMEvenPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") &  Site == "Punta Morro" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_Mex)) +
  geom_point() +
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank())

# Campo Kennedy
CKEvenPlot_Mex <- ggplot(data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Site == "Campo Kennedy" & Weeks_Deployed == 8), aes(x = RunAvg_Whelks, y = Change_Evenness_Mex)) +
  geom_point() +
  labs(x = expression(atop(Whelk~Density, paste(("Avg. # Ind. × Plot"^-1)))), y = " ") +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(-0.05, 24)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))


## Create a figure that is grouped within site (Shannon, Rich, Evenness in a column), with all 4 sites as a square
# Will need to create 4 arranges, 1 per site, 1 to combine 

DanaDivMex <- ggarrange(DanaShannonPlot_Mex, DanaRichPlot_Mex, DanaEvenPlot_Mex,
                        nrow = 1, ncol = 3, align = "h")

ScrippsDivMex <- ggarrange(ScrippsShannonPlot_Mex, ScrippsRichPlot_Mex, ScrippsEvenPlot_Mex,
                           nrow = 1, ncol = 3, align = "h")

PMDivMex <- ggarrange(PMShannonPlot_Mex, PMRichPlot_Mex, PMEvenPlot_Mex,
                      nrow = 1, ncol = 3, align = "h")

CKDivMex <- ggarrange(CKShannonPlot_Mex, CKRichPlot_Mex, CKEvenPlot_Mex,
                      nrow = 1, ncol = 3, align = "h")

# Combine all four site-diversity groups into main figure
AllDiversity_Mex <- ggarrange(DanaDivMex, ScrippsDivMex, PMDivMex, CKDivMex, 
                              nrow = 4, ncol = 1, 
                              widths = c(0.8, 1, 1), heights = c(1, 1, 1, 1.2))

# Print the figure
AllDiversity_Mex

## Clear the working space
rm(DanaShannonPlot_Mex, DanaRichPlot_Mex, DanaEvenPlot_Mex)
rm(ScrippsShannonPlot_Mex, ScrippsRichPlot_Mex, ScrippsEvenPlot_Mex)
rm(PMShannonPlot_Mex, PMRichPlot_Mex, PMEvenPlot_Mex)
rm(CKShannonPlot_Mex, CKRichPlot_Mex, CKEvenPlot_Mex)
rm(DanaDivMex, ScrippsDivMex, PMDivMex, CKDivMex, AllDiversity_Mex)



################################################################################
## Diversity Impact Effect Size (Figure 4) -------------------------------------
################################################################################
# This will be accomplished in several steps:
# (1) Redefine all models
# (2) Extracting coefficients 
# (3) Calculating 95% Confidence Intervals 
# (4) Recombining into a graphable dataframe
# (5) Graphing per Whelk-Prey pair and then combining

## Recreate the LMs from which the slope coefficients and SEs will be extracted

# Impacts of of Acanthinucella on Diversity Metrics ----------------------------

## Acanthinucella --> Shannon Diversity
CM_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))
CMS_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))
Dana_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_AsShannon <- lm(formula = Change_Shannon_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Acanthinucella --> Richness
CM_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))
CMS_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))
Dana_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_AsRichness <- lm(formula = Change_Richness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Models for Acanthinucella --> Evenness
CM_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino North"))
CMS_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Mendocino South"))
Dana_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_AsEvenness <- lm(formula = Change_Evenness_As ~ RunAvg_Whelks, data = filter(DiversityAs, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))

# Impacts of Mexacanthina on Diversity Metrics ---------------------------------

# Models for Mexacanthina --> Shannon Diversity
Dana_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))
PM_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))
CK_MexShannon <- lm(formula = Change_Shannon_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

# Models for Mexacanthina --> Richness
Dana_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))
PM_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))
CK_MexRichness <- lm(formula = Change_Richness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

# Models for Mexacanthina --> Evenness
Dana_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Dana Point"))
Scripps_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Scripps"))
PM_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Punta Morro"))
CK_MexEvenness <- lm(formula = Change_Evenness_Mex ~ RunAvg_Whelks, data = filter(DiversityMex, Plot_ID %!in% c("No Cage", "Partial") & Weeks_Deployed == 8 & Site == "Campo Kennedy"))

## Extract coefficients from regressions, caluclate CIs, and combine into new dataframe 
regressions_list <- list(Dana_MexShannon, Scripps_MexShannon, PM_MexShannon, CK_MexShannon,
                         Dana_MexRichness, Scripps_MexRichness, PM_MexRichness, CK_MexRichness,
                         Dana_MexEvenness, Scripps_MexEvenness, PM_MexEvenness, CK_MexEvenness,
                         CM_AsShannon, CMS_AsShannon, Dana_AsShannon, Scripps_AsShannon,
                         CM_AsRichness, CMS_AsRichness, Dana_AsRichness, Scripps_AsRichness,
                         CM_AsEvenness, CMS_AsEvenness, Dana_AsEvenness, Scripps_AsEvenness)

# Extract coefficients and confidence intervals using broom's tidy function
Coeffs_list <- lapply(regressions_list, broom::tidy)

# Combine coefficients and confidence intervals into a dataframe
Coeffs_df <- do.call(rbind, Coeffs_list)

# Remove unnecessary rows from Coeffs_df
Coeffs_df <- Coeffs_df[-grep("Intercept", Coeffs_df$term),]

# Add Identifying columns to the Coeffs_df
Coeffs_df$Species <- c("M", "M", "M", "M",
                       "M", "M", "M", "M",
                       "M", "M", "M", "M",
                       "As", "As", "As", "As",
                       "As", "As", "As", "As",
                       "As", "As", "As", "As")

Coeffs_df$Site <- c("Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                    "Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                    "Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                    "Mendocino North", "Mendocino South", "Dana Point", "Scripps",
                    "Mendocino North", "Mendocino South", "Dana Point", "Scripps",
                    "Mendocino North", "Mendocino South", "Dana Point", "Scripps")

Coeffs_df$Metric <- c("Shannon", "Shannon", "Shannon", "Shannon",
                      "Richness", "Richness", "Richness", "Richness",
                      "Evenness", "Evenness", "Evenness", "Evenness",
                      "Shannon", "Shannon", "Shannon", "Shannon",
                      "Richness", "Richness", "Richness", "Richness",
                      "Evenness", "Evenness", "Evenness", "Evenness")


# Calculate 95% CI
CI_list <- lapply(regressions_list, confint)

# Combine confidence intervals into a matrix and convert to a dataframe
CI_df <- do.call(rbind, CI_list)  
  # Creates the matrix

CI_df <- as.data.frame(CI_df)     
  # Converts to list 

CI_df <- CI_df[-grep("X.Intercept.", rownames(CI_df)),]  
  # Removes unnecessary rows

# Rename the CI rows
CI_df <- CI_df %>% 
  rename(CILow = "2.5 %", CIHigh = "97.5 %")

# Add the same identifying rows to the CI_df
CI_df$Species <- c("M", "M", "M", "M",
                   "M", "M", "M", "M",
                   "M", "M", "M", "M",
                   "As", "As", "As", "As",
                   "As", "As", "As", "As",
                   "As", "As", "As", "As")

CI_df$Site <- c("Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                "Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                "Dana Point", "Scripps", "Punta Morro", "Campo Kennedy",
                "Mendocino North", "Mendocino South", "Dana Point", "Scripps",
                "Mendocino North", "Mendocino South", "Dana Point", "Scripps",
                "Mendocino North", "Mendocino South", "Dana Point", "Scripps")

CI_df$Metric <- c("Shannon", "Shannon", "Shannon", "Shannon",
                  "Richness", "Richness", "Richness", "Richness",
                  "Evenness", "Evenness", "Evenness", "Evenness",
                  "Shannon", "Shannon", "Shannon", "Shannon",
                  "Richness", "Richness", "Richness", "Richness",
                  "Evenness", "Evenness", "Evenness", "Evenness")

# Combine the two datasets into one functional set 
Slopesdf <- full_join(Coeffs_df, CI_df, by = c("Species", "Site", "Metric"))

# Remove Mex models to clean workspace
rm(Dana_MexShannon, Scripps_MexShannon, PM_MexShannon, CK_MexShannon)
rm(Dana_MexRichness, Scripps_MexRichness, PM_MexRichness, CK_MexRichness)
rm(Dana_MexEvenness, Scripps_MexEvenness, PM_MexEvenness, CK_MexEvenness)

# Remove As models to clean workspace
rm(CM_AsShannon, CMS_AsShannon, Dana_AsShannon, Scripps_AsShannon)
rm(CM_AsRichness, CMS_AsRichness, Dana_AsRichness, Scripps_AsRichness)
rm(CM_AsEvenness, CMS_AsEvenness, Dana_AsEvenness, Scripps_AsEvenness)

# Remove intermediary lists/dfs
rm(regressions_list, Coeffs_list, Coeffs_df, CI_list, CI_df)

# Add new columns for graphing purposes
Slopesdf$Range <- c("Expanded", "Expanded", "Historic", "Historic", 
                    "Expanded", "Expanded", "Historic", "Historic", 
                    "Expanded", "Expanded", "Historic", "Historic",
                    "Expanded", "Expanded", "Historic", "Historic",
                    "Expanded", "Expanded", "Historic", "Historic",
                    "Expanded", "Expanded", "Historic", "Historic")

# Define factor order for plotting purposes 
Slopesdf$Site <- factor(Slopesdf$Site, levels = c("Campo Kennedy", "Punta Morro", "Scripps", "Dana Point", "Mendocino South", "Mendocino North"))

### Create the figures ---------------------------------------------------------

## First create  all the Mexacanthina slope figs 

MexShannon_Est <- ggplot(data = filter(Slopesdf, Species == "M" & Metric == "Shannon"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), linewidth = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = "", y = expression(atop("Effect Size", (Delta~"Shannon * Whelk"^-1)))) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("d.") + 
  scale_y_continuous(breaks = seq(-0.20, 0.20, by = 0.10), limits = c(-0.20, 0.20)) +  
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12))

MexRichness_Est <- ggplot(data = filter(Slopesdf, Species == "M" & Metric == "Richness"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), size = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = "", y = expression(atop("Effect Size", (Delta~"Richness * Whelk"^-1)))) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("e.") + 
  scale_y_continuous(breaks = seq(-1.0, 1.5, by =0.5), limits = c(-1.0, 1.5)) + 
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank())

MexEvenness_Est <- 
  ggplot(data = filter(Slopesdf, Species == "M" & Metric == "Evenness"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), size = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = "", y = expression(atop("Effect Size", (Delta~"Evenness * Whelk"^-1)))) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("f.") + 
  scale_y_continuous(breaks = seq(-0.8, 0.8, by = 0.04), limits = c(-0.09, 0.09)) + 
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank())

## First create  all the Acanthinucella slope figs 

AsShannon_Est <- ggplot(data = filter(Slopesdf, Species == "As" & Metric == "Shannon"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), size = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = " ", y = " ") +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("a.") + 
  scale_y_continuous(breaks = seq(-0.20, 0.20, by = 0.10), limits = c(-0.20, 0.20)) +  
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12))

AsRichness_Est <- ggplot(data = filter(Slopesdf, Species == "As" & Metric == "Richness"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), size = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = "", y = " ") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("b.") + 
  scale_y_continuous(breaks = seq(-1.0, 1.5, by =0.5), limits = c(-1.0, 1.5)) +
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank())


AsEvenness_Est <-
  ggplot(data = filter(Slopesdf, Species == "As" & Metric == "Evenness"), aes(x = Site, y = estimate)) +
  geom_point(aes(shape = Range, color = Range), size = 1.5, show.legend = F) + 
  geom_errorbar(aes(ymin = CILow, ymax = CIHigh, color = Range), size = 0.75, width = 0, show.legend = F) +
  scale_color_manual(values = c("Historic" = "#000000", "Expanded" = "#999999")) +
  labs(x = "", y = " ") +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("c.") + 
  scale_y_continuous(breaks = seq(-0.8, 0.8, by = 0.04), limits = c(-0.09, 0.09)) + 
  coord_flip() +
  theme_classic2() +
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank())


# Now combine together
DiversityEstimatesAll <- ggarrange(AsShannon_Est, AsRichness_Est, AsEvenness_Est,
                                   MexShannon_Est, MexRichness_Est, MexEvenness_Est,
                                   nrow = 2, ncol = 3,
                                   align = "hv")


# Print the Figure
DiversityEstimatesAll
# When saving, can use ratio of 900 x 450


# Remove clutter
rm(MexShannon_Est, MexRichness_Est, MexEvenness_Est)
rm(AsShannon_Est, AsRichness_Est, AsEvenness_Est)
rm(Slopesdf, DiversityEstimatesAll)

