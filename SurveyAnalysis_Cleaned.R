# Created by [name redacted for peer review]
# Last edited: 01/08/2025

## Load Packages ---------------------------------------------------------------

library(tidyverse)    #loads tidyr, readr, and dplyr
library(vegan)        #community analysis functions
library(lme4)         #for running linear mixed model statistics
library(lmerTest)     #for determining the statistical significance
library(performance)  #can use to check model assumptions with check_model()
library(viridis)      #for colorblind-friendly plots
library(ggpubr)       #for visualizing qq plots, enabling ggarrange


# Check Working Directory ------------------------------------------------------
getwd()

# Clear the Workspace ----------------------------------------------------------
rm(list=ls())

# Code basic functions that I like to use --------------------------------------
# Code the "not in" symbol 
"%!in%" <- Negate("%in%")

# Read in all of the cleaned data ----------------------------------------------

# Start by loading all sheets
CommunityData <- read.csv("SurveyCommunityData.csv")
TransectData <- read.csv("SurveyTransectData.csv")
WhelkLengthData <- read.csv("SurveyWhelkLengthData.csv")

# Edit the CommunityData sheet  ------------------------------------------------

## Fill all NA's in the CommunityData frame as zeros
CommunityData[is.na(CommunityData)] = 0

## Remove rows with no community data (i.e., data may have been lost from these locations)
CommunityData <- CommunityData[-c(80),] 
  # VB T12, Q2.3, quad washed away
  # Retained in the dataset to show quad was lost, not forgotten

## Remove Rancho Marino from community data
CommunityData <- filter(CommunityData, Site != "Rancho Marino")
  # Community data is not from field surveys due to time constraints
  # Still valuable but lacks the biological "layering" available in field surveys 

## Add the whelk counts (# individuals) to the CommunityData sheet 
  # This will give two different estimates of # individuals
QuadWhelkCount_DF <- subset(WhelkLengthData, In_Timed_Count == "N" & Species %in% c("As", "M")) %>%
  # Subsets data not in TC and confines to only the focal whelks
  group_by(Region, Site, Survey_Date, Vert_Transect_Dist_m, Quad_Dist_m, Species) %>% 
  # Group together
  count()
  # Count the number of whelks found outside of timed counts

# Pivot QuadWhelkCount wider and rename the column
QuadWhelkCount_DF <- QuadWhelkCount_DF %>%
  # pull in the data
  pivot_wider(names_from = Species, values_from = n) %>% 
  # Pivot wider
  rename(Mex_Count = M, As_Count = As)
  # Rename columns 

# Ensure QuadWhelkCount_DF is a dataframe and define column types
QuadWhelkCount_DF <- as.data.frame(QuadWhelkCount_DF)
  # Convert to dataframe 

QuadWhelkCount_DF$Vert_Transect_Dist_m <- as.numeric(QuadWhelkCount_DF$Vert_Transect_Dist_m)   
  # Change Vertical Distance to numeric

QuadWhelkCount_DF$Quad_Dist_m <- as.numeric(QuadWhelkCount_DF$Quad_Dist_m)
  # Change Whelk Count (# ind) to numeric

# Join this with the initial CommunityData sheet 
  # This joins the whelk counts from quads with whelks to the quads without whelks
CommunityData <- left_join(x = CommunityData, y = QuadWhelkCount_DF, by = c("Region", "Site", "Survey_Date", "Vert_Transect_Dist_m", "Quad_Dist_m"))

# Fill NAs in the new count columns with zeros
CommunityData$As_Count[is.na(CommunityData$As_Count)] <- 0
CommunityData$Mex_Count[is.na(CommunityData$Mex_Count)] <- 0

## Assign "historic" or "expanded" to each site by the species
# Can do this in one long pipe 
CommunityData <- CommunityData %>% 
  # Pull in the data
  mutate(Mex_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", 
                                            "Punta Morro", "San Miguel", "Saldamando") ~ "Historic", 
                                Site %in% c("Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", 
                                            "Victoria Beach", "Heisler Park", "Venice Breakwater", "Santa Monica") ~ "Expanded",
                                Site %in% c("Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", 
                                            "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", 
                                            "Mendocino North") ~ "NA")) %>% 
    # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel",
                                           "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", 
                                           "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove",
                                           "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                                   Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                                   Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ "NA"))
    # This chunk assigns "Status" for Acanthinucella at all sites
    # Sites where species is not observed are counted as NA
 
## Calculate a "present/absent" column for both focal species to be included later models
CommunityData$MexPresence <- as.numeric(CommunityData$Dark_Unicorn_Whelk_M._lugubris != 0)
CommunityData$AcanthPresence <- as.numeric(CommunityData$Angled_Unicorn_Whelk_A._spirata != 0)

## Extract longitude from TransectData and assign to each row
# Start by dropping Coal Oil from the TransectData (no community survey data)
TransectData <- filter(TransectData, Site != "Coal Oil")

# Make the GPS_Lat column numeric and the date column as date
TransectData$GPS_Lat <- as.numeric(TransectData$GPS_Lat)
TransectData$Survey_Date <- as.Date(TransectData$Survey_Date, tryFormats = "%m/%d/%Y")

# Calculate the latitude of the survey site (as the average of all transects) 
Avg_Lat <- TransectData %>% 
  # Pull in the data
  group_by(Region, Site, Survey_Date) %>%    
  # Group by site and survey date
  summarize(Avg_Lat = mean(GPS_Lat))
  # Calculate mean latitude (within each survey date)

# Set Survey Date to "date" and merge the latitude data with the CommunityData
CommunityData$Survey_Date <- as.Date(CommunityData$Survey_Date, tryFormats = "%m/%d/%Y")
CommunityData <- left_join(CommunityData, Avg_Lat, by = c("Region", "Site", "Survey_Date"))

# Remove the TransectData, Avg_Latitude Dfs after the merge
rm(TransectData)


################################################################################
## Whelk Density Models and Relative Abundance by Site -------------------------
################################################################################

# Determine change in density with latitude ------------------------------------

## First, we will calculate the density of whelks observed in all quadrats during surveys 
  # To calculate estimated density, we need: quad size, total area surveyed, and average # per quad
  # Quad size = constant = 0.0625m^2

# Start by extracting only whelk abundance data for both focal species 
WhelkDensity_Df <- CommunityData %>% 
  select(Region:Season, Vert_Transect_Dist_m:Quad_TH_m, As_Count:Mex_Count) 
  # Select the necessary data from the CommunityData column
  # Make sure that the last two numbers are the count data

# Calculate density (per m^2) in two new columns
WhelkDensity_Df <- WhelkDensity_Df %>%
  # Pull in the data
  mutate(As_Density = (As_Count*16), 
         M_Density = (Mex_Count*16))
  # Quads are 1/16 of a square meter, so multiply by 16 to get per m2

# Then pivot longer
WhelkDensity_Df <- WhelkDensity_Df %>%
  # Pull in the data
  pivot_longer(cols = c(As_Density, M_Density), names_to = "Species", values_to = "Density")
  # Conduct the pivot

# Rename WhelkDensity_DF's species column to make sense 
WhelkDensity_Df$Species[WhelkDensity_Df$Species == "As_Density"] <- "As"
WhelkDensity_Df$Species[WhelkDensity_Df$Species == "M_Density"] <- "M"

# Assign whelk status to WhelkDensity_DF
WhelkDensity_Df <- WhelkDensity_Df %>% 
  #Pull in the data
  mutate(Mex_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", 
                                          "Punta Morro", "San Miguel", "Saldamando") ~ "Historic", 
                              Site %in% c("Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", 
                                          "Victoria Beach", "Heisler Park", "Venice Breakwater", "Santa Monica") ~ "Expanded",
                              Site %in% c("Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", 
                                          "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", 
                                          "Mendocino North") ~ "NA")) %>% 
  # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel",
                                           "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", 
                                           "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove",
                                           "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                               Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                               Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ "NA"))
  # This chunk assigns "Status" for Acanthinucella at all sites
  # Sites where species is not observed are counted as NA

# Calculate the average latitude of each site to append to this dataframe 
Avg_Lat <- Avg_Lat %>% 
  # Pull in the data
  group_by(Region, Site) %>% 
  # Group
  summarise(Avg_Lat = mean(Avg_Lat))
  # Calculate the site latitude as the global average (all transects across all seasons)

# Then join this with the WhelkDensity_DF for the final dataset prior to modelling
WhelkDensity_Df <- left_join(WhelkDensity_Df, Avg_Lat, by = c("Region", "Site"))



# Define and run LMMs assessing whelk density ----------------------------------

# Build model to test how Acanthinucella's density varies with latitude
As_Density_Lat_Mod <- lmer(Density ~ Avg_Lat + Quad_TH_m + As_Status + (1|Vert_Transect_Dist_m), 
                           data = filter(WhelkDensity_Df, Species == "As" & As_Status != "NA"))
  # Site not included as random effect because it is singular with Avg_Lat
  # NA sites are not included because I don't want to try to predict whelk counts at sites where we 
  # never found whelks 

# Check model assumptions 
check_model(As_Density_Lat_Mod)

# Test significance 
summary(As_Density_Lat_Mod)


# Build model to test how Acanthinucella's density varies with latitude 
Mex_Density_Lat_Mod <- lmer(Density ~ Avg_Lat + Quad_TH_m + Mex_Status + (1|Vert_Transect_Dist_m), 
                            data = filter(WhelkDensity_Df, Species == "M" & Mex_Status != "NA"))
  # Site not included as random effect because it is singular with Avg_Lat
  # NA sites are not included because I don't want to try to predict whelk counts at sites where we 
  # never found whelks 

# Check model assumptions
check_model(Mex_Density_Lat_Mod)

# Test the significance
summary(Mex_Density_Lat_Mod)

# Clean the Workspace
rm(As_Density_Lat_Mod, Mex_Density_Lat_Mod)

# Graph relative abundance with latitude ---------------------------------------

# To estimate relative abundance, calculate the number of whelks found in timed counts
Whelks_Timed_Count <- subset(WhelkLengthData, In_Timed_Count == "Y" & Species %in% c("As", "M")) %>%   
  # Pull in data only from only data from timed counts 
  group_by(Region, Species, Site, Survey_Date) %>%     
  # Grouping by site and date
  count() %>%                                          
  # Count all records that occur in new column "n"
  rename(Num_Whelks = n)
  # Rename column "n" to "Num_Whelks"

# Drop all Coal Oil data (no surveys conducted)
Whelks_Timed_Count <- subset(Whelks_Timed_Count, Site != "Coal Oil")

# Then average the # whelks found across sampling seasons 
Whelks_Timed_Count <- Whelks_Timed_Count %>%
  # Pull in Time count dataframe
  group_by(Site, Region, Species) %>% 
  # Group
  summarise(Avg_Mex = mean(Num_Whelks), SE_Mex = (sd(Num_Whelks)/sqrt(n())),
            Avg_Acanth = mean(Num_Whelks), SE_Acanth = (sd(Num_Whelks)/sqrt(n()))) %>% 
  # Calculate the mean number of each whelk observed at each site across all surveys
  mutate(Upper_Mex = (Avg_Mex+SE_Mex), Lower_Mex = (Avg_Mex-SE_Mex),
         Upper_Acanth = (Avg_Acanth+SE_Acanth), Lower_Acanth = (Avg_Acanth-SE_Acanth))
  # Calculate the SE, upper SE, and lower SE values 

## To get these into a graphable format, we will need to split, rename, and then recombine
# Start with the Acanthinucella data
As_Num <- subset(Whelks_Timed_Count, Species == "As") %>% 
  # Pull in As timed count data
  select(Site, Species, Avg_Acanth, SE_Acanth, Upper_Acanth, Lower_Acanth) %>% 
  # Determine the order
  rename(Num_Whelks = Avg_Acanth, SE = SE_Acanth, Upper = Upper_Acanth, Lower = Lower_Acanth)
  # Rename columns 

# Do the same with Mexacanthina data
Mex_Num <- subset(Whelks_Timed_Count, Species == "M") %>% 
  select(Site, Species, Avg_Mex, SE_Mex, Upper_Mex, Lower_Mex) %>% 
  rename(Num_Whelks = Avg_Mex, SE = SE_Mex, Upper = Upper_Mex, Lower = Lower_Mex)

# Recombine both datasets into one longer dataframe that will be nice for graphing
Whelks_Timed_Count <- rbind(As_Num, Mex_Num)

# Again, assign "historic" or "expanded" to each site by the species
Whelks_Timed_Count <- Whelks_Timed_Count %>% 
  # Pull in the data
  mutate(Mex_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", 
                                            "Punta Morro", "San Miguel", "Saldamando") ~ "Historic", 
                                Site %in% c("Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", 
                                            "Victoria Beach", "Heisler Park", "Venice Breakwater", "Santa Monica") ~ "Expanded",
                                Site %in% c("Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", 
                                            "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", 
                                            "Mendocino North") ~ "NA")) %>% 
  # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel",
                                           "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", 
                                           "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove",
                                           "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                               Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                               Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ "NA"))
  # This chunk assigns "Status" for Acanthinucella at all sites
  # Sites where species is not observed are counted as NA


# Calculate the average latitude of each site to append to this dataframe 
Avg_Lat <- Avg_Lat %>% 
  # Pull in the data
  group_by(Region, Site) %>% 
  # Group
  summarise(Avg_Lat = mean(Avg_Lat))
  # Calculate the site latitude as the global average (all transects across all seasons)

# Then join to create the final, graphable dataset
Whelks_Timed_Count <- left_join(Whelks_Timed_Count, Avg_Lat, by = c("Region", "Site"))

# Graph the densities and relative abundances ----------------------------------

# Create the graph for Acanthinucella 
AsRelAbundPlot <- 
  ggplot(data = filter(Whelks_Timed_Count, Site != "Rancho Marino" & Species == "As"), aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Venice Breakwater", "Santa Monica", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North"), exclude = "NA"), y = Num_Whelks)) +
    geom_col(aes(fill = As_Status, color = As_Status), width = 0.8) +
    annotate("text", y = 40, x = c(20, 19), label = c(".", "."), size = 12) +
    coord_flip() +
    scale_color_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000"), guide = F) +
    scale_fill_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000")) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge2()) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 105)) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Site", y = bquote(atop(Relative~Abundance, ("# Obs. in 1 Hr Timed Count"))), fill = "") +
    ggtitle((expression(italic("Acanthinucella")))) +
    theme_classic2() +
    theme(legend.position = c(0.85, 0.15))

MexRelAbundPlot <- 
  ggplot(data = filter(Whelks_Timed_Count, Site != "Rancho Marino" & Species == "M"), aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Venice Breakwater", "Santa Monica", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North"), exclude = "NA"), y = Num_Whelks)) +
  geom_col(aes(fill = Mex_Status, color = Mex_Status), width = 0.8) + 
  annotate("text", y = 60, x = c(20, 19), label = c(".", "."), size = 12) +
    # This adds periods to Venice Breakwater and Santa Monica
  coord_flip() +
  scale_color_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000"), guide = F) +
  scale_fill_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000"), guide = F) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge2()) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 105)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Site", y = bquote(atop(Relative~Abundance, ("# Obs. in 1 Hr Timed Count")))) +
  ggtitle((expression(italic("Mexacanthina")))) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.80))

# Combine the single species relative abundance figures into a single plot
Graph_Timed_Count_Abund <- ggarrange(AsRelAbundPlot, MexRelAbundPlot,
                                            labels = c("a.", "b."),
                                            nrow = 1)

# Print the new figure (Figure 2 in Manuscript)
Graph_Timed_Count_Abund
    # Print as 850 x 450

# Remove the graphs and dataframes from the working space
rm(Avg_Lat, As_Density, Mex_Density, As_Num, Mex_Num, WhelkDensity_Df, QuadWhelkCount_DF)
rm(Graph_Site_Density, Graph_Timed_Count_Abund, Graph_Combined_Whelk_SiteAbund)
rm(AsRelAbundPlot, MexRelAbundPlot)
rm(Site_Whelk_Density, SiteQuadCount_DF, SiteWhelkCount_DF, Whelks_Timed_Count)


################################################################################
## Testing and Graphing Impacts on Prey (Barnacles and Mussels) ----------------
################################################################################
# We will focus on Acorn barnacles (Chthamalus/Balanus) and California mussels (Mytilus californianus)
# as preferred prey, a-la Wallingford and Sorte 2022 (and references therein). We will look at the 
# proportion of dead barnacles near whelks as a measure of mortality

# Start by selecting only relevant prey and whelk columns 
CommunityData_Prey <- CommunityData %>% 
  # Pull in the data
  select(Region:Season, Vert_Transect_Dist_m:Quad_Dist_m, Quad_TH_m, Dead_Barnacles:Acorn_Barnacle_Balanus_Chthamalus, 
         California_Mussel_M._californianus, Angled_Unicorn_Whelk_A._spirata:Dark_Unicorn_Whelk_M._lugubris, 
         As_Count:Avg_Lat)
  # Select the columns


# Now calculate the proportion of living barnacles at each site
CommunityData_Prey <- CommunityData_Prey %>% 
  # Pull in the data
  mutate(Prop_LiveBarn = (Acorn_Barnacle_Balanus_Chthamalus/(Acorn_Barnacle_Balanus_Chthamalus + Dead_Barnacles)))

# Replace NaN with NA
CommunityData_Prey$Prop_LiveBarn[is.nan(CommunityData_Prey$Prop_LiveBarn)] <- NA
  # NaNs result from quadrats with 0 Balanus/Chthamalus
  # Converting to NA will omit them from analysis and figures


# Build and Run LMMs Predator-Prey Abundance Relationships ---------------------

# There will be a total of 4 main tests run: 
  # Acanthinucella vs (1) proportion living barnacles and (2) California mussels
  # Mexacanthina vs (3) proportion living barnacles and (4) California mussels
  # For each test, we will define two models (interactive and reduced) and pick one to interpret using AIC

# (1) Acanthinucella v  Prop. Living Barnacles (2 Models, pick 1 with AIC)
LiveBarn_AsGlobal <- lmer(Prop_LiveBarn ~ As_Count*As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, As_Status != "NA"))

LiveBarn_AsGlobal2 <- lmer(Prop_LiveBarn ~ As_Count + As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, As_Status != "NA"))

# Check AICs
AIC(LiveBarn_AsGlobal, LiveBarn_AsGlobal2)
  # Second (non-interactive) model fits better

# Check model assumptions
check_model(LiveBarn_AsGlobal2)

# Check significance
As_LiveBarn_Mod <- summary(LiveBarn_AsGlobal2)


# (2) Acanthinucella v Mussels (2 Models, pick 1 with AIC)
Mussels_AsGlobal <- lmer(California_Mussel_M._californianus ~ As_Count*As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, As_Status != "NA"))

Mussels_AsGlobal2 <- lmer(California_Mussel_M._californianus ~ As_Count + As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, As_Status != "NA"))

# Check AICs
AIC(Mussels_AsGlobal, Mussels_AsGlobal2)
  # First (interactive) model fits better

# Check model assumptions
check_model(Mussels_AsGlobal)

# Check significance
As_Mussel_Mod <- summary(Mussels_AsGlobal)


# (3) Mexacanthina vs Living Barnacles (2 Models, pick 1 with AIC)
LiveBarn_MexGlobal <- lmer(Prop_LiveBarn ~  Mex_Count*Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, Mex_Status != "NA"))

LiveBarn_MexGlobal2 <- lmer(Prop_LiveBarn ~  Mex_Count + Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, Mex_Status != "NA"))

# Check AICs
AIC(LiveBarn_MexGlobal, LiveBarn_MexGlobal2)
  # Second (non-interactive) model fits better

# Check model assumptions
check_model(LiveBarn_MexGlobal2)

# Check significances
M_LiveBarn_Mod <- summary(LiveBarn_MexGlobal2)


# (4) Mexacanthina v Mussels (2 Models, pick 1 with AIC) 
Mussels_MexGlobal <- lmer(California_Mussel_M._californianus ~  Mex_Count*Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, Mex_Status != "NA"))

Mussels_MexGlobal2 <- lmer(California_Mussel_M._californianus ~  Mex_Count + Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData_Prey, Mex_Status != "NA"))

# Check AICs
AIC(Mussels_MexGlobal, Mussels_MexGlobal2)
  # First (interactive) model fits better

# Check model assumptions
check_model(Mussels_MexGlobal)

# Check significance
M_Mussels_Mod <- summary(Mussels_MexGlobal)


# Print each model below to see significance of observational Pred-Prey relationship tests
As_LiveBarn_Mod
As_Mussel_Mod
M_LiveBarn_Mod
M_Mussels_Mod


## Clean workspace 

# Removing Acanthinucella models 
rm(LiveBarn_AsGlobal, LiveBarn_AsGlobal2)
rm(Mussels_AsGlobal, Mussels_AsGlobal2)

# Removing Mexacanthina models 
rm(LiveBarn_MexGlobal, LiveBarn_MexGlobal2)
rm(Mussels_MexGlobal, Mussels_MexGlobal2)

# Remove the remaining items from this section 
rm(As_LiveBarn_Mod, As_Mussel_Mod, M_LiveBarn_Mod, M_Mussels_Mod)



# Graph Whelk Abundances vs Prey Abundances (Supp Fig 1) -----------------------
# To denote significance, geom_smooths will be added to plots significant relationships from above

# Acanthinucella abundance vs proportion living barnacles
Plot_As_LiveBarn <- 
ggplot(data = filter(CommunityData_Prey, As_Status != "NA" & Prop_LiveBarn != "NA"), 
       aes(x = As_Count, y = Prop_LiveBarn)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = "Proportion Living Barnacles") +
  ggtitle("a.") +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(0, 1.0)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(-0.5, 12)) + 
  theme_classic2(base_size = 10) 

# Acanthinucella abundance vs California mussels 
Plot_As_Mussels <- 
ggplot(data = filter(CommunityData_Prey, As_Status != "NA"), 
       aes(x = As_Count, y = California_Mussel_M._californianus)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(-0.5, 12)) + 
  labs(x = expression(italic("Acanthinucella")~"Abundance (# ind.)"), 
       y = "California Mussel Abundance\n(% cover)") +
  ggtitle("c.") +
  theme_classic2(base_size = 10) +
  theme(legend.position = c(0.80,0.85))

# Combine the prey plots for Acanthinucella
Graph_Combined_Prey_As <- ggarrange(Plot_As_LiveBarn, Plot_As_Mussels,
                                        ncol = 1, nrow = 2)

# Print the figure 
Graph_Combined_Prey_As


# Mexacanthina abundance vs the proportion of living barnacles 
Plot_Mex_LiveBarn <- 
  ggplot(data = subset(CommunityData_Prey, Mex_Status != "NA" & Prop_LiveBarn != "NA"), 
         aes(x = Mex_Count, y = Prop_LiveBarn)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +
  geom_smooth(color = "black", method = lm, se = FALSE, show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(0, 1.0)) +
  scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-5, 125)) + 
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("b.") +
  theme_classic2(base_size = 10) 

# Mexacanthina abundance vs California mussels 
Plot_Mex_Mussels <- 
  ggplot(data = subset(CommunityData_Prey, Mex_Status != "NA"), 
         aes(x = Mex_Count, y = California_Mussel_M._californianus)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = F) +
  geom_smooth(aes(color = Mex_Status), method = lm, show.legend = F, se = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(-0.5, 100)) +
  scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-5, 125)) + 
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(italic("Mexacanthina")~"Abundance (# ind.)"), y = " ", color = "Status") +
  ggtitle("d.") +
  theme_classic2(base_size = 10) +
  theme(legend.position = c(0.80,0.85))

# Combine the prey plots for Mexacanthina
Graph_Combined_Prey_Mex <- ggarrange(Plot_Mex_LiveBarn, Plot_Mex_Mussels,
                                     ncol = 1, nrow = 2, align = "hv")

# Print the figure
Graph_Combined_Prey_Mex

# Combine Mexacanthina and Acanthinucella plots
PreyPlotsAll <- ggarrange(Graph_Combined_Prey_As, Graph_Combined_Prey_Mex,
                          ncol = 2, nrow = 1, align = 'hv')

# Print final figure (export as 600 x 500)
# This will be Fig Supplemental Figure 1
PreyPlotsAll


# Clean workspace
rm(Plot_As_LiveBarn, Plot_As_Mussels)
rm(Plot_Mex_LiveBarn, Plot_Mex_Mussels)
rm(Graph_Combined_Prey_As, Graph_Combined_Prey_Mex, PreyPlotsAll)
rm(CommunityData_Prey)




################################################################################
## Whelk Abundances and Community Diversity  -----------------------------------
################################################################################

# Start by calculating community diversity metrics (for each species)
  # Diversity metrics do not include the focal whelk itself 

# Calculate SW diversity for each focal species (i.e., the site diversity without that species included)
CommunityData$SW_Diversity_Mex <- diversity(x = CommunityData[c(18:93, 95:129)], index = "shannon")
CommunityData$SW_Diversity_As <- diversity(x = CommunityData[c(18:92, 94:129)], index = "shannon")

# Calculate Richness for each focal species (i.e., the site diversity without that species included)
CommunityData$Richness_Mex <- specnumber(x = CommunityData[c(18:93, 95:129)])
CommunityData$Richness_As <- specnumber(x = CommunityData[c(18:92, 94:129)])

# Calculate Evenness for each focal species (i.e., the site diversity without that species included)
CommunityData <- CommunityData %>% 
  # Pull in the data
  mutate(Evenness_Mex = SW_Diversity_Mex/log(Richness_Mex),
         Evenness_As = SW_Diversity_As/log(Richness_As))
  # Calculate evvenness for each species 


# Build and Run LMMs Assessing Predator-Diversity Relationships ----------------
# There will be 6 total analyses run: 
  # Acanthinucella vs (1) Shannon Diversity, (2) Richness, (3) Evenness
  # Mexacanthina vs (4) Shannon Diversity, (5) Richness, (6) Evenness

# Acanthinucella - Shannon Diversity (2 Models, pick 1 with AIC)
Diversity_AsGlobal <- lmer(SW_Diversity_As ~ As_Count*As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

Diversity_AsGlobal2 <- lmer(SW_Diversity_As ~ As_Count + As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check AICs
AIC(Diversity_AsGlobal, Diversity_AsGlobal2)
  # Second (non-interactive) model fits better

# Check model assumptions
check_model(Diversity_AsGlobal2)

# Check significance
As_Diversity_Mod <- summary(Diversity_AsGlobal2)


# Acanthinucella - Richness (2 Models, pick 1 with AIC)
Richness_AsGlobal <- lmer(Richness_As ~ As_Count*As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

Richness_AsGlobal2 <- lmer(Richness_As ~ As_Count + As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check AICs
AIC(Richness_AsGlobal, Richness_AsGlobal2)
  # First (interactive) model fits better

# Check model assumptions
check_model(Richness_AsGlobal2)

# Check significance
As_Richness_Mod <- summary(Richness_AsGlobal)


# Acanthinucella - Evenness (2 Models, pick 1 with AIC)
Evenness_AsGlobal <- lmer(Evenness_As ~ As_Count*As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

Evenness_AsGlobal2 <- lmer(Evenness_As ~ As_Count + As_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check AICs
AIC(Evenness_AsGlobal, Evenness_AsGlobal2)
# Second (non-interactive) model fits better

# Check model assumptions
check_model(Evenness_AsGlobal2)

# Check significance
As_Evenness_Mod <- summary(Evenness_AsGlobal2)


# Mexacanthina - Shannon Diversity (2 Models, pick 1 with AIC)
Diversity_MexGlobal <- lmer(SW_Diversity_Mex ~ Mex_Count*Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

Diversity_MexGlobal2 <- lmer(SW_Diversity_Mex ~ Mex_Count + Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check AICs
AIC(Diversity_MexGlobal, Diversity_MexGlobal2)
  # Effectively identical, so both models are equally likely
  # Will retain interaction because it answers the overall question

# Check model assumptions
check_model(Diversity_MexGlobal)

# Check significance
Mex_Diversity_Mod <- summary(Diversity_MexGlobal)

# Mexacanthina - Richness (2 Models, pick 1 with AIC)
Richness_MexGlobal <- lmer(Richness_Mex ~ Mex_Count*Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

Richness_MexGlobal2 <- lmer(Richness_Mex ~ Mex_Count + Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check AICs
AIC(Richness_MexGlobal, Richness_MexGlobal2)
  # First (interactive) model clearly better

# Check model assumptions
check_model(Richness_MexGlobal)

# Check significance
Mex_Richness_Mod <- summary(Richness_MexGlobal)

## Mexacanthina v Evenness (2 Models, pick 1 with AIC)
Evenness_MexGlobal <- lmer(Evenness_Mex ~ Mex_Count*Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

Evenness_MexGlobal2 <- lmer(Evenness_Mex ~ Mex_Count + Mex_Status + Quad_TH_m + (1|Site/Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check AICs
AIC(Evenness_MexGlobal, Evenness_MexGlobal2)
  # Second (non-interactive) model clearly better

# Check model assumptions
check_model(Evenness_MexGlobal2)

# Check significance
Mex_Evenness_Mod <-summary(Evenness_MexGlobal2)


# Print Out Model Summaries
As_Diversity_Mod
As_Richness_Mod
As_Evenness_Mod
Mex_Diversity_Mod
Mex_Richness_Mod
Mex_Evenness_Mod


## Clean the workspace

# Clean all models
rm(Diversity_AsGlobal, Diversity_AsGlobal2, 
   Richness_AsGlobal, Richness_AsGlobal2, 
   Evenness_AsGlobal, Evenness_AsGlobal2,
   Diversity_MexGlobal, Diversity_MexGlobal2, 
   Richness_MexGlobal, Richness_MexGlobal2, 
   Evenness_MexGlobal, Evenness_MexGlobal2)

# Clean all model summaries 
rm(As_Diversity_Mod,
   As_Richness_Mod,
   As_Evenness_Mod,
   Mex_Diversity_Mod,
   Mex_Richness_Mod,
   Mex_Evenness_Mod)

# Clean all figures
rm(Graph_DiversityvMex, Graph_DiversityvAcanth, 
   Graph_RichnessvMex, Graph_RichnessvAcanth,
   Graph_EvennessvMex, Graph_EvennessvAcanth)

# Clean 
 
# Graph Whelk Abundances vs Diversity Metrics (Supp Fig 2) ---------------------
# To denote significance, geom_smooths will be added to plots significant relationships from above

# Plots Acanthinucella abundance v Shannon diversity 
Graph_DiversityvAcanth <- 
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = Angled_Unicorn_Whelk_A._spirata, y = SW_Diversity_As)) +
  geom_point(aes(color = As_Status), position = "jitter",  show.legend = F) +
  geom_smooth(color = "black", method = lm, fullrange = T, se = F, show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = "", y = "Shannon", color = "Status") + 
  ggtitle("a.") +
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(-0.5, 14)) + 
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(-0.05, 2.5)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v Shannon diversity 
Graph_DiversityvMex <- 
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Count, y = SW_Diversity_Mex)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +
  geom_smooth(aes(color = Mex_Status), method = lm, fullrange = T, se = F, show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("b.") +
  scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-1, 125)) + 
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(-0.05, 2.5)) +
  theme_classic2()
    

# Plots Acanthinucella abundance v richness
Graph_RichnessvAcanth <-
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = As_Count, y = Richness_As)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  #geom_smooth(method = lm) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = "Richness", color = "Status") +
  ggtitle("c.") +
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(-0.5, 14)) + 
  scale_y_continuous(breaks = seq(0, 20, by = 4), limits = c(0, 20)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v richness
Graph_RichnessvMex <-
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Count, y = Richness_Mex)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +
  geom_smooth(aes(color = Mex_Status), method = lm, fullrange = T, se = F, show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("d.") +
  scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-1, 125)) + 
  scale_y_continuous(breaks = seq(0, 20, by = 4), limits = c(-0.05, 20)) +
  theme_classic2()


# Plots Acanthinucella abundance v evenness
Graph_EvennessvAcanth <-
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = As_Count, y = Evenness_As)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  geom_smooth(color = "black", method = lm, fullrange = T, se = F, show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(atop(italic("Acanthinucella")~"Abundance",~"(# ind.)")), y = "Evenness", color = "Status") +
  ggtitle("e.") +
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(-0.5, 14)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(-0.05, 1.05)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v evenness
Graph_EvennessvMex <-
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Count, y = Evenness_Mex)) +
  geom_point(aes(color = Mex_Status), show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(atop(italic("Mexacanthina")~"Abundance",~"(# ind.)")), y = " ") +
  ggtitle("f.") +
  scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-1, 125)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(-0.05, 1.05)) +
  theme_classic2()

# Create the combined figure with all diverity metrics
Graph_Diversity_Full <- ggarrange(Graph_DiversityvAcanth, Graph_DiversityvMex,
                                  Graph_RichnessvAcanth, Graph_RichnessvMex,
                                  Graph_EvennessvAcanth, Graph_EvennessvMex,
                                  ncol = 2, nrow = 3, align = "hv")

# Notes on dropped values - 
  # Most are NAN that result from "communities" with only one sp observed.
  # This occurs because the formula equates to 0/log(1), which is undefined.
  # Some are not visualized because there is not enough space to jitter 
  # or cleanly draw on top of other points. 


# Print the figure (Supp Fig. 2)
Graph_Diversity_Full

# Remove figures
rm(Graph_DiversityvMex, Graph_DiversityvAcanth, Graph_Combined_WhelksvDiversity)
rm(Graph_RichnessvMex, Graph_RichnessvAcanth, Graph_Combined_WhelksvRichness)
rm(Graph_EvennessvMex, Graph_EvennessvAcanth, Graph_Combined_WhelksvEvenness)



################################################################################
# Additional Supplemental Figures and Tables (Based on Survey Data) ------------
################################################################################

# Supplementary Figure 3 (Prey Availability) -----------------------------------
# This figure focuses on main prey (Mytilus mussels and acorn barnacles)

# First calculate the mean abundance of mussels at all sites 
MusselAbundanceData <- CommunityData %>% 
  group_by(Site) %>% 
  summarise(Mean_Mussels = mean(California_Mussel_M._californianus),
            Error = sd(California_Mussel_M._californianus)/sqrt(n())) %>% 
  mutate(SE_Upper = Mean_Mussels + Error,
         SE_Lower = Mean_Mussels - Error)

# Next, calculate the mean abundance of mussels in southern California (n = 21.03 % cover)
filter(CommunityData, Region == "SoCal") %>% 
  summarise(SoCalMean = mean(California_Mussel_M._californianus))

# Create the plot 
MusselAbundancePlot <- 
ggplot(data = MusselAbundanceData, aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North")), y = Mean_Mussels)) +
  geom_col(fill = "black", width = 0.80) +
  geom_errorbar(aes(ymin = Mean_Mussels, ymax = SE_Upper), linewidth = 0.75, width = 0) +
  coord_flip() +
  geom_hline(yintercept = 21.03, linewidth = 0.75, linetype = "dashed") + 
  labs(x = "Site", y = "Mussel Abundance\n(Mean % cover ± SE)") +
  scale_y_continuous(limits = c(0, 80)) +
  theme_classic2(base_size = 12)

# Then calculate the mean abundance of barnacles at all sites 
BarnacleAbundanceData <- CommunityData %>% 
  group_by(Site) %>% 
  summarise(Mean_Barnacle = mean(Acorn_Barnacle_Balanus_Chthamalus),
            Error = sd(Acorn_Barnacle_Balanus_Chthamalus)/sqrt(n())) %>% 
  mutate(SE_Upper = Mean_Barnacle + Error,
         SE_Lower = Mean_Barnacle - Error)

# Next, calculate the mean abundance of barnacles in southern California  (n = 11.63 % cover)
filter(CommunityData, Region == "SoCal") %>% 
  summarise(SoCalMean = mean(Acorn_Barnacle_Balanus_Chthamalus))

# Create the plot 
BarnacleAbundancePlot <- 
  ggplot(data = BarnacleAbundanceData, aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North")), y = Mean_Barnacle)) +
  geom_col(fill = "black", width = 0.80) +
  geom_errorbar(aes(ymin = Mean_Barnacle, ymax = SE_Upper), linewidth = 0.75, width = 0) +
  coord_flip() +
  geom_hline(yintercept = 11.63, linewidth = 0.75, linetype = "dashed") + 
  labs(x = "Site", y = "Acorn Barnacle Abundance\n(Mean % cover ± SE)") +
  scale_y_continuous(limits = c(0, 50)) +
  theme_classic2(base_size = 12) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Combine the plots 
SuppFig3 <- ggarrange(MusselAbundancePlot, BarnacleAbundancePlot,
                                nrow = 1, ncol = 2,
                                widths = c(1.0, 0.75))

# Print (700 x 450)
SuppFig3

# Clean the workspace
rm(MusselAbundanceData, MusselAbundancePlot, BarnacleAbundanceData, BarnacleAbundancePlot)
rm(SuppFig3)



# Supplementary Figure 4 (Whelk Abundances) ------------------------------------

# First count all whelks occurring in Timed Counts
QuadWhelkCount_DF <- filter(WhelkLengthData, In_Timed_Count == "N" & Species != "Ol") %>%
  # Subsets data not in TC and confines to only the focal whelks
  group_by(Region, Site, Survey_Date, Vert_Transect_Dist_m, Quad_Dist_m, Species) %>% 
  # Group together
  count() %>% 
  # Count the number of whelks found outside of timed counts
  rename(Number_Whelks = n)

# Make sure dates are date type 
QuadWhelkCount_DF$Survey_Date <- as.Date(QuadWhelkCount_DF$Survey_Date, format = "%m/%d/%Y")

# Extract whelk information from the CommunityData dataframe
WhelkDensity_Df <- CommunityData[, c(3:6, 10:13, 93:98, 101)]
  # Pulls all quadrat metadata (inc. season, which will be used in a moment)

# Merge the two dataframes to get seasonal information 
QuadWhelkCount_DF <- QuadWhelkCount_DF %>% 
  # Pull in the data
  left_join(select(WhelkDensity_Df, Region, Site, Season, Survey_Date, Vert_Transect_Dist_m, Quad_Dist_m),
            # Transect information from the Whelks_DF 
            by = c("Region", "Site", "Survey_Date", "Vert_Transect_Dist_m", "Quad_Dist_m"))
            # Use those columns to join

# Select only spring data
QuadWhelkCount_DF <- filter(QuadWhelkCount_DF, Season == "Spring")
  # This accounts for the fact that SoCal has 2x the data of NorCal and Baja (due to fall and spring surveys)

# Count the total number of them across all surveys
QuadWhelkCount_DF <- QuadWhelkCount_DF %>% 
  # Pull in the data
  group_by(Region, Site, Species) %>% 
  # Group
  summarise(Total_Ind = sum(Number_Whelks))
  # Calculate number of whelks in each group

# Enumerate whelk names 
QuadWhelkCount_DF <- QuadWhelkCount_DF %>% 
  ungroup %>% 
  mutate(Species = case_when(Species == "M" ~ "M. lugubris",
            Species == "As" ~ "A. spirata",
            Species == "Ap" ~ "A. punctulata",
            Species == "Nl" ~ "N. lamellosa",
            Species == "Nc" ~ "N. canaliculata",
            Species == "Ne" ~ "N. emarginata / ostrina",
            Species == "Rp" ~ "R. poulsoni"))

# Order the species 
QuadWhelkCount_DF$Species <- factor(QuadWhelkCount_DF$Species, order = TRUE, levels = c("M. lugubris", "N. emarginata / ostrina", "A. spirata", "N. lamellosa", "R. poulsoni", "A. punctulata", "N. canaliculata"))

# Now graph
ggplot(data = QuadWhelkCount_DF, aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North")), y = Total_Ind)) + 
  geom_col(aes(fill = Species), width = 0.80) +
  scale_fill_viridis(option = "C", discrete = TRUE) +
  coord_flip() +
  labs(x = "Site", y = "Whelk Abundance\n(# Whelks Observed in Spring Surveys)") +
  scale_y_continuous(limits = c(0, 300)) +
  theme_classic2(base_size = 12) + 
  theme(legend.position = c(0.85, 0.70),  legend.text = element_text(face = "italic"))



# Supplementary Table 6 (Whelk Lengths by Region) ------------------------------

# Determine the average lengths of whelks in each region
WhelkLengthData %>% 
  group_by(Species, Region) %>% 
  summarise(Mean_Length = mean(Total_Length, na.rm = TRUE),
            Error = sd(Total_Length, na.rm = TRUE)/sqrt(n()))

# Data was copied and pasted into a Word Table
# Note that Ol (Paciocinebrina lurida) is rare and not a competitor with either species for food  




