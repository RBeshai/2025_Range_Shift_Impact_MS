# Created by Ryan A. Beshai
# Last edited: 11/13/2025

## Load Packages ---------------------------------------------------------------

library(tidyverse)    #loads tidyr, readr, and dplyr
library(vegan)        #community analysis functions
library(lme4)         #for running linear mixed model statistics
library(glmmTMB)      #for running linear mixed model statistics
library(lmerTest)     #for determining the statistical significance
library(cowplot)      #for combining plots
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
CommunityData <- read.csv("SurveyCommunityData.csv")        # Community survey data
TransectData <- read.csv("SurveyTransectData.csv")          # Transect level information 
WhelkLengthData <- read.csv("SurveyWhelkLengthData.csv")    # Timed count and whelk length info

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
QuadWhelkCount_DF <- filter(WhelkLengthData, In_Timed_Count == "N" & Species %in% c("As", "M")) %>%
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
                                            "Mendocino North") ~ NA)) %>% 
    # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel",
                                           "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", 
                                           "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove",
                                           "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                                   Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                                   Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ NA))
    # This chunk assigns "Status" for Acanthinucella at all sites
    # Sites where species is not observed are counted as NA
 
## Calculate a "present/absent" column for both focal species to be included later models
CommunityData$Mex_Present <- as.numeric(CommunityData$Dark_Unicorn_Whelk_M._lugubris != 0)
CommunityData$As_Present <- as.numeric(CommunityData$Angled_Unicorn_Whelk_A._spirata != 0)

## Calculate a log-transformed (log[x + 1]) whelk density 
CommunityData <- CommunityData %>% 
  # Pull in the CommunityData frame
  mutate(As_Log_Trans = log(As_Count + 1),
         Mex_Log_Trans = log(Mex_Count + 1))
  # Log transform whelk densities (b/c some are much higher), adding 1 to account for 0 values

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
                                          "Mendocino North") ~ NA)) %>% 
  # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel",
                                           "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", 
                                           "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove",
                                           "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                               Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                               Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ NA))
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
Whelks_Timed_Count <- filter(Whelks_Timed_Count, Site != "Coal Oil")

# Again, assign "historic" or "expanded" to each site by the species
Whelks_Timed_Count <- Whelks_Timed_Count %>% 
  # Pull in the data
  mutate(Mex_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando") ~ "Historic", 
                                Site %in% c("Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Venice Breakwater", "Santa Monica") ~ "Expanded",
                                Site %in% c("Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North") ~ NA)) %>% 
  # This chunk assigns "Status" for Mexacanthina at all sites
  mutate(As_Status = case_when(Site %in% c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach") ~ "Historic", 
                               Site %in% c("Mendocino South", "Mendocino North") ~ "Expanded",
                               Site %in% c("Sea Ranch", "Moat Creek", "Mussel Rock") ~ NA))
  # This chunk assigns "Status" for Acanthinucella at all sites
  # Sites where species is not observed are counted as NA


# Then join to create the final, graphable dataset
Whelkdata <- left_join(Whelks_Timed_Count, Avg_Lat, by = c("Region", "Site"))

# Graph the densities and relative abundances ----------------------------------

# Remove Rancho Marino from the dataset (no in-person survey, only from photos)
Whelkdata <- filter(Whelkdata, Site != "Rancho Marino") 

# Define site as an ordered factor (important for plotting)
Whelkdata$Site <- factor(Whelkdata$Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Venice Breakwater", "Santa Monica", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North"))

# Calculate means for each site to plot
Whelkdata_means <- Whelkdata %>% 
  # Pull in the Whelkdata
  group_by(Species, As_Status, Mex_Status, Site) %>% 
  # Group data by sites
  summarize(Mean_Whelks = mean(Num_Whelks),
Error = sd(Num_Whelks)/sqrt(n())) %>% 
  mutate(SE_Upper = Mean_Whelks + Error,
         SE_Lower = Mean_Whelks - Error)

# Create the plot for Acanthinucella
AsRelAbundPlot <- 
ggplot() +
  # Open a new ggplot element
  geom_col(data = filter(Whelkdata_means, Species == "As" & Site %!in%  c("Santa Monica", "Venice Breakwater")), 
           aes(x = Site, y = Mean_Whelks, color = As_Status, fill = As_Status), show.legend = F) + 
  # Add a column to show the means 
  geom_errorbar(data = filter(Whelkdata_means, Species == "As" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
               aes(x = Site, ymin = Mean_Whelks, ymax = SE_Upper, color = As_Status), width = 0, alpha = 0.5, 
               show.legend = F) + 
  # Add errorbars to the means 
  geom_point(data = filter(Whelkdata, Species == "As" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks), color = "white", size = 2.25, show.legend = F) +
  # Add a white point underneath each datapoint
  geom_point(data = filter(Whelkdata, Species == "As" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks, color = As_Status), show.legend = F) +
  # Add all of the points I want to show, excluding Santa Monica and Venice Breakwater
  geom_point(data = filter(Whelkdata, Species == "As" & Site %in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks), shape = 23, fill = "steelblue", color = "steelblue") +
  # Add in Santa Monica and Venice Breakwater 
  coord_flip() +
  # Flip the x and y so that sites run north to south 
  scale_color_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000")) +
  # Add a custom point color for historic and expanded regions 
  scale_fill_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000")) +
  # Add a custom fill for historic and expanded regions
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  # Add a scale for the whelk abundance
  scale_x_discrete(drop = FALSE) +
  # Add the x axis, but don't drop sites with no whelks
  labs(x = "Site", y = bquote(atop(Relative~Abundance, ("# Obs. in 1 Hr Timed Count")))) +
  # Change the axis labels
  ggtitle((expression(italic("Acanthinucella")))) +
  # Add a title
  theme_classic2()
  # Give it a nice theme

# Create the plot for Mexacanthina
MexRelAbundPlot <- 
ggplot() +
  # Open new ggplot
  geom_col(data = filter(Whelkdata_means, Species == "M" & Site %!in%  c("Santa Monica", "Venice Breakwater")), 
           aes(x = Site, y = Mean_Whelks, fill = Mex_Status), show.legend = F) + 
  # Add a column to show the means 
  geom_errorbar(data = filter(Whelkdata_means, Species == "M" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
                aes(x = Site, ymin = Mean_Whelks, ymax = SE_Upper, color = Mex_Status), width = 0, alpha = 0.5,
                show.legend = F) + 
  # Add errorbars to the means 
  geom_point(data = filter(Whelkdata, Species == "M" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks), color = "white", size = 2.25, show.legend = F) +
  # Add a white point underneath each datapoint
  geom_point(data = filter(Whelkdata, Species == "M" & Site %!in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks, color = Mex_Status), show.legend = F) +
  # Add all of the points I want to show, excluding Santa Monica and Venice Breakwater
  geom_point(data = filter(Whelkdata, Species == "M" & Site %in% c("Santa Monica", "Venice Breakwater")), 
             aes(x = Site, y = Num_Whelks), shape = 23, fill = "steelblue", color = "steelblue") +
  # Add in Santa Monica and Venice Breakwater 
  coord_flip() +
  # Flip the x and y so that sites run north to south 
  scale_color_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000")) +
  # Add a custom point color for historic and expanded regions 
  scale_fill_manual(breaks = c("Expanded", "Historic"), values = c("#999999", "#000000")) +
  # Add a custom fill for historic and expanded regions
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  # Add a scale for the whelk abundance
  scale_x_discrete(drop = FALSE) +
  # Add the x axis, but don't drop sites with no whelks
  labs(x = " ", y = bquote(atop(Relative~Abundance, ("# Obs. in 1 Hr Timed Count")))) +
  # Change the axis labels
  ggtitle((expression(italic("Mexacanthina")))) +
  # Add a title
  theme_classic2() +
  # Give it a nice theme
  theme(legend.position = c(0.80, 1))


# Combine the single species relative abundance figures into a single plot
Figure2 <- ggarrange(AsRelAbundPlot, MexRelAbundPlot,
                                            labels = c("a.", "b."),
                                            nrow = 1)

# Print the new figure (Figure 2 in Manuscript)
Figure2
    # Print as 750 x 450


## Statistical Analysis of relative abundance patterns

# Build a model testing Acanthinucella abundance over space
As_Density_Lat_Mod <- lmer(Num_Whelks ~ Avg_Lat + As_Status + (1|Site), data = filter(Whelkdata, Species == "As" & !is.na(As_Status) & Site %!in% c("Santa Monica", "Venice")))
  # Site not included as random effect because it is singular with Avg_Lat
  # NA sites are not included because I don't want to try to predict whelk counts at sites where we 
  # never found whelks 

# Check model assumptions
check_model(As_Density_Lat_Mod)

# Test the significance
summary(As_Density_Lat_Mod)


# Build a model testing Mexacanthina abundance over space
Mex_Density_Lat_Mod <- lmer(Num_Whelks ~ Avg_Lat + Mex_Status + (1|Site), data = filter(Whelkdata, Species == "M" & !is.na(Mex_Status) & Site %!in% c("Santa Monica", "Venice")))
  # Site not included as random effect because it is singular with Avg_Lat
  # NA sites are not included because I don't want to try to predict whelk counts at sites where we 
  # never found whelks 

# Check model assumptions
check_model(Mex_Density_Lat_Mod)

# Test the significance
summary(Mex_Density_Lat_Mod)

# Remove statistical models
rm(As_Density_Lat_Mod, Mex_Density_Lat_Mod)

# Remove the graphs and dataframes from the working space
rm(Whelkdata, Avg_Lat, WhelkDensity_Df, QuadWhelkCount_DF)
rm(Figure2, AsRelAbundPlot, MexRelAbundPlot)



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


# Now calculate the proportion of living barnacles at each site and log-transform mussel abundance
CommunityData_Prey <- CommunityData_Prey %>% 
  # Pull in the data
  mutate(Prop_LiveBarn = (Acorn_Barnacle_Balanus_Chthamalus/(Acorn_Barnacle_Balanus_Chthamalus + Dead_Barnacles)),
         Mussel_Prop = California_Mussel_M._californianus/100)
  # Calculating the log of the mussel abundance here for reproducability 

# Replace NaN with NA
CommunityData_Prey$Prop_LiveBarn[is.nan(CommunityData_Prey$Prop_LiveBarn)] <- NA
  # NaNs result from quadrats with 0 Balanus/Chthamalus
  # Converting to NA will omit them from analysis and figures


# Build and Run LMMs Predator-Prey Abundance Relationships ---------------------

# There will be a total of 4 main tests run: 
  # Acanthinucella vs (1) proportion living barnacles and (2) California mussels
  # Mexacanthina vs (3) proportion living barnacles and (4) California mussels
  # For each test, we will define two models (interactive and reduced) and pick one to interpret using AIC

# (1) Acanthinucella v  Prop. Living Barnacles 
# Build model to test how Acanthinucella's density varies with latitude
LiveBarn_AsGlobal <- glmmTMB(Prop_LiveBarn ~ As_Log_Trans*As_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), family = ordbeta(), data = filter(CommunityData_Prey, !is.na(As_Status)))

# Check model assumptions
check_model(LiveBarn_AsGlobal)
  # VIF are high in interactive model (which is expected), but fine when not looking at interactions 

# Check significance
As_LiveBarn_Mod <- summary(LiveBarn_AsGlobal)


# (2) Acanthinucella v Mussels 
Mussels_AsGlobal <- glmmTMB(Mussel_Prop ~ As_Log_Trans + As_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m),
                             dispformula = ~ Quad_TH_m,
                             zi = ~ Quad_TH_m + Site, 
                             family = ordbeta(link = "probit"), 
                             data = filter(CommunityData_Prey, As_Status != "NA"))
  # Because of exceptionally high number of predicted zeros, we are now including 
  # a zero-inflated component and a dispersion factor and switching to a probit link. 
  # Dispersion varies with tide height --> becomes less variable
  # NOTE: this interactive model is singular (likely due to the very low number of quads
  # with whelks and mussels in the expanded region). We are interpreting cautiously because
  # the interpretations qualitatively align with the non-interactive model. 

# Check model assumptions
check_model(Mussels_AsGlobal)
  # VIF good
  # Variance is poor (because of very few expanded quads with mussels and Acanthinucella

# Check significance
As_Mussel_Mod <- summary(Mussels_AsGlobal)


# (3) Mexacanthina vs Living Barnacles 
LiveBarn_MexGlobal <- glmmTMB(Prop_LiveBarn ~ Mex_Log_Trans + Mex_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), family = ordbeta(), data = filter(CommunityData_Prey, Mex_Status != "NA"))

# Check model assumptions
check_model(LiveBarn_MexGlobal)

# Check significances
M_LiveBarn_Mod <- summary(LiveBarn_MexGlobal)


# (4) Mexacanthina v Mussels 
Mussels_MexGlobal <- glmmTMB(Mussel_Prop ~ Mex_Log_Trans*Mex_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m),
                             dispformula = ~ Quad_TH_m,
                             zi = ~ Quad_TH_m + Site, 
                             family = ordbeta(link = "probit"), 
                             data = filter(CommunityData_Prey, !is.na(Mex_Status)))

# Check model assumptions
check_model(Mussels_MexGlobal)

# Check significance
M_Mussels_Mod <- summary(Mussels_MexGlobal)


# Print each model below to see significance of observational Pred-Prey relationship tests
As_LiveBarn_Mod
As_Mussel_Mod
M_LiveBarn_Mod
M_Mussels_Mod


# Graph Whelk Abundances vs Prey Abundances (Supp Fig 1) -----------------------

# Section 1: predict from the data
  # Here we will create new dataframes for each whelk species 
  # Then, we will predict the model outputs (averaged across sites)
  # Model predictions can then be plotted over the data later

# Create a dataframe to predict Acanthinucella models 
filtered_data <- CommunityData_Prey %>% 
  filter(., !is.na(As_Status))

as_prediction_data <- expand.grid(
  As_Log_Trans = seq(min(filtered_data$As_Count, na.rm = TRUE),
                      max(filtered_data$As_Count, na.rm = TRUE),
                      length.out = 100),
  As_Status = unique(filtered_data$As_Status),
  Quad_TH_m = mean(filtered_data$Quad_TH_m, na.rm = TRUE),
  Site = unique(filtered_data$Site),
  Vert_Transect_Dist_m = mean(filtered_data$Vert_Transect_Dist_m, na.rm = TRUE)
)

# Predict values from the Acanthinucella - Barnacle Model
as_prediction_data$predicted_barn <- predict(LiveBarn_AsGlobal, newdata = as_prediction_data, type = "response", allow.new.levels = TRUE)

# Predict values from the Acanthinucella - Mussel Model
as_prediction_data$predicted_mussel <- predict(Mussels_AsGlobal, newdata = as_prediction_data, type = "response", allow.new.levels = TRUE)

# Then average across all sites
as_avg_predictions <- as_prediction_data %>%
  group_by(As_Log_Trans, As_Status) %>%
  summarise(mean_barn = mean(predicted_barn, na.rm = TRUE),
            mean_mussel = mean(predicted_mussel, na.rm = TRUE), 
            .groups = "drop")

# Average across sites
as_avg_predictions <- as_prediction_data %>%
  group_by(As_Log_Trans, As_Status) %>%
  summarise(mean_barn = mean(predicted_barn, na.rm = TRUE),
            mean_mussel = mean(predicted_mussel, na.rm = TRUE),
            .groups = "drop")

# Determine the maximum whelk abundance values to trim prediction plots
max_x_by_status <- CommunityData_Prey %>% 
  filter(!is.na(As_Status)) %>% 
  group_by(As_Status) %>% 
  summarize(max_x = max(As_Log_Trans, na.rm = TRUE))

# Join these to the prediction frame and trim
as_avg_predictions_trimmed <- as_avg_predictions %>% 
  left_join(max_x_by_status, by = "As_Status") %>% 
  filter(As_Log_Trans <= max_x)


## Repeat the process for Mexacanthina models

# Create a dataframe to predict Mexacanthina models 
filtered_data <- CommunityData_Prey %>% 
  filter(., Mex_Status != "NA")

mex_prediction_data <- expand.grid(
  Mex_Log_Trans = seq(min(filtered_data$Mex_Count, na.rm = TRUE),
                     max(filtered_data$Mex_Count, na.rm = TRUE),
                     length.out = 100),
  Mex_Status = unique(filtered_data$Mex_Status),
  Quad_TH_m = mean(filtered_data$Quad_TH_m, na.rm = TRUE),
  Site = unique(filtered_data$Site),
  Vert_Transect_Dist_m = mean(filtered_data$Vert_Transect_Dist_m, na.rm = TRUE)
)

# Predict values from the Mexacanthina - Barnacle Model
mex_prediction_data$predicted_barn <- predict(LiveBarn_MexGlobal, newdata = mex_prediction_data, type = "response", allow.new.levels = TRUE)

# Predict values from the Mexacanthina - Mussel Model
mex_prediction_data$predicted_mussel <- predict(Mussels_MexGlobal, newdata = mex_prediction_data, type = "response", allow.new.levels = TRUE)

# Then average across all sites
mex_avg_predictions <- mex_prediction_data %>%
  group_by(Mex_Log_Trans, Mex_Status) %>%
  summarise(mean_barn = mean(predicted_barn, na.rm = TRUE),
            mean_mussel = mean(predicted_mussel, na.rm = TRUE), 
            .groups = "drop")

# Determine the maximum whelk abundance values to trim prediction plots
max_x_by_status <- CommunityData_Prey %>% 
  filter(!is.na(Mex_Status)) %>% 
  group_by(Mex_Status) %>% 
  summarize(max_x = max(Mex_Log_Trans, na.rm = TRUE))

# Join these to the prediction frame and trim
mex_avg_predictions_trimmed <- mex_avg_predictions %>% 
  left_join(max_x_by_status, by = "Mex_Status") %>% 
  filter(Mex_Log_Trans <= max_x)


# Section 2: Plot raw data and model prediction line

# Acanthinucella abundance vs proportion living barnacles
Plot_As_LiveBarn <- 
ggplot(data = filter(CommunityData_Prey, As_Status != "NA"), 
       aes(x = As_Log_Trans, y = Prop_LiveBarn)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = FALSE) +
  geom_line(data = as_avg_predictions_trimmed, aes(x = As_Log_Trans, y = mean_barn, color = As_Status), show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = "Acorn Barnacle Abundance\n(Proportion Live Cover)") +
  ggtitle("a.") +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(-0.01, 1.01)) +
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(-0.1, 3)) + 
  theme_classic2(base_size = 10) 

# Acanthinucella abundance vs California mussels 
Plot_As_Mussels <- 
ggplot(data = filter(CommunityData_Prey, As_Status != "NA"), 
       aes(x = As_Log_Trans, y = Mussel_Prop)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(-0.01, 1.01)) +
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(-0.1, 3)) +
  labs(x = expression(italic("Acanthinucella")~"Abundance (ln["*"#" ~ ind*"])"), 
       y = "California Mussel Abundance\n(Proportion Cover)") +
  ggtitle("c.") +
  theme_classic2(base_size = 10) +
  theme(legend.position = c(0.80,0.85))


 # Mexacanthina abundance vs the proportion of living barnacles 
Plot_Mex_LiveBarn <- 
ggplot(data = filter(CommunityData_Prey, Mex_Status != "NA"), 
         aes(x = Mex_Log_Trans, y = Prop_LiveBarn)) +
  geom_line(data = mex_avg_predictions_trimmed, 
            aes(x = Mex_Log_Trans, y = mean_barn, color = Mex_Status), show.legend = F,) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +

  scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(-0.01, 1.0)) +
  scale_x_continuous(breaks = seq(0, 5, by = 1), limits = c(-0.1, 5)) + 
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("b.") +
  theme_classic2(base_size = 10) 

# Mexacanthina abundance vs California mussels 
Plot_Mex_Mussels <- 
  ggplot(data = filter(CommunityData_Prey, Mex_Status != "NA"), 
         aes(x = Mex_Log_Trans, y = Mussel_Prop)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = F) +
    geom_line(data = mex_avg_predictions_trimmed, 
           aes(x = Mex_Log_Trans, y = mean_mussel, color = Mex_Status), show.legend = F) +
    scale_y_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(-0.01, 1.0)) +
  scale_x_continuous(breaks = seq(0, 5, by = 1), limits = c(-0.1, 5)) + 
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(italic("Mexacanthina")~"Abundance (ln["*"#" ~ ind*"])"), y = " ", color = "Status") +
  ggtitle("d.") +
  theme_classic2(base_size = 10) +
  theme(legend.position = c(0.80,0.85))
  
# Combine the plots 
FigureS1 <- plot_grid(Plot_As_LiveBarn, Plot_Mex_LiveBarn, 
                                 Plot_As_Mussels, Plot_Mex_Mussels,
                                     ncol = 2, nrow = 2, align = "hv")

# Print the figure (export as 600 x 500)
FigureS1



## Clean workspace of Statistical models

# Removing Acanthinucella models 
rm(LiveBarn_AsGlobal, Mussels_AsGlobal)

# Removing Mexacanthina models 
rm(LiveBarn_MexGlobal, Mussels_MexGlobal)

# Remove model summaries 
rm(As_LiveBarn_Mod, As_Mussel_Mod, M_LiveBarn_Mod, M_Mussels_Mod)



# Clean workspace of plotting objects ------------------------------------------
rm(as_prediction_data, as_avg_predictions, as_avg_predictions_trimmed, as_groups)
rm(mex_prediction_data, mex_avg_predictions, mex_avg_predictions_trimmed, filtered_data, max_x_by_status)
rm(FigureS1, Plot_As_LiveBarn, Plot_As_Mussels, Plot_Mex_LiveBarn, Plot_Mex_Mussels)



# Alternatively, Analyze with Quantile Regression Analysis ---------------------

# Load the necessary package
library(quantreg)

# Define the quantiles of interest
taus <- c(0.25, 0.50, 0.75)

##  Define the predator - prey combinations and run the models
# Acanthinucella - Barnacles
As_Barn_Quantreg <- lapply(taus, function(t) {
  rq(Prop_LiveBarn ~ poly(As_Log_Trans, 2, raw = TRUE) * As_Status + Quad_TH_m, 
     tau = t,
     data = filter(CommunityData_Prey, !is.na(As_Status)))
})

# Mexacanthina - Barnacles
Mex_Barn_Quantreg <- lapply(taus, function(t) {
  rq(Prop_LiveBarn ~ poly(Mex_Log_Trans, 2, raw = TRUE) * Mex_Status + Quad_TH_m, 
     tau = t,
     data = filter(CommunityData_Prey, !is.na(Mex_Status)))
})

# Mexacanthina - Mussels
Mex_Mussel_Quantreg <- lapply(taus, function(t) {
  rq(Mussel_Prop ~ poly(Mex_Log_Trans, 2, raw = TRUE) * Mex_Status + Quad_TH_m, 
     tau = t,
     data = filter(CommunityData_Prey, !is.na(Mex_Status)))
})


## Summarize with bootstrapped SEs

# Set the seed
set.seed(123)

# Run the bootstrap for Acanthinucella - Barnacles
As_Barn_Quant_Sum <- lapply(As_Barn_Quantreg, function(m) {
  summary(m, se = "boot", R = 1000)
})

# Run the bootstrap for Mexacanthina - Barnacles
Mex_Barn_Quant_Sum <- lapply(Mex_Barn_Quantreg, function(m) {
  summary(m, se = "boot", R = 1000)
})

# Run the bootstrap for Mexacanthina - Mussels
Mex_Mus_Quant_Sum <- lapply(Mex_Mussel_Quantreg, function(m) {
  summary(m, se = "boot", R = 1000)
})

# Get the summaries
As_Barn_Quant_Sum
Mex_Barn_Quant_Sum
Mex_Mus_Quant_Sum


# Clean the Workspace ----------------------------------------------------------
rm(taus)
rm(As_Barn_Quantreg, Mex_Barn_Quantreg, Mex_Mussel_Quantreg)
rm(As_Barn_Quant_Sum, Mex_Barn_Quant_Sum, Mex_Mus_Quant_Sum)




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
  # Calculate evenness for each species 


# Build and Run LMMs Assessing Predator-Diversity Relationships ----------------
# There will be 6 total analyses run: 
  # Acanthinucella vs (1) Shannon Diversity, (2) Richness, (3) Evenness
  # Mexacanthina vs (4) Shannon Diversity, (5) Richness, (6) Evenness

# Acanthinucella - Shannon Diversity 
Diversity_AsGlobal <- lmer(SW_Diversity_As ~ As_Log_Trans*As_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check model assumptions
check_model(Diversity_AsGlobal)

# Check significance
As_Diversity_Mod <- summary(Diversity_AsGlobal)


# Acanthinucella - Richness 
Richness_AsGlobal <- lmer(Richness_As ~ As_Log_Trans*As_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check model assumptions
check_model(Richness_AsGlobal)

# Check significance
As_Richness_Mod <- summary(Richness_AsGlobal)


# Acanthinucella - Evenness 
Evenness_AsGlobal <- lmer(Evenness_As ~ As_Log_Trans*As_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, As_Status != "NA"))

# Check model assumptions
check_model(Evenness_AsGlobal)

# Check significance
As_Evenness_Mod <- summary(Evenness_AsGlobal)


# Mexacanthina - Shannon Diversity 
Diversity_MexGlobal <- lmer(SW_Diversity_Mex ~ Mex_Log_Trans*Mex_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check model assumptions
check_model(Diversity_MexGlobal)

# Check significance
Mex_Diversity_Mod <- summary(Diversity_MexGlobal)



# Mexacanthina - Richness 
Richness_MexGlobal <- lmer(Richness_Mex ~ Mex_Log_Trans*Mex_Status + Quad_TH_m + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check model assumptions
check_model(Richness_MexGlobal)

# Check significance
Mex_Richness_Mod <- summary(Richness_MexGlobal)

## Mexacanthina v Evenness (2 Models, pick 1 with AIC)
Evenness_MexGlobal <- lmer(Evenness_Mex ~ Mex_Log_Trans*Mex_Status + Quad_TH_m + + (1|Site) + (1|Vert_Transect_Dist_m), data = filter(CommunityData, Mex_Status != "NA"))

# Check model assumptions
check_model(Evenness_MexGlobal)

# Check significance
Mex_Evenness_Mod <- summary(Evenness_MexGlobal)


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


 
# Graph Whelk Abundances vs Diversity Metrics (Supp Fig 4) ---------------------

# To denote significance, geom_smooths will be added to plots significant relationships from above

# Plots Acanthinucella abundance v Shannon diversity 
Graph_DiversityvAcanth <- 
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = As_Log_Trans, y = SW_Diversity_As)) +
  geom_point(aes(color = As_Status), position = "jitter",  show.legend = F) +
  #geom_smooth(color = "black", method = lm, fullrange = T, se = F, show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = "", y = "Shannon", color = "Status") + 
  ggtitle("a.") +
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(-0.1, 3)) + 
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(-0.05, 2.5)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v Shannon diversity 
Graph_DiversityvMex <- 
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Log_Trans, y = SW_Diversity_Mex)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +
  geom_smooth(aes(color = Mex_Status), method = lm, fullrange = T, se = F, show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("b.") +
  scale_x_continuous(breaks = seq(0, 5, by = 1), limits = c(-0.1, 5)) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(-0.05, 2.5)) +
  theme_classic2()
    

# Plots Acanthinucella abundance v richness
Graph_RichnessvAcanth <-
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = As_Log_Trans, y = Richness_As)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  #geom_smooth(method = lm) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = "Richness", color = "Status") +
  ggtitle("c.") +
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(-0.1, 3)) + 
  scale_y_continuous(breaks = seq(0, 20, by = 4), limits = c(0, 20)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v richness
Graph_RichnessvMex <-
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Log_Trans, y = Richness_Mex)) +
  geom_point(aes(color = Mex_Status), position = "jitter", show.legend = FALSE) +
  geom_smooth(aes(color = Mex_Status), method = lm, fullrange = T, se = F, show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = " ", y = " ") +
  ggtitle("d.") +
  #scale_x_continuous(breaks = seq(0, 125, by = 25), limits = c(-1, 125)) + 
  scale_y_continuous(breaks = seq(0, 20, by = 4), limits = c(-0.05, 20)) +
  theme_classic2()


# Plots Acanthinucella abundance v evenness
Graph_EvennessvAcanth <-
ggplot(data = filter(CommunityData, As_Status != "NA"), aes(x = As_Log_Trans, y = Evenness_As)) +
  geom_point(aes(color = As_Status), position = "jitter", show.legend = F) +
  #geom_smooth(color = "black", method = lm, fullrange = T, se = F, show.legend = F) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(atop(italic("Acanthinucella")~"Abundance (ln["*"#" ~ ind*"])")), y = "Evenness", color = "Status") +
  ggtitle("e.") +
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(-0.1, 3)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(-0.05, 1.05)) +
  theme_classic2() +
  theme(legend.position = c(0.85, 0.15))


# Plots Mexacanthina abundance v evenness
Graph_EvennessvMex <-
ggplot(data = filter(CommunityData, Mex_Status != "NA"), 
       aes(x = Mex_Log_Trans, y = Evenness_Mex)) +
  geom_point(aes(color = Mex_Status), show.legend = FALSE) +
  #geom_smooth(aes(color = Mex_Status), method = lm, fullrange = T, se = F, show.legend = FALSE) +
  scale_color_manual(breaks = c("Historic", "Expanded"), values = c("#000000", "#999999")) +
  labs(x = expression(atop(italic("Mexacanthina")~"Abundance (ln["*"#" ~ ind*"])")), y = " ") +
  ggtitle("f.") +
  scale_x_continuous(breaks = seq(0, 5, by = 1), limits = c(-0.1, 5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(-0.05, 1.05)) +
  theme_classic2()

# Create the combined figure with all diverity metrics
FigureS4 <- plot_grid(Graph_DiversityvAcanth, Graph_DiversityvMex,
                                  Graph_RichnessvAcanth, Graph_RichnessvMex,
                                  Graph_EvennessvAcanth, Graph_EvennessvMex,
                                  ncol = 2, nrow = 3, align = "hv")

# Notes on dropped values - 
  # Most are NAN that result from "communities" with only one sp observed.
  # This occurs because the formula equates to 0/log(1), which is undefined.
  # Some are not visualized because there is not enough space to jitter 
  # or cleanly draw on top of other points. 


# Print the figure (Supp Fig.4, export as 500 x 800)
FigureS4

# Remove figures
rm(Graph_DiversityvMex, Graph_DiversityvAcanth, Graph_Combined_WhelksvDiversity)
rm(Graph_RichnessvMex, Graph_RichnessvAcanth, Graph_Combined_WhelksvRichness)
rm(Graph_EvennessvMex, Graph_EvennessvAcanth, Graph_Combined_WhelksvEvenness)



################################################################################
# Additional Supplemental Figures and Tables (Based on Survey Data) ------------
################################################################################

# Supplementary Figure 5 (Prey Availability) -----------------------------------
# This figure focuses on main prey (Mytilus mussels and acorn barnacles)

# Create new dataframe
musbarn_data <- CommunityData

# Define the site order
musbarn_data$Site <- factor(musbarn_data$Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North"))

# Create a mussel plot (show raw data and distributions)
MusselAbundancePlot <- ggplot(data = musbarn_data, aes(x = Site, y = California_Mussel_M._californianus)) +
  geom_point(aes(color = Region), show.legend = F) +
  stat_summary(fun.y = mean, geom = "crossbar", color = "black") +
  geom_hline(yintercept = 21.03, linewidth = 0.75, color = "grey40", linetype = "dashed") + 
  coord_flip() +
  labs(title = "a.", x = "Site", y = "Mussel Abundance\n(% Cover)") +
  scale_color_viridis(discrete = T) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  theme_classic2(base_size = 12)

BarnacleAbundancePlot <- ggplot(data = musbarn_data, aes(x = Site, y = Acorn_Barnacle_Balanus_Chthamalus)) +
  geom_point(aes(color = Region), show.legend = F) +
  stat_summary(fun.y = mean, geom = "crossbar", color = "black") +
  geom_hline(yintercept = 11.63, linewidth = 0.75, color = "grey40", linetype = "dashed") + 
  coord_flip() +
  labs(title = "b.", x = "Site", y = "Acorn Barnacle Abundance\n(% Cover)") +
  scale_color_viridis(discrete = T) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  theme_classic2(base_size = 12) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Combine the plots 
FigureS5 <- ggarrange(MusselAbundancePlot, BarnacleAbundancePlot,
                                nrow = 1, ncol = 2,
                                widths = c(1.0, 0.75))

# Print (700 x 450)
FigureS5

# Clean the workspace
rm(musbarn_data, BarnacleAbundanceData, BarnacleAbundancePlot)
rm(FigureS5)



# Supplementary Figure 6 (Whelk Abundances) ------------------------------------

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
FigureS6 <- ggplot(data = QuadWhelkCount_DF, aes(x = factor(Site, levels = c("La Chorera", "La Chorera Norte", "Campo Kennedy", "Punta Morro", "San Miguel", "Saldamando", "Cabrillo", "Scripps", "Cardiff", "Swami's", "Dana Point", "Goff Island", "Victoria Beach", "Heisler Park", "Shaw's Cove", "Crystal Cove", "Little Corona", "Rancho Marino", "Dillon Beach", "Sea Ranch", "Moat Creek", "Mussel Rock", "Mendocino South", "Mendocino North")), y = Total_Ind)) + 
  geom_col(aes(fill = Species), width = 0.80) +
  scale_fill_viridis(option = "C", discrete = TRUE) +
  coord_flip() +
  labs(x = "Site", y = "Whelk Abundance\n(# Whelks Observed in Spring Surveys)") +
  scale_y_continuous(limits = c(0, 300)) +
  theme_classic2(base_size = 12) + 
  theme(legend.position = c(0.85, 0.70),  legend.text = element_text(face = "italic"))

# Print the figure
FigureS6

# Clean workspace
rm(QuadWhelkCount_DF, FigureS6)


# Supplementary Table 6 (Whelk Lengths by Region) ------------------------------

# Determine the average lengths of whelks in each region
WhelkLengthData %>% 
  group_by(Species, Region) %>% 
  summarise(Mean_Length = mean(Total_Length, na.rm = TRUE),
            Error = sd(Total_Length, na.rm = TRUE)/sqrt(n()))

# Data was copied and pasted into a Word Table
# Note that Ol (Paciocinebrina lurida, previously Ocinebra lurida) is rare and not a competitor with either species for food  




