# Whelk Analysis R Scraps 



################################################################################
## Estimate Whelk Counts in Quadrats -------------------------------------------
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
As_Density_Lat_Mod <- glmmTMB(As_Count ~ Avg_Lat + Quad_TH_m + As_Status + (1|Site/Vert_Transect_Dist_m),
                              ziformula = ~ Quad_TH_m,
                              family = nbinom2, 
                              data = dplyr::filter(WhelkDensity_Df, Species == "As" & As_Status != "NA"))
# Predicts Acanthinucella counts using Lat, Quad_TH, and Status as fixed effects
# Random effects of Transect nested within Site
# Uses a zero inflated formula with Quad TH as the best predictor
# Per model checking, negative binomial yields better fits than Poisson 
# NA sites are not included because I don't want to try to predict whelk counts at sites where we 
# never found whelks 

# Check model assumptions 
check_model(As_Density_Lat_Mod)

# Test significance 
summary(As_Density_Lat_Mod)


# Build model to test how Mexacanthina's density varies with latitude 
Mex_Density_Lat_Mod <- glmmTMB(As_Count ~ Avg_Lat + Quad_TH_m + Mex_Status + (1|Site/Vert_Transect_Dist_m),
                               ziformula = ~ Quad_TH_m,
                               family = nbinom1, 
                               data = dplyr::filter(WhelkDensity_Df, Species == "M" & Mex_Status != "NA"))
# Predicts Mexacanthina counts using Lat, Quad_TH, and Status as fixed effects
# Random effects of Transect nested within Site
# Uses a zero inflated formula with Quad TH as the best predictor
# Per model checking, negative binomial yields better fits than Poisson 
# NA sites are not included because I don't want to try to predict whelk counts at sites where we 
# never found whelks 


# Check model assumptions
check_model(Mex_Density_Lat_Mod)

# Test the significance
summary(Mex_Density_Lat_Mod)

# Clean the Workspace
rm(As_Density_Lat_Mod, Mex_Density_Lat_Mod)

################################################################################
################################################################################













