# Applying invasion biology frameworks to predict impacts of range-expanding predators 
This README file was generated on 2024-01-26 

## Data Collection Methods
*Overview:*  
There are two main data sources in this dataset: (1) observational survey data collected from 2022-2023, spanning 23 intertidal sites from Baja California, Mexico, through northern California, USA, and (2) manipulative density data from 6 sites in Baja California, Mexico, and northern California USA. These geographic ranges compromise the historic and expanded ranges of the two range-expanding intertidal whelks, Acanthinucella spirata and Mexacanthina lugubris. 
 
*Observational Data:*
To assess the current distributions and abundances of Acanthinucella and Mexacanthina, as well as in situ associations with abundances of prey species and community diversity metrics, we conducted field surveys at 23 sites across California, USA and Baja California, Mexico between April 2022 and October 2023. Sites in northern California, USA, and Mexico were surveyed once per year in the spring (2 surveys per site), and sites in southern California, USA, were surveyed in the spring and fall each year (4 surveys per site). In addition, to ensure we included the most recent range information for these species, we conducted a single, partial survey (timed count only, using methods as described below) at the Santa Monica Pier and Venice Breakwater in March 2024, following recent reports posted on iNaturalist.com. 
At each site, we laid a 25 m horizontal transect at the top of the intertidal (i.e., parallel to the water at the top of the barnacle zone) and 5 vertical transects (perpendicular to and extending toward the waterline) at random locations along the horizontal transect. Along each vertical transect, we placed 0.25 m × 0.25 m (0.0625 m2) quadrats every 0.25 m drop in vertical elevation. Within each quadrat, we conducted community diversity surveys where we assessed the relative abundance of all species, both mobile and sessile, as percent cover. To estimate the relative abundance of each expanding whelk species, we conducted an additional, 1-hour timed search for Acanthinucella and Mexacanthina individuals following the diversity surveys; any individual whelks found during timed counts were associated with a shore height using an Apache Tools Mesa SL101 laser level and local tide predictions (https://tidesandcurrents.noaa.gov/tide_predictions.html).
 
*Manipulative Density Experiment: *
To test the relationship between the abundance of range-expanding whelks and impacts on (a) prey abundance and (b) community diversity, as well as compare impacts across historic and expanded regions, we conducted a manipulative caging experiment from May to September in 2023 at a subset of our survey sites. We identified two Acanthinucella-expanded sites in northern California (Mendocino North and Mendocino South) and two Acanthinucella-historic sites in southern California (Dana Point and Scripps Reserve). Similarly, we identified two Mexacanthina-expanded sites in southern California (Dana Point and Scripps Reserve; at a different elevation than Acanthinucella-historic plots) and two Mexacanthina-historic sites in Baja California, Mexico (Punta Morro and Campo Kennedy).

At each site, we identified seven 33 × 33 cm (0.1 m2) plots parallel to the water line at shore heights of greatest focal whelk abundance, as determined by survey data from 2022 (~0.5 m above MLLW for Acanthinucella plots [mean plot shore height = 0.56 m, range = 0.28-0.81 m] and ~1.0 m above MLLW for Mexacanthina plots [mean plot shore height = 0.87 m, range = 0.65-1.13 m]). Mexacanthina plots contained at least two of these three target prey groups: acorn barnacles (Chthamalus & Balanus spp.), gooseneck barnacles (Pollicipes polymerus), and/or mussels (Mytilus spp.). Gooseneck barnacles and mussels were rare at Acanthinucella’s most common shore height, where plots contained primarily acorn barnacles, shown by previous studies to be their main prey. We randomly assigned a treatment to each of the seven plots: one of two controls—no cage or partial cage (where a cage was deployed but two of the sides are left open)—or one the following density treatments: 0, 1, 3, 6, or 12 whelks added per cage. Thus, the whelk-addition plots represented densities of 0, 10, 30, 60, or 120 individuals per m2, spanning a range of densities encountered in field surveys. We occasionally observed higher densities of Mexacanthina when they were aggregated at sites in Mexico (commonly 100-300 whelks per m2), so we established an additional 24 whelk-addition plot (representing 240 individuals per m2) in both Mexacanthina-historic sites. To control for potential differences in feeding rates between whelks of various sizes within a site, we put equal numbers of small, medium, and large whelks in each plot (e.g., four small, four medium, and four large whelks in the 12-whelk plot), except the 1-whelk plot, which received a single medium individual. Size classes were defined at the site level based on the population structure of the whelks we encountered at each site. We removed all other whelk species from the whelk-addition plots, so cages contained only one focal whelk species. All cages were constructed of welded stainless steel mesh cloth (Master-Carr part ID: 9322T64).   

After we established the plots, we monitored (via community surveys, focal whelk measurements, and archival photographs) and maintained cages (cleaned and whelks replaced as needed to maintain density treatments) every 2 weeks for 8 weeks. Any focal whelks lost from the cages were replaced with medium-sized individuals. We initially conducted all community counts in the field and later verified or corrected counts from archival photographs. Mexacanthina cages were deployed from May to August 2023, and Acanthinucella cages were deployed from late June through September 2023.

## Data Collection Period
Spring (Apr-June) and Fall (Sept-early Nov), 2022-2023.

## Files and Description

### Data and Metadata Files 
**ExptCommunityData (.csv)**
Provides community data from the manipulative experiment. Data were collected in approx. 2 week intervals (see Methods for additional details). Rows correlate to individual observations of each plot. Columns provide plot identifying information and community data. 

**ExptData_Metadata (.csv)** 
Provides the metadata for all columns in the ExptCommunityData file, including all abbreviations and codes. 

**SurveyCommunityData (.csv)**
Contains all observational community survey data. Rows correspond to individual quadrats, and columns provide quadrat identifying data and community information. 

**SurveyTransectData (.csv)**
Contains information about each vertical transect from the observational data. Rows correspond to individual transects, and columns are used to calculate shore height at the top of each transect (relative to mean lower low water, MLLW). 

**SurveyWhelkLengthData (.csv)**
Contains information about the size of all whelk individuals captured in the observational data. Rows correspond to individual whelks. 

**SurveyData_Metadata (.csv)**
Provides the metadata for all observational survey ("Survey") sheets, including all abbreviations and codes. 

### Script Files
**Cage_Predation_Cleaned.R**
This script can be used to generate the statistical analysis of all predation rates and predation figures derived from the manipulative cage experiment. Figures are named according to their placement in the manuscript. 

**Cage_Diversity_Cleaned.R**
This script can be used to generate the statistical analysis of all changes in diversity and diversity figures derived from the manipulative cage experiment. Figures are named according to their placement in the manuscript. 

**Survey_Analysis_Cleaned.R**
This script can be used to generate the statistical analysis of all data and figures generated from the observational study.   Figures are named according to their placement in the manuscript. 

## Recommended Workflow
Download all datafiles, R scripts, and the R project. In addition to statistical analysis and figure creation, data cleaning and manipulation is performed within each script. 

After downloading, begin with the "Cage" R files (as the manuscript focuses on the analysis of the manipulative experiment data). Each script will call the necessary packages and, if open in R studio, will prompt the installation of said packages. Each script is broken into sections of analysis (delineated by "## [Description] ------------"). Each analysis and associated figures can be generated by running all code within a section. Both "Cage" scripts rely on the same data (ExptCommunityData.csv). 

Proceed to the observational survey data files. All three "Survey" data files are used in the SurveyAnalysis_Cleaned.R script. Again, this script is broken into sections, which can be run independently, and approximately follows the order of presentation in the manuscript. 