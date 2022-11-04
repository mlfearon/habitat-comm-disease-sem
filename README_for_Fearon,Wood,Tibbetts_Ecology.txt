Metadata for Fearon, Wood, and Tibbetts. 2023. "Habitat Quality Influences Pollinator Pathogen Prevalence Through Both Habitat–Disease and Biodiversity–Disease Pathways". Ecology.

Data collected by: Michelle L. Fearon
Contact: mlfearon@umich.edu


The R code that generates the path analysis and statistics for the figures included in this manuscript and uses all of these data sets 
can be found in BeeVirusHabitat+CommunityPathModels_29Sept2022.R


File name:
Virus_Positive_Community+Habitat_Covariates_1Mar2022.csv

Description: 
Records of which individual bees from four focal host species (Apis mellifera, Bombus impatiens, Eucera pruinosa, and Lasioglossum spp.) 
were positive or negative when tested for the viral positive strand for deformed wing virus (DWV), black queen cell virus (BQCV), and 
sacbrood virus (SBV), and additional pollinator community and local and landscape habtiat quality co-variates data from each site. All bees were collected in or adjacent to winter squash fields in 2015 or 2016. 
All samples were randomly selected from larger pollinator collections detailed in https://doi.org/10.5061/dryad.zpc866t7g.

Headings:
Sample			Unique sample number for each individual in the study. These sample numbers match those found in PollinatorComm2015_2016_publish.csv
Transect.ID		Unique transect ID that indicates the site, year of collection, transect (A B C or D), visit number to the site, and net or pan trap collected. All individuals collected from the same transect share a transect ID.
Individual.ID 		Includes the transect ID with a sequential number for each sample. These unique IDs match those found in PollinatorComm2015_2016_publish.csv
Year 			Year of sample collection, either 2015 or 2016.
Site			Unique site code for each of the 14 sites (see Appendix S1: Table S1 for site names, dates visited, location, etc). BB site is excluded from the analysies in this mansucript due to missing local habitat data for one visit.
VisitNumber		Visit number to the site, each site was visited twice for this study. The first and second visit are indicated by "2" and "3" in 2015 and “1” and “2” in 2016, respectively. A pilot visit was conducted in 2015 but yielded few pollinators since the squash flowers were not blooming at that time, and therefore were excluded from all analyses.
Lat			Latitude coordinate for the Site
Long			Longitude coordinate for the Site
Collected		Categorical variable for the type of flowers that the bee was collected from, either “squash” flowers or a mixture of native and weed “flowers” along the hedgerows of the squash field.
Transect		Pollinators were sampled along 4 randomly placed transects at each site, and labeled “A”, “B”, “C”, or “D”. Transects were sampled via hand nets in alphabetical order. Transects A, B, and C were always within the squash field, and transect D was always located within a nearby hedgerow that often contained native and weedy plant species.
NetorPan		Categorical variable indicating the method of collection for each bee sample.
			Net = hand net
			Pan = pan trap
DateCollected		Date the sample was collected.
Type			Either “APIS” to indicate a managed, Apis mellifera sample, or “NON” to indicate a native bee and non-Apis sample.
Species			Species identity code (code determined by the first two letters of the genus and species names)
			APME = Apis mellifera
			BOIM = Bombus impatiens
			LASI = Lasioglossum spp.
			PEPR = Eucera pruinosa (species was previously named Peponapis pruinosa)
Genus			Genus identity code (abbreviated code from genus name)
			APIS = Apis mellifera
			BOMB = Bombus impatiens
			LASI = Lasioglossum spp.
			PEPO = Eucera pruinosa (species was previously named Peponapis pruinosa)
Sex			Sex of the bee sampled (most samples were female, but male E. pruinosa were present on flowers. Prioritized females for testing for viral presence)
			F = female
			M = male
RNA conc.		RNA concentration (ng/µl) obtained from the Broad Range Qubit kit after RNA extraction for each sample. Samples with “too high” had > 600 ng/µl, and samples with “Sample too low” had < 5 ng/µl, both of which were outside of the range of detection for the kit used. All samples included in the data set amplified the 18S rRNA gene used as a control, indicating sufficient RNA for detection of viral RNA.
DWV			Samples found to be positive (1) or negative (0) for the presence of the deformed wing virus (DWV).
BQCV			Samples found to be positive (1) or negative (0) for the presence of the black queen cell virus (BQCV).
SBV			Samples found to be positive (1) or negative (0) for the presence of the sacbrood virus (SBV).
Richness		Species richness per site (both visits combined) 
EstRich_338		Estimated species richness for each site at the average number of individuals (338 pollinators) collected at any site based on the rarefaction curve.
Abundance		Total number of pollinators collected at each site (both visits combined)
Diversity		Simpson diversity index (1-D) at each site (both visits combined)
APME_ABUND		Number of Apis mellifera samples collected at each site (both visits combined)
BOIM_ABUND		Number of Bombus impatiens samples collected at each site (both visits combined)
LASI_ABUND		Number of Lasioglossum spp. samples collected at each site (both visits combined)
PEPR_ABUND		Number of Eucera pruinosa samples collected at each site (both visits combined). This species was previously named Peponapis pruinosa.
Nat.1000		The proportion of natural area, including forest, wetland, and grassland lancover types, within a 1000 m radius of the center of each squash field where pollinator collection occurred. Landcover types were taken from the USDA cropland datalayer from 2015 and 2016, respecitively (https://nassgeodata.gmu.edu/CropScape/). Data is per field site.
Land_Richness_1000	The number of different landcover types found within a 1000 m radius of the center of the squash field sampled. Landcover types were taken from the USDA cropland datalayer from 2015 and 2016, respecitively (https://nassgeodata.gmu.edu/CropScape/). Data is per field site.
floralrichness		The total number of florwering plant species found per site. Floral richness was sampled in 1 m-squared plots every 10 meters along each of the four transects per site (n=48 plots).
totaldensity		The total number of open flowers of all plant species found per site and divided by the total area sampled per site (per m-squared). Total floral density was sampled in 1 m-squared plots every 10 meters along each of the four transects per site (n=48 plots).

*For raw Nat.1000, Landscape_Richness_1000, floralrichness, and totaldensity data and calculations, contact mlfearon@umich.edu


R script for path models 
The code used to run the main and species-specific path analyses presented in the manuscript. This code produces the information presented in Figures 1-4, and all tables in Appendix S2 snd Appendix S3.
BeeVirusHabitat+CommunityPathModels_29Sept2022.R


File name: Fig 4_indirect and total effect calculations from main path model.xlsx
Calculations of relative strength of direct, indirect, and total effects from the significant scale standardized pathway coefficients in the main path model.
These calculations were done in Microsoft Excel based on the scale standardized coefficients produced from the main path model (Appendix S2 Table S3).