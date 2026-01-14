########################################################################################################################
#reading in functions and libraries
library(data.table)
library(ncdf4)
#source('C:\\Users\\arik\\Documents\\GitHub\\TCFD_ISIMIP_ingestion_and_reanalysis\\TCFD_intersections_code_functions_v5.r')
source('C:\\Users\\arikt\\Documents\\GitHub\\TCFD_ISIMIP_ingestion_and_reanalysis\\TCFD_intersections_code_functions_v5.r')
#source('~\\TCFD_ISIMIP_ingestion_and_reanalysis\\TCFD_intersections_code_functions_v5.r')
#thisWD = 'J:\\Cai_data\\TCFD\\'	# '~//'
thisWD = 'C:\\Cai_data\\TCFD\\'	# '~//'

########################################################################################################################
# paths and values that are unlikely to change between runs
thisDate = Sys.Date()
hazardFolder = paste0(thisWD, 'ProcessedNCs\\')
ncpathRiverFloods = paste0(thisWD, 'Flash Floods\\')
floodMapTiffLoc = paste0(thisWD, 'CurrentFloodHazard')
waterMaskLoc = paste0(thisWD, 'CurrentFloodHazard\\LandMaskArray.rds')
ncpathDEM = paste0(thisWD, 'SeaLevelRise\\globalDEM\\srtm30plus_v11_land.nc')
ncpathSeaLevelRise = paste0(thisWD, 'SeaLevelRise\\sea_level_change')
locationFootprint = 3		# how big is the footprint of the location of interest? in number of 'boxes' to search to the left and right (so 0 is equal to 1 km^2, 1 is 3x3=9 km^2, 2 is 5x5=25 km^2, 3 is 7x7=49, 4 is 9x9=81, etc.
appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
aggScoreExceptionsValues = c(0, 38, 30)
#startTime = proc.time() only used for testing run times

########################################################################################################################
# new customer data
# !!!!!!!!!!!!!!!!!!!!!!!!
# this is the only section that needs to be updated between runs
# !!!!!!!!!!!!!!!!!!!!!!!!

# 	# Nuveen - Aus_Viticulture
userName = 'Nuveen - Aus_Viticulture'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenAusViticulture_Oct2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenAusViticulture_Oct2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\NuveenAusViticulture_Oct2023\\NuveenAusViticulture_temp_precip_hazards_oct_17.csv'
waterOnly = FALSE

# 	# Nuveen - LATAM Table Grapes
userName = 'Nuveen - LATAM Table Grapes'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenLATAMTableGrapes_Oct2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenLATAMTableGrapes_Oct2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# Nuveen - LATAM Cherries
userName = 'Nuveen - LATAM Cherries'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenLATAMCherries_Oct2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenLATAMCherries_Oct2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# Nuveen - Spanish Almonds
userName = 'Nuveen - Spanish Almonds'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenSpanishAlmonds_Oct2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenSpanishAlmonds_Oct2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# NuveenLATAM
userName = 'NuveenLATAM'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenLATAM_Oct2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenLATAM_Oct2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# Nuveen_Almond
userName = 'Nuveen - Almond'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenAlmonds_Aug2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_Almond_Aug2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# Nuveen_Pistachio
userName = 'Nuveen - Pistachio'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenPistachio_Jun2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_Pistachio_Jun2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# ASR
userName = 'ASR'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-ASR-Jun2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - Hazard_Weights.csv'))				
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_sep_28.csv'
waterOnly = FALSE


# 	# Nuveen FL
userName = 'Nuveen'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Nuveen_May2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_May2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# ITC for June 2023
userName = 'ITC_Rice'	#ITC_Buffalo, ITC_VegOil, ITC_Wood, ITC_Spice, ITC_FruitAndVeg, ITC_Shrimp, ITC_Potato, ITC_Rice
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Jun2023\\')
customerTable_redundant = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-', userName, '_June2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
customerTable = customerTable_redundant[!duplicated(customerTable_redundant$Location),]

hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# EQCap for July 2023
userName = "EQ-Capital"
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Jul2023\\')
customerTable_redundant = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-', userName, '_Jul2023.csv')) #'HMClause_locations_allCucurbit.csv'
customerTable = customerTable_redundant[!duplicated(customerTable_redundant$Location),]

hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
appendedHazardFileLoc = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Jul2023\\EQ-Cap_temp_precip_hazards_july_12.csv')
waterOnly = FALSE


# 	# Nuveen FL
userName = 'Nuveen'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Nuveen_Jul2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen-HendryCoFlorida.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE


# 	# EQ Cap - new format
userName = 'EQ-Capital'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Sep2023\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQ-Capital_Aug2023 - Correct Data.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQ-Cap_temp_precip_hazards_sep_6.csv')
waterOnly = FALSE


# 	# Corbion - new format
userName = 'Corbion'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\'

customerTable = fread(paste0(customerFolder, 'Corbion_Locations_Sept2023 - New Format.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - Hazard_Weights.csv'))				
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\Corbion_temp_precip_hazards_sep_28.csv'
waterOnly = FALSE

# 	# Corbion - old format
userName = 'Corbion'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\'

customerTable = fread(paste0(customerFolder, 'Corbion_Locations_Nov2022 - Sheet1_newOldVersion.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_Corbion_Nov2022 - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_Corbion_Nov2022 - Hazard Scores.csv'))				
appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\Corbion_temp_precip_hazards_sep_28_old_method.csv'
aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
aggScoreExceptionsValues = c(0, 38, 30)
waterOnly = FALSE


# 	# ITC all locs for November 2023
userName = 'ITC_allLocs'	#ITC_Buffalo, ITC_VegOil, ITC_Wood, ITC_Spice, ITC_FruitAndVeg, ITC_Shrimp, ITC_Potato, ITC_Rice
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Nov2023\\')
customerTable_redundant = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-', userName, '_Nov2023.csv')) #'HMClause_locations_allCucurbit.csv'
customerTable = customerTable_redundant[!duplicated(customerTable_redundant$Location),]

hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Weights.csv'))				
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ASR_May2023\\ASR_temp_precip_hazards_may_17.csv'
waterOnly = TRUE

# 	# Accenture - new format
userName = 'Accenture'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Nov2023\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Accenture-Nov2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'Cimpress_temp_precip_hazards_nov_14.csv')
waterOnly = TRUE

# 	# DIU - new format
userName = 'DIU'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Nov2023\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-DIU-Nov2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'DIU_temp_precip_hazards_nov_30.csv')
waterOnly = FALSE


# 	# EQ Cap - new format
userName = 'EQ-Capital'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Dec2023\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQCap-Dec2023.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQ_CAP_DEC_2023_temp_precip_hazards_dec_11.csv')
waterOnly = FALSE

# 	# EQ Cap - new format
userName = 'EQ-Capital'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '_Jan2024\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQCap-JAN2024.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQ_CAP_JAN_2024_temp_precip_hazards_jan_2.csv')
waterOnly = FALSE

# 	# EQ Cap - new format
userName = 'EQ-CAP-JAN-2024-2-OFFICES'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\', userName, '\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQCap-JAN2024.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQ_CAP_JAN_2024_2_OFFICES_temp_precip_hazards_jan_24.csv')
waterOnly = FALSE

# 	# EQ Cap - new format
userName = 'EQCap-31JAN2024'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\EQ-CAP-31JAN-2024\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQCap-31JAN2024.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQCap-31JAN2024_temp_precip_hazards_feb_6.csv')
waterOnly = FALSE

# 	# EQ Cap - new format
userName = 'EQ-CAP-JUN-2024'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\EQ-CAP-Jun2024\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-EQ-CAP-Jun2024.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Definitions - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard Scores - Hazard_Tables - Hazard Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights - hazard_weights.csv'))				
appendedHazardFileLoc =  paste0(customerFolder, 'EQCap-13JUN2024_temp_precip_hazards.csv')
waterOnly = FALSE


 	# Nuveen Timberland - new format
userName = 'Nuveen-Timberland'	
customerFolder = paste0('C:\\Cai_data\\TCFD\\locations\\NuveenTimberland_Apr2025\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenTimberland_Apr2025.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights.csv'))
#appendedHazardFileLoc =  paste0(customerFolder, 'EQCap-13JUN2024_temp_precip_hazards.csv')
addFloods = FALSE
addThirdPartyHazards = FALSE


 	# ASR - data updated 2025
userName = 'ASR-2025'	
customerFolder = paste0('J:\\Cai_data\\TCFD\\locations\\ASR_Mar2025\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-ASR_Mar2025.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights.csv'))
appendedHazardFileLoc =  paste0(customerFolder, 'ASR-APR-25-2025_precip_temp_hazards.csv')
addFloods = TRUE
addThirdPartyHazards = TRUE


 	# Nuveen Timberland - add local focus
userName = 'Nuveen-Timberland'	
customerFolder = paste0('C:\\Cai_data\\TCFD\\locations\\NuveenTimberland_Apr2025\\')

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-NuveenTimberland_Apr2025.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Scores.csv'))				
hazardWeighting = fread(paste0(customerFolder, 'Hazard_Weights.csv'))
#appendedHazardFileLoc =  paste0(customerFolder, 'EQCap-13JUN2024_temp_precip_hazards.csv')
addFloods = FALSE
addThirdPartyHazards = FALSE



########################################################################################################################
# functions to run analysis


f_mainExposureData(
	customerTable = customerTable,
	customerFolder = customerFolder,
	hazardFolder = hazardFolder,
	hazardTable = hazardTable,
	relHazScores = relHazScores,
	appendedHazardNames = appendedHazardNames,
	thisDate = thisDate)

if(addFloods)	{# previously waterOnly
	f_riverFloods(
		customerTable = customerTable,
		customerFolder = customerFolder,
		userName = userName,
		locationFootprint = locationFootprint,
		ncpath = ncpathRiverFloods,
		fileLoc = floodMapTiffLoc,
		waterMaskLoc = waterMaskLoc,
		thisDate = thisDate)

	f_seaLevelRise(
		customerTable = customerTable,
		customerFolder = customerFolder,
		userName = userName,
		ncpathDEM = ncpathDEM,
		ncpathSeaLevelRise = ncpathSeaLevelRise,
		thisDate = thisDate)
}

f_hazardAggregation(
	customerTable = customerTable,
	customerFolder = customerFolder,
	thisDate = thisDate,
	userName = userName,
	addThirdPartyHazards = addThirdPartyHazards,
	appendedHazardFileLoc = appendedHazardFileLoc,
	locationFootprint = locationFootprint,
	hazardWeighting = hazardWeighting)



###################################
# end of exposure; subsequent lines are not necessary to run new exposure data
###################################








rootmoistWRIbasedDroughtSeverityv2_processed
rootmoist_WRIbasedDroughtSeverityv2_processed
rootmoist_WRIbasedDroughtSeverityv2_processed.nc



################################ 
#### just checkin data
dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_', userName, '_', thisDate, '.csv'))
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_Nuveen_2023-05-30.csv'))
#customerFolder = "J:\\Cai_data\\TCFD\\locations\\ITC_Buffalo_Jun2023\\"
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_Nuveen - Romania_2023-04-30.csv'))
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Rice_2023-06-13.csv'))
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Potato_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Potato_2023-06-13.csv'))
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Spice_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Spice_2023-06-13.csv'))
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_FruitAndVeg_Jun2023\\'
dataOutput_old = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_FruitAndVeg_2023-06-13.csv'))
dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_FruitAndVeg_2023-06-13_fixedName.csv'))

customerFolder = "J:\\Cai_data\\TCFD\\locations\\NuveenLATAMCherries_Oct2023\\"
dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_Nuveen - LATAM Cherries_2023-10-17.csv'))

#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Wood_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Wood_2023-06-13.csv'))
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_VegOil_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_VegOil_2023-06-13.csv'))
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Shrimp_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Shrimp_2023-06-13.csv'))
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Buffalo_Jun2023\\'
#dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_ITC_Buffalo_2023-06-13.csv'))
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\NuveenPistachio_Jun2023\\'
dataOutput = fread(paste0(customerFolder, 'processedOutputForAllHazards_Nuveen - Pistachio_2023-06-26.csv'))
dataOutput = fread(paste0(customerFolder,"processedOutputForAllHazards_EQ-CAP-Jun2024_2024-06-13.csv"))


summary(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value - subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_25th)
summary(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value - subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_75th)
summary(subset(dataOutput, Hazard != 'Coastal Flood')$Raw_Hazard_Value - subset(dataOutput, Hazard != 'Coastal Flood')$Raw_Hazard_Value_25th)
summary(subset(dataOutput, Hazard != 'Coastal Flood')$Raw_Hazard_Value - subset(dataOutput, Hazard != 'Coastal Flood')$Raw_Hazard_Value_75th)
#dataOutput[which(dataOutput$Hazard == 'Coastal Flood'), ]$Raw_Hazard_Value_25th[which(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value > subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_25th)] = 
#	subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value[which(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value > subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_25th)]
#dataOutput[which(dataOutput$Hazard == 'Coastal Flood'), ]$Raw_Hazard_Value_75th[which(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value < subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_75th)] = 
#	subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value[which(subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value < subset(dataOutput, Hazard == 'Coastal Flood')$Raw_Hazard_Value_75th)]


summary(subset(dataOutput, Decade == 2090 & Scenario == "2. Middle of the Road" & Hazard == 'Water Stress'))
summary(subset(dataOutput, Decade == 2090 & Scenario == "3. High Emissions" & Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score'))
summary(subset(dataOutput, Decade == 2090 & Scenario == "3. High Emissions" & Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Weighted Aggregate Score'))
summary(subset(dataOutput, Decade == 2090 & Scenario == "2. Middle of the Road" & Hazard == 'Aggregate Climate Score'))
summary(subset(dataOutput, Decade == 2090 & Scenario == "1. Low Emissions" & Hazard == 'Aggregate Climate Score'))
par(mfrow = c(3,3))
thisLoc = 2
thisScen = unique(dataOutput$Scenario)[3]
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Water Stress' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Hurricanes' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Drought' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'River Flood (Local)' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Coastal Flood' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Water Availability' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Intense Precipitation' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Extreme Heat' & Location == customerTable$Location[thisLoc])$Percentile)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Extreme Cold' & Location == customerTable$Location[thisLoc])$Percentile)	
#plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Water Demand' & Location == customerTable$Location[thisLoc])$Percentile)	

plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Water Stress' & Location == customerTable$Location[thisLoc])$Raw_Hazard_Value)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Water Availability' & Location == customerTable$Location[thisLoc])$Raw_Hazard_Value)	
plot(subset(dataOutput, Scenario == 'Middle of the Road' & Hazard == 'Coastal Flood' & Location == customerTable$Location[1])$Raw_Hazard_Value)	
plot(subset(dataOutput, Scenario == thisScen & Hazard == 'Aggregate Climate Score' & Location %in% customerTable$Location)$Percentile)	
	
proc.time() - startTime




# identifying the top n worst / best performing
totalTrend = NULL
scarcityTrend = NULL
availabilityTrend = NULL
dataTrends = data.frame(User = NA, Location = NA, Hazard = NA, Hazard_Measure = NA, Raw_Hazard_Current_Value = NA, Raw_Hazard_Trend_to_2050s = NA, Percentile_Current_Value = NA, Percentile_Trend_to_2050s = NA)
for(thisLoc in unique(dataOutput$Location))	{
	for(thisHazard in unique(dataOutput$Hazard))	{
		dataSub = subset(dataOutput, Location == thisLoc & Scenario == "2. Middle of the Road" & Hazard == thisHazard)
	
		for(thisHazardMeasure in unique(dataSub$Hazard_Measure))	{
			if(!(thisHazardMeasure %in% c("Aggregate Score", "Weighted Aggregate Score"))){
				dataTrends = rbind(dataTrends,
							data.frame(
								User = userName,
								Location = thisLoc,
								Hazard = thisHazard,
								Hazard_Measure = thisHazardMeasure,
								Raw_Hazard_Current_Value = mean(subset(dataSub, Decade %in% c('2010','2020') & Hazard_Measure == thisHazardMeasure)$Raw_Hazard_Value),
								Raw_Hazard_Trend_to_2050s = mean(subset(dataSub, Decade %in% c('2040','2050','2060') & Hazard_Measure == thisHazardMeasure)$Raw_Hazard_Value) - 
									mean(subset(dataSub, Decade %in% c('2010', '2020') & Hazard_Measure == thisHazardMeasure)$Raw_Hazard_Value), 
								Percentile_Current_Value = mean(subset(dataSub, Decade %in% c('2010','2020') & Hazard_Measure == thisHazardMeasure)$Percentile), 
								Percentile_Trend_to_2050s = mean(subset(dataSub, Decade %in% c('2040','2050','2060') & Hazard_Measure == thisHazardMeasure)$Percentile) - 
									mean(subset(dataSub, Decade %in% c('2010', '2020') & Hazard_Measure == thisHazardMeasure)$Percentile)
								)
							)
			}
		}
	}
	totalTrend = c(totalTrend, mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Aggregate Climate Score" & Decade %in% c("2010", "2020"))$Percentile) - 
		mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Aggregate Climate Score" & Decade %in% c("2040", "2050", "2060"))$Percentile))
	scarcityTrend = c(scarcityTrend, mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Water Availability" & Hazard_Measure == "Aggregate Score" & Decade %in% c("2010", "2020"))$Percentile) - 
		mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Water Availability" & Hazard_Measure == "Aggregate Score" & Decade %in% c("2040", "2050", "2060"))$Percentile))
	availabilityTrend = c(availabilityTrend, mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Water Scarcity" & Hazard_Measure == "Aggregate Score" & Decade %in% c("2010", "2020"))$Percentile) - 
		mean(subset(dataOutput, Location == thisLoc & Scenario == "Middle of the Road" & Hazard == "Water Scarcity" & Hazard_Measure == "Aggregate Score" & Decade %in% c("2040", "2050", "2060"))$Percentile))
}	
fwrite(dataTrends, paste0(customerFolder, 'currentAnd2050sAllLocs.csv'))



## reformatting for a table view
dataOutputReformatted = data.table(Location = NA, Region = NA, Subregion = NA, Scenario = NA, Variable = NA, 
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend = NA)
for(thisLoc in unique(dataOutput$Location))	{
	for(thisScen in unique(dataOutput$Scenario))	{
		for(thisVariable in unique(dataOutput$Hazard_Measure))	{
			if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
				theseOutputs = subset(dataOutput, Location == thisLoc & Scenario == thisScen & Hazard_Measure == thisVariable)
				dataOutputReformatted = rbind(dataOutputReformatted, 
					data.table(
						Location = thisLoc,
						Region = theseOutputs$Region[1],
						Subregion = theseOutputs$Subregion[1],
						Scenario = thisScen,
						Variable = thisVariable,
						D2010s= subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value, 
						D2020s= subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value, 
						D2030s= subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value, 
						D2040s= subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value, 
						D2050s= subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value, 
						D2060s= subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value, 
						D2070s= subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value, 
						D2080s= subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value, 
						D2090s= subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value,
						Trend = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength)
					)
			}
		}
	}
}
fwrite(dataOutputReformatted, paste0(customerFolder, 'allLocsRawValues.csv'))
				
				
dataOutputSuperSimple = data.table(Variable = NA, 
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend = NA)
thisScen = '2. Middle of the Road'
for(thisVariable in unique(dataOutputReformatted$Variable))	{
	theseOutputs = subset(dataOutputReformatted, Scenario == thisScen & Variable == thisVariable)
	dataOutputSuperSimple = rbind(dataOutputSuperSimple, 
		data.table(
			Variable = thisVariable,
			D2010s= mean(theseOutputs$D2010), 
			D2020s= mean(theseOutputs$D2020), 
			D2030s= mean(theseOutputs$D2030), 
			D2040s= mean(theseOutputs$D2040), 
			D2050s= mean(theseOutputs$D2050), 
			D2060s= mean(theseOutputs$D2060), 
			D2070s= mean(theseOutputs$D2070), 
			D2080s= mean(theseOutputs$D2080), 
			D2090s= mean(theseOutputs$D2090),
			Trend = mean(theseOutputs$Trend))
	)
}
fwrite(dataOutputSuperSimple, paste0(customerFolder, 'portfolioAvgValues.csv'))

dataOutputSuperSimple = data.table(Variable = NA, Region = NA,
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend = NA)
thisScen = '2. Middle of the Road'
for(thisVariable in unique(dataOutputReformatted$Variable))	{
	for(thisRegion in unique(dataOutputReformatted$Region))	{
		if(!is.na(thisRegion))	{
			theseOutputs = subset(dataOutputReformatted, Scenario == thisScen & Variable == thisVariable & Region == thisRegion)
			dataOutputSuperSimple = rbind(dataOutputSuperSimple, 
				data.table(
					Variable = thisVariable,
					Region = thisRegion,
					D2010s= mean(theseOutputs$D2010), 
					D2020s= mean(theseOutputs$D2020), 
					D2030s= mean(theseOutputs$D2030), 
					D2040s= mean(theseOutputs$D2040), 
					D2050s= mean(theseOutputs$D2050), 
					D2060s= mean(theseOutputs$D2060), 
					D2070s= mean(theseOutputs$D2070), 
					D2080s= mean(theseOutputs$D2080), 
					D2090s= mean(theseOutputs$D2090),
					Trend = mean(theseOutputs$Trend))
			)
		}
	}
}
fwrite(dataOutputSuperSimple, paste0(customerFolder, 'portfolioRegionAvgValues.csv'))

dataOutputSuperSimple = data.table(Variable = NA, Subregion = NA,
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend = NA)
thisScen = '2. Middle of the Road'
for(thisVariable in unique(dataOutputReformatted$Variable))	{
	for(thisSubregion in unique(dataOutputReformatted$Subregion))	{
		if(!is.na(thisSubregion))	{
			theseOutputs = subset(dataOutputReformatted, Scenario == thisScen & Variable == thisVariable & Subregion == thisSubregion)
			dataOutputSuperSimple = rbind(dataOutputSuperSimple, 
				data.table(
					Variable = thisVariable,
					Subregion = thisSubregion,
					D2010s= mean(theseOutputs$D2010), 
					D2020s= mean(theseOutputs$D2020), 
					D2030s= mean(theseOutputs$D2030), 
					D2040s= mean(theseOutputs$D2040), 
					D2050s= mean(theseOutputs$D2050), 
					D2060s= mean(theseOutputs$D2060), 
					D2070s= mean(theseOutputs$D2070), 
					D2080s= mean(theseOutputs$D2080), 
					D2090s= mean(theseOutputs$D2090),
					Trend = mean(theseOutputs$Trend))
			)
		}
	}
}
fwrite(dataOutputSuperSimple, paste0(customerFolder, 'portfolioSubregionAvgValues.csv'))





dataSummaryRegion =    data.frame(Region = NA, Subregion = NA, Measure = NA, values_2010s = NA, changeBy_2020s = NA, changeBy_2030s = NA, changeBy_2040s = NA, changeBy_2050s = NA, changeBy_2060s = NA, changeBy_2070s = NA, changeBy_2080s = NA, changeBy_2090s = NA)
dataSummarySubregion = data.frame(Region = NA, Subregion = NA, Measure = NA, values_2010s = NA, changeBy_2020s = NA, changeBy_2030s = NA, changeBy_2040s = NA, changeBy_2050s = NA, changeBy_2060s = NA, changeBy_2070s = NA, changeBy_2080s = NA, changeBy_2090s = NA)
for(thisMeasure in unique(dataOutput$Hazard_Measure)){
	if(!(thisMeasure %in% c('Weighted Aggregate Score', 'Aggregate Score')))	{
		measureOutput = subset(dataOutput, Scenario == "2. Middle of the Road" & Hazard_Measure == thisMeasure & Long_Term_Trend_Significance < 0.1)
		for(cc in unique(measureOutput$Region)){
			countryOutput = subset(measureOutput, Region == cc)
			dataSummaryRegion = rbind(dataSummaryRegion,
				data.frame(
					Region = cc,
					Subregion = "All",
					Measure = thisMeasure,
					values_2010s = mean(subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2020s = mean(subset(countryOutput, Decade == 2020)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2030s = mean(subset(countryOutput, Decade == 2030)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2040s = mean(subset(countryOutput, Decade == 2040)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2050s = mean(subset(countryOutput, Decade == 2050)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2060s = mean(subset(countryOutput, Decade == 2060)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2070s = mean(subset(countryOutput, Decade == 2070)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2080s = mean(subset(countryOutput, Decade == 2080)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2090s = mean(subset(countryOutput, Decade == 2090)$Raw_Hazard_Value - subset(countryOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE)
			))

			for(ss in unique(countryOutput$Subregion))	{
				subregOutput = subset(countryOutput, Subregion == ss)	
				dataSummarySubregion = rbind(dataSummarySubregion,
					data.frame(
						Region = cc,
						Subregion = ss,
						Measure = thisMeasure,
					values_2010s = mean(subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2020s = mean(subset(subregOutput, Decade == 2020)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2030s = mean(subset(subregOutput, Decade == 2030)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2040s = mean(subset(subregOutput, Decade == 2040)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2050s = mean(subset(subregOutput, Decade == 2050)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2060s = mean(subset(subregOutput, Decade == 2060)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2070s = mean(subset(subregOutput, Decade == 2070)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2080s = mean(subset(subregOutput, Decade == 2080)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE),
					changeBy_2090s = mean(subset(subregOutput, Decade == 2090)$Raw_Hazard_Value - subset(subregOutput, Decade == 2010)$Raw_Hazard_Value, na.rm=TRUE)
					))
			}
		}
	}
}
#fwrite(dataSummaryRegion,  paste0(customerFolder, 'regionAllDifferences.csv'))
#fwrite(dataSummarySubregion,  paste0(customerFolder, 'subregionAllDifferences.csv'))

fwrite(dataSummaryRegion,  paste0(customerFolder, 'regionSignifDifferences.csv'))
fwrite(dataSummarySubregion,  paste0(customerFolder, 'subregionSignifDifferences.csv'))
#fwrite(dataSummary[which(abs(dataSummary[,-c(1:3)] , 'Differences.csv')




## reformatting for a table view but with all regional changes
dataOutputReformatted = data.table(Location = NA, Country = NA, Region = NA, Scenario = NA, Variable = NA, 
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend_Strength = NA, Trend_Significance = NA)
for(thisLoc in unique(dataOutput$Location))	{
	for(thisScen in unique(dataOutput$Scenario))	{
		for(varNum in 1:length(unique(dataOutput$Hazard_Measure)))	{
			thisVariable = unique(dataOutput$Hazard_Measure)[varNum]
			if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
				theseOutputs = subset(dataOutput, Location == thisLoc & Scenario == thisScen & Hazard_Measure == thisVariable)
				if(thisScen == "2. Middle of the Road") {# &
#					theseOutputs$Long_Term_Trend_Significance[1] < 0.1 &
#					abs(theseOutputs$Long_Term_Trend_Strength[1]) < impactfulDiffs[varNum])
#					abs(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value - subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value) > impactfulDiffs[varNum])	{
				
					dataOutputReformatted = rbind(dataOutputReformatted, 
						data.table(
							Location = thisLoc,
							Country = theseOutputs$Country,
							Region = theseOutputs$Region,
							Scenario = thisScen,
							Variable = thisVariable,
							D2010s= subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value, 
							D2020s= subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value, 
							D2030s= subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value, 
							D2040s= subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value, 
							D2050s= subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value, 
							D2060s= subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value, 
							D2070s= subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value, 
							D2080s= subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value, 
							D2090s= subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value,
							Trend_Strength = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength,
							Trend_Significance = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Significance)
						)
				}
			}
		}
	}
}
fwrite(dataOutputReformatted, paste0(customerFolder, 'allLocsRawValues.csv'))
#fwrite(dataOutputReformatted, paste0(customerFolder, 'signifLocsRawValues.csv'))








########### end of standard reports #################










# bullshit fucking people not able to communicate what they want then waiting weeks to respond then giving thoughtless responses for fucks sake such a goddam waste of time

# fucking splitting on countries and regions bc they were too dumb to know how to when they passed the fucking bs db over

fixedData = data.frame(do.call("rbind", strsplit(as.character(dataOutput$Region), ", ", fixed = TRUE, useBytes = TRUE)))
colnames(fixedData) = c("Region_fixed", "Country_fixed")
fixedOutput = cbind(dataOutput, fixedData)
fixedOutput = subset(fixedOutput, Hazard == 'Water Availability')

impactfulDiffs = c(20, 0.5, 0.3, 1000, 1000, 1000, 1000, 1000, 1000, 1000)
#impactfulDiffs = c(1, 0.1, (1.5191648 * 10^-6) / 20, 1000, 1000, 1000, 1000, 1000, 1000, 1000)

## reformatting for a table view but with all regional changes
dataOutputReformatted = data.table(Location = NA, Country = NA, Region = NA, Scenario = NA, Variable = NA, 
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend_Strength = NA, Trend_Significance = NA)
for(thisLoc in unique(fixedOutput$Location))	{
	for(thisScen in unique(fixedOutput$Scenario))	{
		for(varNum in 1:length(unique(fixedOutput$Hazard_Measure)))	{
			thisVariable = unique(fixedOutput$Hazard_Measure)[varNum]
			if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
				theseOutputs = subset(fixedOutput, Location == thisLoc & Scenario == thisScen & Hazard_Measure == thisVariable)
				if(thisScen == "Middle of the Road") {# &
#					theseOutputs$Long_Term_Trend_Significance[1] < 0.1 &
#					abs(theseOutputs$Long_Term_Trend_Strength[1]) < impactfulDiffs[varNum])
#					abs(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value - subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value) > impactfulDiffs[varNum])	{
				
					dataOutputReformatted = rbind(dataOutputReformatted, 
						data.table(
							Location = thisLoc,
							Country = theseOutputs$Country_fixed[1],
							Region = theseOutputs$Region_fixed[1],
							Scenario = thisScen,
							Variable = thisVariable,
							D2010s= subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value, 
							D2020s= subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value, 
							D2030s= subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value, 
							D2040s= subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value, 
							D2050s= subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value, 
							D2060s= subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value, 
							D2070s= subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value, 
							D2080s= subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value, 
							D2090s= subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value,
							Trend_Strength = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength,
							Trend_Significance = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Significance)
						)
				}
			}
		}
	}
}
fwrite(dataOutputReformatted, paste0(customerFolder, 'allLocsRawValues.csv'))
#fwrite(dataOutputReformatted, paste0(customerFolder, 'signifLocsRawValues.csv'))




## reformatting for a table view but with only statistically significant changes
dataOutputReformatted = data.table(Country = NA, Region = NA, Scenario = NA, Variable = NA, 
	avg2010s= NA, avg2020s= NA, avg2030s= NA, avg2040s= NA, avg2050s= NA, avg2060s= NA, avg2070s= NA, avg2080s= NA, avg2090s= NA, 
	diff2020s= NA, diff2030s= NA, diff2040s= NA, diff2050s= NA, diff2060s= NA, diff2070s= NA, diff2080s= NA, diff2090s= NA, 
	Trend_Strength = NA, Trend_Significance = NA)
for(thisCountry in unique(fixedOutput$Country_fixed))	{
		for(thisScen in "Middle of the Road")	{ #unique(fixedOutput$Scenario))	{
			for(varNum in 1:length(unique(fixedOutput$Hazard_Measure)))	{
				thisVariable = unique(fixedOutput$Hazard_Measure)[varNum]
				if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
					theseOutputs = subset(fixedOutput, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable)
			
					dataOutputReformatted = rbind(dataOutputReformatted, 
						data.table(
							Country = thisCountry,
							Region = "Avg of All Locations in Country",
							Scenario = thisScen,
							Variable = thisVariable,
							avg2010s= mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							avg2020s= mean(subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value), 
							avg2030s= mean(subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value), 
							avg2040s= mean(subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value), 
							avg2050s= mean(subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value), 
							avg2060s= mean(subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value), 
							avg2070s= mean(subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value), 
							avg2080s= mean(subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value), 
							avg2090s= mean(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value),
							diff2020s= mean(subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2030s= mean(subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2040s= mean(subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2050s= mean(subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2060s= mean(subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2070s= mean(subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2080s= mean(subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2090s= mean(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							Trend_Strength = mean(subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength),
							Trend_Significance = mean(subset(theseOutputs, Decade == 2090)$Decadal_Trend_Significance))
						)

					if(TRUE)	# &
#						any(theseOutputs$Long_Term_Trend_Significance < 0.1) &
#						abs(theseOutputs$Long_Term_Trend_Strength[1]) < impactfulDiffs[varNum])
#						any(abs(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value - subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value) > impactfulDiffs[varNum]))
					{
						
					theseOutputsSignif = subset(theseOutputs, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable)# & Long_Term_Trend_Significance < 0.1)


						for(thisRegion in unique(theseOutputsSignif$Region_fixed))	{
							theseOutputsSignif = subset(fixedOutput, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable & Region_fixed == thisRegion)
							dataOutputReformatted = rbind(dataOutputReformatted, 
								data.table(
									Country = theseOutputsSignif$Country_fixed[1],
									Region = theseOutputsSignif$Region_fixed[1],
									Scenario = thisScen,
									Variable = thisVariable,
									avg2010s= mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									avg2020s= mean(subset(theseOutputsSignif, Decade == 2020)$Raw_Hazard_Value), 
									avg2030s= mean(subset(theseOutputsSignif, Decade == 2030)$Raw_Hazard_Value), 
									avg2040s= mean(subset(theseOutputsSignif, Decade == 2040)$Raw_Hazard_Value), 
									avg2050s= mean(subset(theseOutputsSignif, Decade == 2050)$Raw_Hazard_Value), 
									avg2060s= mean(subset(theseOutputsSignif, Decade == 2060)$Raw_Hazard_Value), 
									avg2070s= mean(subset(theseOutputsSignif, Decade == 2070)$Raw_Hazard_Value), 
									avg2080s= mean(subset(theseOutputsSignif, Decade == 2080)$Raw_Hazard_Value), 
									avg2090s= mean(subset(theseOutputsSignif, Decade == 2090)$Raw_Hazard_Value),
									diff2020s= mean(subset(theseOutputsSignif, Decade == 2020)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2030s= mean(subset(theseOutputsSignif, Decade == 2030)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2040s= mean(subset(theseOutputsSignif, Decade == 2040)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2050s= mean(subset(theseOutputsSignif, Decade == 2050)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2060s= mean(subset(theseOutputsSignif, Decade == 2060)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2070s= mean(subset(theseOutputsSignif, Decade == 2070)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2080s= mean(subset(theseOutputsSignif, Decade == 2080)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2090s= mean(subset(theseOutputsSignif, Decade == 2090)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									Trend_Strength = median(subset(theseOutputsSignif, Decade == 2090)$Decadal_Trend_Strength),
									Trend_Significance = median(subset(theseOutputsSignif, Decade == 2090)$Decadal_Trend_Significance))
								)
						}
					}
				}
			}
		}
}
fwrite(dataOutputReformatted, paste0(customerFolder, 'allRegsRawValues.csv'))
#fwrite(dataOutputReformatted, paste0(customerFolder, 'signifRegsRawValues.csv'))




	# p2
## reformatting for a table view but with only statistically significant changes
dataOutputReformatted = data.table(Location = NA, Country = NA, Region = NA, Scenario = NA, Variable = NA, 
	D2010s= NA, D2020s= NA, D2030s= NA, D2040s= NA, D2050s= NA, D2060s= NA, D2070s= NA, D2080s= NA, D2090s= NA, Trend_Strength = NA, Trend_Significance = NA)
for(thisLoc in unique(fixedOutput$Location))	{
	for(thisScen in unique(fixedOutput$Scenario))	{
		for(varNum in 1:length(unique(fixedOutput$Hazard_Measure)))	{
			thisVariable = unique(fixedOutput$Hazard_Measure)[varNum]
			if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
				theseOutputs = subset(fixedOutput, Location == thisLoc & Scenario == thisScen & Hazard_Measure == thisVariable)
				if(thisScen == "Middle of the Road" &
					theseOutputs$Long_Term_Trend_Significance[1] < 0.1 &
#					abs(theseOutputs$Long_Term_Trend_Strength[1]) < impactfulDiffs[varNum])
					abs(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value - subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value) > impactfulDiffs[varNum])	{
				
					dataOutputReformatted = rbind(dataOutputReformatted, 
						data.table(
							Location = thisLoc,
							Country = theseOutputs$Country_fixed[1],
							Region = theseOutputs$Region_fixed[1],
							Scenario = thisScen,
							Variable = thisVariable,
							D2010s= subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value, 
							D2020s= subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value, 
							D2030s= subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value, 
							D2040s= subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value, 
							D2050s= subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value, 
							D2060s= subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value, 
							D2070s= subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value, 
							D2080s= subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value, 
							D2090s= subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value,
							Trend_Strength = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength,
							Trend_Significance = subset(theseOutputs, Decade == 2090)$Decadal_Trend_Significance)
						)
				}
			}
		}
	}
}
#fwrite(dataOutputReformatted, paste0(customerFolder, 'allLocsRawValues.csv'))
fwrite(dataOutputReformatted, paste0(customerFolder, 'signifLocsRawValues.csv'))




## reformatting for a table view but with only statistically significant changes
dataOutputReformatted = data.table(Country = NA, Region = NA, Scenario = NA, Variable = NA, 
	avg2010s= NA, avg2020s= NA, avg2030s= NA, avg2040s= NA, avg2050s= NA, avg2060s= NA, avg2070s= NA, avg2080s= NA, avg2090s= NA, 
	diff2020s= NA, diff2030s= NA, diff2040s= NA, diff2050s= NA, diff2060s= NA, diff2070s= NA, diff2080s= NA, diff2090s= NA, 
	Trend_Strength = NA, Trend_Significance = NA)
for(thisCountry in unique(fixedOutput$Country_fixed))	{
		for(thisScen in "Middle of the Road")	{ #unique(fixedOutput$Scenario))	{
			for(varNum in 1:length(unique(fixedOutput$Hazard_Measure)))	{
				thisVariable = unique(fixedOutput$Hazard_Measure)[varNum]
				if(!(thisVariable %in% c('Aggregate Score', 'Weighted Aggregate Score')))	{
					theseOutputs = subset(fixedOutput, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable)
			
					dataOutputReformatted = rbind(dataOutputReformatted, 
						data.table(
							Country = thisCountry,
							Region = "Avg of All Locations in Country",
							Scenario = thisScen,
							Variable = thisVariable,
							avg2010s= mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							avg2020s= mean(subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value), 
							avg2030s= mean(subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value), 
							avg2040s= mean(subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value), 
							avg2050s= mean(subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value), 
							avg2060s= mean(subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value), 
							avg2070s= mean(subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value), 
							avg2080s= mean(subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value), 
							avg2090s= mean(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value),
							diff2020s= mean(subset(theseOutputs, Decade == 2020)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2030s= mean(subset(theseOutputs, Decade == 2030)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2040s= mean(subset(theseOutputs, Decade == 2040)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2050s= mean(subset(theseOutputs, Decade == 2050)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2060s= mean(subset(theseOutputs, Decade == 2060)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value),  
							diff2070s= mean(subset(theseOutputs, Decade == 2070)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2080s= mean(subset(theseOutputs, Decade == 2080)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							diff2090s= mean(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value) - mean(subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value), 
							Trend_Strength = mean(subset(theseOutputs, Decade == 2090)$Decadal_Trend_Strength),
							Trend_Significance = mean(subset(theseOutputs, Decade == 2090)$Decadal_Trend_Significance))
						)

					if(TRUE &
						any(theseOutputs$Long_Term_Trend_Significance < 0.1) &
#						abs(theseOutputs$Long_Term_Trend_Strength[1]) < impactfulDiffs[varNum])
						any(abs(subset(theseOutputs, Decade == 2090)$Raw_Hazard_Value - subset(theseOutputs, Decade == 2010)$Raw_Hazard_Value) > impactfulDiffs[varNum]))
					{
						
					theseOutputsSignif = subset(theseOutputs, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable)# & Long_Term_Trend_Significance < 0.1)


						for(thisRegion in unique(theseOutputsSignif$Region_fixed))	{
							theseOutputsSignif = subset(fixedOutput, Country_fixed == thisCountry & Scenario == thisScen & Hazard_Measure == thisVariable & Region_fixed == thisRegion)
							dataOutputReformatted = rbind(dataOutputReformatted, 
								data.table(
									Country = theseOutputsSignif$Country_fixed[1],
									Region = theseOutputsSignif$Region_fixed[1],
									Scenario = thisScen,
									Variable = thisVariable,
									avg2010s= mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									avg2020s= mean(subset(theseOutputsSignif, Decade == 2020)$Raw_Hazard_Value), 
									avg2030s= mean(subset(theseOutputsSignif, Decade == 2030)$Raw_Hazard_Value), 
									avg2040s= mean(subset(theseOutputsSignif, Decade == 2040)$Raw_Hazard_Value), 
									avg2050s= mean(subset(theseOutputsSignif, Decade == 2050)$Raw_Hazard_Value), 
									avg2060s= mean(subset(theseOutputsSignif, Decade == 2060)$Raw_Hazard_Value), 
									avg2070s= mean(subset(theseOutputsSignif, Decade == 2070)$Raw_Hazard_Value), 
									avg2080s= mean(subset(theseOutputsSignif, Decade == 2080)$Raw_Hazard_Value), 
									avg2090s= mean(subset(theseOutputsSignif, Decade == 2090)$Raw_Hazard_Value),
									diff2020s= mean(subset(theseOutputsSignif, Decade == 2020)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2030s= mean(subset(theseOutputsSignif, Decade == 2030)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2040s= mean(subset(theseOutputsSignif, Decade == 2040)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2050s= mean(subset(theseOutputsSignif, Decade == 2050)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2060s= mean(subset(theseOutputsSignif, Decade == 2060)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value),  
									diff2070s= mean(subset(theseOutputsSignif, Decade == 2070)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2080s= mean(subset(theseOutputsSignif, Decade == 2080)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									diff2090s= mean(subset(theseOutputsSignif, Decade == 2090)$Raw_Hazard_Value) - mean(subset(theseOutputsSignif, Decade == 2010)$Raw_Hazard_Value), 
									Trend_Strength = median(subset(theseOutputsSignif, Decade == 2090)$Decadal_Trend_Strength),
									Trend_Significance = median(subset(theseOutputsSignif, Decade == 2090)$Decadal_Trend_Significance))
								)
						}
					}
				}
			}
		}
}
#fwrite(dataOutputReformatted, paste0(customerFolder, 'allRegsRawValues.csv'))
fwrite(dataOutputReformatted, paste0(customerFolder, 'signifRegsRawValues.csv'))























dataSub = subset(dataOutput, Scenario == "Middle of the Road" & Hazard_Measure =='Groundwater Recharge (mm per yr)' & Region =="China")
#dataSub = subset(dataOutput, Scenario == "Middle of the Road" & Hazard_Measure =='Root Zone Soil Moisture (mm per yr)' & Region =="China")
plot(subset(chinaSM, Location == unique(chinaSM$Location)[1])$Trend_Agg)
for(i in unique(chinaSM$Location)){ lines(subset(chinaSM, Location == i)$Trend_Agg)}

trend_avg = matrix(NA, nrow=9, ncol=6)
trend_sum = matrix(NA, nrow=9, ncol=6)
trend_diff = matrix(NA, nrow=9, ncol=6)
trend_signif_diff = matrix(NA, nrow=9, ncol=6)
iter = 0
for(i in unique(dataSub$Decade))	{
	iter = iter + 1
	trend_avg[iter,] = summary(subset(dataSub, Decade == i)$Trend_Agg)
	trend_sum[iter,] = summary(subset(dataSub, Decade == i)$Trend_Agg) * iter
	trend_diff[iter,] = summary(subset(dataSub, Decade == i)$Raw_Hazard_Value - subset(dataSub, Decade == 2010)$Raw_Hazard_Value)
	trend_signif_diff[iter,] = summary(subset(dataSub, Decade == i & Long_Term_Trend_Significance < 0.05)$Raw_Hazard_Value - subset(dataSub, Decade == 2010 & Long_Term_Trend_Significance < 0.05)$Raw_Hazard_Value)
}
trend_avg=as.data.frame(trend_avg)
trend_sum=as.data.frame(trend_sum)
trend_diff=as.data.frame(trend_diff)
trend_signif_diff=as.data.frame(trend_signif_diff)
names(trend_avg) = c('Min','Q1','Med','Mean','Q3','Max')
names(trend_sum) = c('Min','Q1','Med','Mean','Q3','Max')
names(trend_diff) = c('Min','Q1','Med','Mean','Q3','Max')
names(trend_signif_diff) = c('Min','Q1','Med','Mean','Q3','Max')
row.names(trend_avg) = seq(2010,2090,10)
row.names(trend_sum) = seq(2010,2090,10)
row.names(trend_diff) = seq(2010,2090,10)
row.names(trend_signif_diff) = seq(2010,2090,10)
fwrite(trend_sum, 'total_rate_of_change_by_decade_GW.csv')
fwrite(trend_avg, 'avg_rate_of_change_by_decade_GW.csv')
fwrite(trend_diff, 'median_difference_by_decade_GW.csv')
fwrite(trend_signif_diff, 'median_difference_of_significant_differences_by_decade_GW.csv')



dataSummary = data.frame(Country = NA, Region = NA, Measure = NA, the_2010s = NA, the_2020s = NA, the_2030s = NA, the_2040s = NA, the_2050s = NA, the_2060s = NA, the_2070s = NA, the_2080s = NA, the_2090s = NA)
for(thisMeasure in c("Groundwater Recharge (mm per yr)", "Root Zone Soil Moisture (mm per yr)", "Surface Water (annual total in cubic km)")){
	measureOutput = subset(dataOutput, Hazard_Measure == thisMeasure)
	for(cc in unique(measureOutput$Region)){
		countryOutput = subset(measureOutput, Region == cc)
		dataSummary = rbind(dataSummary,
			data.frame(
				Country = cc,
				Region = "All",
				Measure = thisMeasure,
				the_2010s = mean(subset(countryOutput, Decade == 2010)$Trend_Aggregated_For_Looker),
				the_2020s = mean(subset(countryOutput, Decade == 2020)$Trend_Aggregated_For_Looker),
				the_2030s = mean(subset(countryOutput, Decade == 2030)$Trend_Aggregated_For_Looker),
				the_2040s = mean(subset(countryOutput, Decade == 2040)$Trend_Aggregated_For_Looker),
				the_2050s = mean(subset(countryOutput, Decade == 2050)$Trend_Aggregated_For_Looker),
				the_2060s = mean(subset(countryOutput, Decade == 2060)$Trend_Aggregated_For_Looker),
				the_2070s = mean(subset(countryOutput, Decade == 2070)$Trend_Aggregated_For_Looker),
				the_2080s = mean(subset(countryOutput, Decade == 2080)$Trend_Aggregated_For_Looker),
				the_2090s = mean(subset(countryOutput, Decade == 2090)$Trend_Aggregated_For_Looker)
		))

		for(ss in unique(countryOutput$Subregion))	{
			subregOutput = subset(countryOutput, Subregion == ss)	
			dataSummary = rbind(dataSummary,
				data.frame(
					Country = cc,
					Region = ss,
					Measure = thisMeasure,
					the_2010s = mean(subset(subregOutput, Decade == 2010)$Trend_Aggregated_For_Looker),
					the_2020s = mean(subset(subregOutput, Decade == 2020)$Trend_Aggregated_For_Looker),
					the_2030s = mean(subset(subregOutput, Decade == 2030)$Trend_Aggregated_For_Looker),
					the_2040s = mean(subset(subregOutput, Decade == 2040)$Trend_Aggregated_For_Looker),
					the_2050s = mean(subset(subregOutput, Decade == 2050)$Trend_Aggregated_For_Looker),
					the_2060s = mean(subset(subregOutput, Decade == 2060)$Trend_Aggregated_For_Looker),
					the_2070s = mean(subset(subregOutput, Decade == 2070)$Trend_Aggregated_For_Looker),
					the_2080s = mean(subset(subregOutput, Decade == 2080)$Trend_Aggregated_For_Looker),
					the_2090s = mean(subset(subregOutput, Decade == 2090)$Trend_Aggregated_For_Looker)
			))
		}
	}
}
fwrite(dataSummary, 'Trends.csv')










gg = fread('J:\\Cai_data\\TCFD\\locations\\ITC_Dec2022\\processedOutputForAllHazards_ITC - Water_2022-12-07.csv')
for(i in unique(gg$Decade))	{
	print(sum(subset(gg, Decade == i & Scenario == 'Low Emissions' & Hazard_Measure == 'Surface Water (annual total in cubic km)' & Region == 'India')$Trend_Aggregated_For_Looker))
}	

	
	
	
	
# some sanity checking
		
subset(dataOutput, Hazard != 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'Middle of the Road')[1:99,6:12]
subset(dataOutput, Hazard_Measure == 'Aggregate Score' & Scenario == 'Middle of the Road')[1:99,c(2,7:12)]
subset(dataOutput, Hazard_Measure == 'Aggregate Score' & Scenario == 'Low Emissions')[1:99,c(2,7:12)]
subset(dataOutput, Hazard_Measure == 'Aggregate Score' & Scenario == 'High Emissions')[1:99,c(2,7:12)]

summary(subset(dataOutput, Hazard != 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'Low Emissions')$Percentile)
summary(subset(dataOutput, Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'Low Emissions')$Percentile)

summary(subset(dataOutput, Hazard != 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'Middle of the Road')$Percentile)
summary(subset(dataOutput, Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'Middle of the Road')$Percentile)

summary(subset(dataOutput, Hazard != 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'High Emissions')$Percentile)
summary(subset(dataOutput, Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score' & Scenario == 'High Emissions')$Percentile)

for(thisScen in unique(dataOutput$Scenario))	{
	for(thisLoc in unique(dataOutput$Location))	{
		dataSub = subset(dataOutput, Scenario == thisScen & Location == thisLoc & Hazard == 'Aggregate Climate Score' & Hazard_Measure == 'Aggregate Score')
		png(filename=paste0('J:\\Downloads\\', substr(thisLoc, 1, 10), " - ", thisScen, '.png'))
		plot(dataSub$Decade, dataSub$Percentile, ylim = c(1,100), lwd = 3, type = 'l',
			ylab = 'Percentile', xlab = 'Decade', main = paste0(substr(thisLoc, 1, 10), " - ", thisScen))
		for(thisHazard in unique(dataOutput$Hazard))	{
			
			dataSub = subset(dataOutput, Scenario == thisScen & Location == thisLoc & Hazard == thisHazard & Hazard_Measure == 'Aggregate Score')
			lines(dataSub$Decade, dataSub$Percentile, col = 'grey30', lwd = 1.5)
		}
		dev.off()
	}
}

	
	
	
	
	
	










################################
# old runs

	# corbion
#userName = 'Corbion'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\'

#customerTable = fread(paste0(customerFolder, 'Corbion_Locations_Nov2022 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_Corbion_Nov2022 - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Corbion_Nov2022\\Corbion_temp_precip_hazards_mar_16.csv'
#aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#aggScoreExceptionsValues = c(0, 38, 30)

	# ITC
#userName = 'ITC'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Dec2022\\'

#customerTable = fread(paste0(customerFolder, 'ITC_Locations_Dec2022 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_ITC_Dec2022 - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_ITC_Dec2022 - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ITC_Dec2022\\ITC_temp_precip_hazards_dec_6 - ITC_temp_precip_hazards_dec_6.csv'

	# ITC - Water
#userName = 'ITC - Water'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ITC_Dec2022\\'

#customerTable = fread(paste0(customerFolder, 'ITC_Locations_Dec2022 - Water Exposure.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_ITC_Dec2022 - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_ITC_Dec2022 - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ITC_Dec2022\\ITC_temp_precip_hazards_dec_6 - ITC_temp_precip_hazards_dec_6.csv'

	# Driscolls
#userName = 'Driscolls'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Driscolls_Dec2022\\'

#customerTable = fread(paste0(customerFolder, 'Driscolls_Locations_Dec2022 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_Driscolls_Dec2022 - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Driscolls_Dec2022\\Driscolls_temp_precip_hazards_dec_6 - Driscolls_temp_precip_hazards_dec_6.csv'
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Driscolls_Dec2022\\Sales_davesData.csv'


	# HMClause
userName = 'HMClause'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\HMClause_Dec2022\\'

customerTable = fread(paste0(customerFolder, 'HMClause_Locations_Nov2022 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables_HMClause_Dec2022 - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables_HMClause_Dec2022 - Hazard Scores.csv'))				
appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\HMClause_Dec2022\\HMClause_temp_precip_hazards_dec_6 - HMClause_temp_precip_hazards_dec_6.csv'


	# Nuveen - Spain
#userName = 'Nuveen - Spain'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Nuveen_Feb2023\\'

#customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_Jan2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\HMClause_Dec2022\\HMClause_temp_precip_hazards_dec_6 - HMClause_temp_precip_hazards_dec_6.csv'


	# EQ Cap
#userName = 'Sales'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\EQCap_Feb2023\\'
#thisDate = Sys.Date()

#customerTable = fread(paste0(customerFolder, 'EQCap_Locations_Dec2022 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\EQCap_Feb2023\\EQ_cap_temp_precip_hazards_feb_11.csv'


	# suntory - Japan Demo
#userName = 'Kurita'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Sunotory_Feb2023\\'

#customerTable = fread(paste0(customerFolder, 'Suntory - Hops Supplier Locations - For_exposure.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Sunotory_Feb2023\\Kurita_temp_precip_hazards_feb_18.csv'

	# Indigo - South Africa
#userName = 'Indigo'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Indigo_Feb2023\\'

#customerTable = fread(paste0(customerFolder, 'SA - Exposure Locations - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\Indigo_Feb2023\\Indigo_temp_precip_hazards_feb_22.csv'

	# UBS Demo
#userName = 'UBS'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\UBS_Feb2023\\'

#customerTable = fread(paste0(customerFolder, 'UBS Demo - Exposure Locations - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\UBS_Feb2023\\UBS_temp_precip_hazards_feb_26.csv'

	# ICL Demo
#userName = 'ICL'	
#customerFolder = 'J:\\Cai_data\\TCFD\\locations\\ICL-Demo_Mar2023\\'

#customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-ICL_Demo - Sheet1.csv'))
#hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
#relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
#appendedHazardNames = c("River Flood (Local)", 'Coastal Flood', 'Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#appendedHazardFileLoc = 'J:\\Cai_data\\TCFD\\locations\\ICL-Demo_Mar2023\\ICL_temp_precip_hazards_mar_21.csv'
#aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
#aggScoreExceptionsValues = c(0, 38, 30)

	# Nuveen - Georgia
userName = 'Nuveen - Georgia'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Nubeen_Apr2023\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_Apr2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
appendedHazardNames = NA
appendedHazardFileLoc = NA
aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
aggScoreExceptionsValues = c(0, 38, 30)
waterOnly = TRUE

	# Nuveen - Romania
userName = 'Nuveen - Romania'	
customerFolder = 'J:\\Cai_data\\TCFD\\locations\\Nuveen_Apr2023b\\'

customerTable = fread(paste0(customerFolder, 'Customer_Hazards_and_Locations-Nuveen_Apr2023 - Sheet1.csv')) #'HMClause_locations_allCucurbit.csv'
hazardTable = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Definitions.csv'))							# 
relHazScores = fread(paste0(customerFolder, 'Hazard_Tables - Hazard Scores.csv'))				
appendedHazardNames = NA
appendedHazardFileLoc = NA
aggScoreExceptions = c('Extreme Cold', 'Extreme Heat', 'Intense Precipitation')
aggScoreExceptionsValues = c(0, 38, 30)

