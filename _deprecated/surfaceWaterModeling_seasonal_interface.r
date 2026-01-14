


	# example input data
basinATLAS_locationAndFile = 'C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10\\BasinATLAS_v10.gdb'
dataOut_location = 'J:\\Cai_data\\Nuveen\\surfaceWaterData_and_Output\\'
basinName = 'SJF_atOutlet'
gageLonLat = c(-119.697,37.004) 
pathToBasinBoundaryGPKG = paste0(dataOut_location, "watershedBoundaries_", basinName, ".gpkg")
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")
minTempNCDF = 'C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\CA_era5_tmin.nc'
maxTempNCDF = 'C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\CA_era5_tmax.nc'
pptNCDF = 'C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\CA_era5_tp.nc'
historicStreamflowFileLoc = 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=sjf&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-06-16'
climateInputsFileLoc = paste0(dataOut_location, basinName, '_processedClimateData.RData')
pathToWatershedsGPKG =  paste0(dataOut_location, "HydroBASINSdata_", basinName, ".gpkg")




#########################################################################################################
	# step 1
	# delineate a new basin
basinDelineation_f(gageLonLat,
	basinATLAS_locationAndFile,
	dataOut_location,
	basinName)
					# the function being called
				basinDelineation_f = function(
					gageLonLat = c(-79,36),
					basinATLAS_locationAndFile = 'basinATLAS_location.gdb',
					dataOut_location = 'save_file_location',
					basinName = 'inster_basin_or_outlet_name')
#########################################################################################################
	# step 2
	# import and convert historical climate data
climateDataframe =	climateInputConversion_f(
		pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
		pathToWatershedsGPKG = pathToWatershedsGPKG,
		basinName = basinName,
		minTempNCDF = minTempNCDF,
		maxTempNCDF = maxTempNCDF,
		pptNCDF = pptNCDF,
		dataOut_location = dataOut_location)
saveRDS(climateDataframe, climateInputsFileLoc)
					# the function being called
				climateInputConversion_f = function(
					pathToBasinBoundaryGPKG = 'file_location_and_name.gpkg',
					pathToWatershedsGPKG = 'file_location_and_name.gpkg',
					basinName = 'inster_basin_or_outlet_name',
					minTempNCDF = 'file_location_and_name.nc',
					maxTempNCDF = 'file_location_and_name.nc',
					pptNCDF = 'file_location_and_name.nc',
					tempConversionFactor = NA,
					pptConversionFactor = NA,
					avgTempGiven = FALSE, 
					multipleModels = FALSE,	# are there multiple models that need to be stored?
					startDate = as.Date("1990-01-01"), 	# when does the clock of the netcdf start?
					timeToDaysConversion = 1,	# convert time increments to days if necessary
					dataOut_location = 'save_file_location',
					optionForPET = 1, 	# 1 = PET_fromTemp modified Pen-Mon, 
					variableOrderOption = 1)	# 1 = [longitude,latitude,time]
#########################################################################################################
	# step 3 
	# calibrate model
					# nested function for running HBV
				runHBV_f = function(
					climateInput = climateInput,
					sfcf = 1,	#snowfall correction factor [-]
					tr   = 0,	#solid and liquid precipitation threshold temperature [C]
					tt   = 0,	#melt temperature [C]
					fm   = 2,	#snowmelt factor [mm/C]
					fi   = 2,	#icemelt factor [mm/C]
					fic  = 3,	#debris-covered icemelt factor [mm/C]
					fc   = 200, # field capacity 
					lp   = 0.75, # parameter for PET --> AET
					beta_soils = 1.5, #soil moisture drainage exponent
					k0   = 0.5,	# storage constant of top bucket
					k1   = 0.1,	# storage constant of middle bucket
					k2   = 0.01,	# storage constant of bottom bucket	
					uz1  = 5, #max flux rate from STZ to SUZ in mm/d
					perc = 2)  # max flux rate from SUZ to SLZ in mm/d

modelCalibration_f(climateInputsFileLoc = climateInputsFileLoc,
	historicStreamflowFileLoc = historicStreamflowFileLoc, 
	pathToWatershedsGPKG = pathToWatershedsGPKG, 
	dataOut_location = dataOut_location, 
	numberOfRuns = 500000)
					# the function being called
				modelCalibration_f = function(
					climateInputsFileLoc = 'file_location_and_name.RData',
					historicStreamflowFileLoc = 'https://someplace.gov',
					pathToWatershedsGPKG = 'file_location_and_name.gpkg',
					dataOut_location = 'save_file_location',
					dataSource = 1,											# 1 for FNF from cal.gov, 
					numberOfRuns = 100000,
					targetMetric = 1, 										# 1 = KGE, 2 = NSE, 3 = MAE, 4 = RMSE, 5 = bias
					sfcf = c(runif(500, .2, 1), runif(500, 1, 3)),			#snowfall correction factor [-]
					tr   = runif(1000, -6, 5),								#solid and liquid precipitation threshold temperature [C]
					tt   = runif(1000, -5, 6),								#melt temperature [C]
					fm   = c(runif(500, .2, 1.5), (runif(500, 1.5, 8))),	#snowmelt factor [mm/C]
					fi   = c(runif(500, .2, 1.5), (runif(500, 1.5, 10))),	#icemelt factor [mm/C]
					fic  = runif(1000, 2, 10),								#debris-covered icemelt factor [mm/C]
					fc   = c(runif(500, 25, 150), (runif(500, 150, 1200))),	#field capacity
					lp   = runif(1000, .2, 1),								#parameter to actual ET
					beta_soils = runif(1000, 1, 3),							#beta - exponential value for nonlinear relations between soil storage and runoff
					k0   = c(runif(500, .05, .5), (runif(500, .5, 1))),		#top bucket drainage
					k1   = c(runif(500, .005, .09), (runif(500, .09, .5))),	#middle bucket drainage
					k2   = c(runif(500, .0001, .01), (runif(500, .01, .1))),#bottom bucket drainage	
					uz1  = c(runif(500, .22, 10), (runif(500, 10, 40))),	#max flux rate from STZ to SUZ in mm/d
					perc = c(runif(500, .1, .5), (runif(500, .5, 20))))		#max flux rate from SUZ to SLZ in mm/d
#########################################################################################################
	# step 4
	# import and convert projection climate data
climateDataframe =	climateInputConversion_f(
		pathToBasinBoundaryGPKG = pathToBasinBoundaryGPKG,
		pathToWatershedsGPKG = pathToWatershedsGPKG,
		basinName = basinName,
		minTempNCDF = minTempNCDF,
		maxTempNCDF = maxTempNCDF,
		pptNCDF = pptNCDF,
		dataOut_location = dataOut_location)


	# step 4 
	# run the model with forecasting data



