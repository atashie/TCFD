	#load libraries
library(lubridate)
library(ncdf4)
library(data.table)


####################################################################################################
##	define all names, file locations, and variables
	# names and variables
growYearStartDate = as.Date('2022-08-15')
dataOrigination = c('ERA5', 'CFS', 'SEAS5')
dataPath = 'J:\\Cai_data\\Advanta\\GDD\\'
era5RecentDataName = 'testing-recent-era.nc'
era5ClimatologyDataName = 'testing-climatology-era.nc'
cfsDataName = 'testing-cfs.nc'
seas5DataName = 'testing-seas5.nc'
startDateEra5 = '2022-06-01'
startDateCfs = '2022-02-28'
startDateSeas5 = '2022-07-02'

	# define the 'sortable' variables
theseQuantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
soilTempThresholds = c(10, 13, 16)
consecDaysThresholds = c(5, 10)
baseTempGDDs = c(10, 12)	# 10C for corn, 12C for sorghum
GDDvalueThresholds = seq(1350,1650,50)
userName = 'Advanta'
forecastDate = ncvar_get(nc_open(paste0(dataPath, cfsDataName)), 'time')[1]  / 24 + as.Date(startDateCfs) 


#################################################################################################
##	historic data / climatology only needs to be run once, then data are stored for future use
historicOutput = f_historicTempGDD(
	theseQuantiles = theseQuantiles,
	soilTempThresholds = soilTempThresholds,
	consecDaysThresholds = consecDaysThresholds,
	baseTempGDDs = baseTempGDDs,
	GDDvalueThresholds = GDDvalueThresholds, 
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	dataName = era5ClimatologyDataName, 
	userName = userName,
	forecastDate = forecastDate)

fwrite(historicOutput, paste0(dataPath, "climatologySoilTempAndGDD_13JUL2022.csv"))
##################################################################################################


##################################################################################################
##	forecast data, to be run every week
projectedOutput = f_projectedTempGDD(
	theseQuantiles = theseQuantiles,
	soilTempThresholds = soilTempThresholds,
	consecDaysThresholds = consecDaysThresholds,
	baseTempGDDs = baseTempGDDs,
	GDDvalueThresholds = GDDvalueThresholds, 
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	era5RecentDataName = era5RecentDataName,
	cfsDataName = cfsDataName,
	seas5DataName = seas5DataName,
	startDateEra5 = startDateEra5,
	startDateCfs = startDateCfs,
	startDateSeas5 = startDateSeas5,
	whichCfsModels = whichCfsModels,
	whichSeas5Models = whichSeas5Models,
	userName = userName,
	forecastDate = forecastDate)	




allOutput = rbind(projectedOutput, fread(paste0(dataPath, "climatologySoilTempAndGDD_13JUL2022.csv")))
	# identifying water locations using soil moisture dataset
smData = fread("J:\\Cai_data\\Advanta\\SoilMoisture_redux_newsm_20JUL2022.csv")
smNonas = subset(smData, !is.na(Lat))
landOutput =  subset(allOutput, Lat == smNonas$Lat[1] & Lon == smNonas$Lon[1])
for(i in 2:nrow(smNonas)){
	landOutput = rbind(landOutput, subset(allOutput, Lat == smNonas$Lat[i] & Lon == smNonas$Lon[i]))
}
rmDesertOutput = subset(landOutput, Lat < -29 | Lon > 144)

fwrite(rmDesertOutput, paste0(dataPath, 'GDD_and_soilTemp_subset_', userName, "_", forecastDate, '.csv'))
#fwrite(allOutput, paste0(dataPath, 'GDD_and_soilTemp_', userName, "_", forecastDate, '.csv'))
##################################################################################################

















############################################################
##	function for defining climatology GDD / temp thresholds



f_historicTempGDD = function(theseQuantiles = theseQuantiles,
	soilTempThresholds = soilTempThresholds,
	consecDaysThresholds = consecDaysThresholds,
	baseTempGDDs = baseTempGDDs,
	GDDvalueThresholds = GDDvalueThresholds, 
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	dataName = era5ClimatologyDataName, 
	userName = userName,
	forecastDate = forecastDate)	{

	library(lubridate)
	save_iter = 0
	
		# read in ncdf data
	ncin = nc_open(paste0(dataPath, dataName)) # t2m_max[longitude,latitude,time]   (Contiguous storage) 
	nc_tavg = (ncvar_get(ncin, 't2m_min') + ncvar_get(ncin, 't2m_max')) / 2
	nc_tsoil = (ncvar_get(ncin, 'stl1_mean')) - 273.15

	nc_date = ncvar_get(ncin, 'time') + as.Date('2002-06-01')
	nc_doy = yday(nc_date)
	nc_year = year(nc_date)
	
	nc_lat = ncvar_get(ncin, 'latitude')
	nc_lon = ncvar_get(ncin, 'longitude')
	
	if(yday(growYearStartDate) == 1)	{
		growYearDOYs = 1:366
	}	else	{growYearDOYs = c(yday(growYearStartDate):366, 1:(yday(growYearStartDate)-1))}
	

#	if(!is.na(soilTempName))	{										# !!!!!!!!!!!!!!!! REVISIT WHEN WE HAVE SOIL TEMP DATA
#		nc_soilTemp = ncvar_get(nc_open(paste0(dataPath, soilTempName)), 'WHAT IS THE VARIABLE NAME')
#	}	else	{ nc_soilTemp = nc_tavg - 1.5 }


	finalIter = 0	
	summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
		variableName = 'string', projectionOrClimatology = 'string',
		quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)

	# generating .ncs for saving average historical values for backfilling projection data
	tempClimatologyTavg_nc = nc_tavg[ , , 1:366]
	tempClimatologySoils_nc = nc_tsoil[ , , 1:366]
	
	for(this_lon in 1:length(nc_lon))	{
		for(this_lat in 1:length(nc_lat))	{
			
			theseSoilTemp = nc_tsoil[this_lon, this_lat, ]
			if(any(!is.na(theseSoilTemp)))	{
				theseTavg = nc_tavg[this_lon, this_lat, ]
				
				for(gg in 1:366)	{
					tempClimatologyTavg_nc[this_lon, this_lat, gg] = median(theseTavg[nc_doy == gg])
					tempClimatologySoils_nc[this_lon, this_lat, gg] = median(theseSoilTemp[nc_doy == gg])
				}
							
				allYearsOutput_df = data.frame(soilTempThresholds = NA, consecDaysThresholds = NA, daysToWarmSoils = NA,
					baseTempGDDs = NA, GDDvalueThresholds = NA, daysToGDDThreshold = NA)
				iter = 0
				for(hh in 1:length(baseTempGDDs))	{
					theseGDD = ifelse(theseTavg > baseTempGDDs[hh], theseTavg - baseTempGDDs[hh], 0)
					for(ii in 1:length(soilTempThresholds))	{
						for(jj in 1:length(consecDaysThresholds))	{
							uniqueYears = unique(nc_year)
							for(thisYear in uniqueYears[-c(1,length(uniqueYears))])	{
								thisGrowYear = c(which(nc_year == thisYear - 1 & nc_doy >= yday(growYearStartDate)), 
									which(nc_year == thisYear & nc_doy < yday(growYearStartDate))) 
								exceedsSoilThresh = which(theseSoilTemp[thisGrowYear] > (soilTempThresholds[ii])) ## adjusting to soil temp
								threshDiff = diff(exceedsSoilThresh, consecDaysThresholds[jj])
								daysToWarmSoils = ifelse(any(threshDiff == consecDaysThresholds[jj][1]), 
									exceedsSoilThresh[which(threshDiff == consecDaysThresholds[jj])[1] + (consecDaysThresholds[jj] - 1)],
									364)
			
								for(kk in 1:length(GDDvalueThresholds))	{
									iter = iter + 1
									thisGDD = cumsum(theseGDD[thisGrowYear][-(1:daysToWarmSoils)])
									maxGDD = ifelse(daysToWarmSoils < 364, max(thisGDD) , 0)
									allYearsOutput_df[iter, ] = c(soilTempThresholds[ii], consecDaysThresholds[jj], daysToWarmSoils,
										baseTempGDDs[hh], GDDvalueThresholds[kk],
										ifelse(maxGDD > GDDvalueThresholds[kk], daysToWarmSoils + which(thisGDD > GDDvalueThresholds[kk])[1], 400))
								}
							}
						}
					}
				}
				

				thisLat = nc_lat[this_lat]
				thisLon = nc_lon[this_lon]
				for(thisSoilTempThresholds in soilTempThresholds)				{
					for(thisConsecDaysThresholds in consecDaysThresholds)		{
						for(thisBaseTempGDDs in baseTempGDDs)					{
							daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
									soilTempThresholds == thisSoilTempThresholds & 
									consecDaysThresholds == thisConsecDaysThresholds & 
									baseTempGDDs == thisBaseTempGDDs)$daysToWarmSoils,
									theseQuantiles))
	# bc looker timeline cannot handle quantiles with the same value, need to artificially shift quantiles
							if(any(diff(daysSinceDate) == 0))	{
								daysSinceDate = daysSinceDate + seq(from=0, by=1, length.out = length(theseQuantiles))
							}	
							datesByQuantile = growYearStartDate + daysSinceDate

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'climatology',
								'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
								paste0(userName), paste0(forecastDate), daysSinceDate[2])

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'climatology',
								'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
								paste0(userName), paste0(forecastDate), daysSinceDate[3])

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'climatology',
								'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
								paste0(userName), paste0(forecastDate), daysSinceDate[4])
							
							for(thisGDDvalueThresholds in GDDvalueThresholds)	{
								daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
										soilTempThresholds == thisSoilTempThresholds &
										consecDaysThresholds == thisConsecDaysThresholds &
										baseTempGDDs == thisBaseTempGDDs &
										GDDvalueThresholds == thisGDDvalueThresholds)$daysToGDDThreshold,
										theseQuantiles))

	# bc looker timeline cannot handle quantiles with the same value, need to artificially shift quantiles
								if(any(diff(daysSinceDate) == 0))	{
									daysSinceDate = daysSinceDate + seq(from=0, by=1, length.out = length(theseQuantiles))
								}	
								datesByQuantile = growYearStartDate + daysSinceDate
								
								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'climatology',
									'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
									paste0(userName), paste0(forecastDate), daysSinceDate[2])

								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'climatology',
									'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
									paste0(userName), paste0(forecastDate), daysSinceDate[3])

								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'climatology',
									'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
									paste0(userName), paste0(forecastDate), daysSinceDate[4])
							}
						}
					}
				}
				if(nrow(summaryOutput_df) > 20000)	{
					save_iter = save_iter + 1
					fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))
					summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
						variableName = 'string', projectionOrClimatology = 'string',
						quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
					finalIter = 0
					print(c(this_lat, this_lon))
				}
			} 
		}
	}
		# saving workspace that holds climatology .ncs
	save(tempClimatologyTavg_nc, file=paste0(dataPath, 'climatologyTavg.RData'))
	save(tempClimatologySoils_nc, file=paste0(dataPath, 'climatologySoils.RData'))
		# combining saved iterations and returning the dataframe with formatted climatology
	save_iter = save_iter + 1
	fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))

	summaryOutput_df = fread(paste0(dataPath, "temp_out_", 1, ".csv"))
	for(files in 2:save_iter)	{
		summaryOutput_df = rbind(summaryOutput_df, fread(paste0(dataPath, "temp_out_", files, ".csv")))
	}
	return(summaryOutput_df)
}	

	
load(paste0(dataPath, 'climatologyTavg.RData'))	# loads tempClimatologyTavg_nc
load(paste0(dataPath, 'climatologySoils.RData'))# loads tempClimatologySoils_nc






















										


	


############################################################################################################################
# incorporating projections data

f_projectedTempGDD = function(theseQuantiles = theseQuantiles,
	soilTempThresholds = soilTempThresholds,
	consecDaysThresholds = consecDaysThresholds,
	baseTempGDDs = baseTempGDDs,
	GDDvalueThresholds = GDDvalueThresholds, 
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	era5RecentDataName = era5RecentDataName,
	cfsDataName = cfsDataName,
	seas5DataName = seas5DataName,
	startDateEra5 = startDateEra5,
	startDateCfs = startDateCfs,
	startDateSeas5 = startDateSeas5,
	whichCfsModels = whichCfsModels,
	whichSeas5Models = whichSeas5Models,
	userName = userName,
	forecastDate = forecastDate,
	seas5Models = 1:51,
	cfsModels = 1:4)	{

	library(lubridate)
	save_iter = 0

		# loading climatology records of Tavg and soils by day of year
	load(paste0(dataPath, 'climatologyTavg.RData'))	# loads tempClimatologyTavg_nc
	load(paste0(dataPath, 'climatologySoils.RData'))# loads tempClimatologySoils_nc


		# read in ncdf data
	ncin_era5 = nc_open(paste0(dataPath, era5RecentDataName))
	ncin_cfs = nc_open(paste0(dataPath, cfsDataName))
	ncin_seas5 = nc_open(paste0(dataPath, seas5DataName),  return_on_error = TRUE)
	
	nc_dateEra5 = ncvar_get(ncin_era5, 'time') + as.Date(startDateEra5)
	nc_dateCfs = ncvar_get(ncin_cfs, 'time') / 24 +  ncvar_get(ncin_cfs, 'step') + as.Date(startDateCfs)
	nc_dateSeas5 = ncvar_get(ncin_seas5, 'lead_time') + as.Date(startDateSeas5)
	 
	nc_doyEra5 = yday(nc_dateEra5)
	nc_doyCfs = yday(nc_dateCfs)
	nc_doySeas5 = yday(nc_dateSeas5)
	
	nc_yearEra5 = year(nc_dateEra5)
	nc_yearCfs = year(nc_dateCfs)
	nc_yearSeas5 = year(nc_dateSeas5)
	
	nc_latEra5 = ncvar_get(ncin_era5, 'latitude')
	nc_latCfs = ncvar_get(ncin_cfs, 'latitude')
	nc_latSeas5 = ncvar_get(ncin_seas5, 'latitude')
	
	nc_lonEra5 = ncvar_get(ncin_era5, 'longitude')
	nc_lonCfs = ncvar_get(ncin_cfs, 'longitude')
	nc_lonSeas5 = ncvar_get(ncin_seas5, 'longitude')

		# defining the growing year
	growYear = seq(growYearStartDate, (growYearStartDate+364), 1)
	growYearDOYs = yday(growYear)
	
		# prioritizing higher to lower quality data
	whichEra5Days = which(nc_dateEra5 %in% growYear) 
	incHistData = ifelse(length(whichEra5Days >=1), TRUE, FALSE)
	if(incHistData)	{
		whichCfsDays = which(as.character(nc_dateCfs) %in% as.character(growYear))[-c(1:length(whichEra5Days))]
		whichSeas5Days = which(as.character(nc_dateSeas5) %in% as.character(growYear))[-c(1:(length(whichEra5Days) + length(whichCfsDays)))]
		whichClimatologyDays = growYearDOYs[-c(1:(length(incHistData) + length(whichCfsDays) + length(whichSeas5Days)))]
	} else {
		whichCfsDays = which(as.character(nc_dateCfs) %in% as.character(growYear))
		whichSeas5Days = which(as.character(nc_dateSeas5) %in% as.character(growYear))[-c(1:length(whichCfsDays))]
		whichClimatologyDays = growYearDOYs[-c(1:(length(whichCfsDays) + length(whichSeas5Days)))]
	}
	

	nc_tavgEra5 = (ncvar_get(ncin_era5, 't2m_min')[ , , whichEra5Days] + ncvar_get(ncin_era5, 't2m_max')[ , , whichEra5Days]) / 2
	nc_tavgCfs = (ncvar_get(ncin_cfs, 't2m_min')[ , , , whichCfsDays] + ncvar_get(ncin_cfs, 't2m_max')[ , , , whichCfsDays]) / 2
	nc_tavgSeas5 = (ncvar_get(ncin_seas5, 't2m_min')[ , , , whichSeas5Days] + ncvar_get(ncin_seas5, 't2m_max')[ , , , whichSeas5Days]) / 2
	
	nc_tsoil = (ncvar_get(ncin_era5, 'stl1_mean'))[ , , 1:100] - 273.15

		# soil temp not working
#	nc_tsoilEra5 = (ncvar_get(ncin_era5, 'stl1_mean') - 273.15
#	nc_tsoilCfs = (ncvar_get(ncin_cfs, 'stl1') - 273.15
#	nc_tsoilSeas5 = (ncvar_get(ncin_seas5, 'stl1_mean') - 273.15


	finalIter = 0	
	summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
		variableName = 'string', projectionOrClimatology = 'string',
		quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
	for(this_lon in 1:length(nc_lonEra5))	{
		for(this_lat in 1:length(nc_latEra5))	{
	
			theseSoilTemp = nc_tsoil[this_lon, this_lat, ]
			if(any(!is.na(theseSoilTemp)))	{
	
				allYearsOutput_df = data.frame(soilTempThresholds = NA, consecDaysThresholds = NA, daysToWarmSoils = NA,
					baseTempGDDs = NA, GDDvalueThresholds = NA, daysToGDDThreshold = NA)
				iter = 0
				for(hh in 1:length(baseTempGDDs))	{
					for(ii in 1:length(soilTempThresholds))	{
						for(jj in 1:length(consecDaysThresholds))	{
							for(thisCfsModel in cfsModels)	{
								for(thisSeas5Model in seas5Models)	{
									if(incHistData)	{
										theseTavg = c(nc_tavgEra5[this_lon, this_lat, ],
											nc_tavgCfs[this_lon, this_lat, thisCfsModel, ],
											nc_tavgSeas5[this_lon, this_lat, thisSeas5Model, ],
											tempClimatologyTavg_nc[this_lon, this_lat, whichClimatologyDays])
										theseTavg[is.na(theseTavg)] = mean(theseTavg, na.rm=TRUE)
									} else {
										theseTavg = c(nc_tavgCfs[this_lon, this_lat, thisCfsModel, ],
											nc_tavgSeas5[this_lon, this_lat, thisSeas5Model, ],
											tempClimatologyTavg_nc[this_lon, this_lat, whichClimatologyDays])
										theseTavg[is.na(theseTavg)] = mean(theseTavg, na.rm=TRUE)
									}

									theseGDD = ifelse(theseTavg > baseTempGDDs[hh], theseTavg - baseTempGDDs[hh], 0)
		
									exceedsSoilThresh = which((theseTavg - 1.5) > (soilTempThresholds[ii])) 
									threshDiff = diff(exceedsSoilThresh, consecDaysThresholds[jj])
									daysToWarmSoils = ifelse(any(threshDiff == consecDaysThresholds[jj][1]), 
										exceedsSoilThresh[which(threshDiff == consecDaysThresholds[jj])[1] + (consecDaysThresholds[jj] - 1)],
										364)
			
									for(kk in 1:length(GDDvalueThresholds))	{
										iter = iter + 1
										thisGDD = cumsum(theseGDD[-(1:daysToWarmSoils)])
										maxGDD = ifelse(daysToWarmSoils < 364, max(thisGDD) , 0)
										allYearsOutput_df[iter, ] = c(soilTempThresholds[ii], consecDaysThresholds[jj], daysToWarmSoils,
											baseTempGDDs[hh], GDDvalueThresholds[kk],
											ifelse(maxGDD > GDDvalueThresholds[kk], daysToWarmSoils + which(thisGDD > GDDvalueThresholds[kk])[1], 365))
									}
								}
							}	
						}
					}
				}
				

				thisLat = nc_latEra5[this_lat]
				thisLon = nc_lonEra5[this_lon]
				for(thisSoilTempThresholds in soilTempThresholds)				{
					for(thisConsecDaysThresholds in consecDaysThresholds)		{
						for(thisBaseTempGDDs in baseTempGDDs)					{
							daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
									soilTempThresholds == thisSoilTempThresholds & 
									consecDaysThresholds == thisConsecDaysThresholds &
									baseTempGDDs == thisBaseTempGDDs)$daysToWarmSoils,
									theseQuantiles))
	# bc looker timeline cannot handle quantiles with the same value, need to artificially shift quantiles
							if(any(diff(daysSinceDate) == 0))	{
								daysSinceDate = daysSinceDate + seq(from=0, by=1, length.out = length(theseQuantiles))
							}	
							datesByQuantile = growYearStartDate + daysSinceDate

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'projection',
								'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
								paste0(userName), paste0(forecastDate), daysSinceDate[2])

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'projection',
								'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
								paste0(userName), paste0(forecastDate), daysSinceDate[3])

							finalIter = finalIter + 1
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToSoilTemp'), 'projection',
								'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
								paste0(userName), paste0(forecastDate), daysSinceDate[4])
								
								
							for(thisGDDvalueThresholds in GDDvalueThresholds)	{
								daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
										soilTempThresholds == thisSoilTempThresholds &
										consecDaysThresholds == thisConsecDaysThresholds &
										baseTempGDDs == thisBaseTempGDDs &
										GDDvalueThresholds == thisGDDvalueThresholds)$daysToGDDThreshold,
										theseQuantiles))
	# bc looker timeline cannot handle quantiles with the same value, need to artificially shift quantiles
								if(any(diff(daysSinceDate) == 0))	{
									daysSinceDate = daysSinceDate + seq(from=0, by=1, length.out = length(theseQuantiles))
								}	
								datesByQuantile = growYearStartDate + daysSinceDate
								
								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'projection',
									'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
									paste0(userName), paste0(forecastDate), daysSinceDate[2])

								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'projection',
									'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
									paste0(userName), paste0(forecastDate), daysSinceDate[3])

								finalIter = finalIter + 1
								summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds), 'projection',
									'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
									paste0(userName), paste0(forecastDate), daysSinceDate[4])
		

								}
						}
					}
					if(nrow(summaryOutput_df) > 20000)	{
						save_iter = save_iter + 1
						fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))
						summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
							variableName = 'string', projectionOrClimatology = 'string',
							quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
						finalIter = 0
						print(c(this_lat, this_lon))
					}
				} 
			}
		}
	}
	# combining saved iterations and returning the dataframe with formatted climatology
	save_iter = save_iter + 1
	fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))

	summaryOutput_df = fread(paste0(dataPath, "temp_out_", 1, ".csv"))
	for(files in 2:save_iter)	{
		summaryOutput_df = rbind(summaryOutput_df, fread(paste0(dataPath, "temp_out_", files, ".csv")))
	}
	return(summaryOutput_df)
}



































							
							
			allYearsOutput_df = data.frame(soilTempThresholds = NA, consecDaysThresholds = NA, daysToWarmSoils = NA,
				baseTempGDDs = NA, GDDvalueThresholds = NA, daysToGDDThreshold = NA)


		for(thisGDDvalueThresholds in GDDvalueThresholds)	{
										
					finalIter = finalIter + 1
								summaryOutput_df = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, 
									thisBaseTempGDDs, thisGDDvalueThresholds,
									variableName = 
								
	
	summaryOutput_df = data.frame(Lat = NA, Lon = NA, soilTempThresholds = NA, consecDaysThresholds = NA, baseTempGDDs = NA, 
		variableName = NA, Q05 = NA, Q25=NA, Q50 = NA, Q75 = NA, Q95 = NA)
		
		allHist_df[final_iter, c('Lat','Lon')] = c(nc_lat_hist[this_lat], nc_lon_hist[this_lon])
		firstCol = which(names(allHist_df) == 'DaysToWSls_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$Days_to_Warm_Soils, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1350_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1350, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1400_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1400, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1450_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1450, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1500_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1500, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1550_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1550, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1600_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1600, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1650_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1650, c(0.05,0.25,0.50,0.75,0.95))





					output_df[iter, ] = c(soilTempThresholds[ii], consecDaysThresholds[jj], daysToWarmSoils)

		this_gdd = cumsum(these_gdd[this_grow_year][-c(1:daysToWarmSoils)])
			max_gdd = max(this_gdd)


		hist_df[iter,c('GDD_1350','GDD_1400','GDD_1450','GDD_1500','GDD_1550','GDD_1600','GDD_1650')] = c(
				ifelse(max_gdd > 1350, which(this_gdd > GDD_thresholds[1])[1],400),
				ifelse(max_gdd > 1400, which(this_gdd > GDD_thresholds[2])[1],400),
				ifelse(max_gdd > 1450, which(this_gdd > GDD_thresholds[3])[1],400),
				ifelse(max_gdd > 1500, which(this_gdd > GDD_thresholds[4])[1],400),
				ifelse(max_gdd > 1550, which(this_gdd > GDD_thresholds[5])[1],400),
				ifelse(max_gdd > 1600, which(this_gdd > GDD_thresholds[6])[1],400),
				ifelse(max_gdd > 1650, which(this_gdd > GDD_thresholds[7])[1],400))
			
		}


					}
				}

			}


for(this_year in unique(nc_years)[-1])	{
			iter = iter + 1
			this_grow_year = c(which(nc_years == this_year - 1 & nc_months %in% thisGrowingMonth:12),
				which(nc_years == this_year & nc_months %in% 1:6))
			these_soilThresh = which(these_tavg[this_grow_year] > (threshTemp+2)) ## adjusting to soil temp
			threshDiff = diff(these_soilThresh, 5)
			daysToWarmSoils = these_soilThresh[which(threshDiff == 5)[1] + 4]
			
			
			
			this_gdd = cumsum(these_gdd[this_grow_year][-c(1:daysToWarmSoils)])
			max_gdd = max(this_gdd)
			
			hist_df[iter,'Days_to_Warm_Soils'] = daysToWarmSoils 
		
						

	soilTempThresholds = soilTempThresholds,
	consecDaysThresholds = cosecDaysThresholds,




			nc_gdd = ifelse(nc_tavg_hist > threshTemp, nc_tavg_hist - threshTemp, 0)
	


	
			these_gdd = nc_gdd_hist[this_lon, this_lat, ]
			hist_df = data.frame(Days_to_Warm_Soils = NA,
				GDD_1350 = NA, GDD_1400 = NA, GDD_1450 = NA, GDD_1500 = NA,
				GDD_1550 = NA, GDD_1600 = NA, GDD_1650 = NA)
		



	
		
	nc_lat = ncvar_get(ncin_tmin, 'latitude')
	nc_lon = ncvar_get(ncin_tmax, 'longitude')
	nc_dates = ncvar_get(ncin_tmin, 'time')  + as.Date("1979-01-01") # time var in days since 1979-01-01
	theseGrowingDates = seq(growYearStartDate, growYearStartDate + 365, 1)
	theseGrowingDOYs = unique(yday(theseGrowingDates))	# unique() for removing repeat days of year in non-leap years

	nc_years = year(nc_dates)
	nc_months = month(nc_dates)
	nc_yday_hist = yday(nc_dates)
	theseGrowingDates = seq(growYearStartDate, growYearStartDate + 365, 1)


	
	for(this_lon in 1:length(nc_lon))	{
		for(this_lat in 1:length(nc_lat))	{
			these_tavg = nc_tavg[this_lon, this_lat, ]



			
			nc_gdd = ifelse(nc_tavg_hist > threshTemp, nc_tavg_hist - threshTemp, 0)
			
			these_gdd = nc_gdd_hist[this_lon, this_lat, ]
			hist_df = data.frame(Days_to_Warm_Soils = NA,
				GDD_1350 = NA, GDD_1400 = NA, GDD_1450 = NA, GDD_1500 = NA,
				GDD_1550 = NA, GDD_1600 = NA, GDD_1650 = NA)
		
		
		
		
		
		
		nc_gdd = ifelse(nc_tavg_hist > threshTemp, nc_tavg_hist - threshTemp, 0)

soilThreshLength = 5	# number of days to exceed threshold temp

thisGrowingYear = 2022 #  matching climatology to year of sewing (this year) and year of harvest (next year) for southern hemisphere
thisGrowingMonth = 8 # which month of the year to we bein assessing accumulation?	
	# timeseries of the dates relevant to the coming growing season
thisGrowingDates = seq(as.Date(paste0(thisGrowingYear, "-0", thisGrowingMonth, "-01")),as.Date(paste0(thisGrowingYear+1, "-0", thisGrowingMonth, "-01"))-1,1)


	
		

		iter = 0
		for(this_year in unique(nc_years)[-1])	{
			iter = iter + 1
			this_grow_year = c(which(nc_years == this_year - 1 & nc_months %in% thisGrowingMonth:12),
				which(nc_years == this_year & nc_months %in% 1:6))
			these_soilThresh = which(these_tavg[this_grow_year] > (threshTemp+2)) ## adjusting to soil temp
			threshDiff = diff(these_soilThresh, 5)
			daysToWarmSoils = these_soilThresh[which(threshDiff == 5)[1] + 4]
			
			
			
			this_gdd = cumsum(these_gdd[this_grow_year][-c(1:daysToWarmSoils)])
			max_gdd = max(this_gdd)
			
			hist_df[iter,'Days_to_Warm_Soils'] = daysToWarmSoils 
			
			
			hist_df[iter,c('GDD_1350','GDD_1400','GDD_1450','GDD_1500','GDD_1550','GDD_1600','GDD_1650')] = c(
				ifelse(max_gdd > 1350, which(this_gdd > GDD_thresholds[1])[1],400),
				ifelse(max_gdd > 1400, which(this_gdd > GDD_thresholds[2])[1],400),
				ifelse(max_gdd > 1450, which(this_gdd > GDD_thresholds[3])[1],400),
				ifelse(max_gdd > 1500, which(this_gdd > GDD_thresholds[4])[1],400),
				ifelse(max_gdd > 1550, which(this_gdd > GDD_thresholds[5])[1],400),
				ifelse(max_gdd > 1600, which(this_gdd > GDD_thresholds[6])[1],400),
				ifelse(max_gdd > 1650, which(this_gdd > GDD_thresholds[7])[1],400))
			
		}
		final_iter = final_iter + 1
		allHist_df[final_iter, c('Lat','Lon')] = c(nc_lat_hist[this_lat], nc_lon_hist[this_lon])
		firstCol = which(names(allHist_df) == 'DaysToWSls_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$Days_to_Warm_Soils, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1350_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1350, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1400_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1400, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1450_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1450, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1500_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1500, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1550_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1550, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1600_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1600, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1650_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1650, c(0.05,0.25,0.50,0.75,0.95))
	}
}			
	

	












final_iter = 0
for(this_lon in 1:length(nc_lon_hist))	{
	for(this_lat in 1:length(nc_lat_hist))	{
		these_tavg = nc_tavg_hist[this_lon, this_lat, ]
		these_gdd = nc_gdd_hist[this_lon, this_lat, ]
		hist_df = data.frame(Days_to_Warm_Soils = NA,
			GDD_1350 = NA, GDD_1400 = NA, GDD_1450 = NA, GDD_1500 = NA,
			GDD_1550 = NA, GDD_1600 = NA, GDD_1650 = NA)
		

		iter = 0
		for(this_year in unique(nc_years)[-1])	{
			iter = iter + 1
			this_grow_year = c(which(nc_years == this_year - 1 & nc_months %in% thisGrowingMonth:12),
				which(nc_years == this_year & nc_months %in% 1:6))
			these_soilThresh = which(these_tavg[this_grow_year] > (threshTemp+2)) ## adjusting to soil temp
			threshDiff = diff(these_soilThresh, 5)
			daysToWarmSoils = these_soilThresh[which(threshDiff == 5)[1] + 4]
			
			
			
			this_gdd = cumsum(these_gdd[this_grow_year][-c(1:daysToWarmSoils)])
			max_gdd = max(this_gdd)
			
			hist_df[iter,'Days_to_Warm_Soils'] = daysToWarmSoils 
			
			
			hist_df[iter,c('GDD_1350','GDD_1400','GDD_1450','GDD_1500','GDD_1550','GDD_1600','GDD_1650')] = c(
				ifelse(max_gdd > 1350, which(this_gdd > GDD_thresholds[1])[1],400),
				ifelse(max_gdd > 1400, which(this_gdd > GDD_thresholds[2])[1],400),
				ifelse(max_gdd > 1450, which(this_gdd > GDD_thresholds[3])[1],400),
				ifelse(max_gdd > 1500, which(this_gdd > GDD_thresholds[4])[1],400),
				ifelse(max_gdd > 1550, which(this_gdd > GDD_thresholds[5])[1],400),
				ifelse(max_gdd > 1600, which(this_gdd > GDD_thresholds[6])[1],400),
				ifelse(max_gdd > 1650, which(this_gdd > GDD_thresholds[7])[1],400))
			
		}
		final_iter = final_iter + 1
		allHist_df[final_iter, c('Lat','Lon')] = c(nc_lat_hist[this_lat], nc_lon_hist[this_lon])
		firstCol = which(names(allHist_df) == 'DaysToWSls_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$Days_to_Warm_Soils, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1350_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1350, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1400_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1400, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1450_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1450, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1500_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1500, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1550_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1550, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1600_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1600, c(0.05,0.25,0.50,0.75,0.95))
		firstCol = which(names(allHist_df) == 'GDD_1650_Q05')
		allHist_df[final_iter,firstCol:(firstCol+4)] = quantile(hist_df$GDD_1650, c(0.05,0.25,0.50,0.75,0.95))
	}
}			
	











