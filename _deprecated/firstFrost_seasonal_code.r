	#load libraries
library(lubridate)
library(ncdf4)
library(data.table)


####################################################################################################
##	define all names, file locations, and variables
	# names and variables
growYearStartDate = as.Date('2022-08-15')
dataOrigination = c('ERA5', 'CFS', 'SEAS5')
dataPath = 'J:\\Cai_data\\Simplot\\firstFrost\\'
era5RecentDataName = 'testing-recent-era.nc'
era5ClimatologyDataName = 'testing-climatology-era.nc'
cfsDataName = 'testing-cfs.nc'
seas5DataName = 'testing-seas5.nc'
startDateEra5 = '2022-06-01'
startDateCfs = '2022-02-28'
startDateSeas5 = '2022-07-02'
cfsModels = 1:4	# number of cfs models used
seas5Models = 1:51	# number of seas5 models used

	# define the 'sortable' variables
theseQuantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
firstFrostThresholds = c(-2, -1, 0, 1)	# in C
consecDaysThresholds = c(1, 2, 5)
userName = 'Simplot'
forecastDate = ncvar_get(nc_open(paste0(dataPath, cfsDataName)), 'time')[1]  / 24 + as.Date(startDateCfs) 


#################################################################################################
##	historic data / climatology only needs to be run once, then data are stored for future use
historicOutput = f_historicFirstFrost(
	theseQuantiles = theseQuantiles,
	firstFrostThresholds = firstFrostThresholds,
	consecDaysThresholds = consecDaysThresholds,
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	dataName = era5ClimatologyDataName, 
	userName = userName,
	forecastDate = forecastDate)

fwrite(historicOutput, paste0(dataPath, "climatologyFirstFrost_18JUL2022.csv"))
##################################################################################################


##################################################################################################
##	forecast data, to be run every week
projectedOutput = f_projectedFirstFrost(
	theseQuantiles = theseQuantiles,
	firstFrostThresholds = firstFrostThresholds,
	consecDaysThresholds = consecDaysThresholds,
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	era5RecentDataName = era5RecentDataName,
	cfsDataName = cfsDataName,
	seas5DataName = seas5DataName,
	startDateEra5 = startDateEra5,
	startDateCfs = startDateCfs,
	startDateSeas5 = startDateSeas5,
	cfsModels = cfsModels,
	seas5Models = seas5Models,
	userName = userName,
	forecastDate = forecastDate)	




allOutput = rbind(projectedOutput, fread(paste0(dataPath, "climatologyFirstFrost_18JUL2022.csv")))
fwrite(allOutput, paste0(dataPath, 'firstFrost_', allOutput[1, 'User'], "_", forecastDate, '.csv'))
##################################################################################################

















############################################################
##	function for defining climatology GDD / temp thresholds



f_historicFirstFrost = function(theseQuantiles = theseQuantiles,
	firstFrostThresholds = firstFrostThresholds,
	consecDaysThresholds = consecDaysThresholds,
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	dataName = era5ClimatologyDataName, 
	userName = userName,
	forecastDate = forecastDate)	{

	library(lubridate)
	save_iter = 0
	
		# read in ncdf data
	ncin = nc_open(paste0(dataPath, dataName)) # t2m_max[longitude,latitude,time]   (Contiguous storage) 
	nc_tmin = (ncvar_get(ncin, 't2m_min'))

	nc_date = ncvar_get(ncin, 'time') + as.Date('2002-06-01')
	nc_doy = yday(nc_date)
	nc_year = year(nc_date)
	
	nc_lat = ncvar_get(ncin, 'latitude')
	nc_lon = ncvar_get(ncin, 'longitude')
	
	if(yday(growYearStartDate) == 1)	{
		growYearDOYs = 1:366
	}	else	{growYearDOYs = c(yday(growYearStartDate):366, 1:(yday(growYearStartDate)-1))}
	

	finalIter = 0	

	# generating .ncs for saving average historical values for backfilling projection data
	tempClimatologyTmin_nc = nc_tmin[ , , 1:366]

	summaryOutput_df = data.frame(Lat = 0, Lon = 0, firstFrostThresholds = 0, consecDaysThresholds = 0, 
		variableName = 'string', projectionOrClimatology = 'string',
		quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)

	for(this_lon in 1:length(nc_lon))	{
		for(this_lat in 1:length(nc_lat))	{
			
			theseTmin = nc_tmin[this_lon, this_lat, ]
			
			for(gg in 1:366)	{
					tempClimatologyTmin_nc[this_lon, this_lat, gg] = median(theseTmin[nc_doy == gg])
				}
				
			
			allYearsOutput_df = data.frame(firstFrostThresholds = NA, consecDaysThresholds = NA, daysToFirstFrost = NA)
			iter = 0

			for(ii in 1:length(firstFrostThresholds))	{
				for(jj in 1:length(consecDaysThresholds))	{
					uniqueYears = unique(nc_year)
					for(thisYear in uniqueYears[-c(1,length(uniqueYears))])	{
						iter = iter+1
						thisGrowYear = c(which(nc_year == thisYear - 1 & nc_doy >= yday(growYearStartDate)), 
							which(nc_year == thisYear & nc_doy < yday(growYearStartDate))) 
						
						atFirstFrostTemp = which(theseTmin[thisGrowYear] <= (firstFrostThresholds[ii])) ## adjusting to soil temp

						threshDiff = diff(atFirstFrostTemp, consecDaysThresholds[jj])
						daysToFirstFrost = ifelse(any(threshDiff == consecDaysThresholds[jj][1]), 
							atFirstFrostTemp[which(threshDiff == consecDaysThresholds[jj])[1] + (consecDaysThresholds[jj] - 1)],
							364)
						allYearsOutput_df[iter, ] = c(firstFrostThresholds[ii], consecDaysThresholds[jj], daysToFirstFrost)
					}
				}
			}
		
				

			thisLat = nc_lat[this_lat]
			thisLon = nc_lon[this_lon]
			for(thisFirstFrostThresholds in firstFrostThresholds)				{
				for(thisConsecDaysThresholds in consecDaysThresholds)		{
					daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
							firstFrostThresholds == thisFirstFrostThresholds & 
							consecDaysThresholds == thisConsecDaysThresholds)$daysToFirstFrost,
							theseQuantiles))
					datesByQuantile = growYearStartDate + daysSinceDate

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds, 
						paste0('DaysToFirstFrost'), 'climatology',
						'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
						paste0(userName), paste0(forecastDate), daysSinceDate[2])

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds, 
						paste0('DaysToFirstFrost'), 'climatology',
						'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
						paste0(userName), paste0(forecastDate), daysSinceDate[3])

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds, 
						paste0('DaysToFirstFrost'), 'climatology',
						'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
						paste0(userName), paste0(forecastDate), daysSinceDate[4])
							
				}
			}
			
			if(nrow(summaryOutput_df) > 20000)	{
				save_iter = save_iter + 1
				fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))
				summaryOutput_df = data.frame(Lat = 0, Lon = 0, firstFrostThresholds = 0, consecDaysThresholds = 0, 
					variableName = 'string', projectionOrClimatology = 'string',
					quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
				finalIter = 0
				print(c(this_lat, this_lon))
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

		# saving workspace that holds climatology .ncs
	save(tempClimatologyTmin_nc, file=paste0(dataPath, 'climatologyTmin.RData'))

	return(summaryOutput_df)
}	





















										


	


############################################################################################################################
# incorporating projections data

f_projectedFirstFrost = function(theseQuantiles = theseQuantiles,
	firstFrostThresholds = firstFrostThresholds,
	consecDaysThresholds = consecDaysThresholds,
	growYearStartDate = growYearStartDate,
	dataPath = dataPath,
	era5RecentDataName = era5RecentDataName,
	cfsDataName = cfsDataName,
	seas5DataName = seas5DataName,
	startDateEra5 = startDateEra5,
	startDateCfs = startDateCfs,
	startDateSeas5 = startDateSeas5,
	cfsModels = 1:4,
	seas5Models = 1:51,
	userName = userName,
	forecastDate = forecastDate)
	{



		# loading climatology records of Tmin by day of year
	load(paste0(dataPath, 'climatologyTmin.RData'))	# loads tempClimatologyTmin_nc

	library(lubridate)
	save_iter = 0

		# read in ncdf data
	ncin_era5 = nc_open(paste0(dataPath, era5RecentDataName))
	ncin_cfs = nc_open(paste0(dataPath, cfsDataName))
	ncin_seas5 = nc_open(paste0(dataPath, seas5DataName))
	
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
	

	nc_tminEra5 = ncvar_get(ncin_era5, 't2m_min')[ , , whichEra5Days]
	nc_tminCfs = ncvar_get(ncin_cfs, 't2m_min')[ , , , whichCfsDays]
	nc_tminSeas5 = ncvar_get(ncin_seas5, 't2m_min')[ , , , whichSeas5Days]


	finalIter = 0	
	summaryOutput_df = data.frame(Lat = 0, Lon = 0, firstFrostThresholds = 0, consecDaysThresholds = 0, 
		variableName = 'string', projectionOrClimatology = 'string',
		quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
	for(this_lon in 1:length(nc_lonEra5))	{
		for(this_lat in 1:length(nc_latEra5))	{
	
			allYearsOutput_df = data.frame(firstFrostThresholds = NA, consecDaysThresholds = NA, daysToFirstFrost = NA)
			iter = 0
	
			for(ii in 1:length(firstFrostThresholds))	{
				for(jj in 1:length(consecDaysThresholds))	{
					for(thisCfsModel in cfsModels)	{
						for(thisSeas5Model in seas5Models)	{
							if(incHistData)	{
								theseTmin = c(nc_tminEra5[this_lon, this_lat, ],
									nc_tminCfs[this_lon, this_lat, thisCfsModel, ],
									nc_tminSeas5[this_lon, this_lat, thisSeas5Model, ],
									tempClimatologyTmin_nc[this_lon, this_lat, whichClimatologyDays])
								theseTmin[is.na(theseTmin)] = mean(theseTmin, na.rm=TRUE)
							} else {
								theseTmin = c(nc_tminCfs[this_lon, this_lat, thisCfsModel, ],
									nc_tminSeas5[this_lon, this_lat, thisSeas5Model, ],
									tempClimatologyTmin_nc[this_lon, this_lat, whichClimatologyDays])
								theseTmin[is.na(theseTmin)] = mean(theseTmin, na.rm=TRUE)
							}

							atFirstFrostTemp = which(theseTmin <= (firstFrostThresholds[ii])) ## adjusting to soil temp

							threshDiff = diff(atFirstFrostTemp, consecDaysThresholds[jj])
							daysToFirstFrost = ifelse(any(threshDiff == consecDaysThresholds[jj][1]), 
								atFirstFrostTemp[which(threshDiff == consecDaysThresholds[jj])[1] + (consecDaysThresholds[jj] - 1)],
								364)
							
							iter = iter + 1
							allYearsOutput_df[iter, ] = c(firstFrostThresholds[ii], consecDaysThresholds[jj], daysToFirstFrost)
						}
					}	
				}
			}
			
			

			thisLat = nc_latEra5[this_lat]
			thisLon = nc_lonEra5[this_lon]
			for(thisFirstFrostThresholds in firstFrostThresholds)				{
				for(thisConsecDaysThresholds in consecDaysThresholds)		{
					daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
							firstFrostThresholds == thisFirstFrostThresholds & 
							consecDaysThresholds == thisConsecDaysThresholds)$daysToFirstFrost,
							theseQuantiles))
					datesByQuantile = growYearStartDate + daysSinceDate

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds, 
						paste0('DaysToFirstFrost'), 'projection',
						'Q5_to_Q25', paste0(datesByQuantile[1]), paste0(datesByQuantile[2]),
						paste0(userName), paste0(forecastDate), daysSinceDate[2])

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds,
						paste0('DaysToFirstFrost'), 'projection',
						'Q25_to_Q75', paste0(datesByQuantile[2]), paste0(datesByQuantile[4]),
						paste0(userName), paste0(forecastDate), daysSinceDate[3])

					finalIter = finalIter + 1
					summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisFirstFrostThresholds, thisConsecDaysThresholds,
						paste0('DaysToFirstFrost'), 'projection',
						'Q75_to_Q95', paste0(datesByQuantile[4]), paste0(datesByQuantile[5]),
						paste0(userName), paste0(forecastDate), daysSinceDate[4])
				}
			}

			if(nrow(summaryOutput_df) > 20000)	{
				save_iter = save_iter + 1
				fwrite(summaryOutput_df, paste0(dataPath, "temp_out_", save_iter, ".csv"))
				summaryOutput_df = data.frame(Lat = 0, Lon = 0, firstFrostThresholds = 0, consecDaysThresholds = 0, 
					variableName = 'string', projectionOrClimatology = 'string',
					quantiles = 'string', startDate = 'string', endDate = 'string', User = 'string', forecastDate = 'string', numericValue = 0)
				finalIter = 0
				print(c(this_lat, this_lon))
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
	











