	#load libraries
library(lubridate)
library(ncdf4)



	# define the 'sortable' variables
theseQuantiles_ = c(0.05, 0.25, 0.5, 0.75, 0.95)
soilTempThresholds_ = c(10, 13, 16)
consecDaysThresholds_ = c(5, 10)
baseTempGDDs_ = c(10, 12)	# 10C for corn, 12C for sorghum
GDDvalueThresholds_ = seq(1350,1650,50)

	# define names and variables
growYearStartDate_ = as.Date('2022-08-15')
timeFrame_ = c('Historic', 'Projections')
dataPath_ = 'j:\\Cai_data\\Advanta\\'
tminName_ = 't2m_min_era5_Australia.nc'
tmaxName_ = 't2m_max_era5_Australia.nc'
soilTempName_ = NA

historicOutput = f_historicTempGDD()


##############
# incorporating projections data

# define names and variables
tminNameHist_ = 't2m_min_era5_Australia.nc'
tmaxNameHist_ = 't2m_max_era5_Australia.nc'
tminNameProj_ = 't2m_min_bias_corrected_curr_Australia.nc'
tmaxNameProj_ = 't2m_max_bias_corrected_curr_Australia.nc'
soilTempNameHist_ = NA
soilTempNameProj_ = NA

projectedOutput = f_projectedTempGDD()

projectedOutput[, c('Q05_hist', 'Q25_hist', 'Q50_hist', 'Q75_hist', 'Q95_hist')] = historicOutput[, c('Q05','Q25','Q50','Q75','Q95')]
projectedOutput[, 'ID'] = seq(1,nrow(projectedOutput))
projectedOutput[, 'User'] = 'Advanta'
projectedOutput[, 'forecastDate'] = ncvar_get(nc_open(paste0(dataPath_, tminNameProj_)), 'time')[1]  / 24 + as.Date('1900-01-01') 


fwrite(projectedOutput, paste0(dataPath, 'GDD_and_soilTemp_', projectedOutput[1, 'User'], projectedOutput[1, 'forecastDate'], '.csv'))










f_historicTempGDD = function(theseQuantiles = theseQuantiles_,
	soilTempThresholds = soilTempThresholds_,
	consecDaysThresholds = consecDaysThresholds_,
	baseTempGDDs = baseTempGDDs_,
	GDDvalueThresholds = GDDvalueThresholds_, 
	growYearStartDate = growYearStartDate_,
	dataPath = dataPath_,
	tminName = tminName_,
	tmaxName = tmaxName_,
	soilTempName = soilTempName_)	{

	
		# read in ncdf data
	ncin_tmax = nc_open(paste0(dataPath, tmaxName)) # t2m_max[longitude,latitude,time]   (Contiguous storage) 
	ncin_tmin = nc_open(paste0(dataPath, tminName)) # t2m_max[longitude,latitude,time]   (Contiguous storage) 
	nc_tavg = (ncvar_get(ncin_tmin, 't2m_min') + ncvar_get(ncin_tmax, 't2m_max')) / 2
	nc_date = ncvar_get(ncin_tmin, 'time') + as.Date('1979-01-01')
	nc_doy = yday(nc_date)
	nc_year = year(nc_date)
	nc_lat = ncvar_get(ncin_tmin, 'latitude')
	nc_lon = ncvar_get(ncin_tmin, 'longitude')
	
	if(yday(growYearStartDate) == 1)	{
		growYearDOYs = 1:366
	}	else	{growYearDOYs = c(yday(growYearStartDate):366, 1:(yday(growYearStartDate)-1))}
	

	if(!is.na(soilTempName))	{										# !!!!!!!!!!!!!!!! REVISIT WHEN WE HAVE SOIL TEMP DATA
		nc_soilTemp = ncvar_get(nc_open(paste0(dataPath, soilTempName)), 'WHAT IS THE VARIABLE NAME')
	}	else	{ nc_soilTemp = nc_tavg - 1.5 }


	finalIter = 0	
	summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
		variableName = 'string', Q05 = 'string', Q25='string', Q50 = 'string', Q75 = 'string', Q95 = 'string')
	for(this_lon in 1:length(nc_lon))	{
		for(this_lat in 1:length(nc_lat))	{
			theseTavg = nc_tavg[this_lon, this_lat, ]
			theseSoilTemp = nc_soilTemp[this_lon, this_lat, ]
			
			allYearsOutput_df = data.frame(soilTempThresholds = NA, consecDaysThresholds = NA, daysToWarmSoils = NA,
				baseTempGDDs = NA, GDDvalueThresholds = NA, daysToGDDThreshold = NA)
			iter = 0
			for(hh in 1:length(baseTempGDDs))	{
				theseGDD = ifelse(theseTavg > baseTempGDDs[hh], theseTavg - baseTempGDDs[hh], 0)
				for(ii in 1:length(soilTempThresholds))	{
					for(jj in 1:length(consecDaysThresholds))	{
						for(thisYear in unique(nc_year)[-1])	{
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
						finalIter = finalIter + 1
						daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
								soilTempThresholds == thisSoilTempThresholds & consecDaysThresholds == thisConsecDaysThresholds)$daysToWarmSoils,
								theseQuantiles))
						datesByQuantile = growYearStartDate + daysSinceDate
						summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
							paste0('DaysToSoilTemp'),
							paste0(datesByQuantile))
						for(thisGDDvalueThresholds in GDDvalueThresholds)	{
							finalIter = finalIter + 1
							daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
									soilTempThresholds == thisSoilTempThresholds & consecDaysThresholds == thisConsecDaysThresholds &
									GDDvalueThresholds == thisGDDvalueThresholds)$daysToGDDThreshold,
									theseQuantiles))
							datesByQuantile = growYearStartDate + daysSinceDate
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
								paste0('DaysToGDD_', thisGDDvalueThresholds),
								paste0(datesByQuantile))
						}
					}
				}
			}
		}
	}
	return(summaryOutput_df)
}	







##############
# incorporating projections data

# define names and variables
tminNameHist_ = 't2m_min_era5_Australia.nc'
tmaxNameHist_ = 't2m_max_era5_Australia.nc'
tminNameProj_ = 't2m_min_bias_corrected_curr_Australia.nc'
tmaxNameProj_ = 't2m_max_bias_corrected_curr_Australia.nc'
soilTempNameHist_ = NA
soilTempNameProj_ = NA


f_projectedTempGDD = function(theseQuantiles = theseQuantiles_,
	soilTempThresholds = soilTempThresholds_,
	consecDaysThresholds = consecDaysThresholds_,
	baseTempGDDs = baseTempGDDs_,
	GDDvalueThresholds = GDDvalueThresholds_, 
	growYearStartDate = growYearStartDate_,
	dataPath = dataPath_,
	tminNameHist = tminNameHist_,
	tmaxNameHist = tmaxNameHist_,
	soilTempNameHist = soilTempNameHist_,	
	tminNameProj = tminNameProj_,
	tmaxNameProj = tmaxNameProj_,
	soilTempNameProj = soilTempNameProj_)	{

	
		# read in ncdf data
	ncin_tmaxHist = nc_open(paste0(dataPath, tmaxNameHist)) # [longitude,latitude,time]   (Contiguous storage) 
	ncin_tminHist = nc_open(paste0(dataPath, tminNameHist)) # [longitude,latitude,time]   (Contiguous storage) 
	ncin_tmaxProj = nc_open(paste0(dataPath, tmaxNameProj)) # [longitude,latitude,member,time]   (Contiguous storage)   
	ncin_tminProj = nc_open(paste0(dataPath, tminNameProj)) # [longitude,latitude,member,time]   (Contiguous storage)   
	nc_tavgHist = (ncvar_get(ncin_tminHist, 't2m_min') + ncvar_get(ncin_tmaxHist, 't2m_max')) / 2
	nc_tavgProj = (ncvar_get(ncin_tminProj, 't2m_min') + ncvar_get(ncin_tmaxProj, 't2m_max')) / 2

	nc_dateHist = ncvar_get(ncin_tminHist, 'time') + as.Date('1979-01-01')
	nc_dateProj = ncvar_get(ncin_tminProj, 'time') / 24 + as.Date('1900-01-01') 
	nc_doyHist = yday(nc_dateHist)
	nc_doyProj = yday(nc_dateProj)
	nc_yearHist = year(nc_dateHist)
	nc_yearProj = year(nc_dateProj)
	nc_lat = ncvar_get(ncin_tminProj, 'latitude')
	nc_lon = ncvar_get(ncin_tminProj, 'longitude')
	
	whenBackfill = nc_dateProj[1] > growYearStartDate
	backfillDays = 1:lengthBackfill
	whenFrontfill = (lengthBackfill + length(theseGDDProj[thisModel, ]) + 1) < 366
	frontfillDays = (lengthBackfill + length(theseGDDProj[thisModel, ]) + 1):366
	
	if(yday(growYearStartDate) == 1)	{
		growYearDOYs = 1:366
	}	else	{growYearDOYs = c(yday(growYearStartDate):366, 1:(yday(growYearStartDate)-1))}
	

	if(!is.na(soilTempNameProj))	{										# !!!!!!!!!!!!!!!! REVISIT WHEN WE HAVE SOIL TEMP DATA
		nc_soilTempHist = ncvar_get(nc_open(paste0(dataPath, soilTempNameHist)), 'WHAT IS THE VARIABLE NAME')
		nc_soilTempProj = ncvar_get(nc_open(paste0(dataPath, soilTempNameProj)), 'WHAT IS THE VARIABLE NAME')
	}	else	{ nc_soilTempHist = nc_tavgHist - 1.5 ;   nc_soilTempProj = nc_tavgProj - 1.5}


	finalIter = 0	
	summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilTempThresholds = 0, consecDaysThresholds = 0, baseTempGDDs = 0,
		variableName = 'string', Q05 = 'string', Q25='string', Q50 = 'string', Q75 = 'string', Q95 = 'string')
	for(this_lon in 1:length(nc_lon))	{
		for(this_lat in 1:length(nc_lat))	{
			startProjGrowingSeason = which(nc_doyProj >= growYearDOYs[1])[1]
			theseTavgProj = nc_tavgProj[this_lon, this_lat, , startProjGrowingSeason:length(nc_doyProj)]
			theseSoilTempProj = nc_soilTempProj[this_lon, this_lat, , startProjGrowingSeason:length(nc_doyProj)]
			
			theseTavgHist = NULL
			theseSoilTempHist = NULL
			for(gg in growYearDOYs)	{
				thisDOY = which(nc_doyHist == gg)
				theseTavgHist = c(theseTavgHist,
					median(nc_tavgHist[this_lon, this_lat, thisDOY]))
				theseSoilTempHist = c(theseSoilTempHist,
					median(nc_soilTempHist[this_lon, this_lat, thisDOY]))
			}
			
			allYearsOutput_df = data.frame(soilTempThresholds = NA, consecDaysThresholds = NA, daysToWarmSoils = NA,
				baseTempGDDs = NA, GDDvalueThresholds = NA, daysToGDDThreshold = NA)
			iter = 0
			
			for(hh in 1:length(baseTempGDDs))	{
				theseGDDHist = ifelse(theseTavgHist > baseTempGDDs[hh], theseTavgHist - baseTempGDDs[hh], 0)
				theseGDDProj = ifelse(theseTavgProj > baseTempGDDs[hh], theseTavgProj - baseTempGDDs[hh], 0)
				for(ii in 1:length(soilTempThresholds))	{
					for(jj in 1:length(consecDaysThresholds))	{
						for(thisModel in 1:nrow(theseGDDProj))	{
							theseGDD = theseGDDProj[thisModel, ]
							theseSoilTemp = theseSoilTempProj[thisModel, ]
							
							if(whenBackfill)	{
								theseGDD = c(theseGDDHist[backfillDays], theseGDD)
								theseSoilTemp = c(theseSoilTempHist[backfillDays], theseSoilTemp)
							}
							
							if(whenFrontfill)	{
								theseGDD = c(theseGDD, theseGDDHist[frontfillDays])
								theseSoilTemp = c(theseSoilTemp, theseSoilTempHist[frontfillDays])
							}
									
							exceedsSoilThresh = which(theseSoilTemp > (soilTempThresholds[ii])) 
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
						finalIter = finalIter + 1
						daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
								soilTempThresholds == thisSoilTempThresholds & consecDaysThresholds == thisConsecDaysThresholds)$daysToWarmSoils,
								theseQuantiles))
						datesByQuantile = growYearStartDate + daysSinceDate
						summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
							paste0('DaysToSoilTemp'),
							paste0(datesByQuantile))
						for(thisGDDvalueThresholds in GDDvalueThresholds)	{
							finalIter = finalIter + 1
							daysSinceDate = ceiling(quantile(subset(allYearsOutput_df,
									soilTempThresholds == thisSoilTempThresholds & consecDaysThresholds == thisConsecDaysThresholds &
									GDDvalueThresholds == thisGDDvalueThresholds)$daysToGDDThreshold,
									theseQuantiles))
							datesByQuantile = growYearStartDate + daysSinceDate
							summaryOutput_df[finalIter, ] = c(thisLat, thisLon, thisSoilTempThresholds, thisConsecDaysThresholds, thisBaseTempGDDs,
									paste0('DaysToGDD_', thisGDDvalueThresholds),
									paste0(datesByQuantile))
							}
					}
				}
			}
		}
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
	











