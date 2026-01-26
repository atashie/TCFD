

############################################################################################
## basic exposure data
f_mainExposureData = function(
	customerTable = customerTable,
	customerFolder = customerFolder,
	hazardFolder = hazardFolder,
	hazardTable = hazardTable,
	relHazScores = relHazScores,
	appendedHazardNames = appendedHazardNames,
	thisDate = thisDate)
	{
	
	dataOutput = data.frame(User = NA, Location = NA, Region = NA, Subregion = NA, Lat = NA, Lon = NA,
		Hazard = NA, Hazard_Measure = NA, Decade = NA, Scenario = NA,
		Raw_Hazard_Value = NA, Percentile_Score = NA, Relative_Hazard_Score = NA, Decadal_Trend_Strength = NA, Decadal_Trend_Significance = NA, Long_Term_Trend_Strength = NA, Long_Term_Trend_Significance = NA,
		Relative_Hazard_Score_Number = NA, Trend_Aggregated_For_Looker = NA, Advanced_Data_Measures = NA, Advanced_Data_Measures_Units = NA, Raw_Hazard_Value_25th = NA, Raw_Hazard_Value_75th =  NA,
		Asset_Type = NA, Business_Unit = NA, Country = NA, State = NA, City = NA)

	iter = 0
	for(thisHazard in 11:ncol(customerTable))	{ ####!!!! temp fix, resolve hard coding  !!!!!!##########
		if(any(customerTable[, ..thisHazard]))	{
			if(names(customerTable)[thisHazard] %in% appendedHazardNames)	{	
				print(c('this hazard must be run independently and appended later:     ', names(customerTable)[thisHazard]))
			}	else	{
				hazardMeasures = subset(hazardTable, Hazard == names(customerTable)[thisHazard] & Included_In_Core_Product)
				for(thisHazardMeasure in 1:nrow(hazardMeasures))	{
					hazardMeasureNC = nc_open(paste0(hazardFolder, hazardMeasures$Hazard_Measure[thisHazardMeasure], 'v2_processed.nc'))
					nc_lats = ncvar_get(hazardMeasureNC, 'lat')
					nc_lons = ncvar_get(hazardMeasureNC, 'lon')
					nc_decades = 2000 + ncvar_get(hazardMeasureNC, 'decade')
								
					for(thisLocation in 1:nrow(customerTable))	{

						if(unlist(customerTable[thisLocation, ..thisHazard]))	{
							closeLon = which.min(abs(nc_lons - customerTable$Lon[thisLocation]))
							closeLat = which.min(abs(nc_lats - customerTable$Lat[thisLocation]))

								# temp fix: ensuring we aren't drawing from a coastal water cell
								#### !!!!!!!!!!!!! #######
								# to be replaced with the updated version developed for the water index
								#### !!!!!!!!!!!!! #######
							if(is.na(ncvar_get(hazardMeasureNC, 'tcfdVariable')[closeLon, closeLat, 1, 1, 2]))	{
								closeLat = closeLat + 1
								if(is.na(ncvar_get(hazardMeasureNC, 'tcfdVariable')[closeLon, closeLat, 1, 1, 2]))	{
									closeLat = closeLat - 2
									if(is.na(ncvar_get(hazardMeasureNC, 'tcfdVariable')[closeLon, closeLat, 1, 1, 2]))	{
										closeLat = closeLat + 1
										closeLon = closeLon + 1
										if(is.na(ncvar_get(hazardMeasureNC, 'tcfdVariable')[closeLon, closeLat, 1, 1, 2]))	{
											closeLon = closeLon - 2
										}
									}
								}
							}
			
							theseData = ncvar_get(hazardMeasureNC, 'tcfdVariable')[closeLon, closeLat, , , ]
							
							
							dataOutput = rbind(dataOutput,
								data.frame(
									User = userName,
									Location = customerTable$Location[thisLocation],
									Region = customerTable$Region[thisLocation],
									Subregion = customerTable$Subregion[thisLocation],
									Lat = customerTable$Lat[thisLocation],
									Lon = customerTable$Lon[thisLocation],
									Hazard = hazardMeasures$Hazard[thisHazardMeasure],
									Hazard_Measure = paste0(hazardMeasures$Hazard_Measure_Common_Name[thisHazardMeasure], ' ', hazardMeasures$Hazard_Measure_Units[thisHazardMeasure]),
									Decade = rep(nc_decades, 3),
									Scenario = rep(c('RCP 2.6', 'RCP 6.0', 'RCP 8.5'), each = 9),
									Raw_Hazard_Value = c(theseData[ , 1, 1], theseData[ , 2, 1], theseData[ , 3, 1]),
									Percentile_Score = c(theseData[ , 1, 2], theseData[ , 2, 2], theseData[ , 3, 2]),
									Relative_Hazard_Score = NA,
									Decadal_Trend_Strength = 	 	c(theseData[ , 1, 3], theseData[ , 2, 3], theseData[ , 3, 3]),
									Decadal_Trend_Significance = 	c(theseData[ , 1, 4], theseData[ , 2, 4], theseData[ , 3, 4]),
									Long_Term_Trend_Strength = 		c(rep(theseData[9, 1, 3], 9), rep(theseData[9, 2, 3], 9), rep(theseData[9, 3, 3], 9)),
									Long_Term_Trend_Significance =	c(rep(theseData[9, 1, 4], 9), rep(theseData[9, 2, 4], 9), rep(theseData[9, 3, 4], 9)),
									Relative_Hazard_Score_Number = NA,
									Trend_Aggregated_For_Looker = NA,
									Advanced_Data_Measures = NA,
									Advanced_Data_Measures_Units = NA,
									Raw_Hazard_Value_25th = c(theseData[ , 1, 5], theseData[ , 2, 5], theseData[ , 3, 5]),
									Raw_Hazard_Value_75th = c(theseData[ , 1, 6], theseData[ , 2, 6], theseData[ , 3, 6]),
									Asset_Type = customerTable$AssetType[thisLocation],
									Business_Unit = customerTable$BusinessUnit[thisLocation],
									Country = customerTable$Country[thisLocation],
									State = customerTable$State[thisLocation],
									City = customerTable$City[thisLocation])
							)
						}
						print(c(thisHazard, thisHazardMeasure, thisLocation))
					}
					nc_close(hazardMeasureNC)
					iter = iter + 1
					fwrite(dataOutput[-1,], paste0(customerFolder, 'temp_out_', iter, '.csv'))
					dataOutput = data.frame(User = NA, Location = NA, Region = NA, Subregion = NA, Lat = NA, Lon = NA,
						Hazard = NA, Hazard_Measure = NA, Decade = NA, Scenario = NA,
						Raw_Hazard_Value = NA, Percentile_Score = NA, Relative_Hazard_Score = NA, Decadal_Trend_Strength = NA, Decadal_Trend_Significance = NA, Long_Term_Trend_Strength = NA, Long_Term_Trend_Significance = NA,
						Relative_Hazard_Score_Number = NA, Trend_Aggregated_For_Looker = NA, Advanced_Data_Measures = NA, Advanced_Data_Measures_Units = NA, Raw_Hazard_Value_25th = NA, Raw_Hazard_Value_75th =  NA,
						Asset_Type = NA, Business_Unit = NA, Country = NA, State = NA, City = NA)
				}
			}
		}
	}
	finalOutput = fread(paste0(customerFolder, 'temp_out_', 1, '.csv'))
	for(thisIter in 2:iter)	{
		finalOutput = rbind(finalOutput, fread(paste0(customerFolder, 'temp_out_', thisIter, '.csv')))
	}

	fwrite(finalOutput, paste0(customerFolder, 'processedOutput_', thisDate, '.csv'))
}












###################################################################################################################################################################################################
#### specialized data, for now just localized flooding and coastal flooding
# step 3a intersecting recurrence intervals with flood depth
f_riverFloods = function(
	customerTable = customerTable,
	customerFolder = customerFolder,
	userName = userName,
	locationFootprint = locationFootprint,
	ncpath = ncpathRiverFloods,
	fileLoc = floodMapTiffLoc,
	waterMaskLoc = waterMaskLoc,
	thisDate = thisDate)
	{
		
	library(raster)

		# defining rcp scenarios, recurrence intervals, and decades of interest; this should not change for the foreseeable future
	rcpScenarios = c('RCP 2.6', 'RCP 6.0', 'RCP 8.5')
	recurIntrvls = c(10, 20, 50, 100, 200, 500)
	whichDecades = seq(10,90,10)

		
		# read in flood recurrence data
	fldRcrIntNC = nc_open(paste0(ncpath, 'floodRecurIntervals_v5.nc'))
	nc_lat = ncvar_get(fldRcrIntNC, 'lat')
	nc_lon = ncvar_get(fldRcrIntNC, 'lon')
	fldRcrVals = ncvar_get(fldRcrIntNC, 'fldRecurIntrvl')


		# initializing data output
	dataOutput = data.frame(User = NA, Location = NA, Region = NA, Subregion = NA, Lat = NA, Lon = NA,
		Hazard = NA, Hazard_Measure = NA, Decade = NA, Scenario = NA,
		Raw_Hazard_Value = NA, Percentile_Score = NA, Relative_Hazard_Score = NA, Decadal_Trend_Strength = NA, Decadal_Trend_Significance = NA, Long_Term_Trend_Strength = NA, Long_Term_Trend_Significance = NA,
		Relative_Hazard_Score_Number = NA, Trend_Aggregated_For_Looker = NA, Advanced_Data_Measures = NA, Advanced_Data_Measures_Units = NA, Raw_Hazard_Value_25th = NA, Raw_Hazard_Value_75th =  NA,
		Asset_Type = NA, Business_Unit = NA, Country = NA, State = NA, City = NA)
		
		# read in historic floods data
	fldDepthList = list()
	fldDepthList[[1]] = raster(paste0(fileLoc, '\\floodMapGL_rp10y.tif'))
	fldDepthList[[2]] = raster(paste0(fileLoc, '\\floodMapGL_rp20y.tif'))
	fldDepthList[[3]] = raster(paste0(fileLoc, '\\floodMapGL_rp50y.tif'))
	fldDepthList[[4]] = raster(paste0(fileLoc, '\\floodMapGL_rp100y.tif'))
	fldDepthList[[5]] = raster(paste0(fileLoc, '\\floodMapGL_rp200y.tif'))
	fldDepthList[[6]] = raster(paste0(fileLoc, '\\floodMapGL_rp500y.tif'))
	myWaterMask = readRDS(waterMaskLoc)

		# identifying lat lon coordinates for tiffs
	tif_lat = rev(seq(-54.00824,83.2251,length.out=16468))
	tif_lon = seq(-166.8, 180, length.out=41616)

	hazardDepthName = paste0("Avg Flood Depth (m)")
	hazardLikliName = paste0('Flood Areal Extent (%)')


	for(thisIntrvl in 1:length(recurIntrvls))	{
		for(thisLoc in 1:nrow(customerTable))	{
			closeTiffLons = which.min(abs(customerTable$Lon[thisLoc] - tif_lon))
			closeTiffLats = which.min(abs(customerTable$Lat[thisLoc] - tif_lat))
			
			theseLats = seq(closeTiffLats - locationFootprint, closeTiffLats + locationFootprint, 1)
			theseLons = seq(closeTiffLons - locationFootprint, closeTiffLons + locationFootprint, 1)
				#ensuring lats / lons don't go outside bounding box
			if(any(theseLats < 1))	{theseLats[theseLats < 1] = 1}
			if(any(theseLats > length(tif_lat)))	{theseLats[theseLats > length(tif_lat)] = length(tif_lat)}
			if(any(theseLons < 1))	{theseLons[theseLons < 1] = 1}
			if(any(theseLons > length(tif_lon)))	{theseLons[theseLons > length(tif_lon)] = length(tif_lon)}

				# defining the water mask
			thisWaterMask = myWaterMask[theseLons, theseLats]

			histFloodDepth = mean(fldDepthList[[6]][theseLats, theseLons], na.rm=TRUE)
			histFloodLikli = 100 * length(which(fldDepthList[[6]][theseLats, theseLons] > 0)) / length(theseLats)^2
			
				# check to see if even 500yr floods trigger historically; if not, the skip next analysis
			if(!is.na(histFloodDepth))	{
					# ensuring we are not drawing from a water tile, then searching box around point if so
			#### !!!! todo: implement new smoothing scheme !!!!!

				closeNCLon = which.min(abs(customerTable$Lon[thisLoc] - nc_lon))
				closeNCLat = which.min(abs(customerTable$Lat[thisLoc] - nc_lat))
				if(is.na(fldRcrVals[closeNCLon, closeNCLat, 1, thisIntrvl, 1, 1]))	{
						closeNCLon = closeNCLon + 1
						if(is.na(fldRcrVals[closeNCLon, closeNCLat, 1, thisIntrvl, 1, 1]))	{
							closeNCLon = closeNCLon - 2
							if(is.na(fldRcrVals[closeNCLon, closeNCLat, 1, thisIntrvl, 1, 1]))	{
								closeNCLon = closeNCLon + 1
								closeNCLat = closeNCLat + 1
								if(is.na(fldRcrVals[closeNCLon, closeNCLat, 1, thisIntrvl, 1, 1]))	{
									closeNCLat = closeNCLat - 2
				}}}}
					
				for(thisScenario in 1:length(rcpScenarios))	{
					fldRcrSignif_all = fldRcrVals[closeNCLon, closeNCLat, 9,  thisIntrvl, 2, thisScenario]
	#				fldRcrDrctn_all = fldRcrVals[closeNCLon, closeNCLat, 9,  thisIntrvl, 3, thisScenario]
					
					for(thisDecade in 1:length(whichDecades))	{
						thisFldRcrVal = 		fldRcrVals[closeNCLon, closeNCLat, thisDecade, thisIntrvl, 1, thisScenario]
						thisFldRcrSignif = 		fldRcrVals[closeNCLon, closeNCLat, thisDecade, thisIntrvl, 2, thisScenario]
	#					thisFldRcrDrctn = 		fldRcrVals[closeNCLon, closeNCLat, thisDecade, thisIntrvl, 3, thisScenario]
						
						if(thisFldRcrVal >= recurIntrvls[1])	{
							closestFldRcrIntrvl = 	last(which(thisFldRcrVal >= recurIntrvls))
							
							theseFloodImpacts = fldDepthList[[closestFldRcrIntrvl]][theseLats, theseLons]
							theseFloodImpacts[is.na(theseFloodImpacts)] = 0
							theseFloodImpacts = theseFloodImpacts * thisWaterMask
		
							prevFloodImpacts = fldDepthList[[thisIntrvl]][theseLats, theseLons]
							prevFloodImpacts[is.na(prevFloodImpacts)] = 0
							prevFloodImpacts = prevFloodImpacts * thisWaterMask
		
							
							if(any(!is.na(theseFloodImpacts)))	{
								avgFloodDepth = mean(theseFloodImpacts, na.rm=TRUE)
								avgFloodDepthChng = avgFloodDepth - mean(prevFloodImpacts, na.rm=TRUE) / max(c(1, thisDecade - 1))

								avgFloodLikli = 100 * (length(which(theseFloodImpacts > 0)) / length(!is.na(theseFloodImpacts)))
								avgFloodLikliChng = avgFloodLikli - 100 * (length(which(prevFloodImpacts > 0)) / length(!is.na(prevFloodImpacts))) / max(c(1, thisDecade - 1))
							

								dataOutput = rbind(dataOutput,
			#### !!!! todo: consolidate and simplify !!!!!
									data.table(
										User = userName,
										Location = customerTable$Location[thisLoc],
										Region = customerTable$Region[thisLoc],
										Subregion = customerTable$Subregion[thisLoc],
										Lat = customerTable$Lat[thisLoc],
										Lon = customerTable$Lon[thisLoc],
										Hazard = "River Flood (Local)",
										Hazard_Measure = c(hazardDepthName, hazardLikliName),
										Decade = 2000 + whichDecades[thisDecade],
										Scenario = rcpScenarios[thisScenario],
										Raw_Hazard_Value = c(avgFloodDepth, avgFloodLikli),					# Raw_Hazard_Value
										Percentile_Score = NA,												# Percentile_Score
										Relative_Hazard_Score = NA,											# Relative_Hazard_Score
										Decadal_Trend_Strength = c(avgFloodDepthChng, avgFloodLikliChng),	# Decadal_Trend_Strength
										Decadal_Trend_Significance = thisFldRcrSignif,						# Decadal_Trend_Significance
										Long_Term_Trend_Strength = NA,										# Long_Term_Trend_Strength
										Long_Term_Trend_Significance = fldRcrSignif_all,					# Long_Term_Trend_Significance					
										Relative_Hazard_Score_Number = NA,
										Trend_Aggregated_For_Looker = NA,
										Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
										Advanced_Data_Measures_Units = "Yr Flood",
										Raw_Hazard_Value_25th = NA,
										Raw_Hazard_Value_75th = NA,
										Asset_Type = customerTable$AssetType[thisLoc],
										Business_Unit = customerTable$BusinessUnit[thisLoc],
										Country = customerTable$Country[thisLoc],
										State = customerTable$State[thisLoc],
										City = customerTable$City[thisLoc]))


							} 	else	{ # if all theseFloodImpacts are NAs
								dataOutput = rbind(dataOutput,
									data.table(
										User = userName,
										Location = customerTable$Location[thisLoc],
										Region = customerTable$Region[thisLoc],
										Subregion = customerTable$Subregion[thisLoc],
										Lat = customerTable$Lat[thisLoc],
										Lon = customerTable$Lon[thisLoc],
										Hazard = "River Flood (Local)",
										Hazard_Measure = c(hazardDepthName, hazardLikliName),
										Decade = 2000 + whichDecades[thisDecade],
										Scenario = rcpScenarios[thisScenario],
										Raw_Hazard_Value = 0,								# Raw_Hazard_Value
										Percentile_Score = NA,								# Percentile_Score
										Relative_Hazard_Score = NA,							# Relative_Hazard_Score
										Decadal_Trend_Strength = 0,							# Decadal_Trend_Strength
										Decadal_Trend_Significance = thisFldRcrSignif,		# Decadal_Trend_Significance
										Long_Term_Trend_Strength = NA,						# Long_Term_Trend_Strength
										Long_Term_Trend_Significance = fldRcrSignif_all,	# Long_Term_Trend_Significance					
										Relative_Hazard_Score_Number = NA,
										Trend_Aggregated_For_Looker = NA,
										Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
										Advanced_Data_Measures_Units = "Yr Flood",
										Raw_Hazard_Value_25th = NA,
										Raw_Hazard_Value_75th = NA,
										Asset_Type = customerTable$AssetType[thisLoc],
										Business_Unit = customerTable$BusinessUnit[thisLoc],
										Country = customerTable$Country[thisLoc],
										State = customerTable$State[thisLoc],
										City = customerTable$City[thisLoc]))

							}
						}	else	{ # if thisFldRcrVal < 10 yr recurrence
							dataOutput = rbind(dataOutput,
								data.table(
									User = userName,
									Location = customerTable$Location[thisLoc],
									Region = customerTable$Region[thisLoc],
									Subregion = customerTable$Subregion[thisLoc],
									Lat = customerTable$Lat[thisLoc],
									Lon = customerTable$Lon[thisLoc],
									Hazard = "River Flood (Local)",
									Hazard_Measure = c(hazardDepthName, hazardLikliName),
									Decade = 2000 + whichDecades[thisDecade],
									Scenario = rcpScenarios[thisScenario],
									Raw_Hazard_Value = 0,								# Raw_Hazard_Value
									Percentile_Score = NA,								# Percentile_Score
									Relative_Hazard_Score = NA,							# Relative_Hazard_Score
									Decadal_Trend_Strength = 0,							# Decadal_Trend_Strength
									Decadal_Trend_Significance = thisFldRcrSignif,		# Decadal_Trend_Significance
									Long_Term_Trend_Strength = NA,						# Long_Term_Trend_Strength
									Long_Term_Trend_Significance = fldRcrSignif_all,	# Long_Term_Trend_Significance
									Relative_Hazard_Score_Number = NA,
									Trend_Aggregated_For_Looker = NA,
									Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
									Advanced_Data_Measures_Units = "Yr Flood",
									Raw_Hazard_Value_25th = NA,
									Raw_Hazard_Value_75th = NA,
									Asset_Type = customerTable$AssetType[thisLoc],
									Business_Unit = customerTable$BusinessUnit[thisLoc],
									Country = customerTable$Country[thisLoc],
									State = customerTable$State[thisLoc],
									City = customerTable$City[thisLoc]))

						}
					} # for each decade

				dataOutput[which(dataOutput$Location == customerTable$Location[thisLoc] & dataOutput$Scenario == rcpScenarios[thisScenario]), ]$Long_Term_Trend_Strength = last(dataOutput)$Decadal_Trend_Strength

				} # for each rcp scenario
			}	else	{	# if there are not historical floods out to 500 yrs
			#### !!!! todo: consolidate and simplify !!!!!
				dataOutput = rbind(rbind(rbind(dataOutput,
					data.table(
						User = userName,
						Location = customerTable$Location[thisLoc],
						Region = customerTable$Region[thisLoc],
						Subregion = customerTable$Subregion[thisLoc],
						Lat = customerTable$Lat[thisLoc],
						Lon = customerTable$Lon[thisLoc],
						Hazard = "River Flood (Local)",
						Hazard_Measure = c(hazardDepthName, hazardLikliName),
						Decade = 2000 + rep(whichDecades, each = 2),
						Scenario = rcpScenarios[1],
						Raw_Hazard_Value = 0,				# Raw_Hazard_Value
						Percentile_Score = NA,				# Percentile_Score
						Relative_Hazard_Score = NA,			# Relative_Hazard_Score
						Decadal_Trend_Strength = NA,			# Decadal_Trend_Strength
						Decadal_Trend_Significance = NA,		# Decadal_Trend_Significance
						Long_Term_Trend_Strength = NA,		# Long_Term_Trend_Strength
						Long_Term_Trend_Significance = NA,	# Long_Term_Trend_Significance
						Relative_Hazard_Score_Number = NA,
						Trend_Aggregated_For_Looker = NA,
						Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
						Advanced_Data_Measures_Units = "Yr Flood",
						Raw_Hazard_Value_25th = NA,
						Raw_Hazard_Value_75th = NA,
						Asset_Type = customerTable$AssetType[thisLoc],
						Business_Unit = customerTable$BusinessUnit[thisLoc],
						Country = customerTable$Country[thisLoc],
						State = customerTable$State[thisLoc],
						City = customerTable$City[thisLoc])), 
					data.table(
						User = userName,
						Location = customerTable$Location[thisLoc],
						Region = customerTable$Region[thisLoc],
						Subregion = customerTable$Subregion[thisLoc],
						Lat = customerTable$Lat[thisLoc],
						Lon = customerTable$Lon[thisLoc],
						Hazard = "River Flood (Local)",
						Hazard_Measure = c(hazardDepthName, hazardLikliName),
						Decade = 2000 + rep(whichDecades, each = 2),
						Scenario = rcpScenarios[2],
						Raw_Hazard_Value = 0,				# Raw_Hazard_Value
						Percentile_Score = NA,				# Percentile_Score
						Relative_Hazard_Score = NA,			# Relative_Hazard_Score
						Decadal_Trend_Strength = NA,			# Decadal_Trend_Strength
						Decadal_Trend_Significance = NA,		# Decadal_Trend_Significance
						Long_Term_Trend_Strength = NA,		# Long_Term_Trend_Strength
						Long_Term_Trend_Significance = NA,	# Long_Term_Trend_Significance
						Relative_Hazard_Score_Number = NA,
						Trend_Aggregated_For_Looker = NA,
						Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
						Advanced_Data_Measures_Units = "Yr Flood",
						Raw_Hazard_Value_25th = NA,
						Raw_Hazard_Value_75th = NA,
						Asset_Type = customerTable$AssetType[thisLoc],
						Business_Unit = customerTable$BusinessUnit[thisLoc],
						Country = customerTable$Country[thisLoc],
						State = customerTable$State[thisLoc],
						City = customerTable$City[thisLoc])), 
					data.table(
						User = userName,
						Location = customerTable$Location[thisLoc],
						Region = customerTable$Region[thisLoc],
						Subregion = customerTable$Subregion[thisLoc],
						Lat = customerTable$Lat[thisLoc],
						Lon = customerTable$Lon[thisLoc],
						Hazard = "River Flood (Local)",
						Hazard_Measure = c(hazardDepthName, hazardLikliName),
						Decade = 2000 + rep(whichDecades, each = 2),
						Scenario = rcpScenarios[3],
						Raw_Hazard_Value = 0,				# Raw_Hazard_Value
						Percentile_Score = NA,				# Percentile_Score
						Relative_Hazard_Score = NA,			# Relative_Hazard_Score
						Decadal_Trend_Strength = NA,			# Decadal_Trend_Strength
						Decadal_Trend_Significance = NA,		# Decadal_Trend_Significance
						Long_Term_Trend_Strength = NA,		# Long_Term_Trend_Strength
						Long_Term_Trend_Significance = NA,	# Long_Term_Trend_Significance
						Relative_Hazard_Score_Number = NA,
						Trend_Aggregated_For_Looker = NA,
						Advanced_Data_Measures =  recurIntrvls[thisIntrvl],
						Advanced_Data_Measures_Units = "Yr Flood",
						Raw_Hazard_Value_25th = NA,
						Raw_Hazard_Value_75th = NA,
						Asset_Type = customerTable$AssetType[thisLoc],
						Business_Unit = customerTable$BusinessUnit[thisLoc],
						Country = customerTable$Country[thisLoc],
						State = customerTable$State[thisLoc],
						City = customerTable$City[thisLoc]))

			}
			print(c(thisLoc, thisIntrvl))
		}
	}	

	dataOutput$Raw_Hazard_Value = as.numeric(dataOutput$Raw_Hazard_Value)
	dataOutput = dataOutput[-1,]

	hazardDepthNames =  c("Avg Flood Depth (m)")
	hazardExtentNames = c("Flood Areal Extent (%)")
		# defining relative flood hazard
	basSeqDepth = c(0,seq(0.01, 10, length.out=99))
	basSeqLikli = c(0,seq(1, 100, length.out=99))
	relFloodHazardDepth = (basSeqDepth^2 / max(basSeqDepth^2)) * 10 + min(basSeqDepth)
	relFloodHazardLikli = basSeqLikli
	#relFloodHazard = c(rep(0, 33), seq(0.01,10,length.out=(67)))
	dataOutput$Percentile_Score = 1
	depthRows = which(dataOutput$Hazard_Measure %in% hazardDepthNames)
	likliRows = which(dataOutput$Hazard_Measure %in% hazardExtentNames)
	for(ll in 1:length(relFloodHazardDepth))	{
		dataOutput$Percentile_Score[depthRows][which(dataOutput$Raw_Hazard_Value[depthRows] > relFloodHazardDepth[ll])] = ll + 0
		dataOutput$Percentile_Score[likliRows][which(dataOutput$Raw_Hazard_Value[likliRows] > relFloodHazardLikli[ll])] = ll + 0
	}

	fileName_localFlood = paste0(userName, '_', thisDate, '_highRestFlood_', (locationFootprint*2+1), 'km')
	fwrite(dataOutput, paste0(customerFolder, fileName_localFlood, thisDate, '.csv'))
	}








#########
# step 3b intersecting sea level rise with DEMs
#########################################
# reading in climate netcdf data
f_seaLevelRise = function(
	customerTable = customerTable,
	customerFolder = customerFolder,
	userName = userName,
	ncpathDEM = ncpathDEM,
	ncpathSeaLevelRise = ncpathSeaLevelRise,
	thisDate = thisDate)
	{
	
	ncin = nc_open(paste0(ncpathDEM))
	nc_lon = ncvar_get(ncin, 'lon')
	nc_lat = rev(ncvar_get(ncin, 'lat'))
	nc_close(ncin)
	nc_elev = brick(paste0(ncpathDEM), 'elev')

	seaL_nc = nc_open(paste0(ncpathSeaLevelRise, '_glblGrdv2_processed.nc'))
	seaL_lat = ncvar_get(seaL_nc, 'lat')	# lat is given from high to low
	seaL_lon = ncvar_get(seaL_nc, 'lon')
	seaL_scen = ncvar_get(seaL_nc, 'rcpScen')
	seal_scenRename = c('RCP 2.6', 'RCP 6.0', 'RCP 8.5')
	seaL_elev = ncvar_get(seaL_nc, 'tcfdVariable')
	nc_close(seaL_nc)

	whichDecades = seq(10,90,10)

	# initializing data output
	dataOutput = data.frame(User = NA, Location = NA, Region = NA, Subregion = NA, Lat = NA, Lon = NA,
		Hazard = NA, Hazard_Measure = NA, Decade = NA, Scenario = NA,
		Raw_Hazard_Value = NA, Percentile_Score = NA, Relative_Hazard_Score = NA, Decadal_Trend_Strength = NA, Decadal_Trend_Significance = NA, Long_Term_Trend_Strength = NA, Long_Term_Trend_Significance = NA,
		Relative_Hazard_Score_Number = NA, Trend_Aggregated_For_Looker = NA, Advanced_Data_Measures = NA, Advanced_Data_Measures_Units = NA, Raw_Hazard_Value_25th = NA, Raw_Hazard_Value_75th =  NA,
		Asset_Type = NA, Business_Unit = NA, Country = NA, State = NA, City = NA)

	for(j in 1:nrow(customerTable))	{

		closeSeaLons = order(abs(customerTable$Lon[j] - seaL_lon))[1:3]
		closeSeaLats = order(abs(customerTable$Lat[j] - seaL_lat))[1:3]
		
		closeDemLons = order(abs(customerTable$Lon[j] - nc_lon))[1]
		closeDemLats = order(abs(customerTable$Lat[j] - nc_lat))[1]

		if(any(!is.na(seaL_elev[closeSeaLons, closeSeaLats, , , 1])))	{

			#### !!!! todo: implement new smoothing scheme!!!!!
				# ensuring dem is on a land tile
			avgElev = nc_elev[closeDemLats, closeDemLons]
			if(is.na(avgElev))	{
				avgElev = mean(nc_elev[seq(closeDemLats -1,closeDemLats+1,1), seq(closeDemLons -1, closeDemLons+1)], na.rm = TRUE)
			}
			if(is.na(avgElev))	{
				avgElev = mean(nc_elev[seq(closeDemLats-2,closeDemLats+2,1), seq(closeDemLons-2, closeDemLons+2)], na.rm = TRUE)
			}

			for(thisScen in 1:length(seaL_scen))	{
				seaLvlTrend_all =  mean(seaL_elev[closeSeaLons, closeSeaLats, 9, thisScen, 3], na.rm=TRUE) / 1000
				seaLvlSignif_all = mean(seaL_elev[closeSeaLons, closeSeaLats, 9, thisScen, 4], na.rm=TRUE)
						
				for(thisDecade in 1:length(whichDecades))	{
					seaLvlVal =	   avgElev - mean(seaL_elev[closeSeaLons, closeSeaLats, thisDecade, thisScen, 1], na.rm=TRUE) / 1000
					seaLvlTrend =  mean(seaL_elev[closeSeaLons, closeSeaLats, thisDecade, thisScen, 3], na.rm=TRUE) / 1000
					seaLvlSignif = mean(seaL_elev[closeSeaLons, closeSeaLats, thisDecade, thisScen, 4], na.rm=TRUE)
					seaLvl25th =   avgElev - mean(seaL_elev[closeSeaLons, closeSeaLats, thisDecade, thisScen, 5], na.rm=TRUE) / 1000
					seaLvl75th =   avgElev - mean(seaL_elev[closeSeaLons, closeSeaLats, thisDecade, thisScen, 6], na.rm=TRUE) / 1000
				

					dataOutput = rbind(dataOutput,
								data.table(
									User = userName,
									Location = customerTable$Location[j],
									Region = customerTable$Region[j],
									Subregion = customerTable$Subregion[j],
									Lat = customerTable$Lat[j],
									Lon = customerTable$Lon[j],
									Hazard = "Coastal Flood",
									Hazard_Measure = "Sea Level Rise",
									Decade = 2000 + whichDecades[thisDecade],
									Scenario = seal_scenRename[thisScen],
									Raw_Hazard_Value = as.numeric(seaLvlVal),					# Raw_Hazard_Value
									Percentile_Score = NA,												# Percentile_Score
									Relative_Hazard_Score = NA,											# Relative_Hazard_Score
									Decadal_Trend_Strength = seaLvlTrend,	# Decadal_Trend_Strength
									Decadal_Trend_Significance = seaLvlSignif,						# Decadal_Trend_Significance
									Long_Term_Trend_Strength = seaLvlTrend_all,										# Long_Term_Trend_Strength
									Long_Term_Trend_Significance = seaLvlSignif_all,					# Long_Term_Trend_Significance					
									Relative_Hazard_Score_Number = NA,
									Trend_Aggregated_For_Looker = NA,
									Advanced_Data_Measures =  NA,
									Advanced_Data_Measures_Units = NA,
									Raw_Hazard_Value_25th = as.numeric(seaLvl25th),
									Raw_Hazard_Value_75th = as.numeric(seaLvl75th),
									Asset_Type = customerTable$AssetType[j],
									Business_Unit = customerTable$BusinessUnit[j],
									Country = customerTable$Country[j],
									State = customerTable$State[j],
									City = customerTable$City[j]))
				}
			}
		}	else	{
			dataOutput = rbind(dataOutput,
						data.table(
							User = userName,
							Location = customerTable$Location[j],
							Region = customerTable$Region[j],
							Subregion = customerTable$Subregion[j],
							Lat = customerTable$Lat[j],
							Lon = customerTable$Lon[j],
							Hazard = "Coastal Flood",
							Hazard_Measure = "Sea Level Rise",
							Decade = rep(2000 + whichDecades, length(seaL_scen)),
							Scenario = rep(seal_scenRename, each=length(whichDecades)),
							Raw_Hazard_Value = 1000,								# Raw_Hazard_Value
							Percentile_Score = NA,								# Percentile_Score
							Relative_Hazard_Score = NA,							# Relative_Hazard_Score
							Decadal_Trend_Strength = 0,							# Decadal_Trend_Strength
							Decadal_Trend_Significance = 1,		# Decadal_Trend_Significance
							Long_Term_Trend_Strength = 0,						# Long_Term_Trend_Strength
							Long_Term_Trend_Significance = 1,	# Long_Term_Trend_Significance
							Relative_Hazard_Score_Number = NA,
							Trend_Aggregated_For_Looker = NA,
							Advanced_Data_Measures = NA,
							Advanced_Data_Measures_Units = NA,
							Raw_Hazard_Value_25th = 1000,
							Raw_Hazard_Value_75th = 1000,
							Asset_Type = customerTable$AssetType[j],
							Business_Unit = customerTable$BusinessUnit[j],
							Country = customerTable$Country[j],
							State = customerTable$State[j],
							City = customerTable$City[j]))
		}
	}
		
	dataOutput$Raw_Hazard_Value = as.numeric(dataOutput$Raw_Hazard_Value)
	dataOutput = dataOutput[-1,]


	# in Key West, buildings must be 7 ft above sea level (or on piers of that height): https://www.cityofkeywest-fl.gov/671/Construction-in-Flood-Zones
		# 7 ft ~ 2.1336
	basSeqSeaLevElev = seq(0, 1, length.out=100)
#	relSeaLevElevHaz = -(basSeqSeaLevElev^2 / max(basSeqSeaLevElev^2)) * (2 + 2.1336 * 4) + 2.1336 * 4
	relSeaLevElevHaz = -(basSeqSeaLevElev) * (2 + 2.1336 * 4) + 2.1336 * 4

	dataOutput$Percentile_Score = 1
	for(ll in 1:length(relSeaLevElevHaz))	{
		dataOutput$Percentile_Score[which(dataOutput$Raw_Hazard_Value < relSeaLevElevHaz[ll])] = ll + 0
	}

	fileName_seaLevelRise = paste0(userName, '_', thisDate, '_seaLevelRise_1km')
	fwrite(dataOutput, paste0(customerFolder, fileName_seaLevelRise, thisDate, '.csv'))
	}






#################











###################################################################################################################################################################################################
#### calculating aggregated score and relative hazard values
f_hazardAggregation = function(
	customerTable = customerTable,
	customerFolder = customerFolder,
	thisDate = thisDate,
	userName = userName,
	addThirdPartyHazards = addThirdPartyHazards,
	appendedHazardFileLoc = appendedHazardFileLoc,
	locationFootprint = locationFootprint,
	hazardWeighting = hazardWeighting)
	{
	
	if(!addThirdPartyHazards)	{
		print("only running basic hazards")
		dataOutput = fread(paste0(customerFolder, 'processedOutput_', thisDate, '.csv'))
	} else	{
		mainHazardOutputs = fread(paste0(customerFolder, 'processedOutput_', thisDate, '.csv'))
		appendedOutputs = subset(fread(appendedHazardFileLoc), Scenario %in% c('RCP 2.6', 'RCP 4.5', 'RCP 8.5'))
		appendedOutputs$Relative_Hazard_Score = NA	;	appendedOutputs$Relative_Hazard_Score_Number = NA	;	appendedOutputs$Trend_Aggregated_For_Looker = NA	
		specializedOutput_SLR = fread(paste0(customerFolder, userName, '_', thisDate, '_seaLevelRise_1km', thisDate, '.csv'))
		specializedOutput_LRF = fread(paste0(customerFolder, userName, '_', thisDate, '_highRestFlood_', (locationFootprint*2+1), 'km', thisDate, '.csv'))
		specializedOutputs = merge(specializedOutput_LRF, specializedOutput_SLR, all = TRUE)
		dataOutput = merge(mainHazardOutputs,  specializedOutputs, all = TRUE)
		if(any(names(appendedOutputs) %in% c("BusinessUnit", "AssetType")))	{
			colnames(appendedOutputs)[which(names(appendedOutputs) == "BusinessUnit")] = "Business_Unit"
			colnames(appendedOutputs)[which(names(appendedOutputs) == "AssetType")] = "Asset_Type"
		}
		dataOutput = merge(dataOutput, appendedOutputs, all = TRUE)
	}

	allHazards = unique(dataOutput$Hazard)

	theScenarios = c('RCP 2.6', 'RCP 4.5', 'RCP 6.0', 'RCP 8.5')
	scenarioRename = c('1. Low Emissions', '2. Middle of the Road', '2. Middle of the Road', '3. High Emissions')

	for(thisDecade in unique(dataOutput$Decade))	{
		for(thisScen in c(1,3,4))	{
			for(thisLoc in unique(dataOutput$Location))	{
				locRow = which(customerTable$Location == thisLoc)
				avgOfAllHazards = NULL
				for(thisHazard in allHazards)	{
					newHazard = subset(dataOutput, Decade == thisDecade & Scenario == theScenarios[thisScen] & Location == thisLoc & Hazard == thisHazard)
		# catching exceptions for agg score calculation (where some are calculated using a subset of hazard measures)
						## !! temp fix, need to automate aggregate score exceptions  !!!!#####
					if(thisHazard %in% aggScoreExceptions)	{
						thisException = which(aggScoreExceptions == thisHazard)
						newHazard = subset(dataOutput, Decade == thisDecade & Scenario == theScenarios[thisScen] & Location == thisLoc & Hazard == thisHazard & Advanced_Data_Measures == aggScoreExceptionsValues[thisException])
					}
					
					# catching exceptions for when we do not have 4.5 or 8.5
					if(nrow(newHazard) == 0 | any(is.na(newHazard$Percentile_Score)))	{
						newHazard = subset(dataOutput, Decade == thisDecade & Scenario == theScenarios[thisScen - 1] & Location == thisLoc & Hazard == thisHazard)
						if(thisHazard %in% aggScoreExceptions)	{
							thisException = which(aggScoreExceptions == thisHazard)
							newHazard = subset(dataOutput, Decade == thisDecade & Scenario == theScenarios[thisScen - 1] & Location == thisLoc & Hazard == thisHazard & Advanced_Data_Measures == aggScoreExceptionsValues[thisException])
						}
					}
						## !! end temp fix, need to automate aggregate score exceptions  !!!!#####
						
					dataOutput = rbind(dataOutput,
						data.frame(
							User = userName,
							Location = customerTable$Location[locRow],
							Region = customerTable$Region[locRow],
							Subregion = customerTable$Subregion[locRow],
							Lat = customerTable$Lat[locRow],
							Lon = customerTable$Lon[locRow],
							Hazard = thisHazard,
							Hazard_Measure = 'Aggregate Score',
							Decade = thisDecade,
							Scenario = scenarioRename[thisScen],
							Raw_Hazard_Value = NA,
							Percentile_Score = mean(newHazard$Percentile_Score, na.rm=TRUE),
							Relative_Hazard_Score = NA,
							Decadal_Trend_Strength = NA,
							Decadal_Trend_Significance = NA,
							Long_Term_Trend_Strength = NA,
							Long_Term_Trend_Significance = NA,
							Relative_Hazard_Score_Number = NA,
							Trend_Aggregated_For_Looker = NA,
							Advanced_Data_Measures = NA,
							Advanced_Data_Measures_Units = NA,
							Raw_Hazard_Value_25th = NA,
							Raw_Hazard_Value_75th = NA,
							Asset_Type = customerTable$AssetType[locRow],
							Business_Unit = customerTable$BusinessUnit[locRow],
							Country = customerTable$Country[locRow],
							State = customerTable$State[locRow],
							City = customerTable$City[locRow])
					)

					avgOfAllHazards = c(avgOfAllHazards, mean(newHazard$Percentile_Score, na.rm=TRUE))
				}
				
				hazardWeights = rep(1, length(allHazards))
				for(thisHazWeight in 1:length(allHazards))	{
					newWeightLoc = which(hazardWeighting$Hazard == allHazards[thisHazWeight] & hazardWeighting$Asset_type == customerTable$AssetType[locRow])
					if(identical(newWeightLoc, integer(0)))	{
						print(c('asset not found', allHazards[thisHazWeight]))
						newWeight = 1
					}	else	{
						newWeight = hazardWeighting$Hazard_weight[newWeightLoc]
					}
					hazardWeights[thisHazWeight] = newWeight
				}

				dataOutput = rbind(dataOutput,
					data.frame(
						User = userName,
						Location = customerTable$Location[locRow],
						Region = customerTable$Region[locRow],
						Subregion = customerTable$Subregion[locRow],
						Lat = customerTable$Lat[locRow],
						Lon = customerTable$Lon[locRow],
						Hazard = "Aggregate Climate Score",
						Hazard_Measure = c("Aggregate Score", "Weighted Aggregate Score"),
						Decade = thisDecade,
						Scenario = scenarioRename[thisScen],
						Raw_Hazard_Value = NA,
						Percentile_Score = c(mean(avgOfAllHazards, na.rm=TRUE), weighted.mean(avgOfAllHazards, hazardWeights, na.rm=TRUE)),
						Relative_Hazard_Score = NA,
						Decadal_Trend_Strength = NA,
						Decadal_Trend_Significance = NA,
						Long_Term_Trend_Strength = NA,
						Long_Term_Trend_Significance = NA,
						Relative_Hazard_Score_Number = NA,
						Trend_Aggregated_For_Looker = NA,
						Advanced_Data_Measures = NA,
						Advanced_Data_Measures_Units = NA,
						Raw_Hazard_Value_25th = NA,
						Raw_Hazard_Value_75th = NA,
						Asset_Type = customerTable$AssetType[locRow],
						Business_Unit = customerTable$BusinessUnit[locRow],
						Country = customerTable$Country[locRow],
						State = customerTable$State[locRow],
						City = customerTable$City[locRow])
				)
					
					# ugh these *!@#$$$% ad hoc requests are driving me crazy..... ugh....
					# to be rewritten: here, we are calculating a reweighted score (1-150) for each hazard
				dataOutput = rbind(dataOutput,
					data.frame(
						User = userName,
						Location = customerTable$Location[locRow],
						Region = customerTable$Region[locRow],
						Subregion = customerTable$Subregion[locRow],
						Lat = customerTable$Lat[locRow],
						Lon = customerTable$Lon[locRow],
						Hazard = allHazards,
						Hazard_Measure = "Weighted Aggregate Score",
						Decade = thisDecade,
						Scenario = scenarioRename[thisScen],
						Raw_Hazard_Value = NA,
						Percentile_Score = avgOfAllHazards * hazardWeights,
						Relative_Hazard_Score = NA,
						Decadal_Trend_Strength = NA,
						Decadal_Trend_Significance = NA,
						Long_Term_Trend_Strength = NA,
						Long_Term_Trend_Significance = NA,
						Relative_Hazard_Score_Number = NA,
						Trend_Aggregated_For_Looker = NA,
						Advanced_Data_Measures = NA,
						Advanced_Data_Measures_Units = NA,
						Raw_Hazard_Value_25th = NA,
						Raw_Hazard_Value_75th = NA,
						Asset_Type = customerTable$AssetType[locRow],
						Business_Unit = customerTable$BusinessUnit[locRow],
						Country = customerTable$Country[locRow],
						State = customerTable$State[locRow],
						City = customerTable$City[locRow])
				)
	

			}
		}
	}


		# quick fix to rename scenarios
	scenariosRename = c('1. Low Emissions', '2. Middle of the Road', '2. Middle of the Road', '3. High Emissions')
	for(i in 1:length(scenariosRename))	{
		dataOutput$Scenario[dataOutput$Scenario == theScenarios[i]] = scenariosRename[i]
	}
		
		# quick fix to set NA tends to 0s for looker
	dataOutput$Trend_Aggregated_For_Looker[which(is.na(dataOutput$Trend_Aggregated_For_Looker))] = 0
		
		# identifying relative hazard scores
	for(thisRow in 1:nrow(relHazScores))	{
		dataOutput$Relative_Hazard_Score[which(dataOutput$Percentile_Score > relHazScores$Hazard_Percentile[thisRow])] = paste0(thisRow, '. ', relHazScores$Hazard_Common_Name[thisRow])
		dataOutput$Relative_Hazard_Score_Number[which(dataOutput$Percentile_Score > relHazScores$Hazard_Percentile[thisRow])] = thisRow
	}		
		
		# initializing each scenario as 2010s middle of the road
		# unfortunately, the data structure is all over the place, so only a for loop can work
	for(thisRow in 1:nrow(dataOutput))	{
		if(dataOutput[thisRow,]$Decade == 2010 & dataOutput[thisRow, ]$Scenario %in% c("Low Emissions", 'High Emissions'))	{
			if(!is.na(dataOutput[thisRow, ]$Advanced_Data_Measures))	{
				dataOutput[thisRow, c('Percentile_Score', 'Raw_Hazard_Value', 'Raw_Hazard_Value_25th', 'Raw_Hazard_Value_75th')]	 =
					subset(dataOutput, Location == dataOutput[thisRow, ]$Location & 
						Hazard == dataOutput[thisRow, ]$Hazard &
						Hazard_Measure == dataOutput[thisRow, ]$Hazard_Measure &
						Advanced_Data_Measures == dataOutput[thisRow, ]$Advanced_Data_Measures &
						Scenario == 'Middle of the Road' & 
						Decade == 2010)[ , c('Percentile_Score', 'Raw_Hazard_Value', 'Raw_Hazard_Value_25th', 'Raw_Hazard_Value_75th')]
			}	else	{
				dataOutput[thisRow, c('Percentile_Score', 'Raw_Hazard_Value', 'Raw_Hazard_Value_25th', 'Raw_Hazard_Value_75th')]	 =
				subset(dataOutput, Location == dataOutput[thisRow, ]$Location & 
					Hazard == dataOutput[thisRow, ]$Hazard &
					Hazard_Measure == dataOutput[thisRow, ]$Hazard_Measure &
					Scenario == 'Middle of the Road' & 
					Decade == 2010)[ , c('Percentile_Score', 'Raw_Hazard_Value', 'Raw_Hazard_Value_25th', 'Raw_Hazard_Value_75th')]
			}
		}
	}


		
		# aggregating trends for looker db
	dataOutput$Trend_Aggregated_For_Looker = 0
	signifiDecreaseRows = which(dataOutput$Decadal_Trend_Significance < 0.5 & dataOutput$Long_Term_Trend_Significance < 0.05 & dataOutput$Decadal_Trend_Strength < 0 & dataOutput$Long_Term_Trend_Strength < 0)
	signifiIncreaseRows = which(dataOutput$Decadal_Trend_Significance < 0.5 & dataOutput$Long_Term_Trend_Significance < 0.05 & dataOutput$Decadal_Trend_Strength > 0 & dataOutput$Long_Term_Trend_Strength > 0)
	dataOutput$Trend_Aggregated_For_Looker[signifiDecreaseRows] = dataOutput$Decadal_Trend_Strength[signifiDecreaseRows]
	dataOutput$Trend_Aggregated_For_Looker[signifiIncreaseRows] = dataOutput$Decadal_Trend_Strength[signifiIncreaseRows]

	fwrite(dataOutput, paste0(customerFolder, 'processedOutputForAllHazards_', userName, '_', thisDate, '.csv'))
	}





