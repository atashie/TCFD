########################################################
### library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)
	sf::sf_use_s2(FALSE) # for problem with intersecting spherical w flat
library(lubridate)
library(ncdf4)
library(magrittr)
library(maptools)
library(ncdf4)
library(easyNCDF)	# for ArrayToNc()
#library(robslopes)	# for TheilSen()
#library(trend)		# for sens.slope()
library(mblm)		# for sens slope mlbm()




#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\TCFD\\BurntArea\\"
ncVarFileName = 'burntarea'
ncOutputPath = 'J:\\Cai_data\\TCFD\\ProcessedNCs\\'
saveDate = '17FEB2023'
rcpScenarios = c(26, 60, 85)
whichDecades = seq(10,90,10)
valueType = 1:6

	# initializing start decade
initDates = 1:168
ncname_gfdl = paste0('clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_gfdl = nc_open(paste0(ncpath, ncname_gfdl))
nc_gfdl_init = ncvar_get(ncin_gfdl,ncVarFileName)[ , , initDates]	# lon, lat, time
	# identifying lat lons before closing data
nc_lat = ncvar_get(ncin_gfdl, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(ncin_gfdl, 'lon')
nc_close(ncin_gfdl)

ncname_hadgem = paste0('clm45_hadgem2-es_ewembi_rcp26_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_hadgem-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_hadgem = nc_open(paste0(ncpath, ncname_hadgem))
nc_hadgem_init = ncvar_get(ncin_hadgem,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_hadgem)

ncname_ipsl = paste0('clm45_ipsl-cm5a-lr_ewembi_rcp26_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_ipsl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_ipsl = nc_open(paste0(ncpath, ncname_ipsl))
nc_ipsl_init = ncvar_get(ncin_ipsl,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_ipsl)

ncname_miroc = paste0('clm45_miroc5_ewembi_rcp26_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_miroc-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_miroc = nc_open(paste0(ncpath, ncname_miroc))
nc_miroc_init = ncvar_get(ncin_miroc,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_miroc)

	# array for holding outputs
myMissingData = NA
dataOutArray = array(rep(myMissingData, length(nc_lon) * length(nc_lat) * length(whichDecades) * length(rcpScenarios) * length(valueType)), 
	dim = c(length(nc_lon), length(nc_lat), length(whichDecades), length(rcpScenarios), length(valueType)))

scalar = 100


for(thisScen in 1:length(rcpScenarios))	{
	rcpScenNum = rcpScenarios[thisScen]
	rcpScen = paste0('rcp', rcpScenNum)

	ncname_gfdl = paste0('clm45_gfdl-esm2m_ewembi_',rcpScen,'_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_gfdl = nc_open(paste0(ncpath, ncname_gfdl))
	nc_gfdl = ncvar_get(ncin_gfdl,ncVarFileName)	# lon, lat, time

	ncname_hadgem = paste0('clm45_hadgem2-es_ewembi_',rcpScen,'_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_hadgem = nc_open(paste0(ncpath, ncname_hadgem))
	nc_hadgem = ncvar_get(ncin_hadgem,ncVarFileName)	# lon, lat, time
	
	ncname_ipsl = paste0('clm45_ipsl-cm5a-lr_ewembi_',rcpScen,'_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_ipsl = nc_open(paste0(ncpath, ncname_ipsl))
	nc_ipsl = ncvar_get(ncin_ipsl,ncVarFileName)	# lon, lat, time
	
	ncname_miroc = paste0('clm45_miroc5_ewembi_',rcpScen,'_2005soc_co2_', ncVarFileName, '_global_monthly_2006_2099.nc4')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_miroc = nc_open(paste0(ncpath, ncname_miroc))
	nc_miroc = ncvar_get(ncin_miroc,ncVarFileName)	# lon, lat, time

	nc_date = as.Date("1661-01-01") + ncvar_get(ncin_miroc, 'time') * 30.4375# time is months after 1661-1-1
	nc_years = unique(year(nc_date))
	numYears = length(nc_years)
	missing_data = 1.00000002004088e+20


	for(i in 1:length(nc_lat))	{
		for(j in 1:length(nc_lon))	{
			nc_dummy = nc_gfdl[j,i, ] # reading in one data set to test for nona
			if(any(!is.na(nc_dummy) & any(nc_dummy != missing_data)))	{
				print(c(i, j))
				gfdl_all = c(nc_gfdl_init[j,i, initDates], nc_gfdl[j,i, -initDates]) * scalar
				hadgem_all = c(nc_hadgem_init[j,i, initDates], nc_hadgem[j,i, -initDates]) * scalar
				ipsl_all = c(nc_ipsl_init[j,i, initDates], nc_ipsl[j,i, -initDates]) * scalar
				miroc_all = c(nc_miroc_init[j,i, initDates], nc_miroc[j,i, -initDates]) * scalar

				gfdl_yrly = NULL
				hadgem_yrly = NULL
				ipsl_yrly = NULL
				miroc_yrly = NULL
				for(thisYear in nc_years)	{
					theseDates = which(year(nc_date) == thisYear)
					gfdl_yrly = c(gfdl_yrly, sum(gfdl_all[theseDates]))
					hadgem_yrly = c(hadgem_yrly, sum(hadgem_all[theseDates]))
					ipsl_yrly = c(ipsl_yrly, sum(ipsl_all[theseDates]))
					miroc_yrly = c(miroc_yrly, sum(miroc_all[theseDates]))
				}
				
				gfdl_smth = ksmooth(nc_years, gfdl_yrly, kernel = 'normal', bandwidth = 15, n.points = numYears)$y
				hadgem_smth = ksmooth(nc_years, hadgem_yrly, kernel = 'normal', bandwidth = 15, n.points = numYears)$y
				ipsl_smth = ksmooth(nc_years, ipsl_yrly, kernel = 'normal', bandwidth = 15, n.points = numYears)$y
				miroc_smth = ksmooth(nc_years, miroc_yrly, kernel = 'normal', bandwidth = 15, n.points = numYears)$y


				for(thisDecade in 1:9)	{
					theseYears = (thisDecade - 1) * 10 + 5:14
						# defining absolute values
					dataSmoothMed = median(c(gfdl_smth[theseYears], hadgem_smth[theseYears], ipsl_smth[theseYears], miroc_smth[theseYears]))
					dataOutArray[j, i, thisDecade, thisScen, 1] = dataSmoothMed
					dataQuantDiffs = diff(quantile(c(gfdl_yrly[theseYears], hadgem_yrly[theseYears], ipsl_yrly[theseYears], miroc_yrly[theseYears]), c(0.25, 0.5, 0.75)))
					dataOutArray[j, i, thisDecade, thisScen, 5] = dataSmoothMed - abs(dataQuantDiffs[1])
					dataOutArray[j, i, thisDecade, thisScen, 6] =  dataSmoothMed + abs(dataQuantDiffs[2])
						# calculating decadal trends (sens slope) and 	decadal significance (spearmans)	
					theseYears = 1:((thisDecade - 1) * 10 + 14)
					theseGfdl = gfdl_yrly[theseYears]
					theseHadgem = hadgem_yrly[theseYears]
					theseIpsl = ipsl_yrly[theseYears]
					theseMiroc =  miroc_yrly[theseYears]
					dataOutArray[j, i, thisDecade, thisScen, 3] = median(mblm(theseGfdl ~ theseYears)$coefficients[2],
						mblm(theseHadgem ~ theseYears)$coefficients[2],
						mblm(theseIpsl ~ theseYears)$coefficients[2],
						mblm(theseMiroc ~ theseYears)$coefficients[2])
					dataOutArray[j, i, thisDecade, thisScen, 4] = median(cor.test(theseYears, theseGfdl, method='spearman')$p.value,
						cor.test(theseYears, theseHadgem, method='spearman')$p.value,
						cor.test(theseYears, theseIpsl, method='spearman')$p.value,
						cor.test(theseYears, theseMiroc, method='spearman')$p.value, na.rm=TRUE)
				}	
					
					# calculating long-term trends (sens slope)
	#			dataOutArray[j, i, , thisScen, 5] = dataOutArray[j, i, 9, thisScen, 3]

					# calculating long-term significance (spearmans)
	#			dataOutArray[j, i, , thisScen, 6] = dataOutArray[j, i, 9, thisScen, 4]				
			}
		}
			}
		}
	saveRDS(dataOutArray, file=paste0(ncpath, 'data_out.rds'))
	}
	nc_close(ncin_gfdl)
	nc_close(ncin_hadgem)
	nc_close(ncin_ipsl)
	nc_close(ncin_miroc)
}

dataOutArray = readRDS(file=paste0(ncpath, 'data_out.rds'))


	# removing zeroes from non-impacted regions
maskedLocs26 = which(is.na(dataOutArray[ , , 1, 1, 1]))
histDatSubset26 =  dataOutArray[ , , 1, 1, 1][-maskedLocs26]
maskedLocs60 = which(is.na(dataOutArray[ , , 1, 2, 1]))
histDatSubset60 =  dataOutArray[ , , 1, 2, 1][-maskedLocs60]
maskedLocs85 = which(is.na(dataOutArray[ , , 1, 3, 1]))
histDatSubset85 =  dataOutArray[ , , 1, 3, 1][-maskedLocs85]
histQuants = quantile(c(histDatSubset26, histDatSubset60, histDatSubset85), seq(0.01, 1, length.out=100))
histQuants

minSignif = 0.1
	# removing low values from non-impacted regions
maskedLocs26_zeroes = which(is.na(dataOutArray[ , , 1, 1, 1]) | dataOutArray[ , , 1, 1, 1] < minSignif)
histDatSubset26_zeroes =  dataOutArray[ , , 1, 1, 1][-maskedLocs26_zeroes]
maskedLocs60_zeroes = which(is.na(dataOutArray[ , , 1, 2, 1]) | dataOutArray[ , , 1, 2, 1] < minSignif)
histDatSubset60_zeroes =  dataOutArray[ , , 1, 2, 1][-maskedLocs60_zeroes]
maskedLocs85_zeroes = which(is.na(dataOutArray[ , , 1, 3, 1]) | dataOutArray[ , , 1, 3, 1] < minSignif)
histDatSubset85_zeroes =  dataOutArray[ , , 1, 3, 1][-maskedLocs85_zeroes]
histQuants = quantile(c(histDatSubset26_zeroes, histDatSubset60_zeroes, histDatSubset85_zeroes), seq(0.01, 1, length.out=100))
histQuants



for(i in 1:length(whichDecades))	{
	dataOutArray[ , , i, 1, 2] = 1
	dataOutArray[ , , i, 2, 2] = 1
	dataOutArray[ , , i, 3, 2] = 1
	for(j in 1:(length(histQuants)))	{
		dataOutArray[ , , i, 1, 2][dataOutArray[ , , i, 1, 1] > histQuants[j]] = j + 0
		dataOutArray[ , , i, 2, 2][dataOutArray[ , , i, 2, 1] > histQuants[j]] = j + 0
		dataOutArray[ , , i, 3, 2][dataOutArray[ , , i, 3, 1] > histQuants[j]] = j + 0
	}
	dataOutArray[ , , i, 1, 2][maskedLocs26] = NA
	dataOutArray[ , , i, 2, 2][maskedLocs60] = NA
	dataOutArray[ , , i, 3, 2][maskedLocs85] = NA
}


tcfdVariable = dataOutArray
metadata = list(tcfdVariable = list(units = 'Fire - % Area Burned / yr'))
attr(tcfdVariable, 'variables') = metadata
names(dim(tcfdVariable)) = c('lon', 'lat', 'decade','rcpScen', 'valueClass')

lon = nc_lon
dim(lon) = length(lon)
metadata = list(lon = list(units = 'degrees'))
attr(lon, 'variables') = metadata
names(dim(lon)) = 'lon'

lat = nc_lat
dim(lat) = length(lat)
metadata = list(lat = list(units = 'degrees'))
attr(lat, 'variables') = metadata
names(dim(lat)) = 'lat'

decade = whichDecades
dim(decade) = length(decade)
metadata = list(decade = list(units = 'decades_of_21st_C'))
attr(decade, 'variables') = metadata
names(dim(decade)) = 'decade'

rcpScen = rcpScenarios
dim(rcpScen) = length(rcpScen)
metadata = list(rcpScen = list(units = 'RCP_scenario'))
attr(rcpScen, 'variables') = metadata
names(dim(rcpScen)) = 'rcpScen'

valueClass = 1:6#valueType
dim(valueClass) = length(valueClass)
metadata = list(valueClass = list(units = 'class'))
attr(valueClass, 'variables') = metadata
names(dim(valueClass)) = 'valueClass'

	# saving ncdf
ArrayToNc(list(tcfdVariable, lon, lat, decade, rcpScen, valueClass), file_path = paste0(ncOutputPath, ncVarFileName, 'v2_processed.nc'))

	# testing output, squinty eye test
myNC = nc_open(paste0(ncOutputPath, ncVarFileName, 'v2_processed.nc'))
nc_lat = ncvar_get(myNC, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(myNC, 'lon')
nc_testDat = ncvar_get(myNC, 'tcfdVariable')


image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,2])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,3])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,4])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,5])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,6])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,2,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,2] - nc_testDat[,,1,1,2])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,1,1] - nc_testDat[,,1,1,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,1] - nc_testDat[,,1,1,1])





















##############################################
# US only


dataOutArray = readRDS(file=paste0(ncpath, 'data_out.rds'))
us_lat = c(24.4, 49.6)
us_lon = c(-125.6, -66.3)

nc_us_lat = which(nc_lat > us_lat[1] & nc_lat < us_lat[2])
nc_us_lon = which(nc_lon > us_lon[1] & nc_lon < us_lon[2])

dataOutArray_us = dataOutArray[nc_us_lon, nc_us_lat, , ,]

	# removing zeroes from non-impacted regions
maskedLocs26 = which(is.na(dataOutArray_us[ , , 1, 1, 1]))
histDatSubset26 =  dataOutArray_us[ , , 1, 1, 1][-maskedLocs26]
maskedLocs60 = which(is.na(dataOutArray_us[ , , 1, 2, 1]))
histDatSubset60 =  dataOutArray_us[ , , 1, 2, 1][-maskedLocs60]
maskedLocs85 = which(is.na(dataOutArray_us[ , , 1, 3, 1]))
histDatSubset85 =  dataOutArray_us[ , , 1, 3, 1][-maskedLocs85]
histQuants = quantile(c(histDatSubset26, histDatSubset60, histDatSubset85), seq(0.01, 1, length.out=100))
histQuants

minSignif = 0.01
	# removing low values from non-impacted regions
maskedLocs26_zeroes = which(is.na(dataOutArray_us[ , , 1, 1, 1]) | dataOutArray_us[ , , 1, 1, 1] < minSignif)
histDatSubset26_zeroes =  dataOutArray_us[ , , 1, 1, 1][-maskedLocs26_zeroes]
maskedLocs60_zeroes = which(is.na(dataOutArray_us[ , , 1, 2, 1]) | dataOutArray_us[ , , 1, 2, 1] < minSignif)
histDatSubset60_zeroes =  dataOutArray_us[ , , 1, 2, 1][-maskedLocs60_zeroes]
maskedLocs85_zeroes = which(is.na(dataOutArray_us[ , , 1, 3, 1]) | dataOutArray_us[ , , 1, 3, 1] < minSignif)
histDatSubset85_zeroes =  dataOutArray_us[ , , 1, 3, 1][-maskedLocs85_zeroes]
histQuants = quantile(c(histDatSubset26_zeroes, histDatSubset60_zeroes, histDatSubset85_zeroes), seq(0.01, 1, length.out=100))
histQuants



for(i in 1:length(whichDecades))	{
	dataOutArray_us[ , , i, 1, 2] = 1
	dataOutArray_us[ , , i, 2, 2] = 1
	dataOutArray_us[ , , i, 3, 2] = 1
	for(j in 1:(length(histQuants)))	{
		dataOutArray_us[ , , i, 1, 2][dataOutArray_us[ , , i, 1, 1] > histQuants[j]] = j + 0
		dataOutArray_us[ , , i, 2, 2][dataOutArray_us[ , , i, 2, 1] > histQuants[j]] = j + 0
		dataOutArray_us[ , , i, 3, 2][dataOutArray_us[ , , i, 3, 1] > histQuants[j]] = j + 0
	}
	dataOutArray_us[ , , i, 1, 2][maskedLocs26] = NA
	dataOutArray_us[ , , i, 2, 2][maskedLocs60] = NA
	dataOutArray_us[ , , i, 3, 2][maskedLocs85] = NA
}


tcfdVariable = dataOutArray_us
metadata = list(tcfdVariable = list(units = 'Fire - % Area Burned / yr'))
attr(tcfdVariable, 'variables') = metadata
names(dim(tcfdVariable)) = c('lon', 'lat', 'decade','rcpScen', 'valueClass')

lon = nc_lon[nc_us_lon]
dim(lon) = length(lon)
metadata = list(lon = list(units = 'degrees'))
attr(lon, 'variables') = metadata
names(dim(lon)) = 'lon'

lat = nc_lat[nc_us_lat]
dim(lat) = length(lat)
metadata = list(lat = list(units = 'degrees'))
attr(lat, 'variables') = metadata
names(dim(lat)) = 'lat'

decade = whichDecades
dim(decade) = length(decade)
metadata = list(decade = list(units = 'decades_of_21st_C'))
attr(decade, 'variables') = metadata
names(dim(decade)) = 'decade'

rcpScen = rcpScenarios
dim(rcpScen) = length(rcpScen)
metadata = list(rcpScen = list(units = 'RCP_scenario'))
attr(rcpScen, 'variables') = metadata
names(dim(rcpScen)) = 'rcpScen'

valueClass = 1:6#valueType
dim(valueClass) = length(valueClass)
metadata = list(valueClass = list(units = 'class'))
attr(valueClass, 'variables') = metadata
names(dim(valueClass)) = 'valueClass'

	# saving ncdf
ArrayToNc(list(tcfdVariable, lon, lat, decade, rcpScen, valueClass), file_path = paste0(ncOutputPath, ncVarFileName, 'USonlyv2_processed.nc'))

	# testing output, squinty eye test
myNC = nc_open(paste0(ncOutputPath, ncVarFileName, 'USonlyv2_processed.nc'))
nc_lat = ncvar_get(myNC, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(myNC, 'lon')
nc_testDat = ncvar_get(myNC, 'tcfdVariable')


image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,2])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,3])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,4])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,5])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,6])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,2,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,2] - nc_testDat[,,1,1,2])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,1,1] - nc_testDat[,,1,1,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,1] - nc_testDat[,,1,1,1])
