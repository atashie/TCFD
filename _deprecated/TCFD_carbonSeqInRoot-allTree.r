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
ncpath = "C:\\ClimateData\\Carbon\\Sequestration_Root\\dcdcldbdltr\\future\\rawNcs\\"
ncVarFileName = 'croot-dcdcldbdltr'
ncOutputPath = 'C:\\ClimateData\\Carbon\\Sequestration_Root\\dcdcldbdltr\\future\\reanalysis\\'
saveDate = '27MAY2025'
rcpScenarios = c(126, 370, 585)
whichDecades = seq(10,90,10)
valueType = 1:6

	# reading in dummy data for lat lons
ncname_dummy = paste0('classic_gfdl-esm4_w5e5_ssp126_2015soc_default_croot-dcdcldbdltr_global_annual_2015_2100.nc')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_dummy = nc_open(paste0(ncpath, ncname_dummy))
nc_lat = ncvar_get(ncin_dummy, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(ncin_dummy, 'lon')

	# array for holding outputs
myMissingData = NA
dataOutArray = array(rep(NA, length(nc_lon) * length(nc_lat) * length(whichDecades) * length(rcpScenarios) * length(valueType)), 
	dim = c(length(nc_lon), length(nc_lat), length(whichDecades), length(rcpScenarios), length(valueType)))






for(thisScen in 1:length(rcpScenarios))	{
	rcpScenNum = rcpScenarios[thisScen]
	rcpScen = paste0('ssp', rcpScenNum)

	ncname = paste0('classic_gfdl-esm4_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_annual_2015_2100.nc')
	ncin = nc_open(paste0(ncpath, ncname))
	nc_l = ncvar_get(ncin,ncVarFileName)	# lon, lat, time

	ncname = paste0('classic_gfdl-esm4_w5e5_', rcpScen, '_2015soc-from-histsoc_default_', ncVarFileName, '_global_annual_2015_2100.nc')
	ncin = nc_open(paste0(ncpath, ncname))
	nc_2 = ncvar_get(ncin,ncVarFileName)	# lon, lat, time

	ncname = paste0('classic_ukesm1-0-ll_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_annual_2015_2100.nc')
	ncin = nc_open(paste0(ncpath, ncname))
	nc_3 = ncvar_get(ncin,ncVarFileName)	# lon, lat, time

	ncname = paste0('classic_ukesm1-0-ll_w5e5_', rcpScen, '_2015soc-from-histsoc_default_', ncVarFileName, '_global_annual_2015_2100.nc')
	ncin = nc_open(paste0(ncpath, ncname))
	nc_4 = ncvar_get(ncin,ncVarFileName)	# lon, lat, time



	nc_date = as.Date("1601-01-01") + ncvar_get(ncin, 'time')# time is years after 1601-1-1
	nc_years = unique(year(nc_date))
	missing_data = 1.00000002004088e+20

	dates1019 = which(year(nc_date) == 2015):which(year(nc_date) == 2019)
	dates2029 = which(year(nc_date) == 2020):which(year(nc_date) == 2029)
	dates3039 = which(year(nc_date) == 2030):which(year(nc_date) == 2039)
	dates4049 = which(year(nc_date) == 2040):which(year(nc_date) == 2049)
	dates5059 = which(year(nc_date) == 2050):which(year(nc_date) == 2059)
	dates6069 = which(year(nc_date) == 2060):which(year(nc_date) == 2069)
	dates7079 = which(year(nc_date) == 2070):which(year(nc_date) == 2079)
	dates8089 = which(year(nc_date) == 2080):which(year(nc_date) == 2089)
	dates9099 = which(year(nc_date) == 2090):which(year(nc_date) == 2099)


	for(i in 1:length(nc_lat))	{
		for(j in 1:length(nc_lon))	{
			nc_dummy = nc_l[j,i,dates1019] # reading in one data set to test for nona
			if(any(!is.na(nc_dummy) & any(nc_dummy != missing_data)))	{
				print(c(i, j))
				data1019 = c(nc_l[j,i,dates1019],nc_2[j,i,dates1019],nc_3[j,i,dates1019],nc_4[j,i,dates1019]) 
				data2029 = c(nc_l[j,i,dates2029],nc_2[j,i,dates2029],nc_3[j,i,dates2029],nc_4[j,i,dates2029]) 
				data3039 = c(nc_l[j,i,dates3039],nc_2[j,i,dates3039],nc_3[j,i,dates3039],nc_4[j,i,dates3039]) 
				data4049 = c(nc_l[j,i,dates4049],nc_2[j,i,dates4049],nc_3[j,i,dates4049],nc_4[j,i,dates4049]) 
				data5059 = c(nc_l[j,i,dates5059],nc_2[j,i,dates5059],nc_3[j,i,dates5059],nc_4[j,i,dates5059]) 
				data6069 = c(nc_l[j,i,dates6069],nc_2[j,i,dates6069],nc_3[j,i,dates6069],nc_4[j,i,dates6069]) 			
				data7079 = c(nc_l[j,i,dates7079],nc_2[j,i,dates7079],nc_3[j,i,dates7079],nc_4[j,i,dates7079]) 			
				data8089 = c(nc_l[j,i,dates8089],nc_2[j,i,dates8089],nc_3[j,i,dates8089],nc_4[j,i,dates8089]) 		
				data9099 = c(nc_l[j,i,dates9099],nc_2[j,i,dates9099],nc_3[j,i,dates9099],nc_4[j,i,dates9099]) 			


				datesSeq1019 = rep(seq(first(dates1019), last(dates1019)), 4)
				datesSeq2029 = rep(seq(first(dates2029), last(dates2029)), 4)
				datesSeq3039 = rep(seq(first(dates3039), last(dates3039)), 4)
				datesSeq4049 = rep(seq(first(dates4049), last(dates4049)), 4)
				datesSeq5059 = rep(seq(first(dates5059), last(dates5059)), 4)
				datesSeq6069 = rep(seq(first(dates6069), last(dates6069)), 4)
				datesSeq7079 = rep(seq(first(dates7079), last(dates7079)), 4)
				datesSeq8089 = rep(seq(first(dates8089), last(dates8089)), 4)
				datesSeq9099 = rep(seq(first(dates9099), last(dates9099)), 4)


					# defining absolute values
				dataOutArray[j, i, 1, thisScen, 1] = mean(data1019)# * 100
				dataOutArray[j, i, 2, thisScen, 1] = mean(data2029)# * 100
				dataOutArray[j, i, 3, thisScen, 1] = mean(data3039)# * 100
				dataOutArray[j, i, 4, thisScen, 1] = mean(data4049)# * 100
				dataOutArray[j, i, 5, thisScen, 1] = mean(data5059)# * 100
				dataOutArray[j, i, 6, thisScen, 1] = mean(data6069)# * 100
				dataOutArray[j, i, 7, thisScen, 1] = mean(data7079)# * 100
				dataOutArray[j, i, 8, thisScen, 1] = mean(data8089)# * 100
				dataOutArray[j, i, 9, thisScen, 1] = mean(data9099)# * 100

					# calculating decadal trends (sens slope) and 	decadal significance (spearmans)	
				theDates = datesSeq1019
				dataOutArray[j, i, 1, thisScen, 3] = lm(data1019 ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 1, thisScen, 4] = cor.test(theDates, data1019, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq2029)
				theValues =  c(data1019, data2029)
				dataOutArray[j, i, 2, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 2, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq3039)
				theValues =  c(data1019, data2029, data3039)
				dataOutArray[j, i, 3, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 3, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq4049)
				theValues =  c(data1019, data2029, data3039, data4049)
				dataOutArray[j, i, 4, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 4, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq5059)
				theValues =  c(data1019, data2029, data3039, data4049, data5059)
				dataOutArray[j, i, 5, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 5, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq6069)
				theValues =  c(data1019, data2029, data3039, data4049, data5059, data6069)
				dataOutArray[j, i, 6, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 6, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq7079)
				theValues =  c(data1019, data2029, data3039, data4049, data5059, data6069, data7079)
				dataOutArray[j, i, 7, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 7, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq8089)
				theValues =  c(data1019, data2029, data3039, data4049, data5059, data6069, data7079, data8089)
				dataOutArray[j, i, 8, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 8, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
				
				theDates = c(theDates, datesSeq9099)
				theValues =  c(data1019, data2029, data3039, data4049, data5059, data6069, data7079, data8089, data9099)
				dataOutArray[j, i, 9, thisScen, 3] = lm(theValues ~ theDates)$coefficients[2] * 10
				dataOutArray[j, i, 9, thisScen, 4] = cor.test(theDates, theValues, method='spearman')$p.value
						
					# calculating long-term trends (sens slope)
				dataOutArray[j, i, , thisScen, 5] = dataOutArray[j, i, 9, thisScen, 3]

					# calculating long-term significance (spearmans)
				dataOutArray[j, i, , thisScen, 6] = dataOutArray[j, i, 9, thisScen, 4]				
			}
		}
	saveRDS(dataOutArray, file=paste0(ncpath, 'data_out.rds'))
	}
	nc_close(ncin)
}

dataOutArray = readRDS(file=paste0(ncpath, 'data_out.rds'))

	# defining quantiles 
maskedLocs26 = which(is.na(dataOutArray[ , , 1, 1, 1]))
histDatSubset26 =  dataOutArray[ , , 1, 1, 1][-maskedLocs26]
maskedLocs60 = which(is.na(dataOutArray[ , , 1, 2, 1]))
histDatSubset60 =  dataOutArray[ , , 1, 2, 1][-maskedLocs60]
maskedLocs85 = which(is.na(dataOutArray[ , , 1, 3, 1]))
histDatSubset85 =  dataOutArray[ , , 1, 3, 1][-maskedLocs85]
#histQuants = quantile(c(histDatSubset26, histDatSubset60), seq(0.01, 1, 0.01))

	# removing zeroes from non-impacted regions
maskedLocs26_zeroes = which(is.na(dataOutArray[ , , 1, 1, 1]) | dataOutArray[ , , 1, 1, 1] == 0)
histDatSubset26_zeroes =  dataOutArray[ , , 1, 1, 1][-maskedLocs26_zeroes]
maskedLocs60_zeroes = which(is.na(dataOutArray[ , , 1, 2, 1]) | dataOutArray[ , , 1, 2, 1] == 0)
histDatSubset60_zeroes =  dataOutArray[ , , 1, 2, 1][-maskedLocs60_zeroes]
maskedLocs85_zeroes = which(is.na(dataOutArray[ , , 1, 3, 1]) | dataOutArray[ , , 1, 3, 1] == 0)
histDatSubset85_zeroes =  dataOutArray[ , , 1, 3, 1][-maskedLocs85_zeroes]
histQuants = rev(quantile(c(histDatSubset26_zeroes, histDatSubset60_zeroes), seq(0.01, 1, length.out=100)))
#oldHistQuants

for(i in 1:length(whichDecades))	{
	dataOutArray[ , , i, 1, 2] = 1
	dataOutArray[ , , i, 2, 2] = 1
	dataOutArray[ , , i, 3, 2] = 1
	for(j in 1:(length(histQuants)))	{
		dataOutArray[ , , i, 1, 2][dataOutArray[ , , i, 1, 1] < histQuants[j]] = j# + 20
		dataOutArray[ , , i, 2, 2][dataOutArray[ , , i, 2, 1] < histQuants[j]] = j# + 20
		dataOutArray[ , , i, 3, 2][dataOutArray[ , , i, 3, 1] < histQuants[j]] = j# + 20
	}
	dataOutArray[ , , i, 1, 2][maskedLocs26] = NA
	dataOutArray[ , , i, 2, 2][maskedLocs60] = NA
	dataOutArray[ , , i, 3, 2][maskedLocs85] = NA
}


tcfdVariable = dataOutArray
metadata = list(tcfdVariable = list(units = 'Carbon Mass in Roots kg/m2'))
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
ArrayToNc(list(tcfdVariable, lon, lat, decade, rcpScen, valueClass), file_path = paste0(ncOutputPath, ncVarFileName, '_processed.nc'))

	# testing output, squinty eye test
myNC = nc_open(paste0(ncOutputPath, ncVarFileName, '_processed.nc'))
nc_lat = ncvar_get(myNC, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(myNC, 'lon')
nc_testDat = ncvar_get(myNC, 'tcfdVariable')


image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,2,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,2])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,3])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,4])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,5])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,6])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,2] - nc_testDat[,,1,1,2])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,1,1] - nc_testDat[,,1,1,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,1] - nc_testDat[,,1,1,1])
