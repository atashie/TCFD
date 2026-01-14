########################################################
library(data.table)
library(lubridate)
library(ncdf4)
library(easyNCDF)	# for ArrayToNc()
#library(robslopes)	# for TheilSen()
#library(trend)		# for sens.slope()
library(mblm)		# for sens slope mlbm()



#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\TCFD\\Evapotranspiration\\"
ncOutputPath = 'J:\\Cai_data\\TCFD\\ProcessedNCs\\'
ncVarFileName = 'potevap'
saveDate = '07FEB2022'
rcpScenarios = c(126, 370, 585)
whichDecades = seq(10,90,10)
valueType = 1:6


	# initializing start decade
initDates = 1:168
ncname_gfdl = paste0('h08_gfdl-esm4_w5e5_ssp126_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_gfdl = nc_open(paste0(ncpath, ncname_gfdl))
nc_gfdl_init = ncvar_get(ncin_gfdl,ncVarFileName)[ , , initDates]	# lon, lat, time
	# identifying lat lons before closing data
nc_lat = ncvar_get(ncin_gfdl, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(ncin_gfdl, 'lon')
nc_close(ncin_gfdl)

ncname_mpi = paste0('h08_mpi-esm1-2-hr_w5e5_ssp126_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_hadgem-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_mpi = nc_open(paste0(ncpath, ncname_mpi))
nc_mpi_init = ncvar_get(ncin_mpi,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_mpi)

ncname_ipsl = paste0('h08_ipsl-cm6a-lr_w5e5_ssp126_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_ipsl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_ipsl = nc_open(paste0(ncpath, ncname_ipsl))
nc_ipsl_init = ncvar_get(ncin_ipsl,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_ipsl)

ncname_mri = paste0('h08_mri-esm2-0_w5e5_ssp126_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_miroc-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_mri = nc_open(paste0(ncpath, ncname_mri))
nc_mri_init = ncvar_get(ncin_mri,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_mri)

ncname_ukesm = paste0('h08_ukesm1-0-ll_w5e5_ssp126_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_miroc-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_ukesm = nc_open(paste0(ncpath, ncname_ukesm))
nc_ukesm_init = ncvar_get(ncin_ukesm,ncVarFileName)[ , , initDates]	# lon, lat, time
nc_close(ncin_ukesm)


	# array for holding outputs
myMissingData = NA
dataOutArray = array(rep(myMissingData, length(nc_lon) * length(nc_lat) * length(whichDecades) * length(rcpScenarios) * length(valueType)), 
	dim = c(length(nc_lon), length(nc_lat), length(whichDecades), length(rcpScenarios), length(valueType)))

scalar = 60*60*24*30.4375  #kg / m2 / s to mm / month

for(thisScen in 1:length(rcpScenarios))	{
	rcpScenNum = rcpScenarios[thisScen]
	rcpScen = paste0('ssp', rcpScenNum)

	ncname_gfdl = paste0('h08_gfdl-esm4_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_gfdl = nc_open(paste0(ncpath, ncname_gfdl))

	ncname_mpi = paste0('h08_mpi-esm1-2-hr_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_hadgem-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_mpi = nc_open(paste0(ncpath, ncname_mpi))

	ncname_ipsl = paste0('h08_ipsl-cm6a-lr_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_ipsl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_ipsl = nc_open(paste0(ncpath, ncname_ipsl))

	ncname_mri = paste0('h08_mri-esm2-0_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_miroc-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_mri = nc_open(paste0(ncpath, ncname_mri))

	ncname_ukesm = paste0('h08_ukesm1-0-ll_w5e5_', rcpScen, '_2015soc_default_', ncVarFileName, '_global_monthly_2015_2100.nc')#"clm45_miroc-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
	ncin_ukesm = nc_open(paste0(ncpath, ncname_ukesm))

	nc_date = as.Date("1600-01-15") + ncvar_get(ncin_ukesm, 'time')# time is months after 1601-1-1
	nc_years = unique(year(nc_date))
	numYears = length(nc_years)
	missing_data = 1.00000002004088e+20

	for(i in 1:length(nc_lat))	{
		nc_gfdl = ncvar_get(ncin_gfdl,ncVarFileName)[ ,i , ]	# lon, lat, time
		nc_mpi = ncvar_get(ncin_mpi,ncVarFileName)[ , i, ]	# lon, lat, time
		nc_ipsl = ncvar_get(ncin_ipsl,ncVarFileName)[ , i, ]	# lon, lat, time
		nc_mri = ncvar_get(ncin_mri,ncVarFileName)[ , i, ]	# lon, lat, time
		nc_ukesm = ncvar_get(ncin_ukesm,ncVarFileName)[ , i, ]	# lon, lat, time

		for(j in 1:length(nc_lon))	{
			nc_dummy = nc_gfdl[j, ] # reading in one data set to test for nona
			if(any(!is.na(nc_dummy) & any(nc_dummy != missing_data)))	{
				print(c(i, j))
				gfdl_all = c(nc_gfdl_init[j, i, initDates], nc_gfdl[j, -initDates]) * scalar
				mpi_all = c(nc_mpi_init[j, i, initDates], nc_mpi[j, -initDates]) * scalar
				ipsl_all = c(nc_ipsl_init[j, i, initDates], nc_ipsl[j, -initDates]) * scalar
				mri_all = c(nc_mri_init[j, i, initDates], nc_mri[j, -initDates]) * scalar
				ukesm_all = c(nc_ukesm_init[j, i, initDates], nc_ukesm[j, -initDates]) * scalar

				gfdl_yrly = NULL
				mpi_yrly = NULL
				ipsl_yrly = NULL
				mri_yrly = NULL
				ukesm_yrly = NULL
				for(thisYear in nc_years)	{
					theseDates = which(year(nc_date) == thisYear)
					gfdl_yrly = c(gfdl_yrly, sum(gfdl_all[theseDates]))
					mpi_yrly = c(mpi_yrly, sum(mpi_all[theseDates]))
					ipsl_yrly = c(ipsl_yrly, sum(ipsl_all[theseDates]))
					mri_yrly = c(mri_yrly, sum(mri_all[theseDates]))
					ukesm_yrly = c(ukesm_yrly, sum(ukesm_all[theseDates]))
				}
				
				gfdl_smth = ksmooth(nc_years, gfdl_yrly, kernel = 'normal', bandwidth = 10, n.points = numYears)$y
				mpi_smth = ksmooth(nc_years, mpi_yrly, kernel = 'normal', bandwidth = 10, n.points = numYears)$y
				ipsl_smth = ksmooth(nc_years, ipsl_yrly, kernel = 'normal', bandwidth = 10, n.points = numYears)$y
				mri_smth = ksmooth(nc_years, mri_yrly, kernel = 'normal', bandwidth = 10, n.points = numYears)$y
				ukesm_smth = ksmooth(nc_years, ukesm_yrly, kernel = 'normal', bandwidth = 10, n.points = numYears)$y

					# calculating the 10s as 2014 to 2019
				theseYears = 1:6
				dataSmoothMed = median(c(gfdl_smth[theseYears], mpi_smth[theseYears], ipsl_smth[theseYears], mri_smth[theseYears], ukesm_smth[theseYears]))
				dataOutArray[j, i, 1, thisScen, 1] = dataSmoothMed
				dataQuantDiffs = diff(quantile(c(gfdl_yrly[theseYears], mpi_yrly[theseYears], ipsl_yrly[theseYears], mri_yrly[theseYears], ukesm_yrly[theseYears]), c(0.25, 0.5, 0.75)))
				dataOutArray[j, i, 1, thisScen, 5] = dataSmoothMed - abs(dataQuantDiffs[1])
				dataOutArray[j, i, 1, thisScen, 6] =  dataSmoothMed + abs(dataQuantDiffs[2])
					# calculating decadal trends (sens slope) and 	decadal significance (spearmans)	
				theseGfdl =   gfdl_yrly[theseYears]
				theseMpi =    mpi_yrly[theseYears]
				theseIpsl =   ipsl_yrly[theseYears]
				theseMri =    mri_yrly[theseYears]
				theseUkesm =  ukesm_yrly[theseYears]
				dataOutArray[j, i, 1, thisScen, 3] = median(mblm(theseGfdl ~ theseYears)$coefficients[2],
					mblm(theseMpi ~ theseYears)$coefficients[2],
					mblm(theseIpsl ~ theseYears)$coefficients[2],
					mblm(theseMri ~ theseYears)$coefficients[2],
					mblm(theseUkesm ~ theseYears)$coefficients[2], na.rm=TRUE)
				dataOutArray[j, i, 1, thisScen, 4] = median(cor.test(theseYears, theseGfdl, method='spearman')$p.value,
					cor.test(theseYears, theseMpi, method='spearman')$p.value,
					cor.test(theseYears, theseIpsl, method='spearman')$p.value,
					cor.test(theseYears, theseMri, method='spearman')$p.value,
					cor.test(theseYears, theseUkesm, method='spearman')$p.value, na.rm=TRUE)

				for(thisDecade in 2:9)	{
					theseYears = (thisDecade - 2) * 10 + 7:16
						# defining absolute values
					dataSmoothMed = median(c(gfdl_smth[theseYears], mpi_smth[theseYears], ipsl_smth[theseYears], mri_smth[theseYears], ukesm_smth[theseYears]))
					dataOutArray[j, i, thisDecade, thisScen, 1] = dataSmoothMed
					dataQuantDiffs = diff(quantile(c(gfdl_yrly[theseYears], mpi_yrly[theseYears], ipsl_yrly[theseYears], mri_yrly[theseYears], ukesm_yrly[theseYears]), c(0.25, 0.5, 0.75)))
					dataOutArray[j, i, thisDecade, thisScen, 5] = dataSmoothMed - abs(dataQuantDiffs[1])
					dataOutArray[j, i, thisDecade, thisScen, 6] =  dataSmoothMed + abs(dataQuantDiffs[2])
						# calculating decadal trends (sens slope) and 	decadal significance (spearmans)	
					theseYears = 1:((thisDecade - 2) * 10 + 16)
					theseGfdl =   gfdl_yrly[theseYears]
					theseMpi =    mpi_yrly[theseYears]
					theseIpsl =   ipsl_yrly[theseYears]
					theseMri =    mri_yrly[theseYears]
					theseUkesm =  ukesm_yrly[theseYears]
					dataOutArray[j, i, thisDecade, thisScen, 3] = median(mblm(theseGfdl ~ theseYears)$coefficients[2],
						mblm(theseMpi ~ theseYears)$coefficients[2],
						mblm(theseIpsl ~ theseYears)$coefficients[2],
						mblm(theseMri ~ theseYears)$coefficients[2],
						mblm(theseUkesm ~ theseYears)$coefficients[2], na.rm=TRUE)
					dataOutArray[j, i, thisDecade, thisScen, 4] = median(cor.test(theseYears, theseGfdl, method='spearman')$p.value,
						cor.test(theseYears, theseMpi, method='spearman')$p.value,
						cor.test(theseYears, theseIpsl, method='spearman')$p.value,
						cor.test(theseYears, theseMri, method='spearman')$p.value,
						cor.test(theseYears, theseUkesm, method='spearman')$p.value, na.rm=TRUE)
				}	
					
					# calculating long-term trends (sens slope)
	#			dataOutArray[j, i, , thisScen, 5] = dataOutArray[j, i, 9, thisScen, 3]

					# calculating long-term significance (spearmans)
	#			dataOutArray[j, i, , thisScen, 6] = dataOutArray[j, i, 9, thisScen, 4]				
			}
		}
		saveRDS(dataOutArray, file=paste0(ncpath, 'data_out.rds'))
	}
	nc_close(ncin_gfdl)
	nc_close(ncin_mpi)
	nc_close(ncin_ipsl)
	nc_close(ncin_mri)
	nc_close(ncin_ukesm)
	saveRDS(dataOutArray, file=paste0(ncpath, thisScen, 'data_out.rds'))
}

dataOutArray = readRDS(file=paste0(ncpath, 'data_out.rds'))

	# defining quantiles 
maskedLocs126 = which(is.na(dataOutArray[ , , 1, 1, 1]))
histDatSubset126 =  dataOutArray[ , , 1, 1, 1][-maskedLocs126]
maskedLocs370 = which(is.na(dataOutArray[ , , 1, 2, 1]))
histDatSubset370 =  dataOutArray[ , , 1, 2, 1][-maskedLocs370]
maskedLocs585 = which(is.na(dataOutArray[ , , 1, 3, 1]))
histDatSubset585 =  dataOutArray[ , , 1, 3, 1][-maskedLocs585]
histQuants = quantile(c(histDatSubset126, histDatSubset370, histDatSubset585), seq(0.01, 1, 0.01))
histQuants

for(i in 1:length(whichDecades))	{
	dataOutArray[ , , i, 1, 2] = 1
	for(j in 1:(length(histQuants)))	{
		dataOutArray[ , , i, 1, 2][dataOutArray[ , , i, 1, 1] >= histQuants[j]] = j
		dataOutArray[ , , i, 2, 2][dataOutArray[ , , i, 2, 1] >= histQuants[j]] = j
		dataOutArray[ , , i, 3, 2][dataOutArray[ , , i, 3, 1] >= histQuants[j]] = j
	}
	dataOutArray[ , , i, 1, 2][maskedLocs126] = NA
	dataOutArray[ , , i, 2, 2][maskedLocs370] = NA
	dataOutArray[ , , i, 3, 2][maskedLocs585] = NA
}



tcfdVariable = dataOutArray
metadata = list(tcfdVariable = list(units = 'Potential Evapotranspiration - mm / yr'))
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
metadata = list(rcpScen = list(units = 'SSP_scenario'))
attr(rcpScen, 'variables') = metadata
names(dim(rcpScen)) = 'rcpScen'

valueClass = 1:6#valueType
dim(valueClass) = length(valueClass)
metadata = list(valueClass = list(units = 'class'))
attr(valueClass, 'variables') = metadata
names(dim(valueClass)) = 'valueClass'

	# saving ncdf
ArrayToNc(list(tcfdVariable, lon, lat, decade, rcpScen, valueClass), file_path = paste0(ncOutputPath, ncVarFileName, '_ssp_v2_processed.nc'))

	# testing output, squinty eye test
myNC = nc_open(paste0(ncOutputPath, ncVarFileName, '_ssp_v2_processed.nc'))
nc_lat = ncvar_get(myNC, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(myNC, 'lon')
nc_testDat = ncvar_get(myNC, 'tcfdVariable')


image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,2,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,3,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,2])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,3])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,4])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,5])
image(nc_lon, rev(nc_lat), nc_testDat[,,1,1,6])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,1] - nc_testDat[,,1,3,1])
image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,1] - nc_testDat[,,1,2,1])

image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,2] - nc_testDat[,,1,3,2])
image(nc_lon, rev(nc_lat), nc_testDat[,,9,3,2] - nc_testDat[,,1,2,2])







# fancy plots
library(sf)
library(ggplot2)
library(viridis)
#geomData = st_as_sf(data.frame(lat = rep(nc_lat, each=length(nc_lon)), lon = rep(nc_lon, length(nc_lat)), plotData = as.vector(nc_testDat[ , , 1, 2, 1])), coords = c('lon', 'lat'))
geomData = data.frame(lat = rep(nc_lat, each=length(nc_lon)), lon = rep(nc_lon, length(nc_lat)), plotData = as.vector(nc_testDat[ , , 1, 2, 1]))
df_sub = subset(geomData, lat > -50)
latlonBox = c(25, 65, -10, 60) #(minlat, maxlat, minlon, maxlon)
myTitle = 'Sugar Beet Yield'
mySubtitle = 'nonirrigated'
ggplot(data = df_sub, aes(fill = plotData, y = lat, x = lon)) +
	geom_tile(width = 1, height = 1)	+
#	geom_sf(data = df_sub, aes(fill = plotData, color = plotData)) +
	borders('world', xlim=range(df_sub$lon), ylim=range(df_sub$lat), 
            colour='gray90', size=.2)	+
	theme_light()	+
	theme(panel.ontop=FALSE, panel.background=element_blank()) +
	scale_fill_viridis()	+
#	scale_colour_distiller(palette='Spectral', 
#                         limits=c(quantile(df_sub, .7, na.rm=T), 
#                                  quantile(df_sub, .999, na.rm=T))) +
    coord_quickmap(xlim=c(latlonBox[3], latlonBox[4]), ylim=c(latlonBox[1], latlonBox[2])) +
    labs(title=myTitle, subtitle=mySubtitle, 
         x="Longitude", y="Latitude", 
         fill=expression(tons / hectare / year))
	
	
geomData = data.frame(lat = rep(nc_lat, each=length(nc_lon)), lon = rep(nc_lon, length(nc_lat)), plotData = as.vector(nc_testDat[ , , 9, 2, 1] - nc_testDat[ , , 1, 2, 1]))
df_sub = subset(geomData, lat > -50)
lmt <- ceiling(10*max(abs(range(df_sub$plotData, na.rm=TRUE))))/20
latlonBox = c(-65, 90, -180, 180)#c(25, 65, -10, 60) #(minlat, maxlat, minlon, maxlon)
myTitle = 'Change in Annual Avg Wood Mass'
mySubtitle = 'for temperate needleleaf evergreen trees [kg / m^2]'
ggplot(data = df_sub, aes(fill = plotData, y = lat, x = lon)) +
	geom_tile(width = .5, height = .5)	+
#	geom_sf(data = df_sub, aes(fill = plotData, color = plotData)) +
	borders('world', xlim=range(df_sub$lon), ylim=range(df_sub$lat), 
            colour='gray10', size=.2)	+
	theme_minimal()	+
	theme(panel.ontop=FALSE, panel.background=element_blank()) +
#	scale_fill_viridis()	+
#	scale_colour_distiller(palette='Spectral', 
#                         limits=c(quantile(df_sub, .7, na.rm=T), 
#                                  quantile(df_sub, .999, na.rm=T))) +
	scale_fill_stepsn(colors=c('#b2182b','#ef8a62','#fddbc7',
                             '#f7f7f7','#d1e5f0','#67a9cf','#2166ac'),
             n.breaks=10, limits=c(-lmt,lmt), show.limits=T)  +
	coord_quickmap(xlim=c(latlonBox[3], latlonBox[4]), ylim=c(latlonBox[1], latlonBox[2])) +
    labs(title=myTitle, subtitle=mySubtitle, 
         x="Longitude", y="Latitude", 
         fill=expression(kg / m^2))
	
	




library(rnaturalearth)
library(maps)
sf_use_s2(FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf")
#usa = ne_states(country='usa')

st_as_sf(world)
geomData = data.frame(lat = rep(nc_lat, each=length(nc_lon)), lon = rep(nc_lon, length(nc_lat)), plotData = as.vector(nc_testDat[ , , 5, 2, 3]))
geomData_sf = st_as_sf(geomData, coords=c('lon','lat'), remove=FALSE, crs=4326, agr='constant')

geomJoin = st_join(world, geomData_sf, join = st_intersects)

geomAvg = aggregate(plotData ~ sovereignt, mean, data = geomJoin)#, na.action = na.pass)
geomAvgMrg = merge(geomAvg, world, by = 'sovereignt')
geomAvgMrg_sf = st_as_sf(geomAvgMrg)


ggplot(data = geomAvgMrg_sf) +
#	scale_colour_viridis_c(alpha = 1) +
#	scale_colour_gradientn(colours=rev(c('#F2F3F3', '#F2F3F3', '#FFF2F2', '#FFE5E5', 'orangered1', 'red1','red3','firebrick')),name='drought') +
#	scale_fill_gradient2(low='firebrick', mid='#DBDDDF', high='#196CE1', midpoint = .5, na.value = NA, limits = c(0,1),
#		breaks = c(.25, .5, .75), labels = c('severe drought', 'drought', 'wet conditions'), name='SPEI 3 month') +
#	geom_sf(data = california, fill = '#DBDDDF', color = '#1A232F', size=1.4) +
#	geom_sf(data = states, fill = '#DBDDDF')	+
#	geom_sf(data = counties, fill = NA, color = gray(.5))
#	geom_sf(data = world, fill = NA, color = gray(0.05), size=0.8) +
	geom_sf(data = geomAvgMrg_sf, size=6.4, shape=15, aes(fill=plotData)) +
#	scale_fill_stepsn(colors=c('#b2182b','#ef8a62','#fddbc7',
#                            '#f7f7f7','#d1e5f0','#67a9cf','#2166ac'),
#				             n.breaks=10, limits=c(-lmt,lmt), show.limits=T)  +
	scale_fill_viridis_c(option = 'turbo', trans = scales::pseudo_log_trans(sigma = 0.001), direction = -1)+#, trans = 'log10') +
#	ylim(-2,2) geomAvgMrg_sf
#	geom_sf(data = mask, fill='#F2F3F3', color='#F2F3F3') +
#	geom_sf(data = counties, fill = NA, color = '#666D74') +
#	geom_sf(data = states, fill = NA, color = '#1A232F', size=1.4) +
#	geom_sf(data = cityLocs_sf, size=4.4, stroke = 1.1, shape=21, col='white', fill='grey20')+ #c(rep('#039CE2',8), rep('#FDB600',7), '#23AF41')) +
	theme(panel.grid = element_line(colour='#F2F3F3')) +
	theme(panel.background = element_rect(fill = '#F2F3F3', colour='#F2F3F3')) +
	theme(legend.position = 'bottom',
		legend.background = element_rect(fill='white', colour='grey10')) +
	guides(fill=guide_legend(title="Change in Annual Avg Wood Mass (temperate needleleaf evergreen) [kg / m^2]")) +
	coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE)
#	


	theme_light()	+
	theme(panel.ontop=FALSE, panel.background=element_blank()) +
	scale_fill_viridis()	+
