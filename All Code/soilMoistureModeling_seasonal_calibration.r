library(ncdf4)
library(magrittr)
library(maptools)
library(lubridate)
library(EcoHydRology)
library(hydroGOF)		
library(data.table)

#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\Advanta\\"
	# reading in temp and ppt data 
	# names and variables
dataPath = 'J:\\Cai_data\\Advanta\\GDD\\'
era5ClimatologyDataName = 'testing-climatology-era.nc'

ncin = nc_open(paste0(dataPath, era5ClimatologyDataName)) # t2m_max[longitude,latitude,time]   (Contiguous storage) 
nc_ppt = ncvar_get(ncin, 'tp_sum')	
nc_tmax = ncvar_get(ncin, 't2m_max')
nc_tmin = ncvar_get(ncin, 't2m_min')
	
nc_lat = ncvar_get(ncin, 'latitude')
nc_lon = ncvar_get(ncin, 'longitude')

nc_dates = ncvar_get(ncin, 'time') + as.Date('2002-06-01')
these_dates = which(nc_dates %in% seq(as.Date('2016-01-01'), as.Date('2021-12-31'), by=1))
calDates = which(nc_dates %in% seq(as.Date('2017-01-01'), as.Date('2021-12-31'), by=1))


	# reading in the soil moisture data for calibration
ncpath_sm = "J:\\Cai_data\\Advanta\\SoilMoisture\\"

ncname_sm = paste0('sm_pct_2017.nc')
ncin_sm17 = nc_open(paste0(ncpath_sm, ncname_sm))

nc_lat_sm = ncvar_get(ncin_sm17, 'latitude')
nc_lon_sm = ncvar_get(ncin_sm17, 'longitude')
theseSMlats = which(nc_lat_sm <= max(nc_lat))
theseSMlons = which(nc_lon_sm >= min(nc_lon))
nc_lat_sm = ncvar_get(ncin_sm17, 'latitude')[theseSMlats]
nc_lon_sm = ncvar_get(ncin_sm17, 'longitude')[theseSMlons]

nc_sm17 = ncvar_get(ncin_sm17, 'sm_pct')[theseSMlons, theseSMlats, ]
ncname_sm = paste0('sm_pct_2018.nc')
nc_sm18 = ncvar_get(nc_open(paste0(ncpath_sm, ncname_sm)), 'sm_pct')[theseSMlons, theseSMlats, ]
ncname_sm = paste0('sm_pct_2019.nc')
nc_sm19 = ncvar_get(nc_open(paste0(ncpath_sm, ncname_sm)), 'sm_pct')[theseSMlons, theseSMlats, ]
ncname_sm = paste0('sm_pct_2020.nc')
nc_sm20 = ncvar_get(nc_open(paste0(ncpath_sm, ncname_sm)), 'sm_pct')[theseSMlons, theseSMlats, ]
ncname_sm = paste0('sm_pct_2021.nc')
nc_sm21 = ncvar_get(nc_open(paste0(ncpath_sm, ncname_sm)), 'sm_pct')[theseSMlons, theseSMlats, ]

rootDepth = 1000

calibDF = data.frame(Lat = NA, Lon = NA, KGE = NA, NSE = NA, MAE = NA, Zr = NA,
	rcScl = NA, rcExp = NA, PETexp = NA, n = NA, smhp = NA, smwp = NA, smfc = NA, crn = NA, Ia_scalar = NA)


iter = 0
for(i in 1:length(nc_lat))	{
		#identifying nearest lat in the SM ncdf
	closeSMLat = which.min(abs(nc_lat[i] - nc_lat_sm))

	for(j in 1:length(nc_lon))	{
		iter = iter + 1
			#identify nearest lon in the SM ncdf
		closeSMLon = which.min(abs(nc_lon[j] - nc_lon_sm))

		allSM = c(nc_sm17[closeSMLon, closeSMLat, ], nc_sm18[closeSMLon, closeSMLat, ], nc_sm19[closeSMLon, closeSMLat, ], nc_sm20[closeSMLon, closeSMLat, ], nc_sm21[closeSMLon, closeSMLat, ])
				
		if(any(!is.na(allSM)))	{
		
				# soil moisture routine
			this_ppt = nc_ppt[j,i,these_dates]	; this_ppt[is.na(this_ppt)] = 0

			PET = PET_fromTemp(yday(nc_dates[these_dates]), nc_tmax[j,i,these_dates], nc_tmin[j,i,these_dates],
						lat_radians =  min((nc_lat[i]*pi/180), 1.1)) * 1000 # output in m, converting to mm

			# start of calibration routine
			nruns = 5000

			rcScl = sample(seq(0.01,1,0.01),nruns, replace = TRUE) #0.1,	# water retention curve scalar
			rcExp = sample(seq(0.4,2.2,0.01),nruns, replace = TRUE)#1.3,	# water retention curve exponent
			PETexp = sample(seq(0.1,1.2,0.01),nruns, replace = TRUE) 	# exponential decay of relative SM
			n = sample(seq(0.05,0.7,0.01),nruns, replace = TRUE) #0.5,	# soil porosity
			smhp = sample(seq(0.00,0.00,0.0001),nruns, replace = TRUE) #0.00,	# soil moisture at hydroscopic point
			smwp = sample(seq(0.025,0.22,0.01),nruns, replace = TRUE) #0.10,	# soil moisture at wilting point
			smfc = sample(seq(0.05,0.5,0.01),nruns, replace = TRUE) #0.25,	# soil moisture field capacity
			crn = sample(seq(20,99,1),nruns, replace = TRUE)#50, # Curve number (1-100)
			Ia_scalar = sample(seq(0.01,0.5,0.01),nruns, replace=TRUE)# ppt vol before runoff, used in curve number
			Zr = rootDepth
	
			theseMAE = NULL
			theseNSE = NULL
			theseKGE = NULL
			for(l in 1:nruns)	{
				infil = infiltration_f(PPT = this_ppt,
					crn = crn[l],#50,
					Smax = Zr * n[l],
					Ia_scalar = Ia_scalar[l])#0.2)


				sm_out = SM_routine_f(infiltration = infil, PET = PET,
						rcScl = rcScl[l],#0.1,	# water retention curve scalar
						rcExp = rcExp[l],#1.3,	# water retention curve exponent
						PETexp = PETexp[l], 	# max PET fraction per day
						Zr = Zr,	# root depth
						n = n[l], #0.5,	# soil porosity
						smhp = smhp[l], #0.00,	# soil moisture at hydroscopic point
						smwp = smhp[l], #0.10,	# soil moisture at wilting point
						smfc = smfc[l], #0.25,	# soil moisture field capacity
						s0 = .5)	# initial soil moisture 
				theseMAE = c(theseMAE, mae(allSM, sm_out[-c(1:366, length(sm_out))] / Zr))
				theseNSE = c(theseNSE, NSE(allSM, sm_out[-c(1:366, length(sm_out))] / Zr))
				theseKGE = c(theseKGE, KGE(allSM, sm_out[-c(1:366, length(sm_out))] / Zr))
				if(l %% 1000 == 1)	{
					plot(allSM)
					lines(sm_out[-c(1:366, length(sm_out))] / Zr)
				}
			}
			bestRun = which.max(theseKGE)

			calibDF[iter, ] = c(nc_lat[i], nc_lon[j], theseKGE[bestRun], theseNSE[bestRun], theseMAE[bestRun], Zr,
				rcScl[bestRun], rcExp[bestRun], PETexp[bestRun], n[bestRun],
				smhp[bestRun], smwp[bestRun], smfc[bestRun], crn[bestRun], Ia_scalar[bestRun])

			fwrite(calibDF, "J:\\Cai_data\\Advanta\\SoilMoisture_redux_newsm_20JUL2022.csv")
			print(summary(theseKGE))
			print(calibDF[iter, ])
		}	else	{
			calibDF[iter, ] = NA
		}
	}
}

fwrite(calibDF, "J:\\Cai_data\\Advanta\\SoilMoisture_redux_newsm_20JUL2022.csv")
	
		
plot(rcScl, theseKGE)
summary(lm(theseKGE ~ rcScl))

plot(rcExp, theseKGE)
summary(lm(theseKGE ~ rcExp))

plot(PETexp, theseKGE)	# most imp
summary(lm(theseKGE ~ PETexp))

plot(n, theseKGE)
summary(lm(theseKGE ~ n))

plot(smhp, theseKGE)
summary(lm(theseKGE ~ smhp))

plot(smwp, theseKGE)
summary(lm(theseKGE ~ smwp))

plot(smfc, theseKGE)
summary(lm(theseKGE ~ smfc))

plot(crn, theseKGE)
summary(lm(theseKGE ~ crn))

plot(Ia_scalar, theseKGE)
summary(lm(theseKGE ~ Ia_scalar))























###########################################################
### forecasting SM

#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\Advanta\\"

# SEAS5 data 
ncname = paste0("tp_bias_corrected_curr_Australia.nc")  
ncfname = paste0(ncpath, ncname)
ncin_ppt = nc_open(ncfname) # tp[longitude,latitude,member,time]   (Contiguous storage) 
ncname = paste0("t2m_max_bias_corrected_curr_Australia.nc")  
ncfname = paste0(ncpath, ncname)
ncin_tmax = nc_open(ncfname) # t2m_max[longitude,latitude,member,time]   (Contiguous storage) 
ncname = paste0("t2m_min_bias_corrected_curr_Australia.nc")  
ncfname = paste0(ncpath, ncname)
ncin_tmin = nc_open(ncfname) # t2m_max[longitude,latitude,member,time]   (Contiguous storage) 

nc_lat = ncvar_get(ncin_tmin, 'latitude')
nc_lon = ncvar_get(ncin_tmin, 'longitude') 
nc_time = ncvar_get(ncin_tmin, 'time')  # hours since 1900/01/01
nc_dates = nc_time / 24 + as.Date("1900-01-01")
numModels = length(ncvar_get(ncin_tmin, 'member'))

nc_tmax = ncvar_get(ncin_tmax, 't2m_max')
nc_tmin = ncvar_get(ncin_tmin, 't2m_min')
nc_ppt = ncvar_get(ncin_ppt, 'tp')	



	# reading in the soil moisture data for calibration
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# note that the updated historical soil moisture database must be downloaded prior to running model
ncpath_sm = "J:\\Cai_data\\Advanta\\SoilMoisture\\"

ncname_sm = paste0('sm_pct_2022(1).nc')
ncin_sm = nc_open(paste0(ncpath_sm, ncname_sm))
nc_sm = ncvar_get(ncin_sm, 'sm_pct')

nc_date_sm = as.Date("1900-01-01") + ncvar_get(ncin_sm, 'time') # time is days after jan 1 1900
nc_lat_sm = ncvar_get(ncin_sm, 'latitude')
nc_lon_sm = ncvar_get(ncin_sm, 'longitude')


calibDF = fread("J:\\Cai_data\\Advanta\\SoilMoisture_redux.csv")

summaryOutput_df = data.frame(Lat = 0, Lon = 0, soilDepth = 0, monthsOut = 0, startDate = as.Date('2000-01-01'),
	variableName = 'string', Q05 = 0, Q25 = 0, Q50 = 0, Q75 = 0, Q95 = 0)
soilDepth = 10 # depth in cm

iter = 0
for(i in 1:length(nc_lat))	{
		#identifying nearest lat in the SM ncdf
	closeSMLat = which.min(abs(nc_lat[i] - nc_lat_sm))

	for(j in 1:length(nc_lon))	{
			#identify nearest lon in the SM ncdf
		closeSMLon = which.min(abs(nc_lon[j] - nc_lon_sm))
		recentSM = (ncvar_get(thisNCDF, 'sm_pct')[closeSMLon, closeSMLat, ])

			#identify the relevant row for the calibration database
		closeCalRow = subset(calibDF, Lat == nc_lat[i] & Lon == nc_lon[j])
		
		if(any(!is.na(recentSM)))	{
			finalRealDataDate = last(nc_date_sm)
			initialSM = recentSM[finalRealDataDate]
			skipDates = which(nc_dates <= as.vector(finalRealDataDate))
			
			incPrevSM = which(nc_date_sm == nc_dates[1]):length(nc_date_sm)
			prevoiusSM = recentSM[incPrevSM]
			
			smOutAll = matrix(nrow = length(nc_dates), ncol = numModels)
			for(k in 1:numModels)	{
				# soil moisture routine
				this_ppt = nc_ppt[j,i,k,-skipDates]	; this_ppt[is.na(this_ppt)] = 0

				PET = PET_fromTemp(yday(nc_dates[-skipDates]), nc_tmax[j,i,k,-skipDates], nc_tmin[j,i,k,-skipDates],
							lat_radians =  min((nc_lat[i]*pi/180), 1.1)) * 1000 # output in m, converting to mm

				infil = infiltration_f(PPT = this_ppt,
					crn = closeCalRow$crn,#50,
					Smax = closeCalRow$Zr * closeCalRow$n,
					Ia_scalar = closeCalRow$Ia_scalar)#0.2)


				sm_out = SM_routine_f(infiltration = infil, PET = PET,
						rcScl = closeCalRow$rcScl,#0.1,	# water retention curve scalar
						rcExp = closeCalRow$rcExp,#1.3,	# water retention curve exponent
						PETexp = closeCalRow$PETexp, 	# max PET fraction per day
						Zr = closeCalRow$Zr,	# root depth
						n = closeCalRow$n, #0.5,	# soil porosity
						smhp = closeCalRow$smhp, #0.00,	# soil moisture at hydroscopic point
						smwp = closeCalRow$smhp, #0.10,	# soil moisture at wilting point
						smfc = closeCalRow$smfc, #0.25,	# soil moisture field capacity
						s0 = initialSM)	# initial soil moisture 
				
				smOutAll[,k] = c(previousSM, sm_out)
			}
		
			for(ll in 1:floor(nc_dates / 30))	{
				iter = iter + 1
				thisMonthsOut = ll - 1
				thisTimeFrame = 1:30 + (thisMonthsOut * 30)
				summaryOutput_df[iter, ] = data.frame(nc_lat[i], nc_lon[j], soilDepth,
					monthsOut = thisMonthsOut, startDate = nc_dates[thisTimeFrame[1]],
					variableName = 'AbsoluteValue_Percent',
					quantile(smOutAll[thisTimeFrame, ], c(0.05, 0.25, 0.50, 0.75, 0.95))
				)
				
				#!!!!!!!!!!!!!!! need to include historical averages
				histAvg = 'TODO TODO TODO TODO TODO TODO TODO TODO'
				summaryOutput_df[iter, ] = data.frame(nc_lat[i], nc_lon[j], soilDepth,
					monthsOut = thisMonthsOut, startDate = nc_dates[thisTimeFrame[1]],
					variableName = 'RelativeValue_PercentChange',
					quantile(smOutAll[thisTimeFrame, ], c(0.05, 0.25, 0.50, 0.75, 0.95)) - histAvg
				)

				fwrite(calibDF, "J:\\Cai_data\\Advanta\\SoilMoisture_redux.csv")
			}
		}
	}
}











































plot(nc_lon, nc_lat, nc_sm[,,10])
