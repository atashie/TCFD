library(HBV.IANIGLA) # for running HBV model
library(prism) # for downloading PRISM climate data (to 1981); might replace by just usign mopex data
library(hddtools)		# for mopex data
library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)
	sf::sf_use_s2(FALSE) # for problem with intersecting spherical w flat
library(feather)
library(MazamaSpatialUtils)	# for linking lat lon to HUC codes
library(lubridate)
library(daymetr)
library(ddplyr) # for left_join() merging by date
library(EcoHydRology)	# for PET_fromTemp
library(ncdf4)

### hydrobasins data
basinAt_NorAm_polys = st_read("C:\\Users\\arik\\Documents\\PhD Research\\D4\\basinAt_NorAm_polys.gpkg")
basinAt_NorAm_centroid = st_centroid(basinAt_NorAm_polys)# some centroids put us in the ocean so check for errors


#################################################################
## Reading in and correcting the Cali FNF data

	# data for the Sacramento Valley
BND = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=BND&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")
BND$Date = ymd(unlist(strsplit(BND$DATE.TIME, " "))[seq(1,nrow(BND)*2,2)])
BND_area = 8900 # in sq miles
BND$Date = ymd(unlist(strsplit(BND$DATE.TIME, " "))[seq(1,nrow(BND)*2,2)])
BND$acrefeet_pd = as.numeric(BND$VALUE) * 1.983
	(sum(BND$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
BND_latlon = c(40.288611, -122.185556)

ORO = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")
ORO_area = 3624
ORO$Date = ymd(unlist(strsplit(ORO$DATE.TIME, " "))[seq(1,nrow(ORO)*2,2)])
ORO$acrefeet_pd = as.numeric(ORO$VALUE) * 1.983
	(sum(ORO$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
ORO_latlon = c(39.521667, -121.546667)
 
YRS = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=YRS&SensorNums=8&dur_code=D&Start=1906-01-01&End=2020-12-31")
YRS_area = 1200
YRS$Date = ymd(unlist(strsplit(YRS$DATE.TIME, " "))[seq(1,nrow(YRS)*2,2)])
YRS$acrefeet_pd = as.numeric(YRS$VALUE) * 1.983
	(sum(YRS$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
YRS_latlon = c(39.223611, -121.292500)

FOL = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=FOL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")
FOL_area = 1863
FOL$Date = ymd(unlist(strsplit(FOL$DATE.TIME, " "))[seq(1,nrow(FOL)*2,2)])
FOL$acrefeet_pd = as.numeric(FOL$VALUE) * 1.983
	(sum(FOL$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
FOL_latlon = c(38.704528, -121.164361)

#######################################
	# data for the San Jaoquin Valley
NML = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=NML&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")#Stanislaus River inflow to New Melones Lake
NML$Date = ymd(unlist(strsplit(NML$DATE.TIME, " "))[seq(1,nrow(NML)*2,2)])
NML_area = NA # in sq miles
NML$Date = ymd(unlist(strsplit(NML$DATE.TIME, " "))[seq(1,nrow(NML)*2,2)])
NML$acrefeet_pd = as.numeric(NML$VALUE) * 1.983
	(sum(NML$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
NML_latlon = c(37.95014259160364, -120.52076334699024)

DNP = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=DNP&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")#Tuolumne River inflow to New Don Pedro Reservoir
DNP_area = NA
DNP$Date = ymd(unlist(strsplit(DNP$DATE.TIME, " "))[seq(1,nrow(DNP)*2,2)])
DNP$acrefeet_pd = as.numeric(DNP$VALUE) * 1.983
	(sum(DNP$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
DNP_latlon = c(37.71860423395725, -120.394743108845)
 
EXC = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")#Merced River inflow to Lake McClure
EXC_area = NA
EXC$Date = ymd(unlist(strsplit(EXC$DATE.TIME, " "))[seq(1,nrow(EXC)*2,2)])
EXC$acrefeet_pd = as.numeric(EXC$VALUE) * 1.983
	(sum(EXC$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
EXC_latlon = c(37.60060608855698, -120.2641797970366)

MIL = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=MIL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")#San Joaquin River inflow to Millerton Lake
MIL_area = NA
MIL$Date = ymd(unlist(strsplit(MIL$DATE.TIME, " "))[seq(1,nrow(MIL)*2,2)])
MIL$acrefeet_pd = as.numeric(MIL$VALUE) * 1.983
	(sum(MIL$acrefeet_pd, na.rm=TRUE) / 33) / 1000000
MIL_latlon = c(37.003857840261354, -119.69170419625397)



this_FNF_data = DNP
this_FNF_latlon = DNP_latlon
this_FNF_area = DNP_area


#########################################
# reading in climai netcdf clim data
ncpath = "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\"
	# reading in tmin 
ncname = "CA_era5_tmin.nc"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'latitude')
nc_lon = ncvar_get(ncin, 'longitude')
nc_date = as.Date("1990-01-01") + ncvar_get(ncin, 'time') # time is days after jan 1 1990
nc_tmin = ncvar_get(ncin, "t2m_min")

	# reading in tmin 
ncname = "CA_era5_tmax.nc"  		# check to ensure the lat lons / times are the same for all files 
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_tmax = ncvar_get(ncin, "t2m_max")

	# reading in ppt 
ncname = "CA_era5_tp.nc"  			# check to ensure the lat lons / times are the same for all files 
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_ppt = ncvar_get(ncin, "tp")



#########################################
# linking the gage to the upstream basins


pt1 = st_point(rev(this_FNF_latlon))# gage location
lon_lat = st_sfc(pt1, crs=4326)
which_hydro_basin = st_intersects(lon_lat, basinAt_NorAm_polys)[[1]]# row of intersecting hydrobasins database 
these_hydro_basins = as.character(basinAt_NorAm_polys$HYBAS_ID[which_hydro_basin])# name of intersecting hydrobasins database
if(file.exists(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData")))	{
	gaged_upstreams = readRDS(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
	these_hydro_basins = unlist(gaged_upstreams)
	} else {
		gaged_upstreams = NULL
		if(these_hydro_basins %in% names(gaged_upstreams))	{
			these_hydro_basins = gaged_upstreams[[which(names(gaged_upstreams) == basinAt_NorAm_centroid[which_hydro_basin,"HYBAS_ID"])]]
		}	else	{
			HB_remaining = basinAt_NorAm_centroid[-which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins),]
			print(nrow(HB_remaining))
					
			while(any(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
				for(trueys in which(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
					new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == these_hydro_basins[trueys])
					new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
					these_hydro_basins = c(these_hydro_basins, new_HB_ID)
					print(these_hydro_basins)
					HB_remaining = HB_remaining[-new_HB_ID_rows,]
					print(nrow(HB_remaining))
				}
			}
			new_row = length(gaged_upstreams) + 1
			gaged_upstreams[[new_row]] = these_hydro_basins
			names(gaged_upstreams)[new_row] = these_hydro_basins[1]
			saveRDS(gaged_upstreams,
				file=paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
		}
	}


#########################################################################################
# initializing data for HBV
 
# first bringing in the climate data
tot_hydro_basins = length(these_hydro_basins)
sample_number = 0
meteor_dat_n = data.frame(Date = rep(NA, length(nc_date)), PPT = rep(NA, length(nc_date)), Tmin = rep(NA, length(nc_date)), Tmax = rep(NA, length(nc_date)))
hydroBAS_area = 0
for(num_hydro_basins in 1:tot_hydro_basins)	{
		
	which_hybas = which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins[num_hydro_basins])
	hydroBAS_area = hydroBAS_area + basinAt_NorAm_centroid$SUB_AREA[which_hybas]	# summing area of subbasins from HydroBASINS
	this_lonlat = st_coordinates(basinAt_NorAm_centroid[which_hybas,])
	nearest_lon = which.min(abs(this_lonlat[1] - nc_lon))
	nearest_lat = which.min(abs(this_lonlat[2] - nc_lat))
	
	meteor_dat_n$PPT = nc_ppt[nearest_lon, nearest_lat, ]
	meteor_dat_n$Tmin = nc_tmin[nearest_lon, nearest_lat, ]
	meteor_dat_n$Tmax = nc_tmax[nearest_lon, nearest_lat, ]

	sample_number = sample_number+1
	
	if(sample_number == 1)	{meteor_dat = meteor_dat_n} else {meteor_dat = meteor_dat_n + meteor_dat}
}
all_upst_hydrobasins = st_coordinates(basinAt_NorAm_centroid[which(basinAt_NorAm_centroid$HYBAS_ID %in% as.numeric(unlist(gaged_upstreams))),])
min(all_upst_hydrobasins[,1])
min(all_upst_hydrobasins[,2])
max(all_upst_hydrobasins[,1])
max(all_upst_hydrobasins[,2])

meteor_dat = meteor_dat / sample_number
meteor_dat$Date = nc_date		
meteor_dat$PET = PET_fromTemp(yday(meteor_dat$Date), meteor_dat$Tmax, meteor_dat$Tmin,
			lat_radians =  min((this_lonlat[1,2]*pi/180), 1.1)) * 1000
plot(meteor_dat$Date[1:500], meteor_dat$Tmin[1:500], type='l')
lines(meteor_dat$Date[1:500], meteor_dat_n$Tmin[1:500],col='red')

#if(is.na(this_FNF_area))	{this_FNF_area = hydroBAS_area * 0.386102}	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs
this_FNF_area = hydroBAS_area * 0.386102	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs



################################################
## combining climate data with streamflow
FNF_area_in_sqmm = this_FNF_area * 2.59*10^12
FNF_Q_in_cbmm = this_FNF_data$acrefeet_pd * 1.233*10^12
this_FNF_data$Q_in_mm = FNF_Q_in_cbmm / FNF_area_in_sqmm

mdt = as.data.table(meteor_dat)
fdt = as.data.table(this_FNF_data)
mopey_df = fdt[mdt, on='Date']
mopey_df$Tavg = apply(mopey_df[,c("Tmin","Tmax")],1,mean)



##############################################################
### setting the parameter ranges for the HBV runoff modules
# snow and glacier module
  sfcf = c(runif(500, .2, 1), runif(500, 1, 3)) #1, 2)	#snowfall correction factor [-]
  tr   = runif(1000, -6, 5)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(1000, -5, 6)#0, 3)	#melt temperature [C]
  fm   = c(runif(500, .2, 1.5), (runif(500, 1.5, 8)))#1, 4)	#snowmelt factor [mm/C]
  fi   = c(runif(500, .2, 1.5), (runif(500, 1.5, 10)))#4, 8)	#icemelt factor [mm/C]
  fic  = runif(1000, 2, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

# soil module
  fc   = c(runif(500, 25, 150), (runif(500, 150, 1200)))
  lp   = runif(1000, .2, 1)#0.5, 1)   # parameter to actual ET
  beta_soils = runif(1000, 1, 3)

# routing module
  k0   = c(runif(500, .05, .5), (runif(500, .5, 1)))#runif(1000, .01, 1)#0.09, 0.1)
  k1   = c(runif(500, .005, .09), (runif(500, .09, .5)))#runif(1000, .001, .5)#0.05, 0.07)
  k2   = c(runif(500, .0001, .01), (runif(500, .01, .1)))#runif(1000, .0001, .1)#0.05)	
  uz1  = c(runif(500, .22, 10), (runif(500, 10, 40)))#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
  perc = c(runif(500, .1, .5), (runif(500, .5, 20)))#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d

	
	#########################################################################################
    # initializing data for model run
    n_runs = 100000

    cal_out = data.frame(
      sfcf=rep(NA,n_runs),
      tr=rep(NA,n_runs), 
      tt=rep(NA,n_runs),
      fm=rep(NA,n_runs),
      fi=rep(NA,n_runs),
      fic=rep(NA,n_runs),
      fc=rep(NA,n_runs),
      lp=rep(NA,n_runs),
      beta_soils = rep(NA,n_runs),
      k0 = rep(NA,n_runs),
      k1 = rep(NA,n_runs),
      k2 = rep(NA,n_runs),
      uz1 = rep(NA,n_runs),
      perc = rep(NA,n_runs),
      hbv_nse = rep(NA,n_runs),
      new_nse = rep(NA,n_runs)
    )

# initializing the dataframe for capturing model performance
model_perf = data.frame(kge=rep(NA,n_runs))
 
 
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
        cal_out$sfcf[jj] = sample(sfcf,1)
        cal_out$tr[jj] = sample(tr,1)
        cal_out$tt[jj] = sample(tt,1)
        cal_out$fm[jj] = sample(fm,1)
        cal_out$fi[jj] = sample(fi,1)
        cal_out$fic[jj] = sample(fic,1)
        cal_out$fc[jj] = sample(fc,1)
        cal_out$lp[jj] = sample(lp,1)
        cal_out$beta_soils[jj] = sample(beta_soils,1)
        # since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
        cal_out$k0[jj] = sample(k0,1)
        cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
        cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
        cal_out$uz1[jj] = sample(uz1,1)
        cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)
        
        mopey_df_ppt = cbind(mopey_df,	
                             SnowGlacier_HBV(model = 1,
                                             inputData = cbind(mopey_df$Tavg, mopey_df$PPT),
                                             initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                             param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
                                                        cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                        cal_out$tt[jj],		#Tt - melt temperature [C]
                                                        cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
                                                        cal_out$fi[jj],		#fi - icemelt factor [mm/C]
                                                        cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
                             )	)
        
        # soils model
        mopey_df_rech = cbind(mopey_df_ppt,
                              Soil_HBV(
                                model = 1,
                                inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$PET),
                                initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                                param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
                                           cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
                                           cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
                              )	)
        
        
        mopey_df_disch = cbind(mopey_df_rech,
                               Routing_HBV(
                                 model = 1,	# model=1 gives three stores with routing for each
                                 lake = FALSE,
                                 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
                                 initCond = c(10,10,10),	# initial storage in each reservoir in mm
                                 param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
                                            cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                            cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                            cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                            cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                               )	)
        
    
	
	
	
# saving model performance and calib parameters
    # hbv modules
    model_perf$hbv_sfcf[jj] = cal_out$sfcf[jj]
    model_perf$hbv_tr[jj] = cal_out$tr[jj]
    model_perf$hbv_tt[jj] = cal_out$tt[jj]
    model_perf$hbv_fm[jj] = cal_out$fm[jj]
    model_perf$hbv_fi[jj] = cal_out$fi[jj]
    model_perf$hbv_fic[jj] = cal_out$fic[jj]
    # soil module
    model_perf$hbv_fc[jj] = cal_out$fc[jj]
    model_perf$hbv_lp[jj] = cal_out$lp[jj]
    model_perf$hbv_beta_soils[jj] = cal_out$beta_soils[jj]
    # routing module
    model_perf$hbv_k0[jj] = cal_out$k0[jj]
    model_perf$hbv_k1[jj] = cal_out$k1[jj]
    model_perf$hbv_k2[jj] = cal_out$k2[jj]
    model_perf$hbv_uz1[jj] = cal_out$uz1[jj]
    model_perf$hbv_perc[jj] = cal_out$perc[jj]
        
	model_perf$kge[jj] = KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
	model_perf$nse[jj] = NSE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
    model_perf$bias[jj] = pbias(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)]) * mean(mopey_df_disch$Q_in_mm[-(1:1090)], na.rm=TRUE)/ 100
    model_perf$rmse[jj] = rmse(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
    model_perf$mse[jj] = mse(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
    model_perf$mae[jj] = mae(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
 
	mopey_df_disch$year = year(mopey_df_disch$Date)
	which_years = unique(mopey_df_disch$year)
	mopey_first = subset(mopey_df_disch, year %in% which_years[4:13])
	mopey_mid = subset(mopey_df_disch, year %in% which_years[13:22])
	mopey_last = subset(mopey_df_disch, year %in% tail(which_years, 10))

    model_perf$kge_first[jj] = KGE(mopey_first$Qg, mopey_first$Q_in_mm)
	model_perf$nse_first[jj] = NSE(mopey_first$Qg, mopey_first$Q_in_mm)
    model_perf$bias_first[jj] = pbias(mopey_first$Qg, mopey_first$Q_in_mm) * mean(mopey_first$Q_in_mm[-(1:1090)], na.rm=TRUE)/ 100
    model_perf$rmse_first[jj] = rmse(mopey_first$Qg, mopey_first$Q_in_mm)
    model_perf$mse_first[jj] = mse(mopey_first$Qg, mopey_first$Q_in_mm)
    model_perf$mae_first[jj] = mae(mopey_first$Qg, mopey_first$Q_in_mm)
 
    model_perf$kge_last[jj] = KGE(mopey_last$Qg, mopey_last$Q_in_mm)
	model_perf$nse_last[jj] = NSE(mopey_last$Qg, mopey_last$Q_in_mm)
    model_perf$bias_last[jj] = pbias(mopey_last$Qg, mopey_last$Q_in_mm) * mean(mopey_last$Q_in_mm[-(1:1090)], na.rm=TRUE)/ 100
    model_perf$rmse_last[jj] = rmse(mopey_last$Qg, mopey_last$Q_in_mm)
    model_perf$mse_last[jj] = mse(mopey_last$Qg, mopey_last$Q_in_mm)
    model_perf$mae_last[jj] = mae(mopey_last$Qg, mopey_last$Q_in_mm)
 
    model_perf$kge_mid[jj] = KGE(mopey_mid$Qg, mopey_mid$Q_in_mm)
	model_perf$nse_mid[jj] = NSE(mopey_mid$Qg, mopey_mid$Q_in_mm)
    model_perf$bias_mid[jj] = pbias(mopey_mid$Qg, mopey_mid$Q_in_mm) * mean(mopey_mid$Q_in_mm[-(1:1090)], na.rm=TRUE)/ 100
    model_perf$rmse_mid[jj] = rmse(mopey_mid$Qg, mopey_mid$Q_in_mm)
    model_perf$mse_mid[jj] = mse(mopey_mid$Qg, mopey_mid$Q_in_mm)
    model_perf$mae_mid[jj] = mae(mopey_mid$Qg, mopey_mid$Q_in_mm)
 
 # plotting model performance
		if(jj %% 1000 == 1){ 
			these_days = 1095:(1095+365*3)
			
			plot(mopey_df_disch$Date[these_days], mopey_df_disch$Q_in_mm[these_days],
				 ylim=c(0.001,10), type='l',
	#             log='y',
				 main=paste("KGE =", round(KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)]),2)),
				 ylab = "Q (mm)", xlab = "Date")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
			lines(mopey_df_disch$Date[these_days], mopey_df_disch$Qg[these_days] + .000001, lwd=2, col='red', lty=1)
			write.csv(model_perf, paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", this_FNF_data$STATION_ID[1], ".csv"))
			print(jj)
			summary(model_perf)
		}
}
write.csv(model_perf, paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", this_FNF_data$STATION_ID[1], ".csv"))
    


	# rerunning each model with the best cal to capture streamflow preds
HBV_in_ls = list()
HBV_areas = NA
subbasin_ids = c("BND","ORO","YRS","FOL")
#subbasin_ids = c("NML","DNP","EXC","MIL")
for(i in 1:length(subbasin_ids))	{
	HBV_in_ls[[i]] = read.csv(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", subbasin_ids[i], ".csv"))
}

HBV_out_ls = list()
cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))
#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(i in c(1:length(HBV_in_ls)))	{

	if(i == 1)	{
		this_FNF_data = BND
		this_FNF_latlon = BND_latlon
		this_FNF_area = BND_area
	} else if(i == 2) {
		this_FNF_data = ORO
		this_FNF_latlon = ORO_latlon
		this_FNF_area = ORO_area
	} else if(i == 3) {
		this_FNF_data = YRS
		this_FNF_latlon = YRS_latlon
		this_FNF_area = YRS_area
	} else if(i == 4) {
		this_FNF_data = FOL
		this_FNF_latlon = FOL_latlon
		this_FNF_area = FOL_area
	}

	pt1 = st_point(rev(this_FNF_latlon))# gage location
	lon_lat = st_sfc(pt1, crs=4326)
	which_hydro_basin = st_intersects(lon_lat, basinAt_NorAm_polys)[[1]]# row of intersecting hydrobasins database 
	these_hydro_basins = as.character(basinAt_NorAm_polys$HYBAS_ID[which_hydro_basin])# name of intersecting hydrobasins database
	gaged_upstreams = NULL  # potentially save then read as gaged_upstreams = readRDS("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams.RData")

	if(file.exists(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData")))	{
		gaged_upstreams = readRDS(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
		these_hydro_basins = unlist(gaged_upstreams)
		} else {
			gaged_upstreams = NULL
			if(these_hydro_basins %in% names(gaged_upstreams))	{
				these_hydro_basins = gaged_upstreams[[which(names(gaged_upstreams) == basinAt_NorAm_centroid[which_hydro_basin,"HYBAS_ID"])]]
			}	else	{
				HB_remaining = basinAt_NorAm_centroid[-which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins),]
				print(nrow(HB_remaining))
						
				while(any(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
					for(trueys in which(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
						new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == these_hydro_basins[trueys])
						new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
						these_hydro_basins = c(these_hydro_basins, new_HB_ID)
						print(these_hydro_basins)
						HB_remaining = HB_remaining[-new_HB_ID_rows,]
						print(nrow(HB_remaining))
					}
				}
				new_row = length(gaged_upstreams) + 1
				gaged_upstreams[[new_row]] = these_hydro_basins
				names(gaged_upstreams)[new_row] = these_hydro_basins[1]
				saveRDS(gaged_upstreams,
					file=paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
			}
		}



	#########################################################################################
	# initializing data for HBV
	 

	# first bringing in the climate data
	tot_hydro_basins = length(these_hydro_basins)
	sample_number = 0
	meteor_dat_n = data.frame(Date = rep(NA, length(nc_date)), PPT = rep(NA, length(nc_date)), Tmin = rep(NA, length(nc_date)), Tmax = rep(NA, length(nc_date)))
	hydroBAS_area = 0
	for(num_hydro_basins in 1:tot_hydro_basins)	{
		
		which_hybas = which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins[num_hydro_basins])
		hydroBAS_area = hydroBAS_area + basinAt_NorAm_centroid$SUB_AREA[which_hybas]	# summing area of subbasins from HydroBASINS
		this_lonlat = st_coordinates(basinAt_NorAm_centroid[which_hybas,])
		nearest_lon = which.min(abs(this_lonlat[1] - nc_lon))
		nearest_lat = which.min(abs(this_lonlat[2] - nc_lat))
		
		meteor_dat_n$PPT = nc_ppt[nearest_lon, nearest_lat, ]
		meteor_dat_n$Tmin = nc_tmin[nearest_lon, nearest_lat, ]
		meteor_dat_n$Tmax = nc_tmax[nearest_lon, nearest_lat, ]

		sample_number = sample_number+1
		
		if(sample_number == 1)	{meteor_dat = meteor_dat_n} else {meteor_dat = meteor_dat_n + meteor_dat}
	}
	all_upst_hydrobasins = st_coordinates(basinAt_NorAm_centroid[which(basinAt_NorAm_centroid$HYBAS_ID %in% as.numeric(unlist(gaged_upstreams))),])

	meteor_dat = meteor_dat / sample_number
	meteor_dat$Date = nc_date		
	meteor_dat$PET = PET_fromTemp(yday(meteor_dat$Date), meteor_dat$Tmax, meteor_dat$Tmin,
				lat_radians =  min((this_lonlat[1,2]*pi/180), 1.1)) * 1000
	plot(meteor_dat$Date[1:500], meteor_dat$Tmin[1:500], type='l')
	lines(meteor_dat$Date[1:500], meteor_dat_n$Tmin[1:500],col='red')

	this_FNF_area = hydroBAS_area * 0.386102	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs
	HBV_areas[i] = this_FNF_area



	################################################
	## combining climate data with streamflow
	FNF_area_in_sqmm = this_FNF_area * 2.59*10^12
	FNF_Q_in_cbmm = this_FNF_data$acrefeet_pd * 1.233*10^12
	this_FNF_data$Q_in_mm = FNF_Q_in_cbmm / FNF_area_in_sqmm

	mdt = as.data.table(meteor_dat)
	fdt = as.data.table(this_FNF_data)
	mopey_df = fdt[mdt, on='Date']
	mopey_df$Tavg = apply(mopey_df[,c("Tmin","Tmax")],1,mean)



	##############################################################
	### setting the parameter ranges for the HBV runoff modules
	# snow and glacier module
	thisBestCal = which.max(HBV_in_ls[[i]]$kge_mid)
	print(HBV_in_ls[[i]][thisBestCal,])
		
	cal_out$sfcf[i] = HBV_in_ls[[i]]$hbv_sfcf[thisBestCal]	#snowfall correction factor [-]
	  cal_out$tr[i]   = HBV_in_ls[[i]]$hbv_tr[thisBestCal]	#solid and liquid precipitation threshold temperature [C]
	  cal_out$tt[i]   = HBV_in_ls[[i]]$hbv_tt[thisBestCal]	#melt temperature [C]
	  cal_out$fm[i]   = HBV_in_ls[[i]]$hbv_fm[thisBestCal]	#snowmelt factor [mm/C]
	  cal_out$fi[i]   = HBV_in_ls[[i]]$hbv_fi[thisBestCal]	#icemelt factor [mm/C]
	  cal_out$fic[i]  = HBV_in_ls[[i]]$hbv_fic[thisBestCal]	#debris-covered icemelt factor [mm/C]

	# soil module
	  cal_out$fc[i]   = HBV_in_ls[[i]]$hbv_fc[thisBestCal]
	  cal_out$lp[i]   = HBV_in_ls[[i]]$hbv_lp[thisBestCal]   # parameter to actual ET
	  cal_out$beta_soils[i] = HBV_in_ls[[i]]$hbv_beta_soils[thisBestCal]

	# routing module
	  cal_out$k0[i]   = HBV_in_ls[[i]]$hbv_k0[thisBestCal]#runif(1000, .01, 1)#0.09, 0.1)
	  cal_out$k1[i]   = HBV_in_ls[[i]]$hbv_k1[thisBestCal]#runif(1000, .001, .5)#0.05, 0.07)
	  cal_out$k2[i]   = HBV_in_ls[[i]]$hbv_k2[thisBestCal]#runif(1000, .0001, .1)#0.05)	
	  cal_out$uz1[i]  = HBV_in_ls[[i]]$hbv_uz1[thisBestCal]#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
	  cal_out$perc[i] = HBV_in_ls[[i]]$hbv_perc[thisBestCal]#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d

		
		#########################################################################################
		# initializing data for model run
	 
		
		mopey_df_ppt = cbind(mopey_df,	
						 SnowGlacier_HBV(model = 1,
										 inputData = cbind(mopey_df$Tavg, mopey_df$PPT),
										 initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
										 param = c(	cal_out$sfcf[i],		#SFCF - snowfall correction factor [-]
													cal_out$tr[i],		#Tr - solid and liquid precipitation threshold temperature [C]
													cal_out$tt[i],		#Tt - melt temperature [C]
													cal_out$fm[i],		#fm - snowmelt factor [mm/C]
													cal_out$fi[i],		#fi - icemelt factor [mm/C]
													cal_out$fic[i])	 	#fic - debris-covered icemelt factor [mm/C]
						 )	)
			
		# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
							  Soil_HBV(
								model = 1,
								inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$PET),
								initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
								param = c(	cal_out$fc[i],			#FC - soil field capacity [mm]
										   cal_out$lp[i],			#LP - parameter ot get actual ET [-]
										   cal_out$beta_soils[i])	#beta - exponential value for nonlinear relations between soil storage and runoff
							  )	)
			
			
		mopey_df_disch = cbind(mopey_df_rech,
							   Routing_HBV(
								 model = 1,	# model=1 gives three stores with routing for each
								 lake = FALSE,
								 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
								 initCond = c(10,10,10),	# initial storage in each reservoir in mm
								 param = c(	cal_out$k0[i],	#KO - top bucket (STZ) storage constant [1/t]
											cal_out$k1[i],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
											cal_out$k2[i],	#K2 - lower bucket (SLZ) storage constant [1/t]
											cal_out$uz1[i],	#UZL - max flux rate between STZ and SUZ [mm/t] 
											cal_out$perc[i])	#PERC - max flux rate between SUZ and SLZ [mm/t]
							   )	)
		
		HBV_out_ls[[i]] = mopey_df_disch
		
		
		
		
	# plotting model performance
			these_days = 1095:(1095+365*3)+2000
			jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration_results_", subbasin_ids[i], ".jpg"), width=800,height=400, quality=100)
			par(cex=1.3)
			plot(mopey_df_disch$Date[these_days], mopey_df_disch$Q_in_mm[these_days],
				 #ylim=c(0.001,10),
				 type='l', lwd=2, cex=2,
	#             log='y',
				 main=paste("KGE =", round(KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)]),2)),
				 ylab = "streamflow (mm/day)", xlab = "year")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
			lines(mopey_df_disch$Date[these_days], mopey_df_disch$Qg[these_days] , lwd=3.5, col='tomato', lty=1, cex=1.5)
			dev.off()
			
	# saving model performance and calib parameters
			
		cal_out$kge[i] = KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
		cal_out$nse[i] = NSE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
		cal_out$bias[i] = pbias(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)]) * mean(mopey_df_disch$Q_in_mm[-(1:1090)], na.rm=TRUE)/ 100
		cal_out$rmse[i] = rmse(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
		cal_out$mse[i] = mse(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
		cal_out$mae[i] = mae(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q_in_mm[-(1:1095)])
}



	# combining data and converting back to dumbass millions of acre-feet
theSac = HBV_out_ls[[1]][,"Date"]
theSac$Qg1 = HBV_out_ls[[1]][,"Qg"] 
theSac$Qg2 = HBV_out_ls[[2]][,"Qg"]
theSac$Qg3 = HBV_out_ls[[3]][,"Qg"]
theSac$Qg4 = HBV_out_ls[[4]][,"Qg"]

mm_per_foot = 304.8
acr_per_sqml = 640





##################################################################################
	# for the Sacramento River
for(j in 1:length(HBV_out_ls))	{
	if(j == 1)	{
		this_FNF_data = BND
		this_FNF_latlon = BND_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_1 = theSac$Qg1 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 2) {
		this_FNF_data = ORO
		this_FNF_latlon = ORO_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_2 = theSac$Qg2 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 3) {
		this_FNF_data = YRS
		this_FNF_latlon = YRS_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_3 = theSac$Qg3 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 4) {
		this_FNF_data = FOL
		this_FNF_latlon = FOL_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_4 = theSac$Qg4 * this_FNF_area * acr_per_sqml / mm_per_foot
	}
}
	
theSac$Qtot_mm = apply(theSac[,2:5], 1, sum)
theSac$Qtot_ml_acrft = apply(theSac[,6:9], 1, sum) / 1000000

	# summing by season and year
theSac$year = year(theSac$Date)
theSac$month = month(theSac$Date)

ann_sum = NULL
octMar_sum = NULL
aprJul_sum = NULL
augSep_sum = NULL
iter = 0
for(k in unique(theSac$year))	{
	print(k)
	thisYearSac = subset(theSac, year == k & month %in% 1:9)
	lastYearSac = subset(theSac, year == k-1 & month %in% 10:12)
	wySac = rbind(lastYearSac, thisYearSac)
	octMarSac = subset(wySac, month %in% c(10:12,1:3))
	aprJulSac = subset(wySac, month %in% c(4:7))
	augSepSac = subset(wySac, month %in% c(8,9))
	iter=iter+1
	ann_sum[iter] = sum(wySac$Qtot_ml_acrft)
	octMar_sum[iter] = sum(octMarSac$Qtot_ml_acrft)
	aprJul_sum[iter] = sum(aprJulSac$Qtot_ml_acrft)
	augSep_sum[iter] = sum(augSepSac$Qtot_ml_acrft)
}

theSac_summary_dat = cbind(unique(theSac$year), ann_sum, octMar_sum, aprJul_sum, augSep_sum)
colnames(theSac_summary_dat)[1] = 'wateryr'

valid_data = read.csv("C:\\Users\\arik\\Documents\\ClimateAi\\Wonderful\\valid_dat_SVnSJ.csv")


merge_dat = merge(theSac_summary_dat, valid_data, by='wateryr')[-1,]
summary(merge_dat)

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_annual_Sacramento_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$wysum, merge_dat$ann_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("Sacramento Annual Total: r-squared =", 
			round(summary(lm(merge_dat$ann_sum[c(21:nrow(merge_dat))] ~ merge_dat$wysum[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Annual Total: Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$wysum[c(21:nrow(merge_dat))], merge_dat$ann_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1,lwd=1.5,  lty=2, cex=1.5, col='grey20')
dev.off()
	

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_aprJul_Sacramento_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$apr.jul, merge_dat$aprJul_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("Sacramento Apr-Jul: r-squared =", 
		round(summary(lm(merge_dat$aprJul_sum[c(21:nrow(merge_dat))] ~ merge_dat$apr.jul[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$apr.jul[c(21:nrow(merge_dat))], merge_dat$aprJul_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()
	

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_octMar_Sacramento_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$oct.mar, merge_dat$octMar_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("Sacramento Oct-Mar: r-squared =", 
		round(summary(lm(merge_dat$octMar_sum[c(21:nrow(merge_dat))] ~ merge_dat$oct.mar[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$oct.mar[c(21:nrow(merge_dat))], merge_dat$octMar_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()
	

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation-v-calib_Sacramento_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(model_perf$kge_first, model_perf$kge_last,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("San Jaoquin cal-val vars: r-squared =", 
		round(summary(lm(model_perf$kge_last ~ model_perf$kge_first))$adj, 2)),
	ylab = "validation KGEs", xlab = "calibration KGEs")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()











##################################################################################
	# for the San Jaoquin River
for(j in 1:length(HBV_out_ls))	{
	if(j == 1)	{
		this_FNF_data = NML
		this_FNF_latlon = NML_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_1 = theSac$Qg1 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 2) {
		this_FNF_data = DNP
		this_FNF_latlon = DNP_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_2 = theSac$Qg2 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 3) {
		this_FNF_data = EXC
		this_FNF_latlon = EXC_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_3 = theSac$Qg3 * this_FNF_area * acr_per_sqml / mm_per_foot
	} else if(j == 4) {
		this_FNF_data = MIL
		this_FNF_latlon = MIL_latlon
		this_FNF_area = HBV_areas[j]
		theSac$Qacrft_4 = theSac$Qg4 * this_FNF_area * acr_per_sqml / mm_per_foot
	}
}
	
theSac$Qtot_mm = apply(theSac[,2:5], 1, sum)
theSac$Qtot_ml_acrft = apply(theSac[,6:9], 1, sum) / 1000000

	# summing by season and year
theSac$year = year(theSac$Date)
theSac$month = month(theSac$Date)

ann_sum = NULL
octMar_sum = NULL
aprJul_sum = NULL
augSep_sum = NULL
iter = 0
for(k in unique(theSac$year))	{
	print(k)
	thisYearSac = subset(theSac, year == k & month %in% 1:9)
	lastYearSac = subset(theSac, year == k-1 & month %in% 10:12)
	wySac = rbind(lastYearSac, thisYearSac)
	octMarSac = subset(wySac, month %in% c(10:12,1:3))
	aprJulSac = subset(wySac, month %in% c(4:7))
	augSepSac = subset(wySac, month %in% c(8,9))
	iter=iter+1
	ann_sum[iter] = sum(wySac$Qtot_ml_acrft)
	octMar_sum[iter] = sum(octMarSac$Qtot_ml_acrft)
	aprJul_sum[iter] = sum(aprJulSac$Qtot_ml_acrft)
	augSep_sum[iter] = sum(augSepSac$Qtot_ml_acrft)
}

theSac_summary_dat = cbind(unique(theSac$year), ann_sum, octMar_sum, aprJul_sum, augSep_sum)
colnames(theSac_summary_dat)[1] = 'wateryr'

valid_data = read.csv("C:\\Users\\arik\\Documents\\ClimateAi\\Wonderful\\valid_dat_SVnSJ.csv")


merge_dat = merge(theSac_summary_dat, valid_data, by='wateryr')[-1,]
summary(merge_dat)

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_annual_SanJaoquin_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$sj.wysum, merge_dat$ann_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("San Jaoquin Annual Total r-squared =",
		round(summary(lm(merge_dat$ann_sum[c(21:nrow(merge_dat))] ~ merge_dat$sj.wysum[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Annual Total: Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$sj.wysum[c(21:nrow(merge_dat))], merge_dat$ann_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1,lwd=1.5,  lty=2, cex=1.5, col='grey20')
dev.off()
	

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_aprJul_SanJaoquin_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$sj.apr.jul, merge_dat$aprJul_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("San Jaoquin Apr-Jul r-squared =",
		round(summary(lm(merge_dat$aprJul_sum[c(21:nrow(merge_dat))] ~ merge_dat$sj.apr.jul[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$sj.apr.jul[c(21:nrow(merge_dat))], merge_dat$aprJul_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()
	

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation_octMar_SanJaoquin_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(merge_dat$sj.oct.mar, merge_dat$octMar_sum,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("San Jaoquin Oct-Mar r-squared =",
		round(summary(lm(merge_dat$octMar_sum[c(21:nrow(merge_dat))] ~ merge_dat$sj.oct.mar[c(21:nrow(merge_dat))]))$adj, 2)),
	ylab = "Modeled Streamflow (mil acr-ft)", xlab = "Actual Streamflow (mil acr-ft)")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
points(merge_dat$sj.oct.mar[c(21:nrow(merge_dat))], merge_dat$octMar_sum[c(21:nrow(merge_dat))],
	col='tomato', lwd=3, cex=2)
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()

jpeg(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\validation-v-calib_SanJaoquin_.jpg"), width=600,height=600, quality=100)
par(cex=1.3)
plot(model_perf$kge_first, model_perf$kge_last,
	#ylim=c(0.001,10),
	type='p', lwd=2, cex=2,
#             log='y',
	main=paste("San Jaoquin cal-val vars: r-squared =", 
		round(summary(lm(model_perf$kge_last ~ model_perf$kge_first))$adj, 2)),
	ylab = "validation KGEs", xlab = "calibration KGEs")#,
	#             main=paste(mopey_catalogue[allws, "USGS_ID"],
	#                       ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
abline(0,1, lwd=1.5, lty=2, cex=1.5, col='grey20')
dev.off()















######################################################################
### end validation routines
### start predicting future values
######################################################################

#########################################
# reading in climai netcdf clim data
ncpath = "J:\\Cai_data\\Wonderful\\future_projs\\"
	# reading ssp126 climate data 
ncname = "ssp126_CA_cmip6_pp_future_daily.nc"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("2015-01-01") + ncvar_get(ncin, 'time') # time is days after jan 1 2015

ncname_ppt_1 = "ssp126_CA_cmip6_pp_future_2015_2044_daily_tp_0_25_deg.nc"  
ncname_ppt_2 = "ssp126_CA_cmip6_pp_future_2045_2074_daily_tp_0_25_deg.nc"  
ncfname_ppt_1 = paste0(ncpath, ncname_ppt_1)
ncfname_ppt_2 = paste0(ncpath, ncname_ppt_2)
ncin_ppt_1 = nc_open(ncfname_ppt_1)
ncin_ppt_2 = nc_open(ncfname_ppt_2)

nc_lat_ppt = ncvar_get(ncin_ppt_1, 'lat')
nc_lon_ppt = ncvar_get(ncin_ppt_1, 'lon')
nc_date_ppt_1 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_1, 'time') # time is days after jan 1 2015
nc_date_ppt_2 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_2, 'time') # time is days after jan 1 2015

date_range_ppt2 = seq(1,(length(nc_date) - length(nc_date_ppt_1)),1)

	# rerunning each model with the best cal to capture streamflow preds
HBV_in_ls = list()
HBV_areas = NA
subbasin_ids = c("BND","ORO","YRS","FOL")	;	Rsc = 1.04	# 
#subbasin_ids = c("NML","DNP","EXC","MIL")	;	Rsc = 0.98
for(i in 1:length(subbasin_ids))	{
	HBV_in_ls[[i]] = read.csv(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", subbasin_ids[i], ".csv"))
}

HBV_out_ls = list()	# for storing each subbasin
HBV_avg_out = data.frame(Date = nc_date) # for storing each model (combined subbasins)
HBVoutCol = 1
cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))
#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(thisClmMod in c(4,6:15))	{
	HBVoutCol = HBVoutCol + 1

#	ncvar_thismod_ppt =  ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] ugh annoying
	ncvar_thismod_ppt_1 = ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] and now merging two separate netcdfs
	ncvar_thismod_ppt_2 = ncvar_get(ncin_ppt_2, "tp")[date_range_ppt2,thisClmMod , , ]
	ncvar_thismod_tmin = ncvar_get(ncin, "t2m_min")[, , , thisClmMod]
	ncvar_thismod_tmax = ncvar_get(ncin, "t2m_max")[, , , thisClmMod]

	for(i in c(1:length(HBV_in_ls)))	{

		if(i == 1)	{
			this_FNF_data = BND
			this_FNF_latlon = BND_latlon
			this_FNF_area = BND_area
		} else if(i == 2) {
			this_FNF_data = ORO
			this_FNF_latlon = ORO_latlon
			this_FNF_area = ORO_area
		} else if(i == 3) {
			this_FNF_data = YRS
			this_FNF_latlon = YRS_latlon
			this_FNF_area = YRS_area
		} else if(i == 4) {
			this_FNF_data = FOL
			this_FNF_latlon = FOL_latlon
			this_FNF_area = FOL_area
		}

		pt1 = st_point(rev(this_FNF_latlon))# gage location
		lon_lat = st_sfc(pt1, crs=4326)
		which_hydro_basin = st_intersects(lon_lat, basinAt_NorAm_polys)[[1]]# row of intersecting hydrobasins database 
		these_hydro_basins = as.character(basinAt_NorAm_polys$HYBAS_ID[which_hydro_basin])# name of intersecting hydrobasins database
		gaged_upstreams = NULL  # potentially save then read as gaged_upstreams = readRDS("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams.RData")

		if(file.exists(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData")))	{
			gaged_upstreams = readRDS(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
			these_hydro_basins = unlist(gaged_upstreams)
			} else {
				gaged_upstreams = NULL
				if(these_hydro_basins %in% names(gaged_upstreams))	{
					these_hydro_basins = gaged_upstreams[[which(names(gaged_upstreams) == basinAt_NorAm_centroid[which_hydro_basin,"HYBAS_ID"])]]
				}	else	{
					HB_remaining = basinAt_NorAm_centroid[-which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins),]
					print(nrow(HB_remaining))
							
					while(any(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
						for(trueys in which(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
							new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == these_hydro_basins[trueys])
							new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
							these_hydro_basins = c(these_hydro_basins, new_HB_ID)
							print(these_hydro_basins)
							HB_remaining = HB_remaining[-new_HB_ID_rows,]
							print(nrow(HB_remaining))
						}
					}
					new_row = length(gaged_upstreams) + 1
					gaged_upstreams[[new_row]] = these_hydro_basins
					names(gaged_upstreams)[new_row] = these_hydro_basins[1]
					saveRDS(gaged_upstreams,
						file=paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
				}
			}



		#########################################################################################
		# initializing data for HBV
		 

		# first bringing in the climate data
		tot_hydro_basins = length(these_hydro_basins)
		sample_number = 0
		meteor_dat_n = data.frame(Date = rep(0, length(nc_date)), PPT = rep(NA, length(nc_date)), Tmin = rep(NA, length(nc_date)), Tmax = rep(NA, length(nc_date)))
		hydroBAS_area = 0
		for(num_hydro_basins in 1:tot_hydro_basins)	{
			
			which_hybas = which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins[num_hydro_basins])
			hydroBAS_area = hydroBAS_area + basinAt_NorAm_centroid$SUB_AREA[which_hybas]	# summing area of subbasins from HydroBASINS
			this_lonlat = st_coordinates(basinAt_NorAm_centroid[which_hybas,])
			nearest_lon = which.min(abs(this_lonlat[1] - nc_lon))
			nearest_lat = which.min(abs(this_lonlat[2] - nc_lat))
			nearest_lon_ppt = which.min(abs(this_lonlat[1] - nc_lon_ppt))
			nearest_lat_ppt = which.min(abs(this_lonlat[2] - nc_lat_ppt))
			
#			meteor_dat_n$PPT =  ncvar_thismod_ppt[ , nearest_lon_ppt, nearest_lat_ppt]	  * Rsc # is [time, lon, lat]
			meteor_dat_n$PPT =  c(ncvar_thismod_ppt_1[ , nearest_lon_ppt, nearest_lat_ppt], ncvar_thismod_ppt_2[ , nearest_lon_ppt, nearest_lat_ppt])	  * Rsc # is [time, lon, lat]
			meteor_dat_n$Tmin = ncvar_thismod_tmin[nearest_lon, nearest_lat, ]
			meteor_dat_n$Tmax = ncvar_thismod_tmax[nearest_lon, nearest_lat, ]

			sample_number = sample_number+1
			
			if(sample_number == 1)	{meteor_dat = meteor_dat_n} else {meteor_dat = meteor_dat_n + meteor_dat}
		}
		all_upst_hydrobasins = st_coordinates(basinAt_NorAm_centroid[which(basinAt_NorAm_centroid$HYBAS_ID %in% as.numeric(unlist(gaged_upstreams))),])

		meteor_dat = as.data.frame(meteor_dat / sample_number)
		
		meteor_dat$Date = nc_date		
		if(any(is.na(meteor_dat))) {meteor_dat[is.na(meteor_dat)] = 0}
		meteor_dat$PET = PET_fromTemp(yday(meteor_dat$Date), meteor_dat$Tmax, meteor_dat$Tmin,
					lat_radians =  min((this_lonlat[1,2]*pi/180), 1.1)) * 1000
		if(any(is.na(meteor_dat$PET))) {meteor_dat$PET[is.na(meteor_dat$PET)] = 0}
		
#		plot(meteor_dat$Date[1:500], meteor_dat$Tmin[1:500], type='l')
#		lines(meteor_dat$Date[1:500], meteor_dat_n$Tmin[1:500],col='red')


		this_FNF_area = hydroBAS_area * 0.386102	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs
		HBV_areas[i] = this_FNF_area



		################################################
		## convert to dt and est tavg
		mopey_df = as.data.table(meteor_dat)
		mopey_df$Tavg = apply(mopey_df[,c("Tmin","Tmax")],1,mean)



		##############################################################
		### setting the parameter ranges for the HBV runoff modules
		# snow and glacier module
		print(HBV_in_ls[[i]]$kge)
		thisBestCal = which.max(HBV_in_ls[[i]]$kge)
		  cal_out$sfcf[i] = HBV_in_ls[[i]]$hbv_sfcf[thisBestCal]	#snowfall correction factor [-]
		  cal_out$tr[i]   = HBV_in_ls[[i]]$hbv_tr[thisBestCal]	#solid and liquid precipitation threshold temperature [C]
		  cal_out$tt[i]   = HBV_in_ls[[i]]$hbv_tt[thisBestCal]	#melt temperature [C]
		  cal_out$fm[i]   = HBV_in_ls[[i]]$hbv_fm[thisBestCal]	#snowmelt factor [mm/C]
		  cal_out$fi[i]   = HBV_in_ls[[i]]$hbv_fi[thisBestCal]	#icemelt factor [mm/C]
		  cal_out$fic[i]  = HBV_in_ls[[i]]$hbv_fic[thisBestCal]	#debris-covered icemelt factor [mm/C]

		# soil module
		  cal_out$fc[i]   = HBV_in_ls[[i]]$hbv_fc[thisBestCal]
		  cal_out$lp[i]   = HBV_in_ls[[i]]$hbv_lp[thisBestCal]   # parameter to actual ET
		  cal_out$beta_soils[i] = HBV_in_ls[[i]]$hbv_beta_soils[thisBestCal]

		# routing module
		  cal_out$k0[i]   = HBV_in_ls[[i]]$hbv_k0[thisBestCal]#runif(1000, .01, 1)#0.09, 0.1)
		  cal_out$k1[i]   = HBV_in_ls[[i]]$hbv_k1[thisBestCal]#runif(1000, .001, .5)#0.05, 0.07)
		  cal_out$k2[i]   = HBV_in_ls[[i]]$hbv_k2[thisBestCal]#runif(1000, .0001, .1)#0.05)	
		  cal_out$uz1[i]  = HBV_in_ls[[i]]$hbv_uz1[thisBestCal]#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
		  cal_out$perc[i] = HBV_in_ls[[i]]$hbv_perc[thisBestCal]#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d

			
			#########################################################################################
			# initializing data for model run
		 
			
			mopey_df_ppt = cbind(mopey_df,	
							 SnowGlacier_HBV(model = 1,
											 inputData = cbind(mopey_df$Tavg, mopey_df$PPT),
											 initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
											 param = c(	cal_out$sfcf[i],		#SFCF - snowfall correction factor [-]
														cal_out$tr[i],		#Tr - solid and liquid precipitation threshold temperature [C]
														cal_out$tt[i],		#Tt - melt temperature [C]
														cal_out$fm[i],		#fm - snowmelt factor [mm/C]
														cal_out$fi[i],		#fi - icemelt factor [mm/C]
														cal_out$fic[i])	 	#fic - debris-covered icemelt factor [mm/C]
							 )	)
				
			# soils model
			mopey_df_rech = cbind(mopey_df_ppt,
								  Soil_HBV(
									model = 1,
									inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$PET),
									initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
									param = c(	cal_out$fc[i],			#FC - soil field capacity [mm]
											   cal_out$lp[i],			#LP - parameter ot get actual ET [-]
											   cal_out$beta_soils[i])	#beta - exponential value for nonlinear relations between soil storage and runoff
								  )	)
				
				
			mopey_df_disch = cbind(mopey_df_rech,
								   Routing_HBV(
									 model = 1,	# model=1 gives three stores with routing for each
									 lake = FALSE,
									 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
									 initCond = c(10,10,10),	# initial storage in each reservoir in mm
									 param = c(	cal_out$k0[i],	#KO - top bucket (STZ) storage constant [1/t]
												cal_out$k1[i],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
												cal_out$k2[i],	#K2 - lower bucket (SLZ) storage constant [1/t]
												cal_out$uz1[i],	#UZL - max flux rate between STZ and SUZ [mm/t] 
												cal_out$perc[i])	#PERC - max flux rate between SUZ and SLZ [mm/t]
								   )	)
			
			plot(mopey_df_disch$Date[1000:2000], mopey_df_disch$Qg[1000:2000], type='l')
			
			mopey_df_disch$QMilAcrFt = (mopey_df_disch$Qg * this_FNF_area * acr_per_sqml / mm_per_foot) / 1000000
			HBV_out_ls[[i]] = mopey_df_disch
	}
	
	HBV_basin_tot = data.frame(Date = HBV_out_ls[[1]]$Date, QMilAcrFt = 0)
	for(koo in 1:length(HBV_out_ls))	{
		HBV_basin_tot$QMilAcrFt = HBV_basin_tot$QMilAcrFt + HBV_out_ls[[koo]]$QMilAcrFt 
	}

	HBV_avg_out[,(ncol(HBV_avg_out)+1)] = HBV_basin_tot$QMilAcrFt
}
write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp126_sacramento.csv")
#write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp126_sanjaoquin.csv")


allModsAnn = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsOctMar = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAprJul = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAugSep = data.frame(Date = unique(year(HBV_avg_out$Date)))
HBV_avg_out$year = year(HBV_avg_out$Date)
HBV_avg_out$month = month(HBV_avg_out$Date)

for(eachClMod in 2:12)	{	# cols 2-13 have the output from each climate model
		# summing by season and year
	ann_sum = NULL
	octMar_sum = NULL
	aprJul_sum = NULL
	augSep_sum = NULL
	iter = 0
	
	for(k in unique(HBV_avg_out$year))	{
		thisYearSac = subset(HBV_avg_out, year == k & month %in% 1:9)
		lastYearSac = subset(HBV_avg_out, year == k-1 & month %in% 10:12)
		wySac = rbind(lastYearSac, thisYearSac)
		octMarSac = subset(wySac, month %in% c(10:12,1:3))
		aprJulSac = subset(wySac, month %in% c(4:7))
		augSepSac = subset(wySac, month %in% c(8,9))
		iter=iter+1
		ann_sum[iter] = sum(wySac[,eachClMod])
		octMar_sum[iter] = sum(octMarSac[,eachClMod])
		aprJul_sum[iter] = sum(aprJulSac[,eachClMod])
		augSep_sum[iter] = sum(augSepSac[,eachClMod])
	}
	allModsAnn[,eachClMod] = ann_sum
	allModsOctMar[,eachClMod] = octMar_sum
	allModsAprJul[,eachClMod] = aprJul_sum
	allModsAugSep[,eachClMod] = augSep_sum
}

summary(allModsAnn)
summary(allModsOctMar)
summary(allModsAprJul)
summary(allModsAugSep)
















#############################################################################################################################
	# reading ssp245 climate data 
ncname = "ssp245_CA_cmip6_pp_future_daily.nc"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("2015-01-01") + ncvar_get(ncin, 'time') # time is days after jan 1 2015

ncname_ppt_1 = "ssp245_CA_cmip6_pp_future_2015_2044_daily_tp_0_25_deg.nc"  
ncname_ppt_2 = "ssp245_CA_cmip6_pp_future_2045_2074_daily_tp_0_25_deg.nc"  
ncfname_ppt_1 = paste0(ncpath, ncname_ppt_1)
ncfname_ppt_2 = paste0(ncpath, ncname_ppt_2)
ncin_ppt_1 = nc_open(ncfname_ppt_1)
ncin_ppt_2 = nc_open(ncfname_ppt_2)

nc_lat_ppt = ncvar_get(ncin_ppt_1, 'lat')
nc_lon_ppt = ncvar_get(ncin_ppt_1, 'lon')
nc_date_ppt_1 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_1, 'time') # time is days after jan 1 2015
nc_date_ppt_2 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_2, 'time') # time is days after jan 1 2015

date_range_ppt2 = seq(1,(length(nc_date) - length(nc_date_ppt_1)),1)


	# rerunning each model with the best cal to capture streamflow preds
HBV_in_ls = list()
HBV_areas = NA
#subbasin_ids = c("BND","ORO","YRS","FOL")	;	Rsc = 1.04	# 
subbasin_ids = c("NML","DNP","EXC","MIL")	;	Rsc = 0.98
for(i in 1:length(subbasin_ids))	{
	HBV_in_ls[[i]] = read.csv(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", subbasin_ids[i], ".csv"))
}

HBV_out_ls = list()	# for storing each subbasin
HBV_avg_out = data.frame(Date = nc_date) # for storing each model (combined subbasins)
HBVoutCol = 1
cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))
#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(thisClmMod in c(4,6:15))	{
	HBVoutCol = HBVoutCol + 1

#	ncvar_thismod_ppt =  ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] ugh annoying
	ncvar_thismod_ppt_1 = ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] and now merging two separate netcdfs
	ncvar_thismod_ppt_2 = ncvar_get(ncin_ppt_2, "tp")[date_range_ppt2,thisClmMod , , ]
	ncvar_thismod_tmin = ncvar_get(ncin, "t2m_min")[, , , thisClmMod]
	ncvar_thismod_tmax = ncvar_get(ncin, "t2m_max")[, , , thisClmMod]

	for(i in c(1:length(HBV_in_ls)))	{

		if(i == 1)	{
			this_FNF_data = NML
			this_FNF_latlon = NML_latlon
			this_FNF_area = NML_area
		} else if(i == 2) {
			this_FNF_data = DNP
			this_FNF_latlon = DNP_latlon
			this_FNF_area = DNP_area
		} else if(i == 3) {
			this_FNF_data = EXC
			this_FNF_latlon = EXC_latlon
			this_FNF_area = EXC_area
		} else if(i == 4) {
			this_FNF_data = MIL
			this_FNF_latlon = MIL_latlon
			this_FNF_area = MIL_area
		}

		pt1 = st_point(rev(this_FNF_latlon))# gage location
		lon_lat = st_sfc(pt1, crs=4326)
		which_hydro_basin = st_intersects(lon_lat, basinAt_NorAm_polys)[[1]]# row of intersecting hydrobasins database 
		these_hydro_basins = as.character(basinAt_NorAm_polys$HYBAS_ID[which_hydro_basin])# name of intersecting hydrobasins database
		gaged_upstreams = NULL  # potentially save then read as gaged_upstreams = readRDS("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams.RData")

		if(file.exists(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData")))	{
			gaged_upstreams = readRDS(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
			these_hydro_basins = unlist(gaged_upstreams)
			} else {
				gaged_upstreams = NULL
				if(these_hydro_basins %in% names(gaged_upstreams))	{
					these_hydro_basins = gaged_upstreams[[which(names(gaged_upstreams) == basinAt_NorAm_centroid[which_hydro_basin,"HYBAS_ID"])]]
				}	else	{
					HB_remaining = basinAt_NorAm_centroid[-which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins),]
					print(nrow(HB_remaining))
							
					while(any(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
						for(trueys in which(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
							new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == these_hydro_basins[trueys])
							new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
							these_hydro_basins = c(these_hydro_basins, new_HB_ID)
							print(these_hydro_basins)
							HB_remaining = HB_remaining[-new_HB_ID_rows,]
							print(nrow(HB_remaining))
						}
					}
					new_row = length(gaged_upstreams) + 1
					gaged_upstreams[[new_row]] = these_hydro_basins
					names(gaged_upstreams)[new_row] = these_hydro_basins[1]
					saveRDS(gaged_upstreams,
						file=paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
				}
			}



		#########################################################################################
		# initializing data for HBV
		 

		# first bringing in the climate data
		tot_hydro_basins = length(these_hydro_basins)
		sample_number = 0
		meteor_dat_n = data.frame(Date = rep(0, length(nc_date)), PPT = rep(NA, length(nc_date)), Tmin = rep(NA, length(nc_date)), Tmax = rep(NA, length(nc_date)))
		hydroBAS_area = 0
		for(num_hydro_basins in 1:tot_hydro_basins)	{
			
			which_hybas = which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins[num_hydro_basins])
			hydroBAS_area = hydroBAS_area + basinAt_NorAm_centroid$SUB_AREA[which_hybas]	# summing area of subbasins from HydroBASINS
			this_lonlat = st_coordinates(basinAt_NorAm_centroid[which_hybas,])
			nearest_lon = which.min(abs(this_lonlat[1] - nc_lon))
			nearest_lat = which.min(abs(this_lonlat[2] - nc_lat))
			nearest_lon_ppt = which.min(abs(this_lonlat[1] - nc_lon_ppt))
			nearest_lat_ppt = which.min(abs(this_lonlat[2] - nc_lat_ppt))
			
#			meteor_dat_n$PPT =  ncvar_thismod_ppt[ , nearest_lon_ppt, nearest_lat_ppt]	  * Rsc # is [time, lon, lat]
			meteor_dat_n$PPT =  c(ncvar_thismod_ppt_1[ , nearest_lon_ppt, nearest_lat_ppt], ncvar_thismod_ppt_2[ , nearest_lon_ppt, nearest_lat_ppt])	  * Rsc # is [time, lon, lat]
			meteor_dat_n$Tmin = ncvar_thismod_tmin[nearest_lon, nearest_lat, ]
			meteor_dat_n$Tmax = ncvar_thismod_tmax[nearest_lon, nearest_lat, ]

			sample_number = sample_number+1
			
			if(sample_number == 1)	{meteor_dat = meteor_dat_n} else {meteor_dat = meteor_dat_n + meteor_dat}
		}
		all_upst_hydrobasins = st_coordinates(basinAt_NorAm_centroid[which(basinAt_NorAm_centroid$HYBAS_ID %in% as.numeric(unlist(gaged_upstreams))),])

		meteor_dat = as.data.frame(meteor_dat / sample_number)
		
		meteor_dat$Date = nc_date		
		if(any(is.na(meteor_dat))) {meteor_dat[is.na(meteor_dat)] = 0}
		meteor_dat$PET = PET_fromTemp(yday(meteor_dat$Date), meteor_dat$Tmax, meteor_dat$Tmin,
					lat_radians =  min((this_lonlat[1,2]*pi/180), 1.1)) * 1000
		if(any(is.na(meteor_dat$PET))) {meteor_dat$PET[is.na(meteor_dat$PET)] = 0}
#		plot(meteor_dat$Date[1:500], meteor_dat$Tmin[1:500], type='l')
#		lines(meteor_dat$Date[1:500], meteor_dat_n$Tmin[1:500],col='red')

		this_FNF_area = hydroBAS_area * 0.386102	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs
		HBV_areas[i] = this_FNF_area



		################################################
		## convert to dt and est tavg
		mopey_df = as.data.table(meteor_dat)
		mopey_df$Tavg = apply(mopey_df[,c("Tmin","Tmax")],1,mean)



		##############################################################
		### setting the parameter ranges for the HBV runoff modules
		# snow and glacier module
		#print(HBV_in_ls[[i]]$kge)
		thisBestCal = which.max(HBV_in_ls[[i]]$kge)
		  cal_out$sfcf[i] = HBV_in_ls[[i]]$hbv_sfcf[thisBestCal]	#snowfall correction factor [-]
		  cal_out$tr[i]   = HBV_in_ls[[i]]$hbv_tr[thisBestCal]	#solid and liquid precipitation threshold temperature [C]
		  cal_out$tt[i]   = HBV_in_ls[[i]]$hbv_tt[thisBestCal]	#melt temperature [C]
		  cal_out$fm[i]   = HBV_in_ls[[i]]$hbv_fm[thisBestCal]	#snowmelt factor [mm/C]
		  cal_out$fi[i]   = HBV_in_ls[[i]]$hbv_fi[thisBestCal]	#icemelt factor [mm/C]
		  cal_out$fic[i]  = HBV_in_ls[[i]]$hbv_fic[thisBestCal]	#debris-covered icemelt factor [mm/C]

		# soil module
		  cal_out$fc[i]   = HBV_in_ls[[i]]$hbv_fc[thisBestCal]
		  cal_out$lp[i]   = HBV_in_ls[[i]]$hbv_lp[thisBestCal]   # parameter to actual ET
		  cal_out$beta_soils[i] = HBV_in_ls[[i]]$hbv_beta_soils[thisBestCal]

		# routing module
		  cal_out$k0[i]   = HBV_in_ls[[i]]$hbv_k0[thisBestCal]#runif(1000, .01, 1)#0.09, 0.1)
		  cal_out$k1[i]   = HBV_in_ls[[i]]$hbv_k1[thisBestCal]#runif(1000, .001, .5)#0.05, 0.07)
		  cal_out$k2[i]   = HBV_in_ls[[i]]$hbv_k2[thisBestCal]#runif(1000, .0001, .1)#0.05)	
		  cal_out$uz1[i]  = HBV_in_ls[[i]]$hbv_uz1[thisBestCal]#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
		  cal_out$perc[i] = HBV_in_ls[[i]]$hbv_perc[thisBestCal]#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d

			
			#########################################################################################
			# initializing data for model run
		 
			
			mopey_df_ppt = cbind(mopey_df,	
							 SnowGlacier_HBV(model = 1,
											 inputData = cbind(mopey_df$Tavg, mopey_df$PPT),
											 initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
											 param = c(	cal_out$sfcf[i],		#SFCF - snowfall correction factor [-]
														cal_out$tr[i],		#Tr - solid and liquid precipitation threshold temperature [C]
														cal_out$tt[i],		#Tt - melt temperature [C]
														cal_out$fm[i],		#fm - snowmelt factor [mm/C]
														cal_out$fi[i],		#fi - icemelt factor [mm/C]
														cal_out$fic[i])	 	#fic - debris-covered icemelt factor [mm/C]
							 )	)
				
			# soils model
			mopey_df_rech = cbind(mopey_df_ppt,
								  Soil_HBV(
									model = 1,
									inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$PET),
									initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
									param = c(	cal_out$fc[i],			#FC - soil field capacity [mm]
											   cal_out$lp[i],			#LP - parameter ot get actual ET [-]
											   cal_out$beta_soils[i])	#beta - exponential value for nonlinear relations between soil storage and runoff
								  )	)
				
				
			mopey_df_disch = cbind(mopey_df_rech,
								   Routing_HBV(
									 model = 1,	# model=1 gives three stores with routing for each
									 lake = FALSE,
									 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
									 initCond = c(10,10,10),	# initial storage in each reservoir in mm
									 param = c(	cal_out$k0[i],	#KO - top bucket (STZ) storage constant [1/t]
												cal_out$k1[i],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
												cal_out$k2[i],	#K2 - lower bucket (SLZ) storage constant [1/t]
												cal_out$uz1[i],	#UZL - max flux rate between STZ and SUZ [mm/t] 
												cal_out$perc[i])	#PERC - max flux rate between SUZ and SLZ [mm/t]
								   )	)
			
			plot(mopey_df_disch$Date[1000:2000], mopey_df_disch$Qg[1000:2000], type='l')
			mopey_df_disch$QMilAcrFt = (mopey_df_disch$Qg  * this_FNF_area * acr_per_sqml / mm_per_foot) / 1000000
			HBV_out_ls[[i]] = mopey_df_disch
	}
	
	HBV_basin_tot = data.frame(Date = HBV_out_ls[[1]]$Date, QMilAcrFt = 0)
	for(koo in 1:length(HBV_out_ls))	{
		HBV_basin_tot$QMilAcrFt = HBV_basin_tot$QMilAcrFt + HBV_out_ls[[koo]]$QMilAcrFt 
	}

	HBV_avg_out[,(ncol(HBV_avg_out)+1)] = HBV_basin_tot$QMilAcrFt
}
#write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp245_sacramento.csv")
write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp245_sanjaoquin.csv")


allModsAnn = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsOctMar = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAprJul = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAugSep = data.frame(Date = unique(year(HBV_avg_out$Date)))
HBV_avg_out$year = year(HBV_avg_out$Date)
HBV_avg_out$month = month(HBV_avg_out$Date)

for(eachClMod in 2:12)	{	# cols 2-13 have the output from each climate model
		# summing by season and year
	ann_sum = NULL
	octMar_sum = NULL
	aprJul_sum = NULL
	augSep_sum = NULL
	iter = 0
	
	for(k in unique(HBV_avg_out$year))	{
		thisYearSac = subset(HBV_avg_out, year == k & month %in% 1:9)
		lastYearSac = subset(HBV_avg_out, year == k-1 & month %in% 10:12)
		wySac = rbind(lastYearSac, thisYearSac)
		octMarSac = subset(wySac, month %in% c(10:12,1:3))
		aprJulSac = subset(wySac, month %in% c(4:7))
		augSepSac = subset(wySac, month %in% c(8,9))
		iter=iter+1
		ann_sum[iter] = sum(wySac[,eachClMod])
		octMar_sum[iter] = sum(octMarSac[,eachClMod])
		aprJul_sum[iter] = sum(aprJulSac[,eachClMod])
		augSep_sum[iter] = sum(augSepSac[,eachClMod])
	}
	allModsAnn[,eachClMod] = ann_sum
	allModsOctMar[,eachClMod] = octMar_sum
	allModsAprJul[,eachClMod] = aprJul_sum
	allModsAugSep[,eachClMod] = augSep_sum
}

allModsAnn$avg = apply(allModsAnn[,-1], 1, mean)
summary(allModsAnn)
summary(allModsOctMar)
summary(allModsAprJul)
summary(allModsAugSep)





























##############################################################################################################################
	# reading ssp585 climate data 
ncname = "ssp585_CA_cmip6_pp_future_daily.nc"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("2015-01-01") + ncvar_get(ncin, 'time') # time is days after jan 1 2015

ncname_ppt_1 = "ssp585_CA_cmip6_pp_future_2015_2044_daily_tp_0_25_deg.nc"  
ncname_ppt_2 = "ssp585_CA_cmip6_pp_future_2045_2074_daily_tp_0_25_deg.nc"  
ncfname_ppt_1 = paste0(ncpath, ncname_ppt_1)
ncfname_ppt_2 = paste0(ncpath, ncname_ppt_2)
ncin_ppt_1 = nc_open(ncfname_ppt_1)
ncin_ppt_2 = nc_open(ncfname_ppt_2)

nc_lat_ppt = ncvar_get(ncin_ppt_1, 'lat')
nc_lon_ppt = ncvar_get(ncin_ppt_1, 'lon')
nc_date_ppt_1 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_1, 'time') # time is days after jan 1 2015
nc_date_ppt_2 = as.Date("2015-01-01") + ncvar_get(ncin_ppt_2, 'time') # time is days after jan 1 2015

date_range_ppt2 = seq(1,(length(nc_date) - length(nc_date_ppt_1)),1)



	# rerunning each model with the best cal to capture streamflow preds
HBV_in_ls = list()
HBV_areas = NA
#subbasin_ids = c("BND","ORO","YRS","FOL")	;	Rsc = 1.04	 
subbasin_ids = c("NML","DNP","EXC","MIL")	;	Rsc = 0.98
for(i in 1:length(subbasin_ids))	{
	HBV_in_ls[[i]] = read.csv(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\calibration__", subbasin_ids[i], ".csv"))
}

HBV_out_ls = list()	# for storing each subbasin
HBV_avg_out = data.frame(Date = nc_date) # for storing each model (combined subbasins)
HBVoutCol = 1
cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))
#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(thisClmMod in c(4,6:15))	{
	HBVoutCol = HBVoutCol + 1

#	ncvar_thismod_ppt =  ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] ugh annoying
	ncvar_thismod_ppt_1 = ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] and now merging two separate netcdfs
	ncvar_thismod_ppt_2 = ncvar_get(ncin_ppt_2, "tp")[date_range_ppt2,thisClmMod , , ]
	ncvar_thismod_tmin = ncvar_get(ncin, "t2m_min")[, , , thisClmMod]
	ncvar_thismod_tmax = ncvar_get(ncin, "t2m_max")[, , , thisClmMod]

	for(i in c(1:length(HBV_in_ls)))	{

		if(i == 1)	{
			this_FNF_data = NML
			this_FNF_latlon = NML_latlon
			this_FNF_area = NML_area
		} else if(i == 2) {
			this_FNF_data = DNP
			this_FNF_latlon = DNP_latlon
			this_FNF_area = DNP_area
		} else if(i == 3) {
			this_FNF_data = EXC
			this_FNF_latlon = EXC_latlon
			this_FNF_area = EXC_area
		} else if(i == 4) {
			this_FNF_data = MIL
			this_FNF_latlon = MIL_latlon
			this_FNF_area = MIL_area
		}

		pt1 = st_point(rev(this_FNF_latlon))# gage location
		lon_lat = st_sfc(pt1, crs=4326)
		which_hydro_basin = st_intersects(lon_lat, basinAt_NorAm_polys)[[1]]# row of intersecting hydrobasins database 
		these_hydro_basins = as.character(basinAt_NorAm_polys$HYBAS_ID[which_hydro_basin])# name of intersecting hydrobasins database
		gaged_upstreams = NULL  # potentially save then read as gaged_upstreams = readRDS("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams.RData")

		if(file.exists(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData")))	{
			gaged_upstreams = readRDS(paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
			these_hydro_basins = unlist(gaged_upstreams)
			} else {
				gaged_upstreams = NULL
				if(these_hydro_basins %in% names(gaged_upstreams))	{
					these_hydro_basins = gaged_upstreams[[which(names(gaged_upstreams) == basinAt_NorAm_centroid[which_hydro_basin,"HYBAS_ID"])]]
				}	else	{
					HB_remaining = basinAt_NorAm_centroid[-which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins),]
					print(nrow(HB_remaining))
							
					while(any(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
						for(trueys in which(these_hydro_basins %in% HB_remaining$NEXT_DOWN))	{
							new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == these_hydro_basins[trueys])
							new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
							these_hydro_basins = c(these_hydro_basins, new_HB_ID)
							print(these_hydro_basins)
							HB_remaining = HB_remaining[-new_HB_ID_rows,]
							print(nrow(HB_remaining))
						}
					}
					new_row = length(gaged_upstreams) + 1
					gaged_upstreams[[new_row]] = these_hydro_basins
					names(gaged_upstreams)[new_row] = these_hydro_basins[1]
					saveRDS(gaged_upstreams,
						file=paste0("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\HB_upstreams_", this_FNF_data$STATION_ID[1], ".RData"))
				}
			}



		#########################################################################################
		# initializing data for HBV
		 

		# first bringing in the climate data
		tot_hydro_basins = length(these_hydro_basins)
		sample_number = 0
		meteor_dat_n = data.frame(Date = rep(0, length(nc_date)), PPT = rep(NA, length(nc_date)), Tmin = rep(NA, length(nc_date)), Tmax = rep(NA, length(nc_date)))
		hydroBAS_area = 0
		for(num_hydro_basins in 1:tot_hydro_basins)	{
			
			which_hybas = which(basinAt_NorAm_centroid$HYBAS_ID == these_hydro_basins[num_hydro_basins])
			hydroBAS_area = hydroBAS_area + basinAt_NorAm_centroid$SUB_AREA[which_hybas]	# summing area of subbasins from HydroBASINS
			this_lonlat = st_coordinates(basinAt_NorAm_centroid[which_hybas,])
			nearest_lon = which.min(abs(this_lonlat[1] - nc_lon))
			nearest_lat = which.min(abs(this_lonlat[2] - nc_lat))
			nearest_lon_ppt = which.min(abs(this_lonlat[1] - nc_lon_ppt))
			nearest_lat_ppt = which.min(abs(this_lonlat[2] - nc_lat_ppt))
			
#			meteor_dat_n$PPT =  ncvar_thismod_ppt[ , nearest_lon_ppt, nearest_lat_ppt]	  * Rsc # is [time, lon, lat]
			meteor_dat_n$PPT =  c(ncvar_thismod_ppt_1[ , nearest_lon_ppt, nearest_lat_ppt], ncvar_thismod_ppt_2[ , nearest_lon_ppt, nearest_lat_ppt])	  * Rsc # is [time, lon, lat]
			meteor_dat_n$Tmin = ncvar_thismod_tmin[nearest_lon, nearest_lat, ]
			meteor_dat_n$Tmax = ncvar_thismod_tmax[nearest_lon, nearest_lat, ]

			sample_number = sample_number+1
			
			if(sample_number == 1)	{meteor_dat = meteor_dat_n} else {meteor_dat = meteor_dat_n + meteor_dat}
		}
		all_upst_hydrobasins = st_coordinates(basinAt_NorAm_centroid[which(basinAt_NorAm_centroid$HYBAS_ID %in% as.numeric(unlist(gaged_upstreams))),])

		meteor_dat = as.data.frame(meteor_dat / sample_number)
		
		meteor_dat$Date = nc_date		
		if(any(is.na(meteor_dat))) {meteor_dat[is.na(meteor_dat)] = 0}
		meteor_dat$PET = PET_fromTemp(yday(meteor_dat$Date), meteor_dat$Tmax, meteor_dat$Tmin,
					lat_radians =  min((this_lonlat[1,2]*pi/180), 1.1)) * 1000
		if(any(is.na(meteor_dat$PET))) {meteor_dat$PET[is.na(meteor_dat$PET)] = 0}
#		plot(meteor_dat$Date[1:500], meteor_dat$Tmin[1:500], type='l')
#		lines(meteor_dat$Date[1:500], meteor_dat_n$Tmin[1:500],col='red')

		this_FNF_area = hydroBAS_area * 0.386102	# estimating area from hydrobasins if USGS not available; also need to convert to sqmiles bs
		HBV_areas[i] = this_FNF_area



		################################################
		## convert to dt and est tavg
		mopey_df = as.data.table(meteor_dat)
		mopey_df$Tavg = apply(mopey_df[,c("Tmin","Tmax")],1,mean)



		##############################################################
		### setting the parameter ranges for the HBV runoff modules
		# snow and glacier module
		#print(HBV_in_ls[[i]]$kge)
		thisBestCal = which.max(HBV_in_ls[[i]]$kge)
		  cal_out$sfcf[i] = HBV_in_ls[[i]]$hbv_sfcf[thisBestCal]	#snowfall correction factor [-]
		  cal_out$tr[i]   = HBV_in_ls[[i]]$hbv_tr[thisBestCal]	#solid and liquid precipitation threshold temperature [C]
		  cal_out$tt[i]   = HBV_in_ls[[i]]$hbv_tt[thisBestCal]	#melt temperature [C]
		  cal_out$fm[i]   = HBV_in_ls[[i]]$hbv_fm[thisBestCal]	#snowmelt factor [mm/C]
		  cal_out$fi[i]   = HBV_in_ls[[i]]$hbv_fi[thisBestCal]	#icemelt factor [mm/C]
		  cal_out$fic[i]  = HBV_in_ls[[i]]$hbv_fic[thisBestCal]	#debris-covered icemelt factor [mm/C]

		# soil module
		  cal_out$fc[i]   = HBV_in_ls[[i]]$hbv_fc[thisBestCal]
		  cal_out$lp[i]   = HBV_in_ls[[i]]$hbv_lp[thisBestCal]   # parameter to actual ET
		  cal_out$beta_soils[i] = HBV_in_ls[[i]]$hbv_beta_soils[thisBestCal]

		# routing module
		  cal_out$k0[i]   = HBV_in_ls[[i]]$hbv_k0[thisBestCal]#runif(1000, .01, 1)#0.09, 0.1)
		  cal_out$k1[i]   = HBV_in_ls[[i]]$hbv_k1[thisBestCal]#runif(1000, .001, .5)#0.05, 0.07)
		  cal_out$k2[i]   = HBV_in_ls[[i]]$hbv_k2[thisBestCal]#runif(1000, .0001, .1)#0.05)	
		  cal_out$uz1[i]  = HBV_in_ls[[i]]$hbv_uz1[thisBestCal]#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
		  cal_out$perc[i] = HBV_in_ls[[i]]$hbv_perc[thisBestCal]#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d

			
			#########################################################################################
			# initializing data for model run
		 
			
			mopey_df_ppt = cbind(mopey_df,	
							 SnowGlacier_HBV(model = 1,
											 inputData = cbind(mopey_df$Tavg, mopey_df$PPT),
											 initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
											 param = c(	cal_out$sfcf[i],		#SFCF - snowfall correction factor [-]
														cal_out$tr[i],		#Tr - solid and liquid precipitation threshold temperature [C]
														cal_out$tt[i],		#Tt - melt temperature [C]
														cal_out$fm[i],		#fm - snowmelt factor [mm/C]
														cal_out$fi[i],		#fi - icemelt factor [mm/C]
														cal_out$fic[i])	 	#fic - debris-covered icemelt factor [mm/C]
							 )	)
				
			# soils model
			mopey_df_rech = cbind(mopey_df_ppt,
								  Soil_HBV(
									model = 1,
									inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$PET),
									initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
									param = c(	cal_out$fc[i],			#FC - soil field capacity [mm]
											   cal_out$lp[i],			#LP - parameter ot get actual ET [-]
											   cal_out$beta_soils[i])	#beta - exponential value for nonlinear relations between soil storage and runoff
								  )	)
				
				
			mopey_df_disch = cbind(mopey_df_rech,
								   Routing_HBV(
									 model = 1,	# model=1 gives three stores with routing for each
									 lake = FALSE,
									 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
									 initCond = c(10,10,10),	# initial storage in each reservoir in mm
									 param = c(	cal_out$k0[i],	#KO - top bucket (STZ) storage constant [1/t]
												cal_out$k1[i],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
												cal_out$k2[i],	#K2 - lower bucket (SLZ) storage constant [1/t]
												cal_out$uz1[i],	#UZL - max flux rate between STZ and SUZ [mm/t] 
												cal_out$perc[i])	#PERC - max flux rate between SUZ and SLZ [mm/t]
								   )	)
			
			plot(mopey_df_disch$Date[1000:2000], mopey_df_disch$Qg[1000:2000], type='l')
			mopey_df_disch$QMilAcrFt = (mopey_df_disch$Qg  * this_FNF_area * acr_per_sqml / mm_per_foot) / 1000000
			HBV_out_ls[[i]] = mopey_df_disch
	}
	
	HBV_basin_tot = data.frame(Date = HBV_out_ls[[1]]$Date, QMilAcrFt = 0)
	for(koo in 1:length(HBV_out_ls))	{
		HBV_basin_tot$QMilAcrFt = HBV_basin_tot$QMilAcrFt + HBV_out_ls[[koo]]$QMilAcrFt 
	}

	HBV_avg_out[,(ncol(HBV_avg_out)+1)] = HBV_basin_tot$QMilAcrFt
}
#write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp585_sacramento.csv")
write.csv(HBV_avg_out, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp585_sanjaoquin.csv")


allModsAnn = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsOctMar = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAprJul = data.frame(Date = unique(year(HBV_avg_out$Date)))
allModsAugSep = data.frame(Date = unique(year(HBV_avg_out$Date)))
HBV_avg_out$year = year(HBV_avg_out$Date)
HBV_avg_out$month = month(HBV_avg_out$Date)

for(eachClMod in 2:12)	{	# cols 2-13 have the output from each climate model
		# summing by season and year
	ann_sum = NULL
	octMar_sum = NULL
	aprJul_sum = NULL
	augSep_sum = NULL
	iter = 0
	
	for(k in unique(HBV_avg_out$year))	{
		thisYearSac = subset(HBV_avg_out, year == k & month %in% 1:9)
		lastYearSac = subset(HBV_avg_out, year == k-1 & month %in% 10:12)
		wySac = rbind(lastYearSac, thisYearSac)
		octMarSac = subset(wySac, month %in% c(10:12,1:3))
		aprJulSac = subset(wySac, month %in% c(4:7))
		augSepSac = subset(wySac, month %in% c(8,9))
		iter=iter+1
		ann_sum[iter] = sum(wySac[,eachClMod])
		octMar_sum[iter] = sum(octMarSac[,eachClMod])
		aprJul_sum[iter] = sum(aprJulSac[,eachClMod])
		augSep_sum[iter] = sum(augSepSac[,eachClMod])
	}
	allModsAnn[,eachClMod] = ann_sum
	allModsOctMar[,eachClMod] = octMar_sum
	allModsAprJul[,eachClMod] = aprJul_sum
	allModsAugSep[,eachClMod] = augSep_sum
}

summary(allModsAnn)
summary(allModsOctMar)
summary(allModsAprJul)
summary(allModsAugSep)
















###############################################################################################################################
# merging and melting annual avg data from all models

sanJ126 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp126_sanjaoquin.csv")[,-1]
	colnames(sanJ126) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sanJ126) = c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sanJ126) = c("Date",rep('SanJ-126',11))
	sanJ126$Date = as.Date(sanJ126$Date)
	sum(sanJ126[,2:12]) / (35*11)
sac126 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp126_sacramento.csv")[,-1]
	colnames(sac126) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sac126) = c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sac126) = c("Date",rep('Sac-126',11))
	sac126$Date = as.Date(sanJ126$Date)
	sum(sac126[,2:12]) / (35*11)
sanJ245 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp245_sanjaoquin.csv")[,-1]
	colnames(sanJ245) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sanJ245) =c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sanJ245) = c("Date",rep('SanJ-245',11))
	sanJ245$Date = as.Date(sanJ126$Date)
	sum(sanJ245[,2:12]) / (35*11)
sac245 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp245_sacramento.csv")[,-1]
	colnames(sac245) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sac245) = c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sac245) = c("Date",rep('Sac-245',11))
	sac245$Date = as.Date(sanJ126$Date)
	sum(sac245[,2:12]) / (35*11)
sanJ585 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp585_sanjaoquin.csv")[,-1]
	colnames(sanJ585) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sanJ585) = c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sanJ585) = c("Date",rep('SanJ-585',11))
	sanJ585$Date = as.Date(sanJ126$Date)
	sum(sanJ585[,2:12]) / (35*11)
sac585 = read.csv("C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_ssp585_sacramento.csv")[,-1]
	colnames(sac585) = c("Date","m4",'m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
	colnames(sac585) = c('Date', "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3"  )
	#colnames(sac585) = c("Date",rep('Sac-585',11))
	sac585$Date = as.Date(sanJ126$Date)
	sum(sac585[,2:12]) / (35*11)


##### for just mean vals
allModsList = list(sanJ126,sac126,sanJ245,sac245,sanJ585,sac585)
allModsNames = c('sanJ126','sac126','sanJ245','sac245','sanJ585','sac585')
allModsAnn = data.frame(Date = unique(year(sac585$Date)))
allModsOctMar = data.frame(Date = unique(year(sac585$Date)))
allModsAprJul = data.frame(Date = unique(year(sac585$Date)))
allModsAugSep = data.frame(Date = unique(year(sac585$Date)))

for(eachWtrshdClm in 1:length(allModsList))	{
	HBV_avg_out = allModsList[[eachWtrshdClm]]
	HBV_avg_out$mean = apply(HBV_avg_out[c(2:12)],1,mean) 
	whichIsMeanCol = which(names(HBV_avg_out) == 'mean')
	HBV_avg_out$year = year(HBV_avg_out$Date)
	HBV_avg_out$month = month(HBV_avg_out$Date)

#	for(eachClMod in 2:11)	{	# cols 2-13 have the output from each climate model
			# summing by season and year
	ann_sum = NULL
	octMar_sum = NULL
	aprJul_sum = NULL
	augSep_sum = NULL
	iter = 0
		
	for(k in unique(HBV_avg_out$year))	{
		thisYearSac = subset(HBV_avg_out, year == k & month %in% 1:9)
		lastYearSac = subset(HBV_avg_out, year == k-1 & month %in% 10:12)
		wySac = rbind(lastYearSac, thisYearSac)
		octMarSac = subset(wySac, month %in% c(10:12,1:3))
		aprJulSac = subset(wySac, month %in% c(4:7))
		augSepSac = subset(wySac, month %in% c(8,9))
		iter=iter+1
		ann_sum[iter] = sum(wySac[,whichIsMeanCol])
		octMar_sum[iter] = sum(octMarSac[,whichIsMeanCol])
		aprJul_sum[iter] = sum(aprJulSac[,whichIsMeanCol])
		augSep_sum[iter] = sum(augSepSac[,whichIsMeanCol])
	}
	allModsAnn = cbind(allModsAnn, ann_sum)
	allModsOctMar = cbind(allModsOctMar, octMar_sum)
	allModsAprJul = cbind(allModsAprJul, aprJul_sum)
	allModsAugSep = cbind(allModsAugSep, augSep_sum)
}


names(allModsAnn) = c("Date", allModsNames); allModsAnn = allModsAnn[-1,]
names(allModsOctMar) = c("Date", allModsNames); allModsOctMar = allModsOctMar[-1,]
names(allModsAprJul) = c("Date", allModsNames); allModsAprJul = allModsAprJul[-1,]

write.csv(allModsAnn, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_annavg_sanj-n-sac.csv")
write.csv(allModsOctMar, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_octmar_sanj-n-sac.csv")
write.csv(allModsAprJul, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_aprjul_sanj-n-sac.csv")


plot(allModsAnn$Date,allModsAnn$sanJ126, type='l',col='blue1')
lines(allModsAnn$Date,allModsAnn$sanJ245, col='green3')
lines(allModsAnn$Date,allModsAnn$sanJ585, col='red4')

plot(allModsOctMar$Date,allModsOctMar$sanJ126, type='l',col='blue1')
lines(allModsOctMar$Date,allModsOctMar$sanJ245, col='green3')
lines(allModsOctMar$Date,allModsOctMar$sanJ585, col='red4')

plot(allModsAprJul$Date,allModsAprJul$sanJ126, type='l',col='blue1')
lines(allModsAprJul$Date,allModsAprJul$sanJ245, col='green3')
lines(allModsAprJul$Date,allModsAprJul$sanJ585, col='red4')




write.csv(allModsAnn, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_annavg_sanj-n-sac.csv")
write.csv(allModsOctMar, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_octmar_sanj-n-sac.csv")
write.csv(allModsAprJul, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_aprjul_sanj-n-sac.csv")

valid_data = read.csv("C:\\Users\\arik\\Documents\\ClimateAi\\Wonderful\\valid_dat_SVnSJ.csv")



######################
## for all individ clim model vals
allModsList = list(sanJ126,sac126,sanJ245,sac245,sanJ585,sac585)
allModsNames = c('sanJ126','sac126','sanJ245','sac245','sanJ585','sac585')
allModsAnn = data.frame(Date = unique(year(sac585$Date)))
allModsOctMar = data.frame(Date = unique(year(sac585$Date)))
allModsAprJul = data.frame(Date = unique(year(sac585$Date)))
allModsAugSep = data.frame(Date = unique(year(sac585$Date)))
model_names = c( "EC-Earth3-Veg-LR","GFDL-ESM4", # rmoving models 1:3 and 5 "ACCESS-CM2","AWI-CM-1-1-MR","CMCC-CM2-SR5","FGOALS-g3",
							"INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
							"MRI-ESM2-0","NESM3")

for(eachWtrshdClm in 1:length(allModsList))	{
	HBV_avg_out = allModsList[[eachWtrshdClm]]
	HBV_avg_out$mean = apply(HBV_avg_out[c(2:12)],1,mean) 
	whichIsMeanCol = which(names(HBV_avg_out) == 'mean')
	HBV_avg_out$year = year(HBV_avg_out$Date)
	HBV_avg_out$month = month(HBV_avg_out$Date)

#	for(eachClMod in 2:11)	{	# cols 2-13 have the output from each climate model
			# summing by season and year
	for(ll in model_names)	{
		which_col = which(names(HBV_avg_out) == ll)
		thisSacMod = data.frame(Date = HBV_avg_out$Date, year = HBV_avg_out$year, month = HBV_avg_out$month, streamy = HBV_avg_out[,which_col])

		ann_sum = NULL
		octMar_sum = NULL
		aprJul_sum = NULL
		augSep_sum = NULL
		iter = 0
		
		for(k in unique(HBV_avg_out$year))	{
			thisYearSac = subset(thisSacMod, year == k & month %in% 1:9)
			lastYearSac = subset(thisSacMod, year == k-1 & month %in% 10:12)
			wySac = rbind(lastYearSac, thisYearSac)
			octMarSac = subset(wySac, month %in% c(10:12,1:3))
			aprJulSac = subset(wySac, month %in% c(4:7))
			augSepSac = subset(wySac, month %in% c(8,9))
			
			iter=iter+1
			ann_sum[iter] = sum(wySac[,4])
			octMar_sum[iter] = sum(octMarSac[,4])
			aprJul_sum[iter] = sum(aprJulSac[,4])
			augSep_sum[iter] = sum(augSepSac[,4])
		}
		
		allModsAnn = cbind(allModsAnn, ann_sum)
		allModsOctMar = cbind(allModsOctMar, octMar_sum)
		allModsAprJul = cbind(allModsAprJul, aprJul_sum)
		allModsAugSep = cbind(allModsAugSep, augSep_sum)
	}
}

names(allModsAnn) = c('Date',
	paste0('SanJ-126-',model_names),
	paste0('Sac-126-',model_names),
	paste0('SanJ-245-',model_names),
	paste0('Sac-245-',model_names),
	paste0('SanJ-585-',model_names),
	paste0('Sac-585-',model_names))
allModsAnn = allModsAnn[-1,]

names(allModsOctMar) = c('Date',
	paste0('SanJ-126-',model_names),
	paste0('Sac-126-',model_names),
	paste0('SanJ-245-',model_names),
	paste0('Sac-245-',model_names),
	paste0('SanJ-585-',model_names),
	paste0('Sac-585-',model_names))
allModsOctMar = allModsOctMar[-1,]

names(allModsAprJul) = c('Date',
	paste0('SanJ-126-',model_names),
	paste0('Sac-126-',model_names),
	paste0('SanJ-245-',model_names),
	paste0('Sac-245-',model_names),
	paste0('SanJ-585-',model_names),
	paste0('Sac-585-',model_names))
allModsAprJul = allModsAprJul[-1,]


write.csv(allModsAnn, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_annAllMods_sanj-n-sac.csv")
write.csv(allModsOctMar, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_octmarAllMods_sanj-n-sac.csv")
write.csv(allModsAprJul, "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\future_projs_aprjulAllMods_sanj-n-sac.csv")


plot(allModsAnn$Date,allModsAnn$sanJ126, type='l',col='blue1')
lines(allModsAnn$Date,allModsAnn$sanJ245, col='green3')
lines(allModsAnn$Date,allModsAnn$sanJ585, col='red4')

plot(allModsOctMar$Date,allModsOctMar$sanJ126, type='l',col='blue1')
lines(allModsOctMar$Date,allModsOctMar$sanJ245, col='green3')
lines(allModsOctMar$Date,allModsOctMar$sanJ585, col='red4')

plot(allModsAprJul$Date,allModsAprJul$sanJ126, type='l',col='blue1')
lines(allModsAprJul$Date,allModsAprJul$sanJ245, col='green3')
lines(allModsAprJul$Date,allModsAprJul$sanJ585, col='red4')







meltData = melt(verPerf_keep)
#var_lab_names = c('RS_D20', 'PPT', 'RSR', 'SWE', 'Thaw', 'AI', 'Min T', '%Prm_a', 'Area_a', 'AET_a', '%Glc_a', 'Elev_a', 'AI_a')
ggplot(meltData)	+
	geom_boxplot(aes(x=X2, y=value, fill=X2, colour=X2),
		width=0.65, lwd=1, outlier.shape=NA, notch=FALSE)	+
	ylim(0,1)	+
	scale_fill_manual(values=c('royalblue1','royalblue2','royalblue3','royalblue4')) +#, labels = var_lab_names)	+
#	scale_shape_manual(values=rep(24,7))	+
	scale_color_manual(values=rep('grey10',4)) +#, labels = var_lab_names)	+
	labs(x='', y ='Performance')	+
	theme_classic()
ggsave("C:\\Users\\arik\\Documents\\Postdoc Research\\Arctic_baseflow\\figures\\pred_D20_varPerf.png",
  scale = 0.72,
  dpi = 1000)	








  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##
##
##	workflow for using daymet data and Noah-MP output instead of mopex
##
##
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
library(daymetr)
library(dataRetrieval)
library(EcoHydRology)


	# reading in the watershed boundaries of all North American Hydrobasins
basinAt_NorAm_polys = st_read("C:\\Users\\arik\\Documents\\PhD Research\\D4\\basinAt_NorAm_polys.gpkg")
basinAt_NorAm_centroid = st_centroid(basinAt_NorAm_polys)	# some centroids put us in the ocean, so manually correcting those cases
	bad_centroids = c("8120000180", "8120402550")
	pt1 = st_point(c(-139.499660, 59.776254)); pt2 = st_point(c(-118.11, 64.907))
	long_lats = st_sfc(c(pt1,pt2), crs=4326)
	list_of_basins = which(basinAt_NorAm_polys$HYBAS_ID %in% bad_centroids)
	st_geometry(basinAt_NorAm_centroid)[list_of_basins] = long_lats
basinAt_NorAm_strip = basinAt_NorAm_polys	;	st_geometry(basinAt_NorAm_strip) = NULL
	# reading in the (ever growing) list of all hydrobasins that are upstream of a selected hydrobasin
gaged_upstreams = readRDS("C:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_n.RData")
 
USGS_classif = read.csv("C:\\Users\\arik\\Documents\\Postdoc Research\\Arctic_baseflow\\spreadsheets-in-csv-format\\conterm_bas_classif.txt",
	colClasses = c("STAID" = "character", "AGGECOREGION" = "character")) #   read in the file
USGS_basinid = read.csv("C:\\Users\\arik\\Documents\\Postdoc Research\\Arctic_baseflow\\spreadsheets-in-csv-format\\conterm_basinid.txt",
	colClasses = c("STAID" = "character"))
USGS_gages = merge(USGS_classif, USGS_basinid)
  
  
	# Identifying the gages that Jing is using 
		# manually
#jings_gage_names = c("10258500","11015000","11046300","11111500","11124500","1113850","11143000","11148900","11162570","11173200","11266500","11274500","11381500","11383500","02053200","02081500")
#USGS_gage_names = USGS_gages$STAID
#jings_gages = USGS_gages
#jings_gages = subset(USGS_gages, STAID %in% jings_gage_names)

		# based on the files she has provided
file_loc = "C:\\Users\\arik\\Documents\\LSM Research\\Noah-MP-GrUB\\NC_PPT_n_Rech\\"	
file_names = list.files(file_loc)
jings_gage_names = NULL
for(num_files in 1:length(file_names))	{
	#jings_gage_names = strsplit(file_names, "_")[[1]][2]
	jings_gage_names[num_files] = paste0("0",strsplit(file_names, "_")[[num_files]][2])		# seems she dropped the leading zero, so adding it back here
}
jings_gages = subset(USGS_gages, STAID %in% jings_gage_names)	

	# reading in teh HUC12 watershed boundaries
WBD_all = st_read("C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\WBD_National_GDB.gdb")
WBD_centroids = st_centroid(WBD_all)


	# identifying the HUC12 that contains each gage - note these data are saved below as it takes a while to run this loop
#gage_n_huc = data.frame(USGS_ID = jings_gages$STAID, HUC12 = rep(NA,nrow(jings_gages)))
#for(mopey_row in 5506:nrow(jings_gages))	{
#	lat = jings_gages[mopey_row,"LAT_GAGE"]
#	lon = jings_gages[mopey_row,"LNG_GAGE"]
#	latlon = st_point(c(lon,lat))
#	latlon_sfc = st_sfc(latlon, crs = 4326)
#	mopey_point = st_transform(latlon_sfc,4269)

#	print(mopey_row)
#	gage_n_huc[mopey_row, "HUC12"] = as.character(WBD_all[which(st_intersects(mopey_point, WBD_all, sparse=FALSE)), ]$huc12)
#	write.csv(gage_n_huc, "C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\USGS_outlet_hucs_n.csv",row.names=FALSE)
#}
#write.csv(gage_n_huc, "C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\jings_outlet_hucs_n.csv",row.names=FALSE)
gage_n_huc=read.csv("C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\USGS_outlet_hucs_all.csv", colClasses=c("character","character"))


	
	# identifying all HUC12s that feed / flow into the outlet HUC12
WBD_strip = WBD_all	; st_geometry(WBD_strip) = NULL	; WBD_strip = as.data.frame(WBD_strip)	; rm(WBD_all)


#usgsy_upstreams = list()
#for(new_row in 5505:nrow(gage_n_huc))	{
#	mopey_hucs = gage_n_huc$HUC12[new_row]
#	print(mopey_hucs)
#	WBD_remaining = WBD_strip[-which(WBD_strip$huc12 == mopey_hucs),]
#	print(nrow(WBD_remaining))
		
#	while(any(mopey_hucs %in% WBD_remaining$tohuc))	{
#		for(trueys in which(mopey_hucs %in% WBD_remaining$tohuc))	{
#			new_mopey_huc_rows = which(WBD_remaining$tohuc == mopey_hucs[trueys])
#			new_mopey_hucs = as.character(WBD_remaining$huc12[new_mopey_huc_rows])
#			mopey_hucs = c(mopey_hucs, new_mopey_hucs)
#			print(mopey_hucs)
#			WBD_remaining = WBD_remaining[-new_mopey_huc_rows,]
#			print(nrow(WBD_remaining))
#		}
#	}
#	usgsy_upstreams[[new_row]] = mopey_hucs
#	names(usgsy_upstreams)[new_row] = usgsy_upstreams[[new_row]][1]
#	saveRDS(usgsy_upstreams, file="C:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_jing.RData")
#}
usgsy_upstreams = readRDS("C:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_jing.RData")	
	




	# reading in the empirical data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER")
ksat_model_out = as.data.table(read_feather("C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\all_models_output_v1.feather"))
NHD_atts = as.data.table(read_feather("C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\WS_data\\HUC12s_wNHDatts.feather"))
all_aq_chars = merge(NHD_atts, ksat_model_out, by="HUC_12")
		# old data here: #aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")



# initializing the dataframe for capturing model performance
model_perf = data.frame(
  gage_id = rep(NA,nrow(jings_gages))
)
  


### model 1: setting the parameter ranges for the HBV runoff modules
# snow and glacier module
  sfcf = 1.1#c(runif(500, .2, 1), runif(500, 1, 3)) #1, 2)	#snowfall correction factor [-]
  tr   = -2#runif(1000, -6, 5)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = .4#runif(1000, -5, 6)#0, 3)	#melt temperature [C]
  fm   = 1.5#2.7#c(runif(500, .2, 1.5), (runif(500, 1.5, 8)))#1, 4)	#snowmelt factor [mm/C]
  fi   = 1.5#2.8#c(runif(500, .2, 1.5), (runif(500, 1.5, 10)))#4, 8)	#icemelt factor [mm/C]
  fic  = 5.9#runif(1000, 2, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

# soil module
#  fc   = runif(1000, 30, 3000)#100, 200)
  fc   = 235#c(runif(500, 25, 150), (runif(500, 150, 1200)))
  lp   = 0.6#runif(1000, .2, 1)#0.5, 1)   # parameter to actual ET
  beta_soils = 2.2#runif(1000, 1, 3)

# routing module
  k0   = .5#c(runif(500, .05, .5), (runif(500, .5, 1)))#runif(1000, .01, 1)#0.09, 0.1)
  k1   = .2#c(runif(500, .005, .09), (runif(500, .09, .5)))#runif(1000, .001, .5)#0.05, 0.07)
  k2   = .03#c(runif(500, .0001, .01), (runif(500, .01, .1)))#runif(1000, .0001, .1)#0.05)	
  uz1  = 16#c(runif(500, .22, 10), (runif(500, 10, 40)))#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
#  uz1  = runif(1000, 10, 100)
  perc = 1.5#50#1.5#c(runif(500, .1, .5), (runif(500, .5, 20)))#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in mm/d
#  perc = runif(1000, 10, 100)


	#initializing some generic varaiables
f = 0.1
L_scalar = 2
p = 1
k_fof_D = function(D_thickness, conmod_yint, conmod_slp)	{
	10^(conmod_yint + conmod_slp * D_thickness)
	}
  #68, 77,79,84,98,224,257,278,344,375,393296:343,345:374,376:392,394:nrow(gage_n_huc)
for(jings_ws in 1:nrow(jings_gages))){
	#jings_ws = 9
	jings_HUC = gage_n_huc$HUC12[which(gage_n_huc$USGS_ID == jings_gages$STAID[jings_ws])]
	allws = which(names(usgsy_upstreams) == jings_HUC)
  
  #################################################################################
  # initializing data for our gw module

	allhucs = data.frame()
	iter = 0
	for(this_huc in usgsy_upstreams[[allws]])	{
		iter = iter+1
		if(any(all_aq_chars$HUC_12 %in% this_huc)){
			which_huc = which(all_aq_chars$HUC_12 == this_huc)
			allhucs[iter,"D_dry"] = all_aq_chars$Storage_Dry[which_huc] / f						# mm --> mm
			allhucs[iter,"D_wet"] = all_aq_chars$Storage_Wet[which_huc] / f 					# mm--> mm
			allhucs[iter,"K_dry"] = all_aq_chars$Ksat_Dry[which_huc] * (60*60*24) * 10			# cm/s --> mm/d
			allhucs[iter,"K_wet"] = all_aq_chars$Ksat_Wet[which_huc] * (60*60*24) * 10			# cm/s --> mm/d
			allhucs[iter,"Wo"] =	all_aq_chars$CAT_STREAM_LENGTH[which_huc] * 1000000 * L_scalar	# km --> mm ; also rescaled by L_scalar
			allhucs[iter,"i_slp"] = atan(all_aq_chars$CAT_BASIN_SLOPE[which_huc] / 100)			#
			allhucs[iter,"A"] = 	max(all_aq_chars$CAT_BASIN_AREA[which_huc] * 1000000^2, 1000000^2)				# km2 --> mm2
			allhucs[iter,"B"] =		allhucs[iter,"A"] / (2*allhucs[iter,"Wo"])					#
			allhucs[iter,"sin_slp"] = allhucs[iter,"i_slp"]												#
			allhucs[iter,"strm_crs_sxn"] = allhucs[iter,"Wo"] * 2 / allhucs[iter,"A"]			# accounting for a stream having 2 banks
			allhucs[iter,"B_half"] = allhucs[iter,"B"] / 2										#

			con_mod = lm(c(log10(allhucs[iter,"K_dry"]), log10(allhucs[iter,"K_wet"])) ~
				c(allhucs[iter,"D_dry"]/10, allhucs[iter,"D_wet"]/10))
			allhucs[iter,"conmod_yint"] = con_mod$coef[1]
			allhucs[iter,"conmod_slp"]  = con_mod$coef[2]
		}
	}
	

		# using HUC centroids to download daymet data
	WBD_huc = which(WBD_centroids$huc12 == names(usgsy_upstreams)[allws])
	this_lonlat = st_coordinates(WBD_centroids[WBD_huc,])
	daymet_data = download_daymet(
		lat = this_lonlat[2],
		lon = this_lonlat[1],
		start = 1980,
		end = 2020,
		internal = TRUE)$data
	num_sub_ws = length(usgsy_upstreams[[allws]])	
	if(num_sub_ws > 1)	{
		for(that_huc in 2:num_sub_ws)	{
			n_WBD_huc = which(WBD_centroids$huc12 == usgsy_upstreams[[allws]][that_huc])
			n_this_lonlat = st_coordinates(WBD_centroids[WBD_huc,])
			n_daymet_data = download_daymet(
				lat = this_lonlat[2],
				lon = this_lonlat[1],
				start = 1980,
				end = 2020,
				internal = TRUE)$data
			daymet_data[,-c(1,2)] = daymet_data[,-c(1,2)] + n_daymet_data[,-c(1,2)]
		}
		daymet_data[,-c(1,2)] = daymet_data[,-c(1,2)] / num_sub_ws 
	}
	
	daymet_data$Date = as.Date(paste(daymet_data$year, daymet_data$yday, sep="-"),"%Y-%j")
	daymet_data$T_mean = apply(daymet_data[,c("tmax..deg.c.","tmin..deg.c.")],1,mean)
	daymet_data$E = PET_fromTemp(Jday=daymet_data$yday, Tmax_C=daymet_data$tmax..deg.c., Tmin_C=daymet_data$tmin..deg.c., lat_radians=this_lonlat[2]*pi/180)
	
		# generating different time series for comparing Noah vs DayMet input
	hydro_ts = daymet_data[,c("Date","T_mean","E","prcp..mm.day.")]
	
		# bringing in Noah-MP data
	noahMP_data = fread(paste0(file_loc, file_names[jings_ws]))
	noahMP_data$Date = as.Date(as.character(noahMP_data$Date),format='%Y%m%d')
	

		# ensuring the datasets start at the same time
	use_noah_data = TRUE
	use_noah_recharge = TRUE
	use_noah_dates = TRUE
	use_bypass = TRUE
	bypass_frac = 0.1
	start_date = max(daymet_data$Date[1], noahMP_data$Date[1])
	end_date = min(last(daymet_data$Date), last(noahMP_data$Date))

	if(use_noah_data)	{
		
		start_daymet = which(daymet_data$Date == start_date)
		end_daymet = which(daymet_data$Date == end_date)
		daymet_ts = daymet_data[start_daymet:end_daymet,c("Date","T_mean","E")]
		
		start_noahMP = which(noahMP_data$Date == start_date)
		end_noahMP = which(noahMP_data$Date == end_date)
		noah_ts = noahMP_data[start_noahMP:end_noahMP, c("Precipitation(mm/day)", "Recharge(mm/day)")]
		noah_days = format(noahMP_data$Date, "%d")					#uggghhhh, daymet skipping leap days gives me a headache again... wtf would they even do that? 
		noah_months = format(noahMP_data$Date, "%m")
		noah_leaps = which(noah_months == "02" & noah_days == "29")
		noah_ts_no_leaps = noah_ts[-noah_leaps,]
		hydro_ts = cbind(daymet_ts, noah_ts_no_leaps)
		names(hydro_ts)[4] = "prcp..mm.day."
	}
	
	if(use_noah_dates)	{
		dailyQ = readNWISdv(jings_gages$STAID[jings_ws], '00060', start_date, end_date)
	}	else	{
		dailyQ = readNWISdv(jings_gages$STAID[jings_ws], '00060', "1980-01-01", "2020-12-31")
	}


	names(dailyQ) = c(names(dailyQ)[1:3], "Q_cfs","flag")
	cubic_feet_per_cubic_meters = 35.3147
	seconds_per_day = 60 * 60 * 24
	square_meters_per_sq_km = 1000000
	dailyQ$Q = dailyQ$Q_cfs * 1000 * seconds_per_day / (cubic_feet_per_cubic_meters * jings_gages$DRAIN_SQKM[jings_ws] * square_meters_per_sq_km) 
	setDT(dailyQ)
	setDT(hydro_ts)
	hydro_ts = merge(dailyQ, hydro_ts, by="Date")
	hydro_ts %>%		# filling in missing dates and values
		mutate(Date = as.Date(Date)) %>%
		complete(Date = seq.Date(min(Date), max(Date), by="day")) %>%
		fill(c('E','T_mean'))
		
	
	hbv_ts_input = hydro_ts[hydro_ts, on="Date"]
	na_rows = which(is.na(hbv_ts_input$T_mean))
	hbv_ts_input[na_rows, "T_mean"] = 0
	na_rows = which(is.na(hbv_ts_input$prcp..mm.day.))
	hbv_ts_input[na_rows, "prcp..mm.day."] = 0
	na_rows = which(is.na(hbv_ts_input$E))
	hbv_ts_input[na_rows, "E"] = 0

	
	#########################################################################################
    # initializing data for model run
    n_runs = 1

    cal_out = data.frame(
      sfcf=rep(NA,n_runs),
      tr=rep(NA,n_runs), 
      tt=rep(NA,n_runs),
      fm=rep(NA,n_runs),
      fi=rep(NA,n_runs),
      fic=rep(NA,n_runs),
      fc=rep(NA,n_runs),
      lp=rep(NA,n_runs),
      beta_soils = rep(NA,n_runs),
      k0 = rep(NA,n_runs),
      k1 = rep(NA,n_runs),
      k2 = rep(NA,n_runs),
      uz1 = rep(NA,n_runs),
      perc = rep(NA,n_runs),
      hbv_nse = rep(NA,n_runs),
      new_nse = rep(NA,n_runs)
    )
    
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
        cal_out$sfcf[jj] = sample(sfcf,1)
        cal_out$tr[jj] = sample(tr,1)
        cal_out$tt[jj] = sample(tt,1)
        cal_out$fm[jj] = sample(fm,1)
        cal_out$fi[jj] = sample(fi,1)
        cal_out$fic[jj] = sample(fic,1)
        cal_out$fc[jj] = sample(fc,1)
        cal_out$lp[jj] = sample(lp,1)
        cal_out$beta_soils[jj] = sample(beta_soils,1)
        # since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
        cal_out$k0[jj] = sample(k0,1)
        cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
        cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
        cal_out$uz1[jj] = sample(uz1,1)
        cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)
 

 
        mopey_df_ppt = cbind(hbv_ts_input,	
                             SnowGlacier_HBV(model = 1,
                                             inputData = cbind(hbv_ts_input$T_mean, hbv_ts_input$prcp..mm.day.),
                                             initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                             param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
                                                        cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                        cal_out$tt[jj],		#Tt - melt temperature [C]
                                                        cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
                                                        cal_out$fi[jj],		#fi - icemelt factor [mm/C]
                                                        cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
                             )	)
 
		if(use_bypass)	{
			bypass_vol = rep(0, nrow(mopey_df_ppt))
			signif_pptNmlt = which(mopey_df_ppt$Total >= 10)
			bypass_vol[signif_pptNmlt] = mopey_df_ppt$Total[signif_pptNmlt] * bypass_frac
			mopey_df_ppt$Total = mopey_df_ppt$Total - bypass_vol
		}	else	{
			bypass_vol = 0
		}
 
        # soils model
        mopey_df_rech = cbind(mopey_df_ppt,
                              Soil_HBV(
                                model = 1,
                                inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
                                initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                                param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
                                           cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
                                           cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
                              )	)
        
        
       
		if(use_noah_recharge)	{							# recharge time series
			recharge_data = hydro_ts$Rech
		}	else	{
			recharge_data = mopey_df_rech$Rech
		}
		

	   mopey_df_disch = cbind(mopey_df_rech,
                               Routing_HBV(
                                 model = 1,	# model=1 gives three stores with routing for each
                                 lake = FALSE,
								 inputData = cbind(recharge_data),
                                 initCond = c(10,10,10),	# initial storage in each reservoir in mm
                                 param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
                                            cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                            cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                            cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                            cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                               )	)
        

		
		
        
        ####################################################################
        # new gw module introduced
        # calculating static values so I don't have to recalculate at every loop
        completehucs = allhucs[complete.cases(allhucs),]
		tot_area = sum(completehucs$A)
		completehucs$fracA = completehucs$A / tot_area
		mopey_df_disch$deep_store = 0
        mopey_df_disch$deep_rech = 0
        mopey_df_disch$gw_out = 0
		mopey_df_disch$gw_et_out = 0
		
		for(allthehucs in 1:nrow(completehucs))	{
			rm(gw_out)
			rm(gw_et_out)
			rm(et_out)
			gw_out = NULL
			gw_et_out = NULL
			et_out = NULL
		
			this_huc = allthehucs
			D_dry = completehucs[this_huc,"D_dry"] 
			D_wet = completehucs[this_huc,"D_wet"]
			K_dry = completehucs[this_huc,"K_dry"]
			K_wet = completehucs[this_huc,"K_wet"]
			Wo = 	completehucs[this_huc,"Wo"]
			#i_slp =	completehucs[iter,"i_slp"]
			A = completehucs[this_huc,"A"]
			B = completehucs[this_huc,"B"]
			sin_slp = completehucs[this_huc,"sin_slp"]
			strm_crs_sxn = completehucs[this_huc,"strm_crs_sxn"]
			B_half = completehucs[this_huc,"B_half"] 
			conmod_yint = completehucs[this_huc,"conmod_yint"] 
			conmod_slp = completehucs[this_huc,"conmod_slp"]  
			
			deep_store = D_dry * f # * 1000 #as.numeric(mopey_df_disch$SLZ[1])
			depth = D_dry
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2 + bypass_vol)
			min_ridge = 1 / (B/2)				 # setting the minimum "ridge" i.e. dx in dh/dx to xx meter
			accum_rech = 0
			riparian_fraction = 15000 * strm_crs_sxn 	# x-meters of "riparian" for the entire length of the stream
			evap_gw = riparian_fraction * (mopey_df_rech$E - mopey_df_rech$Eac) 
			n0_frac = 0.95
			gw_r_exp = 3
			deep_rech_mn = mean(deep_rech)
			for(gg in 1:(nrow(mopey_df_disch)))	{
			  accum_rech = min(accum_rech * n0_frac + max(deep_rech[gg] / 1, 0), D_wet)
			  k_arb = as.numeric(k_fof_D(conmod_yint = conmod_yint, conmod_slp = conmod_slp, D_thickness = depth))
			  gw_ridge_calc = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^gw_r_exp, min_ridge)
			  gw_ridge = max(gw_ridge_calc, min_ridge)
			  gw_out_iter =
				#           f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
				k_arb *			# 
				#          (depth) / (B * gw_ridge) *	# dh / dx, without accounting for slope of landscape
				(((sin_slp * B * gw_ridge ) + depth ) / (B * gw_ridge )) *	# dh / dx, while accounting for slope of imperm layer
				#            (sin_slp * sqrt(B_half) + depth) / sqrt(B_half) *	# dh / dx, while accounting for slope of imperm layer
				(depth / 100  + 0*accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
				#(accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
				#           (depth / 4) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
				strm_crs_sxn 
				
			  if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
			  
			  new_deep_store = deep_store[gg] +
				deep_rech[gg] -
				gw_out_iter
			  deep_store = c(deep_store, new_deep_store)
			  depth = new_deep_store / f 
			  
			  gw_out = c(gw_out, gw_out_iter)
			  gw_et_out = c(gw_et_out, max(gw_out_iter - evap_gw[gg], 0)) 
			  
			}
			print(summary(deep_store[-(1:1090)] / (f*1000)))
			print(c(D_dry, D_wet) / 1000)
			print(depth / 1000)
			print(accum_rech)
			print(gw_ridge)
			print(c(K_dry, K_wet))
			print(k_arb)
			print(riparian_fraction)

        these_days = which(!is.na(mopey_df_disch$Q))[1095:1825]
        mopey_df_disch$gw_out_one = gw_out
		mopey_df_disch$gw_et_out_one = gw_et_out

       old_nse = KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q[-(1:1095)])
       new_nse = KGE((mopey_df_disch$Q0[-(1:1095)] + mopey_df_disch$Q1[-(1:1095)] + 	mopey_df_disch$gw_out_one[-(1:1095)]), 
                      mopey_df_disch$Q[-(1:1095)])

		plot(mopey_df_disch$Q[these_days] + .000001, lwd=2, type='l',
             ylim=c(0.001,10),
#             ylim=c(0,max(mopey_df_disch$Q[these_days]/2, na.rm=TRUE)),
#             log='y',
             ylab = "Q (mm)", xlab = "Date",
             #!!! TEMPFIX
             #		main=paste(mopey_catalogue[this_ws,"Name"],
             main=paste(jings_gages[jings_ws, "STAID"],
                        #!!! TEMPFIX	           
                        ": \n old KGE =", round(old_nse,2), "new =", round(new_nse,2)))
  #      lines(mopey_df_disch$Qg[these_days] + .000001, lwd=2, col='red', lty=1)
        lines(mopey_df_disch$Qg[these_days] + .000001, lwd=2, col='red', lty=1)
        lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out_one[these_days] + .000001,
              col='blue1', lwd=2.5)
        lines(mopey_df_disch$gw_out_one[these_days]+ .000001,col='purple1',lwd=2)
        lines(mopey_df_disch$gw_et_out_one[these_days]+ .000001,col='yellow4',lwd=1)

			mopey_df_disch$deep_store = mopey_df_disch$deep_store + deep_store[-length(deep_store)] * completehucs[this_huc,"fracA"]
			mopey_df_disch$deep_rech = mopey_df_disch$deep_rech + deep_rech * completehucs[this_huc,"fracA"]
			mopey_df_disch$gw_out = mopey_df_disch$gw_out + gw_out * completehucs[this_huc,"fracA"]
			mopey_df_disch$gw_et_out = mopey_df_disch$gw_et_out + gw_et_out * completehucs[this_huc,"fracA"]
			
		}

        ##############################################################################
        # plotting and calibrating
 
	Q5_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .04, na.rm=TRUE))
	Q5_old_mae = mae(mopey_df_disch$Qg[Q5_days], mopey_df_disch$Q[Q5_days])
	Q5_new_mae = mae(mopey_df_disch$Q0[Q5_days] + mopey_df_disch$Q1[Q5_days] + 	mopey_df_disch$gw_out[Q5_days], 
                      mopey_df_disch$Q[Q5_days])
    
        
        these_days = which(!is.na(mopey_df_disch$Q))[(1366:2095)+365*2]
  #      old_nse = NSE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q[-(1:1095)])
 #       new_nse = NSE((mopey_df_disch$Q0[-(1:1095)] + mopey_df_disch$Q1[-(1:1095)] + 	mopey_df_disch$gw_out[-(1:1095)]), 
#                    mopey_df_disch$Q[-(1:1095)])
#        old_nse = KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q[-(1:1095)])
#        new_nse = KGE((mopey_df_disch$Q0[-(1:1095)] + mopey_df_disch$Q1[-(1:1095)] + 	mopey_df_disch$gw_out[-(1:1095)]), 
#                      mopey_df_disch$Q[-(1:1095)])
        png(filename=paste0("C:\\Users\\arik\\Documents\\LSM Research\\Noah-MP-GrUB\\",jings_gages$STAID[jings_ws],
			use_noah_data,use_noah_recharge,n0_frac,"exp",gw_r_exp,"bypass",bypass_frac,".png"),width=1200,height=300)
		plot(mopey_df_disch$Date[these_days],mopey_df_disch$Q[these_days], lwd=2, type='l',
             ylim=c(0.01,10),
#             ylim=c(0,max(mopey_df_disch$Q[these_days]/2, na.rm=TRUE)),
             log='y',
             ylab = "Q (mm)", xlab = "Date",
             #!!! TEMPFIX
             #		main=paste(mopey_catalogue[this_ws,"Name"],
             main=paste(jings_gages[jings_ws, "STAID"],
                        #!!! TEMPFIX	           
                        ": \n old Q5 MAE =", round(Q5_old_mae,2), "new =", round(Q5_new_mae,2)))
        lines(mopey_df_disch$Date[these_days], mopey_df_disch$Qg[these_days] + .000001, lwd=2, col='red', lty=1)
  #      lines(mopey_df_disch$Q2[these_days] + .000001, lwd=2, col='red', lty=1)
        lines(mopey_df_disch$Date[these_days], mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days] + .000001,
              col='blue1', lwd=2.5)
        lines(mopey_df_disch$Date[these_days], mopey_df_disch$gw_out[these_days]+ .000001,col='purple1',lwd=2)
        lines(mopey_df_disch$Date[these_days], mopey_df_disch$gw_et_out[these_days]+ .000001,col='yellow4',lwd=1)
		dev.off()
		
        # saving objective function data
        cal_out$hbv_nse[jj] = old_nse
        cal_out$new_nse[jj] = new_nse
      # calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
      #  actual_nonas = which(!is.na(mopey_df_disch$Q))[-(1:1090)]
      #  actual_flow = sum(mopey_df_disch$Q[actual_nonas])
      #  old_flow = sum(mopey_df_disch$Qg[actual_nonas])
      #  new_flow = sum(mopey_df_disch$Q0[actual_nonas]) + sum(mopey_df_disch$Q1[actual_nonas]) + 	sum(mopey_df_disch$gw_out[actual_nonas])
      #  old_nse_p = old_nse - .1 * abs(old_flow - actual_flow) / actual_flow
      #  new_nse_p = new_nse - .1 * abs(new_flow - actual_flow) / actual_flow
      #  cal_out$hbv_nse_p[jj] = old_nse_p
      #  cal_out$new_nse_p[jj] = new_nse_p
        
    }
      
    hbv_best = which.max(cal_out$hbv_nse)
    new_best = which.max(cal_out$new_nse)
    #hbv_best = which.max(cal_out$hbv_nse_p)
    #new_best = which.max(cal_out$new_nse_p)
    
    model_perf$hbv_nse[allws] = cal_out$hbv_nse[hbv_best]
    model_perf$new_nse[allws] = cal_out$new_nse[new_best]
       # hbv modules
    model_perf$hbv_sfcf[allws] = cal_out$sfcf[hbv_best]
    model_perf$hbv_tr[allws] = cal_out$tr[hbv_best]
    model_perf$hbv_tt[allws] = cal_out$tt[hbv_best]
    model_perf$hbv_fm[allws] = cal_out$fm[hbv_best]
    model_perf$hbv_fi[allws] = cal_out$fi[hbv_best]
    model_perf$hbv_fic[allws] = cal_out$fic[hbv_best]
    # soil module
    model_perf$hbv_fc[allws] = cal_out$fc[hbv_best]
    model_perf$hbv_lp[allws] = cal_out$lp[hbv_best]
    model_perf$hbv_beta_soils[allws] = cal_out$beta_soils[hbv_best]
    # routing module
    model_perf$hbv_k0[allws] = cal_out$k0[hbv_best]
    model_perf$hbv_k1[allws] = cal_out$k1[hbv_best]
    model_perf$hbv_k2[allws] = cal_out$k2[hbv_best]
    model_perf$hbv_uz1[allws] = cal_out$uz1[hbv_best]
    model_perf$hbv_perc[allws] = cal_out$perc[hbv_best]
    # hbv modules
    model_perf$new_sfcf[allws] = cal_out$sfcf[new_best]
    model_perf$new_tr[allws] = cal_out$tr[new_best]
    model_perf$new_tt[allws] = cal_out$tt[new_best]
    model_perf$new_fm[allws] = cal_out$fm[new_best]
    model_perf$new_fi[allws] = cal_out$fi[new_best]
    model_perf$new_fic[allws] = cal_out$fic[new_best]
    # soil module
    model_perf$new_fc[allws] = cal_out$fc[new_best]
    model_perf$new_lp[allws] = cal_out$lp[new_best]
    model_perf$new_beta_soils[allws] = cal_out$beta_soils[new_best]
    # routing module
    model_perf$new_k0[allws] = cal_out$k0[new_best]
    model_perf$new_k1[allws] = cal_out$k1[new_best]
    model_perf$new_k2[allws] = cal_out$k2[new_best]
    model_perf$new_uz1[allws] = cal_out$uz1[new_best]
    model_perf$new_perc[allws] = cal_out$perc[new_best]
	
	
	Q5_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .04, na.rm=TRUE))
    Q20_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .20, na.rm=TRUE))
    mean_Q = mean(mopey_df_disch$Q, na.rm=TRUE)
    tnt_days = which(mopey_df_disch$Q < (.30 * mean_Q))  #Tennant method, 30% of mean annual flow
      
    old_Q = mopey_df_disch$Qg
    act_Q = mopey_df_disch$Q + 10^(-5)
    old_nse = NSE(old_Q[-(1:1090)], act_Q[-(1:1090)])
    old_kge = KGE(old_Q[-(1:1090)], act_Q[-(1:1090)])
        #        old_nse_log = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)], FUN=log)
    old_mse = mse(old_Q[-(1:1090)], act_Q[-(1:1090)])
    old_mae = mae(old_Q[-(1:1090)], act_Q[-(1:1090)])
    old_pctB = pbias(old_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
    Q5_old_nse = NSE(old_Q[Q5_days], act_Q[Q5_days])
    Q5_old_kge = KGE(old_Q[Q5_days], act_Q[Q5_days])
    Q5_old_mse = mse(old_Q[Q5_days], act_Q[Q5_days])
    Q5_old_mae = mae(old_Q[Q5_days], act_Q[Q5_days])
    Q5_old_pctB = pbias(old_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days])/ 100
 
    new_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out
    act_Q = mopey_df_disch$Q + 10^(-5)
        
    new_nse = NSE(new_Q[-(1:1090)], act_Q[-(1:1090)])
    new_kge = KGE(new_Q[-(1:1090)], act_Q[-(1:1090)])
    new_mse = mse(new_Q[-(1:1090)], act_Q[-(1:1090)])
    new_mae = mae(new_Q[-(1:1090)], act_Q[-(1:1090)])
    new_pctB = pbias(new_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
    Q5_new_nse = NSE(new_Q[Q5_days], act_Q[Q5_days])
    Q5_new_kge = KGE(new_Q[Q5_days], act_Q[Q5_days])
    Q5_new_mse = mse(new_Q[Q5_days], act_Q[Q5_days])
    Q5_new_mae = mae(new_Q[Q5_days], act_Q[Q5_days])
    Q5_new_pctB = pbias(new_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
 
    model_perf[allws,'hbv_nse'] = old_nse
    model_perf[allws,'hbv_kge'] = old_kge
    #        model_perf[allws,'hbv_log_nse'] = old_nse_log
    model_perf[allws,'hbv_mse'] = old_mse
    model_perf[allws,'hbv_mae'] = old_mae
    model_perf[allws,'hbv_pctB'] = old_pctB
    model_perf[allws,'hbv_Q5_nse'] = Q5_old_nse
    model_perf[allws,'hbv_Q5_kge'] = Q5_old_kge
    model_perf[allws,'hbv_Q5_mse'] = Q5_old_mse
    model_perf[allws,'hbv_Q5_mae'] = Q5_old_mae
    model_perf[allws,'hbv_Q5_pctB'] = Q5_old_pctB

	model_perf[allws,'new_nse'] = new_nse
    model_perf[allws,'new_kge'] = new_kge
    #        model_perf[allws,'new_log_nse'] = new_nse_log
    model_perf[allws,'new_mse'] = new_mse
    model_perf[allws,'new_mae'] = new_mae
    model_perf[allws,'new_pctB'] = new_pctB
    model_perf[allws,'new_Q5_nse'] = Q5_new_nse
    model_perf[allws,'new_Q5_kge'] = Q5_new_kge
    model_perf[allws,'new_Q5_mse'] = Q5_new_mse
    model_perf[allws,'new_Q5_mae'] = Q5_new_mae
    model_perf[allws,'new_Q5_pctB'] = Q5_new_pctB
		


    newE_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_et_out
    act_Q = mopey_df_disch$Q + 10^(-5)
        
    newE_nse = NSE(newE_Q[-(1:1090)], act_Q[-(1:1090)])
    newE_kge = KGE(newE_Q[-(1:1090)], act_Q[-(1:1090)])
    Q5_newE_mae = mae(newE_Q[Q5_days], act_Q[Q5_days])
    Q5_newE_pctB = pbias(newE_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
 
   
	model_perf[allws,'newE_nse'] = newE_nse
    model_perf[allws,'newE_kge'] = newE_kge
   model_perf[allws,'newE_Q5_mae'] = Q5_newE_mae
    model_perf[allws,'newE_Q5_pctB'] = Q5_newE_pctB



  write.csv(model_perf, "C:\\Users\\arik\\Documents\\LSM Research\\GrUB-ER\\test2fast_smallaccumrech.csv")
  print(allws)
}

   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
################################################################
##### step 2: rerunning each hbv_ and new_ optimized model for each watershed to capture model performance
  
  
  for(allws in c(1:74,76:length(ws_ihave))){
    
    #################################################################################
    # initializing data for our gw module
    
    #!!! TEMPFIX
    which_ws = allws
    #this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
    model_perf$gage_id[allws] = ws_ihave[which_ws]
    #!!! TEMPFIX
    
    D_dry = aq_chars$B_depth_dry[this_ws_props]  # in m?
    D_wet = aq_chars$B_depth_wet[this_ws_props] # in m?
    K_dry = aq_chars$B_K2_dry[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    K_wet = aq_chars$B_K2_wet[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
    i_slp =	atan(aq_chars$TOT_BASIN_SLOPE[this_ws_props] / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
    p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
    f = 0.1																# set to 0.1 since that's what I did when originally solved
    A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
    #    B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth
    B = (A / (2*Wo)) #/ cos(i_slp)										# aquifer breadth
    sin_slp = sin(i_slp)
    strm_crs_sxn = 2*Wo/A   # accounting for streams having to banks
    f = 0.1 # the drainable porosoity we set for recession analysis
    B_half = B/2
    
    
#    plot(x=c(D_dry, D_wet), y=c(log10(K_dry), log10(K_wet)))
    con_mod = lm(c(log10(K_dry), log10(K_wet)) ~ 	c(D_dry, D_wet))
#    con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))
    
    # funscion for variable conductivity
    k_fof_D = function(D_thickness) {
#      max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here?? oh yeah, the linear model sometimes goes negative
      10^(con_mod$coef[1] + con_mod$coef[2] * D_thickness)	# now trying a log mod??
    }

    
    
    #########################################################################################
    # initializing data for HBV
    
    #!!! TEMPFIX
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    this_ws = ws_ihave[which_ws]
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
    less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
    write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    names(mopey_df) = c('P','E','Q','T_max','T_min')
    for(aa in c(1,2,4,5)){
      if(any(mopey_df[ , ..aa] < -90))
        mopey_df[which(mopey_df[ , ..aa] < -90), aa] = 0
    }
    if(any(mopey_df$Q < -90)){ mopey_df$Q[which(mopey_df$Q < -90)] = NA}
    
    #!!! TEMPFIX
    mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)
    
    
    n_runs = 3  # on run for each optimized model
    
    
     
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
    
      # snow-glacier model
      mopey_df_ppt = cbind(mopey_df,	
                           SnowGlacier_HBV(model = 1,
                                           inputData = cbind(mopey_df$T_mean, mopey_df$P),
                                           initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                          if(jj != 3){
                                           param = c(	model_perf$hbv_sfcf[allws],		#SFCF - snowfall correction factor [-]
                                                      model_perf$hbv_tr[allws],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                      model_perf$hbv_tt[allws],		#Tt - melt temperature [C]
                                                      model_perf$hbv_fm[allws],		#fm - snowmelt factor [mm/C]
                                                      model_perf$hbv_fi[allws],		#fi - icemelt factor [mm/C]
                                                      model_perf$hbv_fic[allws])	 	#fic - debris-covered icemelt factor [mm/C]
                                          } else  {
                                           param = c( model_perf$new_sfcf[allws],		#SFCF - snowfall correction factor [-]
                                                      model_perf$new_tr[allws],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                      model_perf$new_tt[allws],		#Tt - melt temperature [C]
                                                      model_perf$new_fm[allws],		#fm - snowmelt factor [mm/C]
                                                      model_perf$new_fi[allws],		#fi - icemelt factor [mm/C]
                                                      model_perf$new_fic[allws])	 	#fic - debris-covered icemelt factor [mm/C]
                                          }
                           )	)
      
      # soils model
      mopey_df_rech = cbind(mopey_df_ppt,
                            Soil_HBV(
                              model = 1,
                              inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
                              initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                            if(jj != 3){
                              param = c( model_perf$hbv_fc[allws],			#FC - soil field capacity [mm]
                                         model_perf$hbv_lp[allws],			#LP - parameter ot get actual ET [-]
                                         model_perf$hbv_beta_soils[allws])	#beta - exponential value for nonlinear relations between soil storage and runoff
                            } else {
                              param = c( model_perf$new_fc[allws],			#FC - soil field capacity [mm]
                                         model_perf$new_lp[allws],			#LP - parameter ot get actual ET [-]
                                         model_perf$new_beta_soils[allws])	#beta - exponential value for nonlinear relations between soil storage and runoff
                              
                            }
                            )	)
      
      
      mopey_df_disch = cbind(mopey_df_rech,
                             Routing_HBV(
                               model = 1,	# model=1 gives three stores with routing for each
                               lake = FALSE,
                               inputData = cbind(mopey_df_rech$Rech),	# recharge time series
                               initCond = c(10,10,10),	# initial storage in each reservoir in mm
                             if(jj != 3){
                               param = c(	model_perf$hbv_k0[allws],	#KO - top bucket (STZ) storage constant [1/t]
                                          model_perf$hbv_k1[allws],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                          model_perf$hbv_k2[allws],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$hbv_uz1[allws],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                          model_perf$hbv_perc[allws])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                             } else {
                               param = c(	model_perf$new_k0[allws],	#KO - top bucket (STZ) storage constant [1/t]
                                          model_perf$new_k1[allws],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                         # ifelse(model_perf$hbv_k2[allws] < model_perf$new_k1[allws], model_perf$hbv_k2[allws], model_perf$new_k1[allws]-.0001),	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$new_k2[allws],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$new_uz1[allws],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                          model_perf$new_perc[allws])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                             }
                             )	)
      
      
      ####################################################################
      # new gw module introduced
      # calculating static values so I don't have to recalculate at every loop
      gw_out = NULL
      deep_store = D_dry * f * 1000 #as.numeric(mopey_df_disch$SLZ[1])
      depth = D_dry
      deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
      min_ridge = 1 / (B_half/2) # setting the minimum "ridge" i.e. dx in dh/dx to xx meter
      #        recent_rech = seq(1/182,1,length=182) / 100 # * drainable porosity and converted to meters
      accum_rech = 0
      riparian_fraction = 0 * (strm_crs_sxn * 2)
      evap_gw = riparian_fraction * (mopey_df_rech$E - mopey_df_rech$Eac) 
      for(gg in 1:(nrow(mopey_df_disch)))	{
        #          if(gg > length(recent_rech))	{
        #            accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) 
        #          }
        accum_rech = min(accum_rech * 0.95 + max(deep_rech[gg] / (f * 1000), 0), D_wet - depth / 100)
        k_arb = as.numeric(k_fof_D(D_thickness = depth))
        gw_ridge_calc = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^3, min_ridge)
        gw_ridge = max(gw_ridge_calc, min_ridge)
        gw_out_iter =
          #           f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
          k_arb *			# 
          #          (depth) / (B * gw_ridge) *	# dh / dx, without accounting for slope of landscape
          (((sin_slp * B_half * gw_ridge / 2) + depth ) / (B_half * gw_ridge / 2)) *	# dh / dx, while accounting for slope of imperm layer
          #            (sin_slp * sqrt(B_half) + depth) / sqrt(B_half) *	# dh / dx, while accounting for slope of imperm layer
          (depth / 100  + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
          #           (depth / 4) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
          strm_crs_sxn	* 2 *
          1000                    # converting back to mm 
        if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
        
        new_deep_store = deep_store[gg] +
          deep_rech[gg] -
          gw_out_iter
        deep_store = c(deep_store, new_deep_store)
        depth = new_deep_store / (f * 1000) # converting back to m and normalizing for f
        
        gw_stream = max(gw_out_iter - evap_gw[gg], 0) 
        gw_out = c(gw_out, gw_stream)
        
      }
      print(summary(deep_store[-(1:1090)] / (f * 1000)))
      print(c(D_dry, D_wet))
      print(accum_rech)
      print(gw_ridge)
      print(k_arb)
      print(c(K_dry, K_wet))
      print(riparian_fraction)
      if(jj == 2){
        mopey_df_disch$deep_store = deep_store[-length(deep_store)]
        mopey_df_disch$deep_rech = deep_rech
        mopey_df_disch$gw_out = gw_out
      }
      if(jj == 3){
        mopey_df_disch$deep_store_wCalib = deep_store[-length(deep_store)]
        mopey_df_disch$deep_rech_wCalib = deep_rech
        mopey_df_disch$gw_out_wCalib = gw_out
      }
      
      ##############################################################################
      # plotting and calibrating
      
      
      Q5_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .04, na.rm=TRUE))
      Q10_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .10, na.rm=TRUE))
      Q20_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .20, na.rm=TRUE))
      mean_Q = mean(mopey_df_disch$Q, na.rm=TRUE)
      tnt_days = which(mopey_df_disch$Q < (.30 * mean_Q))  #Tennant method, 30% of mean annual flow
      
      if(jj == 1){
        old_Q = mopey_df_disch$Qg
        act_Q = mopey_df_disch$Q + 10^(-5)
        old_nse = NSE(old_Q[-(1:1090)], act_Q[-(1:1090)])
        old_kge = KGE(old_Q[-(1:1090)], act_Q[-(1:1090)])
        #        old_nse_log = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)], FUN=log)
        old_mse = mse(old_Q[-(1:1090)], act_Q[-(1:1090)])
        old_mae = mae(old_Q[-(1:1090)], act_Q[-(1:1090)])
        old_pctB = pbias(old_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
        Q5_old_nse = NSE(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_kge = KGE(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_mse = mse(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_mae = mae(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_pctB = pbias(old_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days])/ 100
        Q10_old_nse = NSE(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_kge = KGE(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_mse = mse(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_mae = mae(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_pctB = pbias(old_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
        Q20_old_nse = NSE(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_kge = KGE(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_mse = mse(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_mae = mae(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_pctB = pbias(old_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
        tnt_old_nse = NSE(old_Q[tnt_days], act_Q[tnt_days])
        tnt_old_kge = KGE(old_Q[tnt_days], act_Q[tnt_days])
        tnt_old_mse = mse(old_Q[tnt_days], act_Q[tnt_days])
        tnt_old_mae = mae(old_Q[tnt_days], act_Q[tnt_days])
        tnt_old_pctB = pbias(old_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
        
        model_perf[allws,'hbv_nse'] = old_nse
        model_perf[allws,'hbv_kge'] = old_kge
        #        model_perf[allws,'hbv_log_nse'] = old_nse_log
        model_perf[allws,'hbv_mse'] = old_mse
        model_perf[allws,'hbv_mae'] = old_mae
        model_perf[allws,'hbv_pctB'] = old_pctB
        model_perf[allws,'hbv_Q5_nse'] = Q5_old_nse
        model_perf[allws,'hbv_Q5_kge'] = Q5_old_kge
        model_perf[allws,'hbv_Q5_mse'] = Q5_old_mse
        model_perf[allws,'hbv_Q5_mae'] = Q5_old_mae
        model_perf[allws,'hbv_Q5_pctB'] = Q5_old_pctB
        model_perf[allws,'hbv_Q10_nse'] = Q10_old_nse
        model_perf[allws,'hbv_Q10_kge'] = Q10_old_kge
        model_perf[allws,'hbv_Q10_mse'] = Q10_old_mse
        model_perf[allws,'hbv_Q10_mae'] = Q10_old_mae
        model_perf[allws,'hbv_Q10_pctB'] = Q10_old_pctB
        model_perf[allws,'hbv_Q20_nse'] = Q20_old_nse
        model_perf[allws,'hbv_Q20_kge'] = Q20_old_kge
        model_perf[allws,'hbv_Q20_mse'] = Q20_old_mse
        model_perf[allws,'hbv_Q20_mae'] = Q20_old_mae
        model_perf[allws,'hbv_Q20_pctB'] = Q20_old_pctB
        model_perf[allws,'hbv_tnt_nse'] = tnt_old_nse
        model_perf[allws,'hbv_tnt_kge'] = tnt_old_kge
        model_perf[allws,'hbv_tnt_mse'] = tnt_old_mse
        model_perf[allws,'hbv_tnt_mae'] = tnt_old_mae
        model_perf[allws,'hbv_tnt_pctB'] = tnt_old_pctB
      } 
      if(jj == 2) {
        new_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out
        act_Q = mopey_df_disch$Q + 10^(-5)
        
        new_nse = NSE(new_Q[-(1:1090)], act_Q[-(1:1090)])
        new_kge = KGE(new_Q[-(1:1090)], act_Q[-(1:1090)])
        #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
        #                      mopey_df_disch$Q[-(1:365)], FUN=log)
        new_mse = mse(new_Q[-(1:1090)], act_Q[-(1:1090)])
        new_mae = mae(new_Q[-(1:1090)], act_Q[-(1:1090)])
        new_pctB = pbias(new_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
        Q5_new_nse = NSE(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_kge = KGE(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_mse = mse(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_mae = mae(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_pctB = pbias(new_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
        Q10_new_nse = NSE(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_kge = KGE(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_mse = mse(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_mae = mae(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_pctB = pbias(new_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
        Q20_new_nse = NSE(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_kge = KGE(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_mse = mse(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_mae = mae(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_pctB = pbias(new_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
        tnt_new_nse = NSE(new_Q[tnt_days], act_Q[tnt_days])
        tnt_new_kge = KGE(new_Q[tnt_days], act_Q[tnt_days])
        tnt_new_mse = mse(new_Q[tnt_days], act_Q[tnt_days])
        tnt_new_mae = mae(new_Q[tnt_days], act_Q[tnt_days])
        tnt_new_pctB = pbias(new_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
        
        model_perf[allws,'new_nse'] = new_nse
        model_perf[allws,'new_kge'] = new_kge
        #        model_perf[allws,'new_log_nse'] = new_nse_log
        model_perf[allws,'new_mse'] = new_mse
        model_perf[allws,'new_mae'] = new_mae
        model_perf[allws,'new_pctB'] = new_pctB
        model_perf[allws,'new_Q5_nse'] = Q5_new_nse
        model_perf[allws,'new_Q5_kge'] = Q5_new_kge
        model_perf[allws,'new_Q5_mse'] = Q5_new_mse
        model_perf[allws,'new_Q5_mae'] = Q5_new_mae
        model_perf[allws,'new_Q5_pctB'] = Q5_new_pctB
        model_perf[allws,'new_Q10_nse'] = Q10_new_nse
        model_perf[allws,'new_Q10_kge'] = Q10_new_kge
        model_perf[allws,'new_Q10_mse'] = Q10_new_mse
        model_perf[allws,'new_Q10_mae'] = Q10_new_mae
        model_perf[allws,'new_Q10_pctB'] = Q10_new_pctB
        model_perf[allws,'new_Q20_nse'] = Q20_new_nse
        model_perf[allws,'new_Q20_kge'] = Q20_new_kge
        model_perf[allws,'new_Q20_mse'] = Q20_new_mse
        model_perf[allws,'new_Q20_mae'] = Q20_new_mae
        model_perf[allws,'new_Q20_pctB'] = Q20_new_pctB
        model_perf[allws,'new_tnt_nse'] = tnt_new_nse
        model_perf[allws,'new_tnt_kge'] = tnt_new_kge
        model_perf[allws,'new_tnt_mse'] = tnt_new_mse
        model_perf[allws,'new_tnt_mae'] = tnt_new_mae
        model_perf[allws,'new_tnt_pctB'] = tnt_new_pctB
      }
      if(jj == 3) {
        newC_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
        act_Q = mopey_df_disch$Q + 10^(-5)
        
        newC_nse = NSE(newC_Q[-(1:1090)], act_Q[-(1:1090)])
        newC_kge = KGE(newC_Q[-(1:1090)], act_Q[-(1:1090)])
        #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
        #                      mopey_df_disch$Q[-(1:365)], FUN=log)
        newC_mse = mse(newC_Q[-(1:1090)], act_Q[-(1:1090)])
        newC_mae = mae(newC_Q[-(1:1090)], act_Q[-(1:1090)])
        newC_pctB = pbias(newC_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
        Q5_newC_nse = NSE(newC_Q[Q5_days], act_Q[Q5_days])
        Q5_newC_kge = KGE(newC_Q[Q5_days], act_Q[Q5_days])
        Q5_newC_mse = mse(newC_Q[Q5_days], act_Q[Q5_days])
        Q5_newC_mae = mae(newC_Q[Q5_days], act_Q[Q5_days])
        Q5_newC_pctB = pbias(newC_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
        Q10_newC_nse = NSE(newC_Q[Q10_days], act_Q[Q10_days])
        Q10_newC_kge = KGE(newC_Q[Q10_days], act_Q[Q10_days])
        Q10_newC_mse = mse(newC_Q[Q10_days], act_Q[Q10_days])
        Q10_newC_mae = mae(newC_Q[Q10_days], act_Q[Q10_days])
        Q10_newC_pctB = pbias(newC_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
        Q20_newC_nse = NSE(newC_Q[Q20_days], act_Q[Q20_days])
        Q20_newC_kge = KGE(newC_Q[Q20_days], act_Q[Q20_days])
        Q20_newC_mse = mse(newC_Q[Q20_days], act_Q[Q20_days])
        Q20_newC_mae = mae(newC_Q[Q20_days], act_Q[Q20_days])
        Q20_newC_pctB = pbias(newC_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
        tnt_newC_nse = NSE(newC_Q[tnt_days], act_Q[tnt_days])
        tnt_newC_kge = KGE(newC_Q[tnt_days], act_Q[tnt_days])
        tnt_newC_mse = mse(newC_Q[tnt_days], act_Q[tnt_days])
        tnt_newC_mae = mae(newC_Q[tnt_days], act_Q[tnt_days])
        tnt_newC_pctB = pbias(newC_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
        
        model_perf[allws,'newC_nse'] = newC_nse
        model_perf[allws,'newC_kge'] = newC_kge
        #        model_perf[allws,'newC_log_nse'] = newC_nse_log
        model_perf[allws,'newC_mse'] = newC_mse
        model_perf[allws,'newC_mae'] = newC_mae
        model_perf[allws,'newC_pctB'] = newC_pctB
        model_perf[allws,'newC_Q5_nse'] = Q5_newC_nse
        model_perf[allws,'newC_Q5_kge'] = Q5_newC_kge
        model_perf[allws,'newC_Q5_mse'] = Q5_newC_mse
        model_perf[allws,'newC_Q5_mae'] = Q5_newC_mae
        model_perf[allws,'newC_Q5_pctB'] = Q5_newC_pctB
        model_perf[allws,'newC_Q10_nse'] = Q10_newC_nse
        model_perf[allws,'newC_Q10_kge'] = Q10_newC_kge
        model_perf[allws,'newC_Q10_mse'] = Q10_newC_mse
        model_perf[allws,'newC_Q10_mae'] = Q10_newC_mae
        model_perf[allws,'newC_Q10_pctB'] = Q10_newC_pctB
        model_perf[allws,'newC_Q20_nse'] = Q20_newC_nse
        model_perf[allws,'newC_Q20_kge'] = Q20_newC_kge
        model_perf[allws,'newC_Q20_mse'] = Q20_newC_mse
        model_perf[allws,'newC_Q20_mae'] = Q20_newC_mae
        model_perf[allws,'newC_Q20_pctB'] = Q20_newC_pctB
        model_perf[allws,'newC_tnt_nse'] = tnt_newC_nse
        model_perf[allws,'newC_tnt_kge'] = tnt_newC_kge
        model_perf[allws,'newC_tnt_mse'] = tnt_newC_mse
        model_perf[allws,'newC_tnt_mae'] = tnt_newC_mae
        model_perf[allws,'newC_tnt_pctB'] = tnt_newC_pctB
      }
      
      these_days = which(!is.na(mopey_df_disch$Q))[2500:2865] # for illustrative plotting
      if(jj==1) {
      plot(mopey_df_disch$Q[these_days], type='l',
           ylim=c(0.01,max(mopey_df_disch$Q[these_days]+8, na.rm=TRUE)),
          log='y',
           col='grey5', lwd=3,
           ylab = "Q (mm)", xlab = "Date",
           #!!! TEMPFIX
           #		main=paste(mopey_catalogue[this_ws,"Name"],
           main=paste(this_ws,
                      #!!! TEMPFIX	           
                      ": \n old KGE =", round(old_kge,2), "new =", round(new_kge,2), "new w.Cal=", round(newC_kge,2)))
      lines(mopey_df_disch$Qg[these_days], lwd=3, col='tomato', lty=1)
      }
      if(jj==2){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
            col='purple3', lwd=3, lty=1)
      }
      if(jj==3){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out_wCalib[these_days],
            col='royalblue1', lwd=3)
      }      
  }
}

  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
basic = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b.csv")
switchy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b_switch.csv")


newy = model_perf

complexy = fread("C:/Users/arik/Documents/LSM Research/model_perf_2apr2021_all.csv")  
newy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_lin_2x_newrech.csv")
oldy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b_switch.csv")


complexy$version = "v.a1"
oldy$version = "v.c1"
newy$version = "v.c2"


allmods = rbind(complexy[,-c(1)], oldy, newy)
  
  # NSE boxplots
boxplot(hbv_nse ~ version, data = allmods,
          boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "NSE")
boxplot(new_nse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  # MAE Q10
boxplot(hbv_Q10_mae ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q10 MAE (mm/d)")
boxplot(new_Q10_mae ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# RMSE Q10
boxplot(hbv_Q10_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q10 RMSE (mm/d)")
boxplot(new_Q10_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# MAE Q5
boxplot(hbv_Q5_mae ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 MAE (mm/d)")
boxplot(new_Q5_mae ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# RMSE Q5
boxplot(hbv_Q5_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 RMSE (mm/d)")
boxplot(new_Q5_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  
boxplot(hbv_Q5_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 RMSE (mm/d)")
boxplot(new_Q5_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')

  
  
boxplot(newy$hbv_Q20_pctB, data = newy,
        boxwex=0.25, at=1 - 0.2,
        xlim = c(0.5, 1.5), ylim = c(-50, 10), yaxs = "i",
        col='tomato',
        ylab = "% bias")
boxplot(newy$new_Q20_pctB , data = newy, add=TRUE,
        boxwex=0.25, at=1 + 0.2,
        col='royalblue3')
abline(h=0)
legend(.5,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  
  ### data to save for the final run
  

  
  









































################################################################
##### step 3: sensitivity analysis of each parameter (w, p, and v) that is estimated without data


for(allws in c(1:74,76:length(ws_ihave))){
  
  #################################################################################
  # initializing data for our gw module
  
  #!!! TEMPFIX
  which_ws = allws
  #this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
  #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
  this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
  model_perf$gage_id[allws] = ws_ihave[which_ws]
  #!!! TEMPFIX
  
  D_dry = aq_chars$B_depth_dry[this_ws_props]  # in m?
  D_wet = aq_chars$B_depth_wet[this_ws_props] # in m?
  K_dry = aq_chars$B_K2_dry[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
  K_wet = aq_chars$B_K2_wet[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
  Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
  i_slp =	atan(aq_chars$TOT_BASIN_SLOPE[this_ws_props] / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
  p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
  f = 0.1																# set to 0.1 since that's what I did when originally solved
  A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
  #    B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth
  B = (A / (2*Wo)) #/ cos(i_slp)										# aquifer breadth
  sin_slp = sin(i_slp)
  strm_crs_sxn = 2*Wo/A   # accounting for streams having to banks
  f = 0.1 # the drainable porosoity we set for recession analysis
  B_half = B/2
  
  
  #    plot(x=c(D_dry, D_wet), y=c(log10(K_dry), log10(K_wet)))
  con_mod = lm(c(log10(K_dry), log10(K_wet)) ~ 	c(D_dry, D_wet))
  #    con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))
  
  # funscion for variable conductivity
  k_fof_D = function(D_thickness) {
    #      max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here?? oh yeah, the linear model sometimes goes negative
    10^(con_mod$coef[1] + con_mod$coef[2] * D_thickness)	# now trying a log mod??
  }
  
  
  
  #########################################################################################
  # initializing data for HBV
  
  #!!! TEMPFIX
  #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
  this_ws = ws_ihave[which_ws]
  #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
  #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
  idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
  less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
  write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
  mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
  names(mopey_df) = c('P','E','Q','T_max','T_min')
  for(aa in c(1,2,4,5)){
    if(any(mopey_df[ , ..aa] < -90))
      mopey_df[which(mopey_df[ , ..aa] < -90), aa] = 0
  }
  if(any(mopey_df$Q < -90)){ mopey_df$Q[which(mopey_df$Q < -90)] = NA}
  
  #!!! TEMPFIX
  mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)
  
  
  n_runs = 5  # on run for each optimized model
  
  
  
   
    # snow-glacier model
    mopey_df_ppt = cbind(mopey_df,	
                         SnowGlacier_HBV(model = 1,
                                         inputData = cbind(mopey_df$T_mean, mopey_df$P),
                                         initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                           param = c(	model_perf$hbv_sfcf[allws],		#SFCF - snowfall correction factor [-]
                                                      model_perf$hbv_tr[allws],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                      model_perf$hbv_tt[allws],		#Tt - melt temperature [C]
                                                      model_perf$hbv_fm[allws],		#fm - snowmelt factor [mm/C]
                                                      model_perf$hbv_fi[allws],		#fi - icemelt factor [mm/C]
                                                      model_perf$hbv_fic[allws])	 	#fic - debris-covered icemelt factor [mm/C]
                          )	)
    
    # soils model
    mopey_df_rech = cbind(mopey_df_ppt,
                          Soil_HBV(
                            model = 1,
                            inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
                            initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                              param = c( model_perf$hbv_fc[allws],			#FC - soil field capacity [mm]
                                         model_perf$hbv_lp[allws],			#LP - parameter ot get actual ET [-]
                                         model_perf$hbv_beta_soils[allws])	#beta - exponential value for nonlinear relations between soil storage and runoff
                          )	)
    
    
    mopey_df_disch = cbind(mopey_df_rech,
                           Routing_HBV(
                             model = 1,	# model=1 gives three stores with routing for each
                             lake = FALSE,
                             inputData = cbind(mopey_df_rech$Rech),	# recharge time series
                             initCond = c(10,10,10),	# initial storage in each reservoir in mm
                               param = c(	model_perf$hbv_k0[allws],	#KO - top bucket (STZ) storage constant [1/t]
                                          model_perf$hbv_k1[allws],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                          model_perf$hbv_k2[allws],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$hbv_uz1[allws],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                          model_perf$hbv_perc[allws])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                           )	)
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
      
    
    ####################################################################
    # new gw module introduced
    # calculating static values so I don't have to recalculate at every loop
    param_w = 3     #c(1,2,3,5,10)[jj]
    param_p = 100   #c(2,10,100,200,1000)[jj]
    param_v = c(0.1,0.52,0.95,0.963,0.9875)[jj]
    
    gw_out = NULL
    deep_store = D_dry * f * 1000 #as.numeric(mopey_df_disch$SLZ[1])
    depth = D_dry
    deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
    min_ridge = 1 / (B_half/2) # setting the minimum "ridge" i.e. dx in dh/dx to xx meter
    #        recent_rech = seq(1/182,1,length=182) / 100 # * drainable porosity and converted to meters
    accum_rech = 0
    riparian_fraction = 0 * (strm_crs_sxn * 2)
    evap_gw = riparian_fraction * (mopey_df_rech$E - mopey_df_rech$Eac) 
    for(gg in 1:(nrow(mopey_df_disch)))	{
      #          if(gg > length(recent_rech))	{
      #            accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) 
      #          }
      accum_rech = min(accum_rech * param_v + max(deep_rech[gg] / (f * 1000), 0), D_wet - depth / 100)
      k_arb = as.numeric(k_fof_D(D_thickness = depth))
      gw_ridge_calc = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^param_w, min_ridge)
      gw_ridge = max(gw_ridge_calc, min_ridge)
      gw_out_iter =
        #           f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
        k_arb *			# 
        #          (depth) / (B * gw_ridge) *	# dh / dx, without accounting for slope of landscape
        (((sin_slp * B_half * gw_ridge / 2) + depth ) / (B_half * gw_ridge / 2)) *	# dh / dx, while accounting for slope of imperm layer
        #            (sin_slp * sqrt(B_half) + depth) / sqrt(B_half) *	# dh / dx, while accounting for slope of imperm layer
        (depth / param_p  + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
        #           (depth / 4) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
        strm_crs_sxn	* 2 *
        1000                    # converting back to mm 
      if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
      
      new_deep_store = deep_store[gg] +
        deep_rech[gg] -
        gw_out_iter
      deep_store = c(deep_store, new_deep_store)
      depth = new_deep_store / (f * 1000) # converting back to m and normalizing for f
      
      gw_stream = max(gw_out_iter - evap_gw[gg], 0) 
      gw_out = c(gw_out, gw_stream)
      
    }
    print(summary(deep_store[-(1:1090)] / (f * 1000)))
    print(c(D_dry, D_wet))
    print(accum_rech)
    print(gw_ridge)
    print(k_arb)
    print(c(K_dry, K_wet))
    print(riparian_fraction)

    mopey_df_disch$deep_store = deep_store[-length(deep_store)]
    mopey_df_disch$deep_rech = deep_rech
    mopey_df_disch$gw_out = gw_out
      mopey_df_disch$deep_store_wCalib = deep_store[-length(deep_store)]
      mopey_df_disch$deep_rech_wCalib = deep_rech
      mopey_df_disch$gw_out_wCalib = gw_out
    
    ##############################################################################
    # plotting and calibrating
    
    
    Q5_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .04, na.rm=TRUE))
    Q10_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .10, na.rm=TRUE))
    Q20_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .20, na.rm=TRUE))
    mean_Q = mean(mopey_df_disch$Q, na.rm=TRUE)
    tnt_days = which(mopey_df_disch$Q < (.30 * mean_Q))  #Tennant method, 30% of mean annual flow
    
    if(jj == 1) {
      param1_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
      act_Q = mopey_df_disch$Q + 10^(-5)
      
      param1_nse = NSE(param1_Q[-(1:1090)], act_Q[-(1:1090)])
      param1_kge = KGE(param1_Q[-(1:1090)], act_Q[-(1:1090)])
      #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
      #                      mopey_df_disch$Q[-(1:365)], FUN=log)
      param1_mse = mse(param1_Q[-(1:1090)], act_Q[-(1:1090)])
      param1_mae = mae(param1_Q[-(1:1090)], act_Q[-(1:1090)])
      param1_pctB = pbias(param1_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
      Q5_param1_nse = NSE(param1_Q[Q5_days], act_Q[Q5_days])
      Q5_param1_kge = KGE(param1_Q[Q5_days], act_Q[Q5_days])
      Q5_param1_mse = mse(param1_Q[Q5_days], act_Q[Q5_days])
      Q5_param1_mae = mae(param1_Q[Q5_days], act_Q[Q5_days])
      Q5_param1_pctB = pbias(param1_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
      Q10_param1_nse = NSE(param1_Q[Q10_days], act_Q[Q10_days])
      Q10_param1_kge = KGE(param1_Q[Q10_days], act_Q[Q10_days])
      Q10_param1_mse = mse(param1_Q[Q10_days], act_Q[Q10_days])
      Q10_param1_mae = mae(param1_Q[Q10_days], act_Q[Q10_days])
      Q10_param1_pctB = pbias(param1_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
      Q20_param1_nse = NSE(param1_Q[Q20_days], act_Q[Q20_days])
      Q20_param1_kge = KGE(param1_Q[Q20_days], act_Q[Q20_days])
      Q20_param1_mse = mse(param1_Q[Q20_days], act_Q[Q20_days])
      Q20_param1_mae = mae(param1_Q[Q20_days], act_Q[Q20_days])
      Q20_param1_pctB = pbias(param1_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
      tnt_param1_nse = NSE(param1_Q[tnt_days], act_Q[tnt_days])
      tnt_param1_kge = KGE(param1_Q[tnt_days], act_Q[tnt_days])
      tnt_param1_mse = mse(param1_Q[tnt_days], act_Q[tnt_days])
      tnt_param1_mae = mae(param1_Q[tnt_days], act_Q[tnt_days])
      tnt_param1_pctB = pbias(param1_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
      
      model_perf[allws,'param1_nse'] = param1_nse
      model_perf[allws,'param1_kge'] = param1_kge
      #        model_perf[allws,'param1_log_nse'] = param1_nse_log
      model_perf[allws,'param1_mse'] = param1_mse
      model_perf[allws,'param1_mae'] = param1_mae
      model_perf[allws,'param1_pctB'] = param1_pctB
      model_perf[allws,'param1_Q5_nse'] = Q5_param1_nse
      model_perf[allws,'param1_Q5_kge'] = Q5_param1_kge
      model_perf[allws,'param1_Q5_mse'] = Q5_param1_mse
      model_perf[allws,'param1_Q5_mae'] = Q5_param1_mae
      model_perf[allws,'param1_Q5_pctB'] = Q5_param1_pctB
      model_perf[allws,'param1_Q10_nse'] = Q10_param1_nse
      model_perf[allws,'param1_Q10_kge'] = Q10_param1_kge
      model_perf[allws,'param1_Q10_mse'] = Q10_param1_mse
      model_perf[allws,'param1_Q10_mae'] = Q10_param1_mae
      model_perf[allws,'param1_Q10_pctB'] = Q10_param1_pctB
      model_perf[allws,'param1_Q20_nse'] = Q20_param1_nse
      model_perf[allws,'param1_Q20_kge'] = Q20_param1_kge
      model_perf[allws,'param1_Q20_mse'] = Q20_param1_mse
      model_perf[allws,'param1_Q20_mae'] = Q20_param1_mae
      model_perf[allws,'param1_Q20_pctB'] = Q20_param1_pctB
      model_perf[allws,'param1_tnt_nse'] = tnt_param1_nse
      model_perf[allws,'param1_tnt_kge'] = tnt_param1_kge
      model_perf[allws,'param1_tnt_mse'] = tnt_param1_mse
      model_perf[allws,'param1_tnt_mae'] = tnt_param1_mae
      model_perf[allws,'param1_tnt_pctB'] = tnt_param1_pctB
    }
    if(jj == 2) {
      param2_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
      act_Q = mopey_df_disch$Q + 10^(-5)
      
      param2_nse = NSE(param2_Q[-(1:1090)], act_Q[-(1:1090)])
      param2_kge = KGE(param2_Q[-(1:1090)], act_Q[-(1:1090)])
      #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
      #                      mopey_df_disch$Q[-(1:365)], FUN=log)
      param2_mse = mse(param2_Q[-(1:1090)], act_Q[-(1:1090)])
      param2_mae = mae(param2_Q[-(1:1090)], act_Q[-(1:1090)])
      param2_pctB = pbias(param2_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
      Q5_param2_nse = NSE(param2_Q[Q5_days], act_Q[Q5_days])
      Q5_param2_kge = KGE(param2_Q[Q5_days], act_Q[Q5_days])
      Q5_param2_mse = mse(param2_Q[Q5_days], act_Q[Q5_days])
      Q5_param2_mae = mae(param2_Q[Q5_days], act_Q[Q5_days])
      Q5_param2_pctB = pbias(param2_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
      Q10_param2_nse = NSE(param2_Q[Q10_days], act_Q[Q10_days])
      Q10_param2_kge = KGE(param2_Q[Q10_days], act_Q[Q10_days])
      Q10_param2_mse = mse(param2_Q[Q10_days], act_Q[Q10_days])
      Q10_param2_mae = mae(param2_Q[Q10_days], act_Q[Q10_days])
      Q10_param2_pctB = pbias(param2_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
      Q20_param2_nse = NSE(param2_Q[Q20_days], act_Q[Q20_days])
      Q20_param2_kge = KGE(param2_Q[Q20_days], act_Q[Q20_days])
      Q20_param2_mse = mse(param2_Q[Q20_days], act_Q[Q20_days])
      Q20_param2_mae = mae(param2_Q[Q20_days], act_Q[Q20_days])
      Q20_param2_pctB = pbias(param2_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
      tnt_param2_nse = NSE(param2_Q[tnt_days], act_Q[tnt_days])
      tnt_param2_kge = KGE(param2_Q[tnt_days], act_Q[tnt_days])
      tnt_param2_mse = mse(param2_Q[tnt_days], act_Q[tnt_days])
      tnt_param2_mae = mae(param2_Q[tnt_days], act_Q[tnt_days])
      tnt_param2_pctB = pbias(param2_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
      
      model_perf[allws,'param2_nse'] = param2_nse
      model_perf[allws,'param2_kge'] = param2_kge
      #        model_perf[allws,'param2_log_nse'] = param2_nse_log
      model_perf[allws,'param2_mse'] = param2_mse
      model_perf[allws,'param2_mae'] = param2_mae
      model_perf[allws,'param2_pctB'] = param2_pctB
      model_perf[allws,'param2_Q5_nse'] = Q5_param2_nse
      model_perf[allws,'param2_Q5_kge'] = Q5_param2_kge
      model_perf[allws,'param2_Q5_mse'] = Q5_param2_mse
      model_perf[allws,'param2_Q5_mae'] = Q5_param2_mae
      model_perf[allws,'param2_Q5_pctB'] = Q5_param2_pctB
      model_perf[allws,'param2_Q10_nse'] = Q10_param2_nse
      model_perf[allws,'param2_Q10_kge'] = Q10_param2_kge
      model_perf[allws,'param2_Q10_mse'] = Q10_param2_mse
      model_perf[allws,'param2_Q10_mae'] = Q10_param2_mae
      model_perf[allws,'param2_Q10_pctB'] = Q10_param2_pctB
      model_perf[allws,'param2_Q20_nse'] = Q20_param2_nse
      model_perf[allws,'param2_Q20_kge'] = Q20_param2_kge
      model_perf[allws,'param2_Q20_mse'] = Q20_param2_mse
      model_perf[allws,'param2_Q20_mae'] = Q20_param2_mae
      model_perf[allws,'param2_Q20_pctB'] = Q20_param2_pctB
      model_perf[allws,'param2_tnt_nse'] = tnt_param2_nse
      model_perf[allws,'param2_tnt_kge'] = tnt_param2_kge
      model_perf[allws,'param2_tnt_mse'] = tnt_param2_mse
      model_perf[allws,'param2_tnt_mae'] = tnt_param2_mae
      model_perf[allws,'param2_tnt_pctB'] = tnt_param2_pctB
    }
    if(jj == 3) {
      param3_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
      act_Q = mopey_df_disch$Q + 10^(-5)
      
      param3_nse = NSE(param3_Q[-(1:1090)], act_Q[-(1:1090)])
      param3_kge = KGE(param3_Q[-(1:1090)], act_Q[-(1:1090)])
      #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
      #                      mopey_df_disch$Q[-(1:365)], FUN=log)
      param3_mse = mse(param3_Q[-(1:1090)], act_Q[-(1:1090)])
      param3_mae = mae(param3_Q[-(1:1090)], act_Q[-(1:1090)])
      param3_pctB = pbias(param3_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
      Q5_param3_nse = NSE(param3_Q[Q5_days], act_Q[Q5_days])
      Q5_param3_kge = KGE(param3_Q[Q5_days], act_Q[Q5_days])
      Q5_param3_mse = mse(param3_Q[Q5_days], act_Q[Q5_days])
      Q5_param3_mae = mae(param3_Q[Q5_days], act_Q[Q5_days])
      Q5_param3_pctB = pbias(param3_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
      Q10_param3_nse = NSE(param3_Q[Q10_days], act_Q[Q10_days])
      Q10_param3_kge = KGE(param3_Q[Q10_days], act_Q[Q10_days])
      Q10_param3_mse = mse(param3_Q[Q10_days], act_Q[Q10_days])
      Q10_param3_mae = mae(param3_Q[Q10_days], act_Q[Q10_days])
      Q10_param3_pctB = pbias(param3_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
      Q20_param3_nse = NSE(param3_Q[Q20_days], act_Q[Q20_days])
      Q20_param3_kge = KGE(param3_Q[Q20_days], act_Q[Q20_days])
      Q20_param3_mse = mse(param3_Q[Q20_days], act_Q[Q20_days])
      Q20_param3_mae = mae(param3_Q[Q20_days], act_Q[Q20_days])
      Q20_param3_pctB = pbias(param3_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
      tnt_param3_nse = NSE(param3_Q[tnt_days], act_Q[tnt_days])
      tnt_param3_kge = KGE(param3_Q[tnt_days], act_Q[tnt_days])
      tnt_param3_mse = mse(param3_Q[tnt_days], act_Q[tnt_days])
      tnt_param3_mae = mae(param3_Q[tnt_days], act_Q[tnt_days])
      tnt_param3_pctB = pbias(param3_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
      
      model_perf[allws,'param3_nse'] = param3_nse
      model_perf[allws,'param3_kge'] = param3_kge
      #        model_perf[allws,'param3_log_nse'] = param3_nse_log
      model_perf[allws,'param3_mse'] = param3_mse
      model_perf[allws,'param3_mae'] = param3_mae
      model_perf[allws,'param3_pctB'] = param3_pctB
      model_perf[allws,'param3_Q5_nse'] = Q5_param3_nse
      model_perf[allws,'param3_Q5_kge'] = Q5_param3_kge
      model_perf[allws,'param3_Q5_mse'] = Q5_param3_mse
      model_perf[allws,'param3_Q5_mae'] = Q5_param3_mae
      model_perf[allws,'param3_Q5_pctB'] = Q5_param3_pctB
      model_perf[allws,'param3_Q10_nse'] = Q10_param3_nse
      model_perf[allws,'param3_Q10_kge'] = Q10_param3_kge
      model_perf[allws,'param3_Q10_mse'] = Q10_param3_mse
      model_perf[allws,'param3_Q10_mae'] = Q10_param3_mae
      model_perf[allws,'param3_Q10_pctB'] = Q10_param3_pctB
      model_perf[allws,'param3_Q20_nse'] = Q20_param3_nse
      model_perf[allws,'param3_Q20_kge'] = Q20_param3_kge
      model_perf[allws,'param3_Q20_mse'] = Q20_param3_mse
      model_perf[allws,'param3_Q20_mae'] = Q20_param3_mae
      model_perf[allws,'param3_Q20_pctB'] = Q20_param3_pctB
      model_perf[allws,'param3_tnt_nse'] = tnt_param3_nse
      model_perf[allws,'param3_tnt_kge'] = tnt_param3_kge
      model_perf[allws,'param3_tnt_mse'] = tnt_param3_mse
      model_perf[allws,'param3_tnt_mae'] = tnt_param3_mae
      model_perf[allws,'param3_tnt_pctB'] = tnt_param3_pctB
    }
    if(jj == 4) {
      param4_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
      act_Q = mopey_df_disch$Q + 10^(-5)
      
      param4_nse = NSE(param4_Q[-(1:1090)], act_Q[-(1:1090)])
      param4_kge = KGE(param4_Q[-(1:1090)], act_Q[-(1:1090)])
      #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
      #                      mopey_df_disch$Q[-(1:365)], FUN=log)
      param4_mse = mse(param4_Q[-(1:1090)], act_Q[-(1:1090)])
      param4_mae = mae(param4_Q[-(1:1090)], act_Q[-(1:1090)])
      param4_pctB = pbias(param4_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
      Q5_param4_nse = NSE(param4_Q[Q5_days], act_Q[Q5_days])
      Q5_param4_kge = KGE(param4_Q[Q5_days], act_Q[Q5_days])
      Q5_param4_mse = mse(param4_Q[Q5_days], act_Q[Q5_days])
      Q5_param4_mae = mae(param4_Q[Q5_days], act_Q[Q5_days])
      Q5_param4_pctB = pbias(param4_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
      Q10_param4_nse = NSE(param4_Q[Q10_days], act_Q[Q10_days])
      Q10_param4_kge = KGE(param4_Q[Q10_days], act_Q[Q10_days])
      Q10_param4_mse = mse(param4_Q[Q10_days], act_Q[Q10_days])
      Q10_param4_mae = mae(param4_Q[Q10_days], act_Q[Q10_days])
      Q10_param4_pctB = pbias(param4_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
      Q20_param4_nse = NSE(param4_Q[Q20_days], act_Q[Q20_days])
      Q20_param4_kge = KGE(param4_Q[Q20_days], act_Q[Q20_days])
      Q20_param4_mse = mse(param4_Q[Q20_days], act_Q[Q20_days])
      Q20_param4_mae = mae(param4_Q[Q20_days], act_Q[Q20_days])
      Q20_param4_pctB = pbias(param4_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
      tnt_param4_nse = NSE(param4_Q[tnt_days], act_Q[tnt_days])
      tnt_param4_kge = KGE(param4_Q[tnt_days], act_Q[tnt_days])
      tnt_param4_mse = mse(param4_Q[tnt_days], act_Q[tnt_days])
      tnt_param4_mae = mae(param4_Q[tnt_days], act_Q[tnt_days])
      tnt_param4_pctB = pbias(param4_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
      
      model_perf[allws,'param4_nse'] = param4_nse
      model_perf[allws,'param4_kge'] = param4_kge
      #        model_perf[allws,'param4_log_nse'] = param4_nse_log
      model_perf[allws,'param4_mse'] = param4_mse
      model_perf[allws,'param4_mae'] = param4_mae
      model_perf[allws,'param4_pctB'] = param4_pctB
      model_perf[allws,'param4_Q5_nse'] = Q5_param4_nse
      model_perf[allws,'param4_Q5_kge'] = Q5_param4_kge
      model_perf[allws,'param4_Q5_mse'] = Q5_param4_mse
      model_perf[allws,'param4_Q5_mae'] = Q5_param4_mae
      model_perf[allws,'param4_Q5_pctB'] = Q5_param4_pctB
      model_perf[allws,'param4_Q10_nse'] = Q10_param4_nse
      model_perf[allws,'param4_Q10_kge'] = Q10_param4_kge
      model_perf[allws,'param4_Q10_mse'] = Q10_param4_mse
      model_perf[allws,'param4_Q10_mae'] = Q10_param4_mae
      model_perf[allws,'param4_Q10_pctB'] = Q10_param4_pctB
      model_perf[allws,'param4_Q20_nse'] = Q20_param4_nse
      model_perf[allws,'param4_Q20_kge'] = Q20_param4_kge
      model_perf[allws,'param4_Q20_mse'] = Q20_param4_mse
      model_perf[allws,'param4_Q20_mae'] = Q20_param4_mae
      model_perf[allws,'param4_Q20_pctB'] = Q20_param4_pctB
      model_perf[allws,'param4_tnt_nse'] = tnt_param4_nse
      model_perf[allws,'param4_tnt_kge'] = tnt_param4_kge
      model_perf[allws,'param4_tnt_mse'] = tnt_param4_mse
      model_perf[allws,'param4_tnt_mae'] = tnt_param4_mae
      model_perf[allws,'param4_tnt_pctB'] = tnt_param4_pctB
    }
    if(jj == 5) {
      param5_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out_wCalib
      act_Q = mopey_df_disch$Q + 10^(-5)
      
      param5_nse = NSE(param5_Q[-(1:1090)], act_Q[-(1:1090)])
      param5_kge = KGE(param5_Q[-(1:1090)], act_Q[-(1:1090)])
      #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
      #                      mopey_df_disch$Q[-(1:365)], FUN=log)
      param5_mse = mse(param5_Q[-(1:1090)], act_Q[-(1:1090)])
      param5_mae = mae(param5_Q[-(1:1090)], act_Q[-(1:1090)])
      param5_pctB = pbias(param5_Q[-(1:1090)], act_Q[-(1:1090)]) * mean(mopey_df_disch$Q[-(1:1090)], na.rm=TRUE)/ 100
      Q5_param5_nse = NSE(param5_Q[Q5_days], act_Q[Q5_days])
      Q5_param5_kge = KGE(param5_Q[Q5_days], act_Q[Q5_days])
      Q5_param5_mse = mse(param5_Q[Q5_days], act_Q[Q5_days])
      Q5_param5_mae = mae(param5_Q[Q5_days], act_Q[Q5_days])
      Q5_param5_pctB = pbias(param5_Q[Q5_days], act_Q[Q5_days]) * mean(mopey_df_disch$Q[Q5_days]) / 100
      Q10_param5_nse = NSE(param5_Q[Q10_days], act_Q[Q10_days])
      Q10_param5_kge = KGE(param5_Q[Q10_days], act_Q[Q10_days])
      Q10_param5_mse = mse(param5_Q[Q10_days], act_Q[Q10_days])
      Q10_param5_mae = mae(param5_Q[Q10_days], act_Q[Q10_days])
      Q10_param5_pctB = pbias(param5_Q[Q10_days], act_Q[Q10_days]) * mean(mopey_df_disch$Q[Q10_days])/ 100
      Q20_param5_nse = NSE(param5_Q[Q20_days], act_Q[Q20_days])
      Q20_param5_kge = KGE(param5_Q[Q20_days], act_Q[Q20_days])
      Q20_param5_mse = mse(param5_Q[Q20_days], act_Q[Q20_days])
      Q20_param5_mae = mae(param5_Q[Q20_days], act_Q[Q20_days])
      Q20_param5_pctB = pbias(param5_Q[Q20_days], act_Q[Q20_days]) * mean(mopey_df_disch$Q[Q20_days])/ 100
      tnt_param5_nse = NSE(param5_Q[tnt_days], act_Q[tnt_days])
      tnt_param5_kge = KGE(param5_Q[tnt_days], act_Q[tnt_days])
      tnt_param5_mse = mse(param5_Q[tnt_days], act_Q[tnt_days])
      tnt_param5_mae = mae(param5_Q[tnt_days], act_Q[tnt_days])
      tnt_param5_pctB = pbias(param5_Q[tnt_days], act_Q[tnt_days]) * mean(mopey_df_disch$Q[tnt_days])/ 100
      
      model_perf[allws,'param5_nse'] = param5_nse
      model_perf[allws,'param5_kge'] = param5_kge
      #        model_perf[allws,'param5_log_nse'] = param5_nse_log
      model_perf[allws,'param5_mse'] = param5_mse
      model_perf[allws,'param5_mae'] = param5_mae
      model_perf[allws,'param5_pctB'] = param5_pctB
      model_perf[allws,'param5_Q5_nse'] = Q5_param5_nse
      model_perf[allws,'param5_Q5_kge'] = Q5_param5_kge
      model_perf[allws,'param5_Q5_mse'] = Q5_param5_mse
      model_perf[allws,'param5_Q5_mae'] = Q5_param5_mae
      model_perf[allws,'param5_Q5_pctB'] = Q5_param5_pctB
      model_perf[allws,'param5_Q10_nse'] = Q10_param5_nse
      model_perf[allws,'param5_Q10_kge'] = Q10_param5_kge
      model_perf[allws,'param5_Q10_mse'] = Q10_param5_mse
      model_perf[allws,'param5_Q10_mae'] = Q10_param5_mae
      model_perf[allws,'param5_Q10_pctB'] = Q10_param5_pctB
      model_perf[allws,'param5_Q20_nse'] = Q20_param5_nse
      model_perf[allws,'param5_Q20_kge'] = Q20_param5_kge
      model_perf[allws,'param5_Q20_mse'] = Q20_param5_mse
      model_perf[allws,'param5_Q20_mae'] = Q20_param5_mae
      model_perf[allws,'param5_Q20_pctB'] = Q20_param5_pctB
      model_perf[allws,'param5_tnt_nse'] = tnt_param5_nse
      model_perf[allws,'param5_tnt_kge'] = tnt_param5_kge
      model_perf[allws,'param5_tnt_mse'] = tnt_param5_mse
      model_perf[allws,'param5_tnt_mae'] = tnt_param5_mae
      model_perf[allws,'param5_tnt_pctB'] = tnt_param5_pctB
    }
    
    these_days = which(!is.na(mopey_df_disch$Q))[2500:2865] # for illustrative plotting
    if(jj==1) {
      plot(mopey_df_disch$Q[these_days], type='l',
           ylim=c(0.01,max(mopey_df_disch$Q[these_days]+8, na.rm=TRUE)),
           log='y',
           col='grey5', lwd=3,
           ylab = "Q (mm)", xlab = "Date",
           #!!! TEMPFIX
           #		main=paste(mopey_catalogue[this_ws,"Name"],
           main=paste(this_ws,
                      #!!! TEMPFIX	           
                      ": \n old KGE =", round(old_kge,2), "new =", round(new_kge,2), "new w.Cal=", round(newC_kge,2)))
      lines(mopey_df_disch$Qg[these_days], lwd=3, col='purple1', lty=1)
    }
    if(jj==2){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
            col='purple2', lwd=3, lty=1)
    }
    if(jj==3){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out_wCalib[these_days],
            col='royalblue1', lwd=3)
    }      
    if(jj==4){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out_wCalib[these_days],
            col='purple3', lwd=3)
    }      
    if(jj==5){
      lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out_wCalib[these_days],
            col='purple4', lwd=3)
    }      
  }
  write.csv(model_perf, "C:\\Users\\arik\\Documents\\LSM Research\\model_perf_26apr2021_exp_2x_100div_3up_95_b_hlfB_0evap_KGE_sens-param_v.csv")
  print(allws)
}

































################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .5, 4)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 4)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -1, 4)#0, 3)	#melt temperature [C]
  fm   = runif(100, .5, 6)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 400)#100, 200)
  lp   = runif(100, .2, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, .9)#0.09, 0.1)
  k1   = runif(100, .001, .5)#0.05, 0.07)
  k2   = runif(100, .0001, .1)#0.05)
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)


this_ws = 8

  #!!!! TEMPFIX
#mopey_df = tsMOPEX(mopey_catalogue$USGS_ID[this_ws])	# identifying the first ts data from the first idnetified mopex gage
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", mopey_catalogue[this_ws], ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", mopey_catalogue[this_ws], ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", mopey_catalogue[this_ws], ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
  mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
  }
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
#!!!! TEMPFIX


#mopey_df_nonas = subset(mopey_df, !is.na(Q) & !is.na(P))
#mopey_df_nonas$T_mean = apply(mopey_df_nonas[,c("T_max","T_min")],1,mean)
mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)


n_loops = 10
n_runs = 1000

cal_out = data.frame(
  sfcf=rep(NA,n_runs),
  tr=rep(NA,n_runs), 
  tt=rep(NA,n_runs),
  fm=rep(NA,n_runs),
  fi=rep(NA,n_runs),
  fic=rep(NA,n_runs),
  fc=rep(NA,n_runs),
  lp=rep(NA,n_runs),
  beta_soils = rep(NA,n_runs),
  k0 = rep(NA,n_runs),
  k1 = rep(NA,n_runs),
  k2 = rep(NA,n_runs),
  uz1 = rep(NA,n_runs),
  perc = rep(NA,n_runs),
  nse = rep(NA,n_runs))


for(ii in 1:n_loops){
	for(jj in 1:n_runs){
		cal_out$sfcf[jj] = sample(sfcf,1)
		cal_out$tr[jj] = sample(tr,1)
		cal_out$tt[jj] = sample(tt,1)
		cal_out$fm[jj] = sample(fm,1)
		cal_out$fi[jj] = sample(fi,1)
		cal_out$fic[jj] = sample(fic,1)
		cal_out$fc[jj] = sample(fc,1)
		cal_out$lp[jj] = sample(lp,1)
		cal_out$beta_soils[jj] = sample(beta_soils,1)
		# since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
		cal_out$k0[jj] = sample(k0,1)
		cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
		cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
		cal_out$uz1[jj] = sample(uz1,1)
		cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)

		mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


		the_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
		plot(mopey_df_disch$Q[365:730], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[365:730],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFFIX
#			main=paste(mopey_catalogue[this_ws,"Name"], ": NSE =", round(the_nse,2)))
      main=paste(mopey_catalogue[this_ws], ": NSE =", round(the_nse,2)))
		#!!! TEMPFFIX
    lines(mopey_df_disch$Qg[365:730], lwd=2, col='blue', lty=1)

		# calibrating by NSE
		cal_out$nse[jj][-(1:365)] = the_nse
		# calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
			# basically NSE - .1 * abs(relative volume error)
		

	}
	print(ii)
	cal_sort = head(cal_out[rev(order(cal_out$nse)),], 100)
	# snow and glacier module
	  sfcf = runif(100, min(cal_sort$sfcf), max(cal_sort$sfcf))#1, 2)	#snowfall correction factor [-]
	  tr   = runif(100, min(cal_sort$tr), max(cal_sort$tr))#-1, 3)	#solid and liquid precipitation threshold temperature [C]
	  tt   = runif(100, min(cal_sort$tt), max(cal_sort$tt))#0, 3)	#melt temperature [C]
	  fm   = runif(100, min(cal_sort$fm), max(cal_sort$fm))#1, 4)	#snowmelt factor [mm/C]
	  fi   = runif(100, min(cal_sort$fi), max(cal_sort$fi))#4, 8)	#icemelt factor [mm/C]
	  fic  = runif(100, min(cal_sort$fic), max(cal_sort$fic))#6, 6)	#debris-covered icemelt factor [mm/C]
	# soil module
	  fc   = runif(100, min(cal_sort$fc), max(cal_sort$fc))#100, 200)
	  lp   = runif(100, min(cal_sort$lp), max(cal_sort$lp))#0.5, 1)
	  beta_soils = runif(100, min(cal_sort$beta_soils), max(cal_sort$beta_soils))
	# routing module
	  k0   = runif(100, min(cal_sort$k0), max(cal_sort$k0))#0.09, 0.1)
	  k1   = runif(100, min(cal_sort$k1), max(cal_sort$k1))#0.05, 0.07)
	  k2   = runif(100, min(cal_sort$k2), max(cal_sort$k2))#0.05)
	  uz1  = runif(100, min(cal_sort$uz1), max(cal_sort$uz1))#1, 5)
	  perc = runif(100, min(cal_sort$perc), max(cal_sort$perc))#.8, 2)
}



#########################################################################################
#### model 2: adding our gw module
library(sf)
	# reading in the data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\PhD Research\\D3")
aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")

	# choose the watersheds
  #!!! TEMPFIX
#for(gg in 1:nrow(mopey_catalogue)){
for(gg in mopey_catalogue){
  #!!! TEMPFIX
  
  this_ws = gg
    #!!! TEMPFIX
  #if(mopey_catalogue[this_ws, "USGS_ID"] %in% aq_chars$gage) 	print(this_ws)
  if(this_ws %in% aq_chars$gage) 	print(this_ws)
    #!!! TEMPFIX
}

  #!!! TEMPFIX
this_ws = "14232500"
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
    mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
}
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
this_ws_props = which(aq_chars$gage == this_ws)
  #!!! TEMPFIX


D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}


mopey_df_disch$deep_rech = c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2
f = 0.1 # the drainable porosoity we set for recession analysis
mopey_df_disch$gw_out = NA
mopey_df_disch$deep_store = NA
mopey_df_disch$deep_store[1] =  mopey_df_disch$SLZ[1]
depth = mopey_df_disch$deep_store[1] / f
B_half = B/2

	# Start the clock
	prm = proc.time()

	# calculating values so I don't have to recalculate at every loop
sin_slp = sin(i_slp)
strm_crs_sxn = 2*Wo/A
		# Start the clock!
		ptm <- proc.time()
for(gg in 1:(nrow(mopey_df_disch)-1))	{
	k_arb = k_fof_D(D_thickness = depth)
	gw_ridge = ifelse(D_wet > depth, ((D_wet - depth) / D_wet)^2, .0001)
	mopey_df_disch$gw_out[gg] = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
		k_arb *			# 2k to account for 2 hillslopes per reach of stream
#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
		(((sin_slp * B_half) + 2 * depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
		(depth / 2) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
		strm_crs_sxn					# 
	mopey_df_disch$deep_store[gg+1] = mopey_df_disch$deep_store[gg] + mopey_df_disch$deep_rech[gg] - mopey_df_disch$gw_out[gg]
	depth = mopey_df_disch$deep_store[gg+1] / f
}
	# Stop the clock
	proc.time() - ptm


		# Start the clock!
		ptm <- proc.time()
f = 0.1 # the drainable porosoity we set for recession analysis
gw_out = NULL
deep_store = as.numeric(mopey_df_disch$SLZ[1])
depth = deep_store / f
deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
for(gg in 1:(nrow(mopey_df_disch)))	{
	k_arb = as.numeric(k_fof_D(D_thickness = depth))
	gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^2, .001)
	gw_out_iter = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
		k_arb *			# 2k to account for 2 hillslopes per reach of stream
#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
		(((sin_slp * B_half) + 2 * depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
		(depth / 2) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
		strm_crs_sxn					# 
	
	new_deep_store = deep_store[gg] +
			deep_rech[gg] -
			gw_out_iter
	deep_store = c(deep_store, new_deep_store)
	depth = new_deep_store / f

	gw_out = c(gw_out, gw_out_iter)
	
}
	
mopey_df_disch$deep_store = deep_store[-length(deep_store)]
mopey_df_disch$deep_rech = deep_rech
mopey_df_disch$gw_out = gw_out
	# Stop the clock
	proc.time() - ptm



these_days = 2600:2965
old_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
new_nse = NSE((mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out), mopey_df_disch$Q)
		plot(mopey_df_disch$Q[these_days], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[2600:2965],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#			main=paste(mopey_catalogue[this_ws,"Name"],
      main=paste(this_ws,
  #!!! TEMPFIX
			": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
		lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
		lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
			col='blue', lwd=2)










#######################################################
### step 3: calibrating HBV with our model included
	# reading in the data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\PhD Research\\D3")
aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")

	# choose the watersheds
ws_ihave = NULL
  #!!! TEMPFIX
#for(gg in 1:nrow(mopey_catalogue)){
  #if(mopey_catalogue[gg, "USGS_ID"] %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
#}
#head(mopey_catalogue[ws_ihave, ])
for(gg in mopey_catalogue){
  if(gg %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
}
head(ws_ihave)
  #!!! TEMPFIX




#################################################################################
	# initializing data for our gw module

  #!!! TEMPFIX
which_ws = all_ws
#this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
  #!!! TEMPFIX

D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}



#########################################################################################
	# initializing data for HBV

  #!!! TEMPFIX
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
this_ws = ws_ihave[which_ws]
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
    mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
}
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
#!!! TEMPFIX




################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .3, 5)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 5)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -2, 5)#0, 3)	#melt temperature [C]
  fm   = runif(100, .2, 8)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 2000)#100, 200)
  lp   = runif(100, .1, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, 1)#0.09, 0.1)
  k1   = runif(100, .0001, .5)#0.05, 0.07)
  k2   = 0.01#runif(100, .0001, .1)#0.05)	# since we won't be using this component of hbv
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)




mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)


n_runs = 100
n_loops = 10


cal_out = data.frame(
	sfcf=rep(NA,n_runs),
	tr=rep(NA,n_runs), 
	tt=rep(NA,n_runs),
	fm=rep(NA,n_runs),
	fi=rep(NA,n_runs),
	fic=rep(NA,n_runs),
	fc=rep(NA,n_runs),
	lp=rep(NA,n_runs),
	beta_soils = rep(NA,n_runs),
	k0 = rep(NA,n_runs),
	k1 = rep(NA,n_runs),
	k2 = rep(NA,n_runs),
	uz1 = rep(NA,n_runs),
	perc = rep(NA,n_runs),
	nse = rep(NA,n_runs))

	# Start the clock
ptm = proc.time()
for(ii in 1:n_loops){
	for(jj in 1:n_runs){
		cal_out$sfcf[jj] = sample(sfcf,1)
		cal_out$tr[jj] = sample(tr,1)
		cal_out$tt[jj] = sample(tt,1)
		cal_out$fm[jj] = sample(fm,1)
		cal_out$fi[jj] = sample(fi,1)
		cal_out$fic[jj] = sample(fic,1)
		cal_out$fc[jj] = sample(fc,1)
		cal_out$lp[jj] = sample(lp,1)
		cal_out$beta_soils[jj] = sample(beta_soils,1)
		# since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
		cal_out$k0[jj] = sample(k0,1)
		cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
		cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
		cal_out$uz1[jj] = sample(uz1,1)
		cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)

		mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


####################################################################
	# new gw module introduced

			f = 0.1 # the drainable porosoity we set for recession analysis
			gw_out = NULL
			deep_store = as.numeric(mopey_df_disch$SLZ[1])
			depth = deep_store / f
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
			recent_rech = seq(1/180,1,length=180)^2
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^2, .002)
				gw_out_iter =
				  f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
					k_arb *			# 2k to account for 2 hillslopes per reach of stream
			#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
					(((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
					(depth / 4 + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
					strm_crs_sxn					# 
				if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
				
				new_deep_store = deep_store[gg] +
						deep_rech[gg] -
						gw_out_iter
				deep_store = c(deep_store, new_deep_store)
				depth = new_deep_store / f

				gw_out = c(gw_out, gw_out_iter)
				
			}
				
			mopey_df_disch$deep_store = deep_store[-length(deep_store)]
			mopey_df_disch$deep_rech = deep_rech
			mopey_df_disch$gw_out = gw_out

##############################################################################
	# plotting and calibrating
	
	
	these_days = 3300:3665
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#		main=paste(mopey_catalogue[this_ws,"Name"],
	  main=paste(this_ws,
	#!!! TEMPFIX	           
		": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
	lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
	lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
		col='blue', lwd=2)

	# calibrating by NSE
	cal_out$nse[jj] = new_nse
	# calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
		# basically NSE - .1 * abs(relative volume error)

	}
  
	print(ii)
	cal_sort = head(cal_out[rev(order(cal_out$nse)),], 10)
	print(head(cal_sort, 20))
	# snow and glacier module
	  sfcf = runif(100, min(cal_sort$sfcf), max(cal_sort$sfcf))#1, 2)	#snowfall correction factor [-]
	  tr   = runif(100, min(cal_sort$tr), max(cal_sort$tr))#-1, 3)	#solid and liquid precipitation threshold temperature [C]
	  tt   = runif(100, min(cal_sort$tt), max(cal_sort$tt))#0, 3)	#melt temperature [C]
	  fm   = runif(100, min(cal_sort$fm), max(cal_sort$fm))#1, 4)	#snowmelt factor [mm/C]
	  fi   = runif(100, min(cal_sort$fi), max(cal_sort$fi))#4, 8)	#icemelt factor [mm/C]
	  fic  = runif(100, min(cal_sort$fic), max(cal_sort$fic))#6, 6)	#debris-covered icemelt factor [mm/C]
	# soil module
	  fc   = runif(100, min(cal_sort$fc), max(cal_sort$fc))#100, 200)
	  lp   = runif(100, min(cal_sort$lp), max(cal_sort$lp))#0.5, 1)
	  beta_soils = runif(100, min(cal_sort$beta_soils), max(cal_sort$beta_soils))
	# routing module
	  k0   = runif(100, min(cal_sort$k0), max(cal_sort$k0))#0.09, 0.1)
	  k1   = runif(100, min(cal_sort$k1), max(cal_sort$k1))#0.05, 0.07)
	  #k2   = runif(100, min(cal_sort$k2), max(cal_sort$k2))#0.05)
	  uz1  = runif(100, min(cal_sort$uz1), max(cal_sort$uz1))#1, 5)
	  perc = runif(100, min(cal_sort$perc), max(cal_sort$perc))#.8, 2)
}
	#Stop clock
proc.time() - ptm








# now assessing variations in our new gw module

mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_sort$sfcf[1],		#SFCF - snowfall correction factor [-]
							cal_sort$tr[1],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_sort$tt[1],		#Tt - melt temperature [C]
							cal_sort$fm[1],		#fm - snowmelt factor [mm/C]
							cal_sort$fi[1],		#fi - icemelt factor [mm/C]
							cal_sort$fic[1])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_sort$fc[1],			#FC - soil field capacity [mm]
							cal_sort$lp[1],			#LP - parameter ot get actual ET [-]
							cal_sort$beta_soils[1])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_sort$k0[1],	#KO - top bucket (STZ) storage constant [1/t]
							cal_sort$k1[1],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_sort$k2[1],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_sort$uz1[1],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_sort$perc[1])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


####################################################################
	# new gw module introduced
			length_rech_dist = 20:190			# opt based on 3 runs
			exp_rech_dist = runif(100,1,3) # opt based on 2 runs: smaller is better
			depth_sclr_dist = runif(100,1,8) # opt based on 2 runs: smaller is better (1)
			ridge_exp_dist = runif(100,2,4)	# opt based on 2 runs: doesn't seem to matter
			ridge_max_dist = runif(100,0.0001,0.05) # optimal based on 2 runs: btw .0001 and .002
			
	
		

n_runs = 1000

cal_out_3 = data.frame(
  length_rech=rep(NA,n_runs),
  exp_rech=rep(NA,n_runs), 
  depth_sclr=rep(NA,n_runs),
  ridge_exp=rep(NA,n_runs),
  ridge_max=rep(NA,n_runs),
  nse=rep(NA,n_runs)
)

	# Stop the clock
ptm = proc.time()
	for(yy in 1:n_runs){
		cal_out_3$length_rech[yy] = sample(length_rech_dist,1)
		cal_out_3$exp_rech[yy] = sample(exp_rech_dist,1)
		cal_out_3$depth_sclr[yy] = sample(depth_sclr_dist,1)
		cal_out_3$ridge_exp[yy] = sample(ridge_exp_dist,1)
		cal_out_3$ridge_max[yy] = sample(ridge_max_dist,1)
		
			f = 0.1 # the drainable porosoity we set for recession analysis
			gw_out = NULL
			deep_store = as.numeric(mopey_df_disch$SLZ[1])
			depth = deep_store / f
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
			recent_rech = seq(.01, 1, length=sample(cal_out_3$length_rech[yy],1))^cal_out_3$exp_rech[yy]
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^cal_out_3$ridge_exp[yy], cal_out_3$ridge_max[yy])
				gw_out_iter =
				  f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
				  k_arb *			# 2k to account for 2 hillslopes per reach of stream
				  #		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
				  (((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
				  (depth / 4 + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
				  strm_crs_sxn					# 
				if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
				
				new_deep_store = deep_store[gg] +
						deep_rech[gg] -
						gw_out_iter
				deep_store = c(deep_store, new_deep_store)
				depth = new_deep_store / f

				gw_out = c(gw_out, gw_out_iter)
				
			}
				
			mopey_df_disch$deep_store = deep_store[-length(deep_store)]
			mopey_df_disch$deep_rech = deep_rech
			mopey_df_disch$gw_out = gw_out
	
	these_days = 3300:3965
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#		main=paste(mopey_catalogue[this_ws,"Name"],
		main=paste(this_ws,
	#!!! TEMPFIX	           
		": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
	lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
	lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
		col='blue', lwd=2)

	cal_out_3$nse[yy] = new_nse
	
	}

	cal_sort_3 = cal_out_3[order(cal_out_3$nse),]
	#cal_sort_ny for ny

 plot(cal_out_3$length_rech, cal_out_3$nse)
 plot(cal_out_3$exp_rech, cal_out_3$nse)
 plot(cal_out_3$depth_sclr, cal_out_3$nse)
 plot(cal_out_3$ridge_exp, cal_out_3$nse)
 plot(cal_out_3$ridge_max, cal_out_3$nse)


##############################################################################
	# plotting and calibrating
	
	
	













	









































this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])

D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}


mopey_df_disch$deep_rech = c(0,diff(mopey_df_disch$SLZ) + mopey_df_disch$Q2)
f = 0.1 # the drainable porosoity we set for recession analysis
mopey_df_disch$gw_out = NA
mopey_df_disch$deep_store = NA
mopey_df_disch$deep_store[1] =  mopey_df_disch$SLZ[1]
depth = mopey_df_disch$deep_store[1] / f

	













































































################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .5, 4)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 4)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -1, 4)#0, 3)	#melt temperature [C]
  fm   = runif(100, .5, 6)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 400)#100, 200)
  lp   = runif(100, .2, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, .9)#0.09, 0.1)
  k1   = runif(100, .001, .5)#0.05, 0.07)
  k2   = runif(100, .0001, .1)#0.05)
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)
















mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


		the_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
		plot(mopey_df_disch$Q[635:1000], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[635:1000,],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
			main=paste(mopey_catalogue[this_ws,"Name"], ": NSE =", round(the_nse,2)))
		lines(mopey_df_disch$Qg[635:1000], lwd=2, col='blue', lty=1)

		cal_out$nse[jj] = the_nse
	}






































# snow and glacier module
glacier_range <- rbind(
  sfcf = c(1, 2),	#snowfall correction factor [-]
  tr   = c(-1, 3),	#solid and liquid precipitation threshold temperature [C]
  tt   = c(0, 3),	#melt temperature [C]
  fm   = c(1, 4),	#snowmelt factor [mm/C]
  fi   = c(4, 8),	#icemelt factor [mm/C]
  fic  = c(6, 6)	#debris-covered icemelt factor [mm/C]
)


  # soil module
soil_range <- rbind(
  fc   = c(100, 200),
  lp   = c(0.5, 1),
  beta = c(1, 3)
)

  # routing module (here I will give you the correct values)
routing_range <- rbind(
  k0   = c(0.09, 0.09),
  k1   = c(0.07, 0.07),
  k2   = c(0.05, 0.05),
  uzl  = c(5, 5),
  perc = c(2, 2)
)

  # transfer function module (I will give the correct value)
tf_range <- rbind(
  bmax = c(2.25, 2.25)
)



    param_snow = c(1.1, 0, 0, 2.5)
    param_soil = c(150, 0.90, 1.5)
    param_route = c(0.09, 0.07, 0.05, 5, 2)
    param_tf = c(3.00)





# paramaters I don't need to calibrate
	
	 # precip model
	precip_range <- rbind(
	  p_grad  = c(5, 25)
	)

	 # air temperature model
	tair_range <- rbind(
	  t_grad  = c(-9.8, -2)
	)



## Case example with the first model
inputMatrix <- cbind(
                     runif(n = 200, max = 100, min = 0),
                     runif(n = 200, max = 50, min = 5),
                     runif(n = 100, max = 3, min = 1)
                     )

routeMod1   <- Routing_HBV(model = 1, lake = TRUE, inputData = inputMatrix,
                     initCond = c(10, 15, 20), param = c(0.1, 0.05, 0.001, 1, 0.8))





############################################################
# part 1
# Precip (from gage and possibly include elevation gradient)

	# we don't need to use this unless we want to introduce a topo gradient for ppt
Precip_model(
       model, # =1 for linear ppt gradient; =2 for linear ppt gradient with upper threshold
       inputData, # ppt (mm/t) from gage (vector)
       zmeteo,	# altitude of ppt gage (masl) (numeric)
       ztopo,	# target height (masl) (numberic)
       param # ppt gradient (%/100m); optional second value (for model=2) that sets max threshold elevation above which ppt doesn't increase
)


############################################################
# part 2
# Timeseries of rainfall and snowmelt (from part 1 and also bring in Temp)

SnowGlacier_HBV(
       model, # =1 for temeprature index model; =2 or =3 for variable snow cover and variable glacier cover
       inputData, # two columns, column_1= air temperature time series, column_2= ppt series
       initCond, # 3 values, 1st is SWE initial condition (mm), 2nd is 1/2/3 for clean ice/soil/dirty ice, 3rd set to 0 since we have no glaciers
	   param # snowfall correction factor (SFCF) [-], solid / liquid ppt threshold temperature (Tr), melt temperature (Tt), snowmelt factor (fm) in mm/T>0, icemelt factor (fi), debris-covered icemelt factor
)
	# outputs 1) rainfall; 2) snowfall; 3) SWE; 4) melted snow; rainfall + melted snow

############################################################
# part 3
# Soil routing model

Soil_HBV(
       model,
       inputData,
       initCond,
       param
       )


############################################################
# part 4
# streamflow partitioning model






# we first consider the SnowGlacier module to take into account 
# precipitation partioning and the snow accumulation/melting. 
snow_module <-
  SnowGlacier_HBV(model = 1, 
                  inputData = as.matrix( lumped_hbv[ , c('T(ºC)', 'P(mm/d)')] ),
                  initCond = c(20, 2), 
                  param = c(1.20, 1.00, 0.00, 2.5) )
