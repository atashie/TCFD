########################################################
### library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)
	sf::sf_use_s2(FALSE) # for problem with intersecting spherical w flat
library(lubridate)
library(ddplyr) # for left_join() merging by date
library(ncdf4)
library(magrittr)
library(maptools)



#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\TCFD\\GWstorage\\"
	# reading ssp126 climate data 
ncname_gfdl = 'watergap2-2c_gfdl-esm2m_ewembi_rcp26_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_gfdl = nc_open(paste0(ncpath, ncname_gfdl))
nc_gfdl = ncvar_get(ncin_gfdl,"groundwstor")	# lon, lat, time

ncname_hadgem = 'watergap2-2c_hadgem2-es_ewembi_rcp26_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_hadgem = nc_open(paste0(ncpath, ncname_hadgem))
nc_hadgem = ncvar_get(ncin_hadgem,"groundwstor")	# lon, lat, time

ncname_ipsl = 'watergap2-2c_ipsl-cm5a-lr_ewembi_rcp26_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_ipsl = nc_open(paste0(ncpath, ncname_ipsl))
nc_ipsl = ncvar_get(ncin_ipsl,"groundwstor")	# lon, lat, time

ncname_miroc = 'watergap2-2c_miroc5_ewembi_rcp26_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncin_miroc = nc_open(paste0(ncpath, ncname_miroc))
nc_miroc = ncvar_get(ncin_miroc,"groundwstor")	# lon, lat, time

nc_lat = ncvar_get(ncin_miroc, 'lat')	# lat is given from high to low
nc_lon = ncvar_get(ncin_miroc, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin_miroc, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
missing_data = 1.00000002004088e+20


dates1019 = which(year(nc_date) == 2010)[1]:which(year(nc_date) == 2019)[12]
dates2029 = which(year(nc_date) == 2020)[1]:which(year(nc_date) == 2029)[12]
dates3039 = which(year(nc_date) == 2030)[1]:which(year(nc_date) == 2039)[12]
dates4049 = which(year(nc_date) == 2040)[1]:which(year(nc_date) == 2049)[12]
dates5059 = which(year(nc_date) == 2050)[1]:which(year(nc_date) == 2059)[12]
dates6069 = which(year(nc_date) == 2060)[1]:which(year(nc_date) == 2069)[12]
dates7079 = which(year(nc_date) == 2070)[1]:which(year(nc_date) == 2079)[12]
dates8089 = which(year(nc_date) == 2080)[1]:which(year(nc_date) == 2089)[12]
dates9099 = which(year(nc_date) == 2090)[1]:which(year(nc_date) == 2099)[12]


quantile_data = data.frame(Lat = NA, Lon = NA, Q_25_abs = NA, Q_50_abs = NA, Q_75_abs = NA, 
	Q_25_rel = NA, Q_50_rel = NA, Q_75_rel = NA,
	Q_25_indx = NA, Q_50_indx = NA, Q_75_indx = NA,
	Likelihood = NA)

iter = 1
save_iter = 0
for(i in 1:length(nc_lat))	{
	for(j in 1:length(nc_lon))	{
		nc_dummy = nc_gfdl[j,i,dates1019] # reading in one data set to test for nona
		if(any(!is.na(nc_dummy) & any(nc_dummy != missing_data)))	{
			nc_1019 = c(nc_gfdl[j,i,dates1019],nc_hadgem[j,i,dates1019],nc_ipsl[j,i,dates1019],nc_miroc[j,i,dates1019]) 
			nc_2029 = c(nc_gfdl[j,i,dates2029],nc_hadgem[j,i,dates2029],nc_ipsl[j,i,dates2029],nc_miroc[j,i,dates2029]) 
			nc_3039 = c(nc_gfdl[j,i,dates3039],nc_hadgem[j,i,dates3039],nc_ipsl[j,i,dates3039],nc_miroc[j,i,dates3039]) 
			nc_4049 = c(nc_gfdl[j,i,dates4049],nc_hadgem[j,i,dates4049],nc_ipsl[j,i,dates4049],nc_miroc[j,i,dates4049]) 
			nc_5059 = c(nc_gfdl[j,i,dates5059],nc_hadgem[j,i,dates5059],nc_ipsl[j,i,dates5059],nc_miroc[j,i,dates5059]) 
			nc_6069 = c(nc_gfdl[j,i,dates6069],nc_hadgem[j,i,dates6069],nc_ipsl[j,i,dates6069],nc_miroc[j,i,dates6069]) 			
			nc_7079 = c(nc_gfdl[j,i,dates7079],nc_hadgem[j,i,dates7079],nc_ipsl[j,i,dates7079],nc_miroc[j,i,dates7079]) 			
			nc_8089 = c(nc_gfdl[j,i,dates8089],nc_hadgem[j,i,dates8089],nc_ipsl[j,i,dates8089],nc_miroc[j,i,dates8089]) 		
			nc_9099 = c(nc_gfdl[j,i,dates9099],nc_hadgem[j,i,dates9099],nc_ipsl[j,i,dates9099],nc_miroc[j,i,dates9099]) 			

			
			data1019 = NULL
			for(kh in 1:(length(nc_1019)/12))	{
				new_dat = median(nc_1019[(1:12)+(12*(kh-1))])
				data1019 = c(data1019,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data2029 = NULL
			for(kh in 1:(length(nc_2029)/12))	{
				new_dat = median(nc_2029[(1:12)+(12*(kh-1))])
				data2029 = c(data2029,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data3039 = NULL
			for(kh in 1:(length(nc_3039)/12))	{
				new_dat = median(nc_3039[(1:12)+(12*(kh-1))])
				data3039 = c(data3039,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data4049 = NULL
			for(kh in 1:(length(nc_4049)/12))	{
				new_dat = median(nc_4049[(1:12)+(12*(kh-1))])
				data4049 = c(data4049,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data5059 = NULL
			for(kh in 1:(length(nc_5059)/12))	{
				new_dat = median(nc_5059[(1:12)+(12*(kh-1))])
				data5059 = c(data5059,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data6069 = NULL
			for(kh in 1:(length(nc_6069)/12))	{
				new_dat = median(nc_6069[(1:12)+(12*(kh-1))])
				data6069 = c(data6069,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data7079 = NULL
			for(kh in 1:(length(nc_7079)/12))	{
				new_dat = median(nc_7079[(1:12)+(12*(kh-1))])
				data7079 = c(data7079,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data8089 = NULL
			for(kh in 1:(length(nc_8089)/12))	{
				new_dat = median(nc_8089[(1:12)+(12*(kh-1))])
				data8089 = c(data8089,
					ifelse(new_dat > 0, new_dat, new_dat))
			}
			data9099 = NULL
			for(kh in 1:(length(nc_9099)/12))	{
				new_dat = median(nc_9099[(1:12)+(12*(kh-1))])
				data9099 = c(data9099,
					ifelse(new_dat > 0, new_dat, new_dat))
			}

			
		the_lat = nc_lat[i]
		the_lon = nc_lon[j]

		medOld = median(data1019)
		maxOld = median(data1019)
		dummy_data = sample(1:100,1)
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data1019, c(0.25,0.50,0.75)),
			quantile(data1019, c(0.25,0.50,0.75)) - medOld,
			quantile(data1019, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data2029, c(0.25,0.50,0.75)),
			quantile(data2029, c(0.25,0.50,0.75)) - medOld,
			quantile(data2029, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data3039, c(0.25,0.50,0.75)),
			quantile(data3039, c(0.25,0.50,0.75)) - medOld,
			quantile(data3039, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data4049, c(0.25,0.50,0.75)),
			quantile(data4049, c(0.25,0.50,0.75)) - medOld,
			quantile(data4049, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data5059, c(0.25,0.50,0.75)),
			quantile(data5059, c(0.25,0.50,0.75)) - medOld,
			quantile(data5059, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data6069, c(0.25,0.50,0.75)),
			quantile(data6069, c(0.25,0.50,0.75)) - medOld,
			quantile(data6069, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data7079, c(0.25,0.50,0.75)),
			quantile(data7079, c(0.25,0.50,0.75)) - medOld,
			quantile(data7079, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data8089, c(0.25,0.50,0.75)),
			quantile(data8089, c(0.25,0.50,0.75)) - medOld,
			quantile(data8089, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		quantile_data[iter,] = c(the_lat, the_lon, quantile(data9099, c(0.25,0.50,0.75)),
			quantile(data9099, c(0.25,0.50,0.75)) - medOld,
			quantile(data9099, c(0.25,0.50,0.75)) / maxOld,
			dummy_data)	; iter = iter+1
		
		if(nrow(quantile_data) > 20000)	{
			save_iter = save_iter + 1
			fwrite(quantile_data, paste0(ncpath, "temp_out_", save_iter, ".csv"))
			quantile_data = data.frame(Lat = NA, Lon = NA, Q_25_abs = NA, Q_50_abs = NA, Q_75_abs = NA, 
				Q_25_rel = NA, Q_50_rel = NA, Q_75_rel = NA,
				Q_25_indx = NA, Q_50_indx = NA, Q_75_indx = NA,
				Likelihood = NA)
			iter = 1
			print(c(i,j))
			print(c(the_lat, the_lon))
			} 
		}
	}
#	print(ftrFldIntrv_df)
}	
save_iter = save_iter + 1
fwrite(quantile_data, paste0(ncpath, "temp_out_", save_iter, ".csv"))


all_data = fread(paste0(ncpath, "temp_out_", 1, ".csv"))
for(files in 2:save_iter)	{
	all_data = rbind(all_data, fread(paste0(ncpath, "temp_out_", files, ".csv")))
}

all_data$Hazard = "Groundwater (mm)"
#all_data$Evaluation_year = 2022
all_data$Decade = rep(c('2010s','2020s','2030s','2040s','2050s','2060s','2070s','2080s','2090s'),
	nrow(all_data) / 9)
#all_data$Variable = rep(c(rep("Burnt Area Fraction (Total Annual)",9), rep("Change in Burnt Area Fraction (Relative to 2010s)",8)),
#	nrow(all_data) / 9)
#all_data$CMIP_Models = "GFDL-esm2m, hadgem2-es, ipsl-cm5a, miroc5"
#all_data$LSM_Models = 'CLM4.5'
#all_data$Correction = "ewembi"
all_data$Climate_Scenario = 'rcp2.6'

	# associating each point with a region
library(sf)
regions_shp = st_read("J:\\Cai_data\\TCFD\\referenceRegions.dbf")[-c(28,29),-c(1,3)] # rows 28 and 29 are the Arctic and Antarctica, which overlap with other regions
allData_sf = st_as_sf(all_data, coords = c("Lon", "Lat"), 	# convert foreign object ot object class sf
	crs = st_crs(4326)) 										# EPSG for polar north is 3995, for north america is 2136

allData_Region = st_join(allData_sf, regions_shp, join = st_intersects)
allData_Region$LAB[is.na(allData_Region$LAB)] = 'ARC' # recreating the arctic region w/o overlaps
summary(allData_Region$LAB)

all_data$Region = allData_Region$LAB

fwrite(all_data, paste0(ncpath, "groundwater_", "rcp2.6", ".csv"))


	# combining the data
all_data = rbind(fread(paste0(ncpath, "groundwater_", "rcp6.0", ".csv")), fread(paste0(ncpath, "groundwater_", "rcp2.6", ".csv")))
IDcol = data.frame(ID = 1:nrow(all_data))
fwrite(cbind(IDcol, all_data), paste0(ncpath, "groundwater_", "allData_30MAY2022", ".csv"))

all_data = fread(paste0(ncpath, "groundwater_", "allData_30MAY2022", ".csv"))

	# indexes were saved as char, converting to numeric
all_data$Q_25_indx = as.numeric(all_data$Q_25_indx)
all_data$Q_50_indx = as.numeric(all_data$Q_50_indx)
all_data$Q_75_indx = as.numeric(all_data$Q_75_indx)

	# limiting the range of some large values
all_data$Q_25_abs = ifelse(all_data$Q_25_abs > 0, all_data$Q_25_abs , 0)
all_data$Q_50_abs = ifelse(all_data$Q_50_abs > 0, all_data$Q_50_abs , 0)
all_data$Q_75_abs = ifelse(all_data$Q_75_abs > 0, all_data$Q_75_abs , 0)

all_data$Q_25_rel = ifelse(all_data$Q_25_rel > 5000, 5000, all_data$Q_25_rel)
all_data$Q_50_rel = ifelse(all_data$Q_50_rel > 5000, 5000, all_data$Q_50_rel)
all_data$Q_75_rel = ifelse(all_data$Q_75_rel > 5000, 5000, all_data$Q_75_rel)
all_data$Q_25_rel = ifelse(all_data$Q_25_rel < -5000, -5000, all_data$Q_25_rel)
all_data$Q_50_rel = ifelse(all_data$Q_50_rel < -5000, -5000, all_data$Q_50_rel)
all_data$Q_75_rel = ifelse(all_data$Q_75_rel < -5000, -5000, all_data$Q_75_rel)


fwrite(all_data, paste0(ncpath, "groundwater_nonegs_", "allData_30MAY2022", ".csv"))


bchars_sf1 = subset(all_data, Decade == "2090s" & Climate_Scenario == 'rcp6.0') 
bchars_sf2 = subset(all_data, Decade == "2090s" & Climate_Scenario == 'rcp2.6') 
bchars_sf1$scenario_diff = bchars_sf1$Q_50_abs - bchars_sf2$Q_50_abs 



































###########################################################################################
##### Plotting Data
###########################################################################################
setwd("J:\\Users\\arik\\Documents\\PhD Research\\D2\\images")
library(viridis)

library("ggplot2")
library("sf")
library("rnaturalearth") # for base map
library("rnaturalearthdata") # for downlading lakes and rivers data
library("ggspatial") # for scale bar

lakes <- ne_download(scale = "medium", type = 'lakes', category = 'physical', returnclass="sf")
rivers <- ne_download(scale = "medium", type = 'rivers_lake_centerlines', category = 'physical', returnclass='sf')
world <- ne_countries(scale = "medium", returnclass = "sf")

allData_sf = cbind(st_as_sf(all_data, coords = c("Lon", "Lat"), 	# convert foreign object ot object class sf
	crs = st_crs(4326)), all_data[,c("Lat","Lon")]) 										# EPSG for polar north is 3995, for north america is 2136

###########################################################################################
##### basic spatial plotter for 'squinty eye' testing the data for errors
###########################################################################################

	#########
	## basic plotter
bchars_sf$pos = 0
bchars_sf$pos[which(bchars_sf$D50_day_slp > 0)] = 2
bchars_sf$pos[which(bchars_sf$D50_day_slp < 0)] = 1
#bchars_sf$pos[which(bchars_sf$rnmlt_by_D5_slp < 0)] = NA
not_signif = unique(c(which(bchars_sf$D50_day_pval > .1), which(is.na(bchars_sf$D50_day_slp))))
bchars_sf$pos[not_signif] = 0
#bchars_sf$pos[bchars_sf$REGULATED] = NA


bchars_sf = subset(allData_sf, Decade == "2050s" & Climate_Scenario == 'rcp6.0')# &
#	Lat > -10 & Lat < 80 & Lon > 0 & Lon < 57) 
bchars_sf$Q_50_rel[which(bchars_sf$Q_50_rel < -50)] = -50
bchars_sf$Q_50_rel[which(bchars_sf$Q_50_rel > 50)] = 50
bchars_sf$Q_50_rel[which(bchars_sf$Q_50_rel <= 1 & bchars_sf$Q_50_rel >= -1)] = 0
bchars_sf$Q_50_rel[which(bchars_sf$Q_50_rel > 1 & bchars_sf$Q_50_rel < 50)] = 10
bchars_sf$Q_50_rel[which(bchars_sf$Q_50_rel < -1 & bchars_sf$Q_50_rel > -50)] = -10

ggplot(data = world) +
geom_sf(data = bchars_sf, size=2.2, shape=15, aes(color=as.factor(Q_50_abs))) +
geom_sf(fill='grey92', color="grey92", alpha=0.1) +
geom_sf(data=rivers, color='white', alpha=0.1) +
geom_sf(data=lakes, fill='white', color='grey92', alpha=0.1) +
#geom_sf(data = bchars_aqs[-1,], shape=16, alpha=.7, aes(stroke=2.5, color=Ksat_Wet)) +
#geom_point(data=legend_points, color='grey30', shape=16, aes(sdex[1], sdwy[1], stroke=(1))) +
#geom_text(data=legend_points, color='grey30', fontface='bold', size=8, aes(sdex[3],sdwy[1]+170000,label="change per decade")) +
#geom_rect(data=legend_points, fill=NA, color='grey20',
#	aes(xmin=sdex[1]-180000, xmax=sdex[3]+119000, ymin=sdwy[1]-225000, ymax=sdwy[3]+118000)) +
#scale_color_manual(name="", aesthetics=c("colour", "fill"),
#	values=c('red3','pink',"grey70",'skyblue1',"blue3")) +
#scale_alpha(aesthetics = "alpha", values = c(1,0.3,0.2,0.3,1)) +
scale_color_viridis(name="GW",aesthetics="colour")	+
#coord_sf(crs = st_crs(2163), xlim = c(-3100000, 3400000), ylim = c(-4000000, 3800000),expand=FALSE) + # for NorthAm
#coord_sf(crs = st_crs(3995), xlim = c(5000000, -6000000), ylim = c(-5000000, 5000000),expand=FALSE) + # for Arctic
#coord_sf(crs = st_crs(3857), xlim = c(2000000, 5000000), ylim = c(0, 4200000),expand=FALSE) + # for Egypt
coord_sf(crs = st_crs(3857), xlim = c(-1200000, 5000000), ylim = c(3000000,6000000),expand=FALSE) + # for Med Basin
annotation_scale(location = "bl", width_hint = 0.25) +
#annotation_north_arrow(location = "bl", which_north = "true", 
#	pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"),
#	style = north_arrow_fancy_orienteering) +
theme(legend.background=element_blank(), 
	legend.position=c(1.2,.5),
	legend.title=element_text(face='bold', size=26),
	legend.text=element_text(size=24, face='bold'),
	panel.background=NULL ,
	panel.grid=NULL , panel.grid.major=element_line(color="grey92"), panel.grid.minor=NULL,
	panel.border=element_rect(fill=NA),
	axis.title=element_blank(),
	axis.text=element_text(size=20)) 




















######################
# old version

########################################################
### library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)
	sf::sf_use_s2(FALSE) # for problem with intersecting spherical w flat
library(lubridate)
library(ddplyr) # for left_join() merging by date
library(ncdf4)
library(magrittr)
library(maptools)



#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\TCFD\\BurntArea\\"
	# reading ssp126 climate data 
ncname = 'clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_burntarea_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_burntArea = ncvar_get(ncin,"burntarea")	# lon, lat, time




old_dates = 1:which(year(nc_date) == 2019)[12]
future_dates2029 = which(year(nc_date) == 2020)[1]:which(year(nc_date) == 2029)[12]
future_dates3039 = which(year(nc_date) == 2030)[1]:which(year(nc_date) == 2039)[12]
future_dates4049 = which(year(nc_date) == 2040)[1]:which(year(nc_date) == 2049)[12]
future_dates5059 = which(year(nc_date) == 2050)[1]:which(year(nc_date) == 2059)[12]
future_dates6069 = which(year(nc_date) == 2060)[1]:which(year(nc_date) == 2069)[12]
future_dates7079 = which(year(nc_date) == 2070)[1]:which(year(nc_date) == 2079)[12]
future_dates8089 = which(year(nc_date) == 2080)[1]:which(year(nc_date) == 2089)[12]
future_dates9099 = which(year(nc_date) == 2090)[1]:which(year(nc_date) == 2099)[12]



burntAreaChng_2029 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_3039 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_4049 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_5059 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_6069 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_7079 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_8089 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
burntAreaChng_9099 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))


for(i in 1:length(nc_lat))	{
	for(j in 1:length(nc_lon))	{
		nc_burntAreatot_old =nc_burntArea[j,i,old_dates] #ncvar_get(ncin,"maxdis")[j,i,]	# in lon,lat,time
		if(any(!is.na(nc_burntAreatot_old)))	{
			nc_burntAreatot_f_2029 = nc_burntArea[j,i,future_dates2029] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_burntAreatot_f_3039 = nc_burntArea[j,i,future_dates3039] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_burntAreatot_f_4049 = nc_burntArea[j,i,future_dates4049] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_burntAreatot_f_5059 = nc_burntArea[j,i,future_dates5059] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_burntAreatot_f_6069 = nc_burntArea[j,i,future_dates6069] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_burntAreatot_f_7079 = nc_burntArea[j,i,future_dates7079] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_burntAreatot_f_8089 = nc_burntArea[j,i,future_dates8089] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_burntAreatot_f_9099 = nc_burntArea[j,i,future_dates9099] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			

			print(c(i,j))

			nc_maxburntArea_y = NULL
			for(kh in 1:(length(nc_burntAreatot_old)/12))	{
				nc_maxburntArea_y = c(nc_maxburntArea_y, sum(nc_burntAreatot_old[(1:12)+(12*(kh-1))]))
			}

			nc_burntAreatot_y_f_2029 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_2029)/12))	{
				nc_burntAreatot_y_f_2029 = c(nc_burntAreatot_y_f_2029, sum(nc_burntAreatot_f_2029[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_3039 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_3039)/12))	{
				nc_burntAreatot_y_f_3039 = c(nc_burntAreatot_y_f_3039, sum(nc_burntAreatot_f_3039[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_4049 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_4049)/12))	{
				nc_burntAreatot_y_f_4049 = c(nc_burntAreatot_y_f_4049, sum(nc_burntAreatot_f_4049[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_5059 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_5059)/12))	{
				nc_burntAreatot_y_f_5059 = c(nc_burntAreatot_y_f_5059, sum(nc_burntAreatot_f_5059[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_6069 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_6069)/12))	{
				nc_burntAreatot_y_f_6069 = c(nc_burntAreatot_y_f_6069, sum(nc_burntAreatot_f_6069[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_7079 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_7079)/12))	{
				nc_burntAreatot_y_f_7079 = c(nc_burntAreatot_y_f_7079, sum(nc_burntAreatot_f_7079[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_8089 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_8089)/12))	{
				nc_burntAreatot_y_f_8089 = c(nc_burntAreatot_y_f_8089, sum(nc_burntAreatot_f_8089[(1:12)+(12*(kh-1))]))
			}
			nc_burntAreatot_y_f_9099 = NULL
			for(kh in 1:(length(nc_burntAreatot_f_9099)/12))	{
				nc_burntAreatot_y_f_9099 = c(nc_burntAreatot_y_f_9099, sum(nc_burntAreatot_f_9099[(1:12)+(12*(kh-1))]))
			}


		burntAreaChng_2029[i,j] = median(nc_burntAreatot_y_f_2029) - median(nc_maxburntArea_y) 
		burntAreaChng_3039[i,j] = median(nc_burntAreatot_y_f_3039) - median(nc_maxburntArea_y) 
		burntAreaChng_4049[i,j] = median(nc_burntAreatot_y_f_4049) - median(nc_maxburntArea_y)  
		burntAreaChng_5059[i,j] = median(nc_burntAreatot_y_f_5059) - median(nc_maxburntArea_y)  
		burntAreaChng_6069[i,j] = median(nc_burntAreatot_y_f_6069) - median(nc_maxburntArea_y)  
		burntAreaChng_7079[i,j] = median(nc_burntAreatot_y_f_7079) - median(nc_maxburntArea_y) 
		burntAreaChng_8089[i,j] = median(nc_burntAreatot_y_f_8089) - median(nc_maxburntArea_y) 
		burntAreaChng_9099[i,j] = median(nc_burntAreatot_y_f_9099) - median(nc_maxburntArea_y)  
		
		
			
		} else { 
		burntAreaChng_2029[i,j] = NA
		burntAreaChng_3039[i,j] = NA
		burntAreaChng_4049[i,j] = NA
		burntAreaChng_5059[i,j] = NA 
		burntAreaChng_6069[i,j] = NA
		burntAreaChng_7079[i,j] = NA
		burntAreaChng_8089[i,j] = NA 
		burntAreaChng_9099[i,j] = NA 
			}
	}
#	print(ftrFldIntrv_df)
}		
write.csv(burntAreaChng_2029, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_2029.csv")
write.csv(burntAreaChng_3039, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_3039.csv")
write.csv(burntAreaChng_4049, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_4049.csv")
write.csv(burntAreaChng_5059, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_5059.csv")
write.csv(burntAreaChng_6069, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_6069.csv")
write.csv(burntAreaChng_7079, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_7079.csv")
write.csv(burntAreaChng_8089, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_8089.csv")
write.csv(burntAreaChng_9099, "J:\\Cai_data\\TCFD\\BurntArea\\BurntAreachng_CLM45_rcp26_9099.csv")


#ftrFldIntrv_df = as.data.frame(read.csv("J:\\Cai_data\\TCFD\\Flash Floods\\test_out.csv"))[,-1]


library(colorRamps)
col5 <- colorRampPalette(c('red4', 'red3', 'red3', 'red2', 'red1', 'grey90', 'blue1', 'blue2', 'blue3',
	'blue3', 'blue4'))  #create color ramp starting from blue to red

par(mar=c(2.1,2.1,2.1,1))
image(nc_lon,nc_lat,t(as.matrix(burntAreaChng_9099[nrow(burntAreaChng_9099):1,])), ylim = c(-57,85),
	col=rev(col5(n=11)), breaks = c(-1,-0.1,-0.01,-0.001,-0.0002,-0.0001,0.0001,0.0002,0.001,0.01,0.1,1),#seq(-.02,.02,length.out=12),
	ylab='Lat', xlab = 'lon', main = "Changes Wildfire Burnt Area (2090s - Prev. 20Yrs) RCP6.0")
data(wrld_simpl)
plot(wrld_simpl,add=TRUE)

# models available in CLM4.5, CARAIB, LPJ-GUESS, LPJmL, ORCHIDEE, ORCHIDEE-DGVM, VISIT


if (any(is.na(as.matrix(ftrFldIntrv_df[,ncol(ftrFldIntrv_df):1])))) {
    mmat <- ifelse(is.na(as.matrix(ftrFldIntrv_df[,ncol(ftrFldIntrv_df):1])), 1, NA)
    image(mmat, axes = FALSE, xlab = "", ylab = "", col = 'white', useRaster=TRUE, add = TRUE)
  }

ftrFldIntrv_mt = as.matrix(unlist(ftrFldIntrv_df))

pt1
ftrFldIntrv_sf = st_as_sf(as.vector(ftrFldIntrv_df), coords = 

















ftrFldIntrv_df[is.na(ftrFldIntrv_df)] = 0

# converting df into array 
lon2 = nc_lon	; lat2 = nc_lat
nlon2 = length(nc_lon)	; nlat2 = length(nc_lat)
#time2 = time; tunits2 = tunits	; nt2 = nt

# convert to 2d array
ftrFldIntrv_ar <- array(ftrFldIntrv_df, dim=c(nlon2,nlat2))
dim(ftrFldIntrv_ar)

# plot to check creation of arrays
library(lattice)
library(RColorBrewer)

grid <- expand.grid(lon=nc_lon, lat=nc_lat)

levelplot(ftrFldIntrv_ar ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
  col.regions=(rev(brewer.pal(10,"RdBu"))), main="MTWA (C)")



	#names of path, file, and dname
ncpath = "J:\\Cai_data\\TCFD\\Flash Floods\\test_out.csv"
ncname = "GLM_ftrFlood"
ncfname = paste(ncpath, ncname, '.nc', sep='')
dname = 'fldInt'


	#define dimensions
londim = ncdim_def('lon','degrees',as.double(nc_lon))
latdim = ncdim_def('lat','degrees',as.double(nc_lat))
#timedim <- ncdim_def("time",tunits3,as.double(time3))

# define variables
fillvalue = NA
dlname = "Flood Interval in Years"
fldInt_def = ncvar_def("fldInt","yrs",list(londim,latdim),fillvalue,dlname,prec="single")
#dlname = "mean_temperture_warmest_month"
#mtwa.def <- ncvar_def("mtwa","deg_C",list(londim,latdim,timedim),fillvalue,dlname,prec="single")


# create netCDF file and put arrays
#ncout <- nc_create(ncfname,list(tmp_def,mtco.def,mtwa.def,mat.def),force_v4=TRUE)
ncout <- nc_create(ncfname,list(fldInt_def),force_v4=FALSE)

# put variables
#ncvar_put(ncout,tmp_def,tmp_array3)
#ncvar_put(ncout,mtwa.def,mtwa_array3)
ncvar_put(ncout,fldInt_def,ftrFldIntrv_df)
ncvar_put(ncout,mat.def,mat_array3)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

# add global attributes
ncatt_put(ncout,0,"title",title$value)
ncatt_put(ncout,0,"institution",institution$value)
ncatt_put(ncout,0,"source",datasource$value)
ncatt_put(ncout,0,"references",references$value)
history <- paste("P.J. Bartlein", date(), sep=", ")
ncatt_put(ncout,0,"history",history)
ncatt_put(ncout,0,"Conventions",Conventions$value)

# Get a summary of the created file:
ncout





dim1 <- dim.def.ncdf( name = "lat", units = "degrees", vals = as.double(nc_lat) )
dim2 <- dim.def.ncdf( name = "lon", units = "degrees", vals = as.double(nc_lon) )
ftrFldIntrv_mt = as.matrix(ftrFldIntrv_df)

# define the EMPTY (elevation) netcdf variable

varz <- var.def.ncdf(name = "floodFreq", units =  "years", 
                    dim = list(dim1, dim2), missval = NA, 
                    longname ="Flood Frequency of what was a 20 year flood")		

nc.ex <- create.ncdf( filename = "J:\\Cai_data\\TCFD\\Flash Floods\\test_out.nc", vars = varz )
put.var.ncdf(nc = nc.ex, varid = varz, vals = ftrFldIntrv_mt)




ggplot(aes(x=X1, y=X2, fill=value), data=mpr) + geom_raster() 
		
myClrs =  colorRamp(c("red", "blue"))
plot(x=nc_lon, y=nc_lat, 



lon1 <- ncdim_def( "Lon", "degrees_east", unique(DB$lon))
lat1 <- ncdim_def( "Lat", "degrees_north", unique(DB$lat))
record <- ncdim_def( "record", "files", list_obs$record_uid,unlim=T)
mv <- -999 # missing value to use
tmpmean_ <- ncvar_def( "tmpmean_", "degrees", list(lon1,lat1,record),mv)

time<-list_obs$record_uid

ncnew <- nc_create( "tmpmean_.nc",tmpmean_)

for( i in 1:length(record))
  ncvar_put( ncnew, tmpmean_, DB$tmpmean_, start=c(1,1,i), 
             count=c(-1,-1,1))
		
		
		
		
		

cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))

#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(thisClmMod in c(4,6:15))	{
	HBVoutCol = HBVoutCol + 1

#	ncvar_thismod_ppt =  ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] ugh annoying
	ncvar_thismod_ppt_1 = ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] and now merging two separate netcdfs
	ncvar_thismod_ppt_2 = ncvar_get(ncin_ppt_2, "tp")[date_range_ppt2,thisClmMod , , ]
	ncvar_thismod_tmin = ncvar_get(ncin, "t2m_min")[, , , thisClmMod]
	ncvar_thismod_tmax = ncvar_get(ncin, "t2m_max")[, , , thisClmMod]



