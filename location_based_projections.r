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

this_lat = which.min(abs(nc_lat - -24.900833013921915))
this_lon = which.min(abs(nc_lon - 152.41212621076588))

this_lat = which.min(abs(nc_lat - -35.7))
this_lon = which.min(abs(nc_lon - 146))

this_lat = which.min(abs(nc_lat - -13.7))
this_lon = which.min(abs(nc_lon - 142.7))

this_lat = which.min(abs(nc_lat - -32.2))
this_lon = which.min(abs(nc_lon - 116.0))

this_lat = which.min(abs(nc_lat - -34.5))
this_lon = which.min(abs(nc_lon - 138.7))


	# reading ssp126 climate data 
ncname = 'watergap2-2c_gfdl-esm2m_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))

	#creating a database for storing all initial data
last_year = last(which(year(nc_date) == 2049))
first_year = which(year(nc_date) == 2010)[1]
nc_gwStor = data.frame(matrix(-999999, nrow=length(first_year:last_year), ncol=4))
names(nc_gwStor) = c("gfdl","hadgem2","ipsl","miroc5")

nc_gwStor$gfdl = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_hadgem2-es_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$hadgem2 = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_ipsl-cm5a-lr_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$ipsl = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_miroc5_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$miroc5 = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year] # lon, lat, time
nc_gwStor$year = rep(2010:2049, each=12)

gwStor_df = data.frame(matrix(-999999, nrow=(length(2010:2049)*4), ncol=5))
names(gwStor_df) = c("year","minStor","maxStor","meanStor","model")
gwStor_df$year = rep(2010:2049, 4)
gwStor_df$decade = rep(rep(c("2010s","2020s",'2030s','2040s'), each=10), 4)#,'2050s'), each=10)
gwStor_df$model = rep(c("gfdl","hadgem2","ipsl","miroc5"), each=40)

iter = 0
for(this_model in unique(gwStor_df$model))	{
	this_gw = nc_gwStor[,which(names(nc_gwStor)==this_model)]
	
	for(kh in unique(gwStor_df$year))	{
		iter=iter+1
		this_year = which(nc_gwStor$year == kh)
		gwStor_df$minStor[iter] = min(this_gw[this_year]) + 100
		gwStor_df$maxStor[iter] = max(this_gw[this_year])+ 100
		gwStor_df$meanStor[iter] = mean(this_gw[this_year])+ 100
	}
}

p = ggplot(gwStor_df, aes(factor(decade),meanStor))
p + geom_boxplot()

plot(subset(gwStor_df,model=='hadgem2')$minStor)

gwStor_df_data = c(100*(gwStor_df$minStor - median(subset(gwStor_df, decade == '2010s')$minStor)) /  median(subset(gwStor_df, decade == '2010s')$minStor),
	100*(gwStor_df$meanStor - median(subset(gwStor_df, decade == '2010s')$meanStor)) /  median(subset(gwStor_df, decade == '2010s')$meanStor),
	100*(gwStor_df$maxStor - median(subset(gwStor_df, decade == '2010s')$maxStor)) /  median(subset(gwStor_df, decade == '2010s')$maxStor))
box_plotter = data.frame(the_data = gwStor_df_data,
	the_decade = rep(gwStor_df$decade,3),
	the_range = rep(c("Annual Min","Annual Avg","Annual Max"), each=40))

ggplot(box_plotter, aes(the_range, the_data))				+
	geom_boxplot(aes(fill = the_decade))					+
	ylim(-50,50)												+
	scale_fill_viridis(discrete=TRUE, begin=0.3)			+
	geom_hline(yintercept = 0)								+
	labs(y="Relative Change in Groundwater Availability (%)",
		x="")												+													
	ggtitle(paste0("Data for: ",
		round(nc_lat[this_lat],1), ", ",
		round(nc_lon[this_lon],1))) +
	theme(plot.title = element_text(hjust = 0.5),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank())
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ###############################
 # same but for entire period of record

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

this_lat = which.min(abs(nc_lat - -24.900833013921915))
this_lon = which.min(abs(nc_lon - 152.41212621076588))

this_lat = which.min(abs(nc_lat - -35.7))
this_lon = which.min(abs(nc_lon - 146))

this_lat = which.min(abs(nc_lat - -13.7))
this_lon = which.min(abs(nc_lon - 142.7))

this_lat = which.min(abs(nc_lat - -32.2))
this_lon = which.min(abs(nc_lon - 116.0))

this_lat = which.min(abs(nc_lat - -34.5))
this_lon = which.min(abs(nc_lon - 138.7))


	# reading ssp126 climate data 
ncname = 'watergap2-2c_gfdl-esm2m_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))

	#creating a database for storing all initial data
last_year = last(which(year(nc_date) == 2099))
first_year = which(year(nc_date) == 2010)[1]
nc_gwStor = data.frame(matrix(-999999, nrow=length(first_year:last_year), ncol=4))
names(nc_gwStor) = c("gfdl","hadgem2","ipsl","miroc5")

nc_gwStor$gfdl = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_hadgem2-es_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$hadgem2 = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_ipsl-cm5a-lr_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$ipsl = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year]	# lon, lat, time

ncname = 'watergap2-2c_miroc5_ewembi_rcp60_2005soc_co2_groundwstor_global_monthly_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_2005co2_thawdepth_global_annual_2006_2099.nc4"  
	# note that watergap treats gw as a relative value (depth below or above a relative datum) so the actual value can go negative (and theoretically can go to any real negative number)
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'lat')#rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_gwStor$miroc5 = ncvar_get(ncin,"groundwstor")[this_lon,this_lat,first_year:last_year] # lon, lat, time

nc_gwStor$year = rep(2010:2099, each=12)

gwStor_df = data.frame(matrix(NA, nrow=(length(2010:2099)*4), ncol=5))
names(gwStor_df) = c("year","minStor","maxStor","meanStor","model")
gwStor_df$year = rep(2010:2099, 4)
gwStor_df$decade = rep(rep(c("2010s","2020s",'2030s','2040s','2050s','2060s','2070s','2080s','2090s'), each=10), 4)#,'2050s'), each=10)
gwStor_df$model = rep(c("gfdl","hadgem2","ipsl","miroc5"), each=90)

iter = 0
for(this_model in unique(gwStor_df$model))	{
	this_gw = nc_gwStor[,which(names(nc_gwStor)==this_model)]
	
	for(kh in unique(gwStor_df$year))	{
		iter=iter+1
		this_year = which(nc_gwStor$year == kh)
		gwStor_df$minStor[iter] = min(this_gw[this_year]) + 100
		gwStor_df$maxStor[iter] = max(this_gw[this_year])+ 100
		gwStor_df$meanStor[iter] = mean(this_gw[this_year])+ 100
	}
}

p = ggplot(gwStor_df, aes(factor(decade),meanStor))
p + geom_boxplot()

plot(subset(gwStor_df,model=='hadgem2')$minStor)

gwStor_df_data = c(100*(gwStor_df$minStor - median(subset(gwStor_df, decade == '2010s')$minStor)) /  median(subset(gwStor_df, decade == '2010s')$minStor),
	100*(gwStor_df$meanStor - median(subset(gwStor_df, decade == '2010s')$meanStor)) /  median(subset(gwStor_df, decade == '2010s')$meanStor),
	100*(gwStor_df$maxStor - median(subset(gwStor_df, decade == '2010s')$maxStor)) /  median(subset(gwStor_df, decade == '2010s')$maxStor))
box_plotter = data.frame(the_data = gwStor_df_data,
	the_decade = rep(gwStor_df$decade,3),
	the_range = rep(c("Annual Min","Annual Avg","Annual Max"), each=90))

ggplot(box_plotter, aes(the_range, the_data))				+
	geom_boxplot(aes(fill = the_decade))					+
	ylim(-50,50)												+
	scale_fill_viridis(discrete=TRUE, begin=0.3)			+
	geom_hline(yintercept = 0)								+
	labs(y="Relative Change in Groundwater Availability (%)",
		x="")												+													
	ggtitle(paste0("Data for: ",
		round(nc_lat[this_lat],1), ", ",
		round(nc_lon[this_lon],1))) +
	theme(plot.title = element_text(hjust = 0.5),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank())
 