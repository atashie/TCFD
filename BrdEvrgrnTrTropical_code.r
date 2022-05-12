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
ncpath = "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\"
	# reading ssp126 climate data 
ncname = 'caraib_gfdl-esm2m_ewembi_rcp60_2005soc_co2_npp-brevtrt_global_annual_2006_2099.nc4'#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_burntarea_global_monthly_2006_2099.nc4"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = rev(ncvar_get(ncin, 'lat'))	# lat is given from high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 365.25# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_npp = ncvar_get(ncin,"npp-brevtrt")	# lon, lat, time




old_dates = 1:which(year(nc_date) == 2019)
future_dates2029 = which(year(nc_date) == 2020):which(year(nc_date) == 2029)
future_dates3039 = which(year(nc_date) == 2030):which(year(nc_date) == 2039)
future_dates4049 = which(year(nc_date) == 2040):which(year(nc_date) == 2049)
future_dates5059 = which(year(nc_date) == 2050):which(year(nc_date) == 2059)
future_dates6069 = which(year(nc_date) == 2060):which(year(nc_date) == 2069)
future_dates7079 = which(year(nc_date) == 2070):which(year(nc_date) == 2079)
future_dates8089 = which(year(nc_date) == 2080):which(year(nc_date) == 2089)
future_dates9099 = which(year(nc_date) == 2090):which(year(nc_date) == 2099)



nppChng_2029 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_3039 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_4049 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_5059 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_6069 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_7079 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_8089 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))
nppChng_9099 = data.frame(matrix(-999999, nrow=length(nc_lat), ncol=length(nc_lon)))


for(i in 1:length(nc_lat))	{
	for(j in 1:length(nc_lon))	{
		nc_npptot_old = median(nc_npp[j,i,old_dates]) #ncvar_get(ncin,"maxdis")[j,i,]	# in lon,lat,time
		if(any(!is.na(nc_npptot_old)))	{
			nc_npptot_f_2029 = nc_npp[j,i,future_dates2029] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_npptot_f_3039 = nc_npp[j,i,future_dates3039] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_npptot_f_4049 = nc_npp[j,i,future_dates4049] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_npptot_f_5059 = nc_npp[j,i,future_dates5059] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_npptot_f_6069 = nc_npp[j,i,future_dates6069] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_npptot_f_7079 = nc_npp[j,i,future_dates7079] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_npptot_f_8089 = nc_npp[j,i,future_dates8089] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			nc_npptot_f_9099 = nc_npp[j,i,future_dates9099] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			

			print(c(i,j))

			nppChng_2029[i,j] = median(nc_npptot_f_2029) - nc_npptot_old 
			nppChng_3039[i,j] = median(nc_npptot_f_3039) - nc_npptot_old 
			nppChng_4049[i,j] = median(nc_npptot_f_4049) - nc_npptot_old  
			nppChng_5059[i,j] = median(nc_npptot_f_5059) - nc_npptot_old  
			nppChng_6069[i,j] = median(nc_npptot_f_6069) - nc_npptot_old  
			nppChng_7079[i,j] = median(nc_npptot_f_7079) - nc_npptot_old 
			nppChng_8089[i,j] = median(nc_npptot_f_8089) - nc_npptot_old 
			nppChng_9099[i,j] = median(nc_npptot_f_9099) - nc_npptot_old  
			
			
				
			} else { 
			nppChng_2029[i,j] = NA
			nppChng_3039[i,j] = NA
			nppChng_4049[i,j] = NA
			nppChng_5059[i,j] = NA 
			nppChng_6069[i,j] = NA
			nppChng_7079[i,j] = NA
			nppChng_8089[i,j] = NA 
			nppChng_9099[i,j] = NA 
				}
		}
#	print(ftrFldIntrv_df)
}		
write.csv(nppChng_2029, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_2029.csv")
write.csv(nppChng_3039, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_3039.csv")
write.csv(nppChng_4049, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_4049.csv")
write.csv(nppChng_5059, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_5059.csv")
write.csv(nppChng_6069, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_6069.csv")
write.csv(nppChng_7079, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_7079.csv")
write.csv(nppChng_8089, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_8089.csv")
write.csv(nppChng_9099, "J:\\Cai_data\\TCFD\\BrdEvrgrnTrTropical\\nppchng_CARAIB_rcp26_9099.csv")


#ftrFldIntrv_df = as.data.frame(read.csv("J:\\Cai_data\\TCFD\\Flash Floods\\test_out.csv"))[,-1]


library(colorRamps)
col5 <- colorRampPalette(c('red4', 'red3', 'red3', 'red2', 'red1', 'grey90', 'green1', 'green2', 'green3',
	'green3', 'green4'))  #create color ramp starting from green to red

par(mar=c(2.1,2.1,2.1,1))
image(nc_lon,nc_lat,t(as.matrix(nppChng_9099[nrow(nppChng_9099):1,]))*10000000, ylim = c(-57,85),
	col=col5(n=11), breaks = c(-1,-0.1,-0.01,-0.001,-0.0002,-0.0001,0.0001,0.0002,0.001,0.01,0.1,1),#seq(-.02,.02,length.out=12),
	ylab='Lat', xlab = 'lon', main = "Changes Productivity of Broadleaf Evergreen Tropical Forests (2090s - Prev. 20Yrs) RCP6.0")
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



