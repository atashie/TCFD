
##################################################
#### reading basinATLAS data

##########################################################################################
##	this section can be skipped after having been run and saved
# bchars_sf is at gages, bchars_cent_sf is centroids
							#setwd("J:\\Users\\arik\\Documents\\PhD Research\\D4\\hybas_af_lev00_v1c")
							#basinAt_2 = st_read("hybas_af_lev00_v1c.shp")
							setwd("C:\\Users\\arik\\Documents\\PhD Research\\D4\\BasinATLAS_Data_v10")
							st_layers("BasinATLAS_v10.gdb")
							basinAt12 = st_read(dsn="BasinATLAS_v10.gdb", layer="BasinATLAS_v10_lev12")
							sf::sf_use_s2(FALSE)

							#basinAt2 = st_read(dsn="BasinATLAS_v10.gdb", layer="BasinATLAS_v10_lev02")#




							# identifying all watersheds that feed / flow into the outlet catchment
############################# takes a long time to run, so file is saved for future use
								NorAm_box = st_as_sf(st_sfc(st_polygon(list(matrix(c(-52,-52,-170,-170,-130,-52,15,90,90,50,15,15),ncol=2))), crs = st_crs(4326)))
								basinAt_NorAm_cent = st_join(basinAt_cent, NorAm_box, join=st_within, left=FALSE)

								basinAt_prj = st_transform(basinAt12, crs = st_crs(4326))
								basinAt_NorAm_polys = st_join(basinAt_prj, NorAm_box, join=st_within, left=FALSE)
								st_write(basinAt_NorAm_polys, 	"C:\\Users\\arik\\Documents\\Postdoc Research\\Arctic_baseflow\\basinAt_NorAm_polys_Meteorology.gpkg")
#############################
							
							basinAt_NorAm_polys = st_read("C:\\Users\\arik\\Documents\\Postdoc Research\\Arctic_baseflow\\basinAt_NorAm_polys_Meteorology.gpkg")
								# begin merging our basin-level trends data with the HydroBasins dataset
							bchars_exdat = st_intersection(bchars_sf, st_buffer(basinAt_NorAm_polys,0)) # keeps points
							st_geometry(bchars_exdat) = NULL											# removing geometry for merging

							bchars_poly = st_transform(st_as_sf(merge(bchars_exdat, basinAt_NorAm_polys, all.x=TRUE)), 4326)



								# freeing memory
								rm(basinAt12)	;	rm(basinAt_prj)



							HB_strip = basinAt_NorAm_polys	; st_geometry(HB_strip) = NULL	; HB_strip = as.data.frame(HB_strip)
							bchars_polys_strip = bchars_poly ; st_geometry(bchars_polys_strip) = NULL ; bchars_polys_strip = as.data.frame(bchars_polys_strip)


							gaged_upstreams =  readRDS("C:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_n.RData")
							for(gaged_row in 1:nrow(bchars_polys_strip))	{
								gaged_ID = bchars_polys_strip$HYBAS_ID[gaged_row]
								if(gaged_ID %in% names(gaged_upstreams))	{
									print(gaged_ID)
								}	else	{
									
									HB_remaining = HB_strip[-which(HB_strip$HYBAS_ID == gaged_ID),]
									print(nrow(HB_remaining))
										
									while(any(gaged_ID %in% HB_remaining$NEXT_DOWN))	{
										for(trueys in which(gaged_ID %in% HB_remaining$NEXT_DOWN))	{
											new_HB_ID_rows = which(HB_remaining$NEXT_DOWN == gaged_ID[trueys])
											new_HB_ID = as.character(HB_remaining$HYBAS_ID[new_HB_ID_rows])
											gaged_ID = c(gaged_ID, new_HB_ID)
											print(gaged_ID)
											HB_remaining = HB_remaining[-new_HB_ID_rows,]
											print(nrow(HB_remaining))
										}
									}
									gaged_upstreams[[gaged_row]] = gaged_ID
									saveRDS(gaged_upstreams, file="C:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_nn.RData")
								}
							}
								#gaged_upstreams = readRDS("J:\\Users\\arik\\Documents\\PhD Research\\D4\\HB_upstreams_nn.RData")
	#						names(gaged_upstreams) = bchars_polys_strip$HYBAS_ID

								# freeing memory
								rm(HB_strip)	
								
								# creating new polygons based on merging all the upstream hydrobasins
							this_basin_ids = which(basinAt_NorAm_polys$HYBAS_ID %in% gaged_upstreams[[1]])
							these_catchments = basinAt_NorAm_polys[this_basin_ids,]
							all_gaged_basins = st_sf(st_union(these_catchments))
							for(pp in 2:length(gaged_upstreams))	{
								this_basin_ids = which(basinAt_NorAm_polys$HYBAS_ID %in% gaged_upstreams[[pp]])
								these_catchments = basinAt_NorAm_polys[this_basin_ids,]
								this_basin = st_sf(st_union(these_catchments))
								all_gaged_basins = rbind(all_gaged_basins, this_basin)
							}
							all_gaged_basins$HYBAS_ID = names(gaged_upstreams)
							bchars_poly_basins = st_sf(merge(bchars_polys_strip, all_gaged_basins))
							#st_write(bchars_poly_basins, "C:\\Users\\arik\\Documents\\Postdoc Research\\bdyko_trends\\bchars_poly_basins_Meteorology.gpkg")
	
##	this section can be skipped after having been run and saved
# bchars_sf is at gages, bchars_cent_sf is centroids
##########################################################################################
