# remotes::install_github("slisovski/migrationSDP")
library(migrationSDP) 
library(dplyr)
library(tibble)
library(tidyr)
library(fGarch)
library(ggplot2)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(geosphere)

source('functions.R')

######DATA######
################

load('site_parameter/eaafMap.rda')

# ggplot() +
#   plot+
#   geom_sf(data = eaafMap$grid) 
#   geom_sf(data = eaafMap$map, fill = adjustcolor('grey40', alpha.f = 0.4)) +
#   geom_sf(data = eaafMap$bbox, fill = NA) +
#   geom_sf(data = crds_sf, aes(geometry = geometry), color = "red") +
#   theme_void()
# 
# crds_sf <- st_as_sf(crds, coords = c("Lon", "Lat"), crs = 4326)

load('site_parameter/mudflatTab.rda')

# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = mudflatTab %>% filter(lake_area>quantile(lake_area, probs = 0.5, na.rm = T)), aes(fill = lake_area), alpha = 0.75) +
#   scale_fill_continuous(type = "viridis",na.value = "transparent")

# siteTab      <- (mudflatTab %>% st_centroid() %>% st_transform(4326) %>%
#                mutate(areaHist = ifelse(is.na(histArea_inner), currArea_outer, histArea_outer),
#                       areaCurr = currArea_outer,
#                       areaMang = mangArea_outer,
#                       lake     = ifelse(!is.na(lake_area) & lake_area > quantile(lake_area, probs = 0.5, na.rm = T), TRUE, FALSE)) %>%
#                dplyr::select(areaHist, areaCurr, areaMang, lake) %>%
#                rownames_to_column(var = "index") %>% mutate(index = as.integer(index))) %>%
#   filter((areaHist>0 | areaCurr>0 | areaMang>0 | lake)) %>%
#   filter(st_coordinates(.)[,1] > 104 | st_coordinates(.)[,1] < -150) %>%
#   relocate(geometry, .after = last_col()) %>%
#   suppressWarnings()
# 
# save(siteTab , file = 'site_parameter/siteTab.rda')
load('site_parameter/siteTab.rda')

# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = siteTab, aes(fill = lake), alpha = 0.75, size = 5, shape = 23) +
#   scale_fill_continuous(type = "viridis", na.value = "transparent")



#### Temp
load('site_parameter/tempTab.rda')

######PLYGONS########
####wint/breed#######


### 2 bboxes
# bbox_crs <-  c(st_bbox(c(xmin =  100, xmax =  180, ymin = -80, ymax = 80), crs = 4326) %>% st_as_sfc(),
#                st_bbox(c(xmin = -180, xmax = -150, ymin = -80, ymax = 80), crs = 4326) %>% st_as_sfc())
# 
# 
# shp <- read_sf('N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/Wader_3.shp') %>%
#   st_intersection(bbox_crs) %>% st_transform(st_crs(eaafMap$map)) %>% st_intersection(eaafMap$bbox)
# # save(shp, file = 'N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/transformed_shp.shp')
# # load('N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/transformed_shp.shp')
# 
spTab <- readxl::read_xlsx('temp/Parameters/spTab2.xlsx') %>% filter(Filter==TRUE)

# 
# start_end_table <- lapplyf(spTab$ScieName, function(sp) {
#   
#   # sp <- "Limosa lapponica"
#   tmp_shp <- shp %>% filter(SCINAME==sp, SEASONA%in%c(3,2))
#   
#   tmp <- mudflatTab %>% rowid_to_column(var = 'id') %>% dplyr::select('id') %>% 
#     mutate(poly = apply(st_intersects(., tmp_shp, sparse = FALSE), 1, function(x) ifelse(any(x), which(x)[1], NA))) %>%
#     filter(!is.na(poly)) %>% left_join(tibble(poly = 1:nrow(tmp_shp), seasonal = tmp_shp$SEASONA), by = 'poly')
#   
#   
#   tibble(sp = sp, id = tmp$id, type = ifelse(tmp$seasonal==3, 'winter', 'breeding'))
#   
# }) %>% Reduce('rbind',.)
# 
# save(start_end_table, file = 'N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/start_end_table.RData')
load('temp/Waders/start_end_table.RData')

# start_end_list <- lapply(unique(start_end_table$sp), function(species) {
# 
#   tmp <- start_end_table %>% filter(sp==species)
#   
#   winterSites <- siteTab[siteTab$index%in%tmp$id[tmp$type=='winter'],] %>% filter(!is.na(index)) %>%
#     mutate(prob = scales:::rescale(log(areaHist + areaMang+0.1), c(0,1))) %>%
#     mutate(prob = ifelse(spTab$Lake_propensity[spTab$ScieName==species] & prob < median(prob, na.rm = T),  median(prob, na.rm = T), prob))
#   
#   # ggplot() +
#   #   geom_sf(data = shp %>% filter(SCINAME== species), mapping = aes(geometry = geometry, fill = as.factor(SEASONAL))) +
#   #   geom_sf(data =  mudflatTab[winterSites$index,] %>% st_centroid() %>% st_geometry(), mapping = aes(size = winterSites$prob), pch = 16, color = 'darkgreen')
#   
#   
#   tibble(breeding_id = sample(tmp %>% filter(type=='breeding') %>% pull(id), spTab$Pop[spTab$ScieName==species], TRUE)) %>%
#       mutate(winter_id = sample(winterSites$index, spTab$Pop[spTab$ScieName==species], TRUE, prob = winterSites$prob)) %>%
#       group_by(winter_id, breeding_id) %>% summarise(n = n())
# 
# })
# names(start_end_list) <- unique(start_end_table$sp)
# save(start_end_list, file = 'temp/Waders/start_end_list.RData')
load('temp/Waders/start_end_list.RData')


#### Sowmelt - Arrival

# breed_polys <- mudflatTab %>% rownames_to_column(var = "index") %>%
#   filter(index %in% (unique(start_end_table %>% filter(type == "breeding") %>% pull(id)))) %>%
#   dplyr::select(index)
# 
# library(stars)
# snow    <- read_stars("site_parameter/nhsce_v01r01_19661004_20220103.nc",
#                    sub = 'snow_cover_extent')
# snowMat <- apply(snow$snow_cover_extent, 3, function(x) c(x))
# 
# dates   <- st_get_dimension_values(snow, which = "time")
# datesD  <- tibble(index = 1:length(dates), year = as.numeric(format(dates, "%Y")),
#                   doy = as.numeric(format(dates, "%j")))
# 
# locs    <- read_stars("site_parameter/nhsce_v01r01_19661004_20220103.nc")
# lon_lat <- tibble(lon = as.numeric(c(locs$longitude)),
#                   lat = as.numeric(c(locs$latitude))) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#            mutate(land = c(read_stars("site_parameter/nhsce_v01r01_19661004_20220103.nc", sub = "land")$land)) %>%
#   st_transform(st_crs(breed_polys))
# 
# snowTab <- parallel::mclapply(breed_polys$index, function(id) {
# 
#  ids <- as.numeric(lon_lat %>% rownames_to_column(var = "index") %>%
#    mutate(dist = as.numeric(st_distance(lon_lat, breed_polys[which(breed_polys$index==id),] %>% st_centroid()))/1000) %>%
#    filter(land==1) %>% arrange(dist) %>% pull(index))[1:4]
# 
# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = lon_lat[ids,])
# 
#  pastD <- tibble(doy  = rep(datesD$doy[datesD$year %in% 1970:1975],  each = length(ids)),
#                  year = rep(datesD$year[datesD$year %in% 1970:1975], each = length(ids)),
#                  snow = c(snowMat[ids, which(datesD$year %in% 1970:1975)])) %>%
#                  mutate(pixel = rep(ids, nrow(.)/length(ids))) %>%
#                  filter(snow<=1)
# 
# 
# 
# # ggplot(pastD, aes(x = doy, y = snow, color = as.factor(year))) +
# #   geom_line() +
# #   facet_wrap(~pixel)
# # 
# # with(pastD %>% filter(pixel == 4526), plot(date, snow, type = "o"))
# 
# 
# res1 <- lapply(pastD %>% group_split(pixel), function(p) gaussMLE(day = p$doy, size = rep(100, nrow(p)), prob = p$snow, thresh = 0.25))
# 
# plot(pastD$doy, pastD$snow)
# lapply(res1, function(l)  lines(l$result$week, l$result$fit))
# lapply(res1, function(l)  points(l$start, 0.25, pch = 16))
# 
#  pastSnow <- median(sapply(group_split(pastD, year), function(s) min(s$doy[s$snow<1]))) #unlist(median(sapply(res1, function(l)  l$start)))
# 
#  currD <- tibble(doy  = rep(datesD$doy[datesD$year %in% 2017:2021],  each = length(ids)),
#                  year = rep(datesD$year[datesD$year %in% 2017:2021], each = length(ids)),
#                  snow = c(snowMat[ids, which(datesD$year %in% 2017:2021)])) %>%
#    mutate(pixel = rep(ids, nrow(.)/length(ids))) %>%
#    filter(snow<=1)
# 
#  # res2 <- lapply(currD %>% group_split(pixel), function(p) gaussMLE(day = p$doy, size = rep(100, nrow(p)), prob = p$snow, thresh = 0.25))
#  #
#  # plot(currD$doy, currD$snow)
#  # lapply(res2, function(l)  lines(l$result$week, l$result$fit))
#  # lapply(res2, function(l)  points(l$start, 0.25, pch = 16))
# 
#  currSnow <- median(sapply(group_split(currD, year), function(s) min(s$doy[s$snow<1]))) # unlist(median(sapply(res2, function(l)  l$start)))
# 
#  tibble(index = id, pastSnow = pastSnow, currSnow = currSnow)
# }, mc.cores = 1) %>% Reduce("rbind",.)
# 
# 
# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = breed_polys %>% mutate(snow = snowTab$currSnow), mapping = aes(fill = snow)) +
#   geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = snowTab$pastSnow), size = 2)
#   # geom_sf(data = breed_polys %>% mutate(snow = ifelse(snowTab$currSnow<100, 1, 2)), mapping = aes(fill = as.factor(snow)))
# 
# 
# plot(snowTab$pastSnow, (mudflatTab[snowTab$index,] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,2])
# 
# 
# 
# 
# plot(tempTab[693,,1])
# abline(v = arrivalTab$arrival_hist[arrivalTab$index==693])
# abline(h = 0)
# 
# 
# 
# 
# 
# ggplot(tibble(snow = snowTab$currSnow, temp = apply(cbind(as.numeric(breed_polys$index), round(snowTab$currSnow)), 1, function(x) tempTab[x[1], x[2], 2])) %>% filter(snow>100),
#        aes(x = snow, y = temp)) +
#   geom_point() +
#   geom_smooth(method='lm')
# 
# lmTab <- tibble(snow_hist = snowTab$pastSnow,
#                 snow_curr = snowTab$currSnow,
#                 temp_hist = apply(cbind(as.numeric(breed_polys$index), round(snowTab$pastSnow)), 1, function(x) tempTab[x[1], x[2], 1]),
#                 temp_curr = apply(cbind(as.numeric(breed_polys$index), round(snowTab$currSnow)), 1, function(x) tempTab[x[1], x[2], 2]),
#                 latitude  = (mudflatTab[snowTab$index,] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,2],
#                 start_hist = snow_hist, start_curr = snow_curr)
# 
# lmTab$start_hist[lmTab$start_hist<100] <- predict(lm(snow_hist~latitude, data = lmTab %>% filter(snow_hist>100)) ,
#                                                   newdata = tibble(latitude = lmTab$latitude))[lmTab$start_hist<100]
# 
# lmTab$start_curr[lmTab$start_curr<100] <- predict(lm(snow_curr~latitude, data = lmTab %>% filter(snow_curr>100)) ,
#                                                   newdata = tibble(latitude = lmTab$latitude))[lmTab$start_curr<100]
# 
# ggplot(lmTab, aes(x = latitude, y = start_hist)) +
#   geom_point() +
#   geom_smooth(method='lm', color = "#EEC900")+
#   labs(x = "latitude" , y = "time of snow melt [day of the year]", title = "Past" ) +
#   theme_bw()
# 
# ggplot(lmTab, aes(x = latitude, y = start_curr)) +
#   geom_point() +
#   geom_smooth(method='lm', color = "#9A32CD")+
#   labs(x = "latitude" , y = "time of snow melt [day of the year]", title = "Present" ) +
#   theme_bw()
# 
# 
# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   # geom_sf(data = breed_polys %>% mutate(snow = snowTab$currSnow), mapping = aes(fill = snow)) +
#   geom_sf(data = breed_polys, mapping = aes(fill = lmTab$start_hist)) +
#   geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = round(lmTab$start_hist, 0)), size = 2)
# 
# 
# arrivalTab <- tibble(index = as.numeric(breed_polys$index), arrival_hist = lmTab$start_hist, arrival_curr = lmTab$start_curr)
# hist(arrivalTab %>% filter(index!=346, index!=441) %>% pull (arrival_hist))
# hist(arrivalTab %>% filter(index!=346, index!=441) %>% pull (arrival_curr))
# save(arrivalTab, file = 'temp/Waders/arrivalTab.RData')
# load('temp/Waders/arrivalTab.RData')

# ggplot()+
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = mudflatTab[arrivalTab$index,], mapping = aes(fill = arrivalTab$arrival_curr)) +
#   geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = round(arrivalTab$arrival_curr, 0)), size = 2)


### plot
# tab <- start_end_list['Arenaria interpres'][[1]]
# 
# sf_tab <- lapply(1:nrow(tab), function(x) {
#   mudflatTab[c(tab$breeding_id[x], tab$winter_id[x]),] %>% st_centroid() %>% st_geometry() %>% st_combine() %>% st_cast('LINESTRING')
# }) %>% Reduce('rbind', .) %>% st_as_sfc() %>% st_set_crs(st_crs(mudflatTab))
# 
# ggplot() +
#   geom_sf(data = shp %>% filter(SCINAME=='Arenaria interpres'), mapping = aes(geometry = geometry, fill = as.factor(SEASONAL))) +
#   geom_sf(data =  mudflatTab[tab$breeding_id,] %>% st_centroid() %>% st_geometry(), pch = 16, color = 'red') +
#   geom_sf(data =  mudflatTab[tab$winter_id,] %>% st_centroid() %>% st_geometry(), pch = 16, color = 'darkgreen') +
#   geom_sf(data = sf_tab, alpha = 0.4)