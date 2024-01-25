####### State Dependent Optimality Model  #########################################################################################################################

########## SET  UP 

remotes::install_github("slisovski/migrationSDP")
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

######Data######


load('sdpModel/site_parameter/eaafMap.rda')

ggplot() +
  plot+
  geom_sf(data = eaafMap$grid)
geom_sf(data = eaafMap$map, fill = adjustcolor('grey40', alpha.f = 0.4)) +
  geom_sf(data = eaafMap$bbox, fill = NA) +
  geom_sf(data = crds_sf, aes(geometry = geometry), color = "red") +
  theme_void()

crds_sf <- st_as_sf(crds, coords = c("Lon", "Lat"), crs = 4326)

load('sdpModel/site_parameter/mudflatTab.rda')

ggplot() +
  geom_sf(data = eaafMap$map) +
  geom_sf(data = mudflatTab %>% filter(lake_area>quantile(lake_area, probs = 0.5, na.rm = T)), aes(fill = lake_area), alpha = 0.75) +
  scale_fill_continuous(type = "viridis",na.value = "transparent")

siteTab      <- (mudflatTab %>% st_centroid() %>% st_transform(4326) %>%
                   mutate(areaHist = ifelse(is.na(histArea_inner), currArea_outer, histArea_outer),
                          areaCurr = currArea_outer,
                          areaMang = mangArea_outer,
                          lake     = ifelse(!is.na(lake_area) & lake_area > quantile(lake_area, probs = 0.5, na.rm = T), TRUE, FALSE)) %>%
                   dplyr::select(areaHist, areaCurr, areaMang, lake) %>%
                   rownames_to_column(var = "index") %>% mutate(index = as.integer(index))) %>%
  filter((areaHist>0 | areaCurr>0 | areaMang>0 | lake)) %>%
  filter(st_coordinates(.)[,1] > 104 | st_coordinates(.)[,1] < -150) %>%
  relocate(geometry, .after = last_col()) %>%
  suppressWarnings()

save(siteTab , file = 'sdpModel/site_parameter/siteTab.rda')
load('sdpModel/site_parameter/siteTab.rda')

# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = siteTab, aes(fill = lake), alpha = 0.75, size = 5, shape = 23) +
#   scale_fill_continuous(type = "viridis", na.value = "transparent")



#### Temp
load('sdpModel/site_parameter/tempTab.rda')

######Polygons########
####wint/breed#######


### 2 bboxes
bbox_crs <-  c(st_bbox(c(xmin =  100, xmax =  180, ymin = -80, ymax = 80), crs = 4326) %>% st_as_sfc(),
               st_bbox(c(xmin = -180, xmax = -150, ymin = -80, ymax = 80), crs = 4326) %>% st_as_sfc())


shp <- read_sf('sdpModel/species_parameters/Wader_3.shp') %>%
  st_intersection(bbox_crs) %>% st_transform(st_crs(eaafMap$map)) %>% st_intersection(eaafMap$bbox)
# save(shp, file = 'N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/transformed_shp.shp')
# load('N:/bioing/data/PathogenTransport/sdp_simulation/temp/Waders/transformed_shp.shp')

spTab <- readxl::read_xlsx('sdpModel/species_parameters/spTab2.xlsx') %>% filter(Filter==TRUE)


start_end_table <- lapplyf(spTab$ScieName, function(sp) {
  
  
  tmp_shp <- shp %>% filter(SCINAME==sp, SEASONA%in%c(3,2))
  
  tmp <- mudflatTab %>% rowid_to_column(var = 'id') %>% dplyr::select('id') %>%
    mutate(poly = apply(st_intersects(., tmp_shp, sparse = FALSE), 1, function(x) ifelse(any(x), which(x)[1], NA))) %>%
    filter(!is.na(poly)) %>% left_join(tibble(poly = 1:nrow(tmp_shp), seasonal = tmp_shp$SEASONA), by = 'poly')
  
  
  tibble(sp = sp, id = tmp$id, type = ifelse(tmp$seasonal==3, 'winter', 'breeding'))
  
}) %>% Reduce('rbind',.)

save(start_end_table, file = 'sdpModel/species_parameters/start_end_table.RData')

start_end_list <- lapply(unique(start_end_table$sp), function(species) {
  
  tmp <- start_end_table %>% filter(sp==species)
  
  winterSites <- siteTab[siteTab$index%in%tmp$id[tmp$type=='winter'],] %>% filter(!is.na(index)) %>%
    mutate(prob = scales:::rescale(log(areaHist + areaMang+0.1), c(0,1))) %>%
    mutate(prob = ifelse(spTab$Lake_propensity[spTab$ScieName==species] & prob < median(prob, na.rm = T),  median(prob, na.rm = T), prob))
  
  # ggplot() +
  #   geom_sf(data = shp %>% filter(SCINAME== species), mapping = aes(geometry = geometry, fill = as.factor(SEASONAL))) +
  #   geom_sf(data =  mudflatTab[winterSites$index,] %>% st_centroid() %>% st_geometry(), mapping = aes(size = winterSites$prob), pch = 16, color = 'darkgreen')
  
  
  tibble(breeding_id = sample(tmp %>% filter(type=='breeding') %>% pull(id), spTab$Pop[spTab$ScieName==species], TRUE)) %>%
    mutate(winter_id = sample(winterSites$index, spTab$Pop[spTab$ScieName==species], TRUE, prob = winterSites$prob)) %>%
    group_by(winter_id, breeding_id) %>% summarise(n = n())
  
})
names(start_end_list) <- unique(start_end_table$sp)
save(start_end_list, file = 'sdpModel/species_parameters/start_end_list.RData')
load('sdpModel/species_parameters/start_end_list.RData')


#### Sowmelt - Arrival #########

breed_polys <- mudflatTab %>% rownames_to_column(var = "index") %>%
  filter(index %in% (unique(start_end_table %>% filter(type == "breeding") %>% pull(id)))) %>%
  dplyr::select(index)

library(stars)
snow    <- read_stars("sdpModel/site_parameter/nhsce_v01r01_19661004_20220103.nc",
                      sub = 'snow_cover_extent')
snowMat <- apply(snow$snow_cover_extent, 3, function(x) c(x))

dates   <- st_get_dimension_values(snow, which = "time")
datesD  <- tibble(index = 1:length(dates), year = as.numeric(format(dates, "%Y")),
                  doy = as.numeric(format(dates, "%j")))

locs    <- read_stars("sdpModel/site_parameter/nhsce_v01r01_19661004_20220103.nc")
lon_lat <- tibble(lon = as.numeric(c(locs$longitude)),
                  lat = as.numeric(c(locs$latitude))) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(land = c(read_stars("sdpModel/site_parameter/nhsce_v01r01_19661004_20220103.nc", sub = "land")$land)) %>%
  st_transform(st_crs(breed_polys))

snowTab <- parallel::mclapply(breed_polys$index, function(id) {
  
  ids <- as.numeric(lon_lat %>% rownames_to_column(var = "index") %>%
                      mutate(dist = as.numeric(st_distance(lon_lat, breed_polys[which(breed_polys$index==id),] %>% st_centroid()))/1000) %>%
                      filter(land==1) %>% arrange(dist) %>% pull(index))[1:4]
  
  ggplot() +
    geom_sf(data = eaafMap$map) +
    geom_sf(data = lon_lat[ids,])
  
  pastD <- tibble(doy  = rep(datesD$doy[datesD$year %in% 1970:1975],  each = length(ids)),
                  year = rep(datesD$year[datesD$year %in% 1970:1975], each = length(ids)),
                  snow = c(snowMat[ids, which(datesD$year %in% 1970:1975)])) %>%
    mutate(pixel = rep(ids, nrow(.)/length(ids))) %>%
    filter(snow<=1)
  
  
  
  # ggplot(pastD, aes(x = doy, y = snow, color = as.factor(year))) +
  #   geom_line() +
  #   facet_wrap(~pixel)
  #
  # with(pastD %>% filter(pixel == 4526), plot(date, snow, type = "o"))
  
  
  res1 <- lapply(pastD %>% group_split(pixel), function(p) gaussMLE(day = p$doy, size = rep(100, nrow(p)), prob = p$snow, thresh = 0.25))
  
  plot(pastD$doy, pastD$snow)
  lapply(res1, function(l)  lines(l$result$week, l$result$fit))
  lapply(res1, function(l)  points(l$start, 0.25, pch = 16))
  
  pastSnow <- median(sapply(group_split(pastD, year), function(s) min(s$doy[s$snow<1]))) #unlist(median(sapply(res1, function(l)  l$start)))
  
  currD <- tibble(doy  = rep(datesD$doy[datesD$year %in% 2017:2021],  each = length(ids)),
                  year = rep(datesD$year[datesD$year %in% 2017:2021], each = length(ids)),
                  snow = c(snowMat[ids, which(datesD$year %in% 2017:2021)])) %>%
    mutate(pixel = rep(ids, nrow(.)/length(ids))) %>%
    filter(snow<=1)
  
  # res2 <- lapply(currD %>% group_split(pixel), function(p) gaussMLE(day = p$doy, size = rep(100, nrow(p)), prob = p$snow, thresh = 0.25))
  #
  # plot(currD$doy, currD$snow)
  # lapply(res2, function(l)  lines(l$result$week, l$result$fit))
  # lapply(res2, function(l)  points(l$start, 0.25, pch = 16))
  
  currSnow <- median(sapply(group_split(currD, year), function(s) min(s$doy[s$snow<1]))) # unlist(median(sapply(res2, function(l)  l$start)))
  
  tibble(index = id, pastSnow = pastSnow, currSnow = currSnow)
}, mc.cores = 1) %>% Reduce("rbind",.)


ggplot() +
  geom_sf(data = eaafMap$map) +
  geom_sf(data = breed_polys%>% filter (index != 346 & index!= 441) %>% mutate(snow = snowTab$currSnow[snowTab$index != 346 & snowTab$index!= 441]), mapping = aes(fill = snow)) +
  geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = snowTab$pastSnow), size = 2)
# geom_sf(data = breed_polys %>% mutate(snow = ifelse(snowTab$currSnow<100, 1, 2)), mapping = aes(fill = as.factor(snow)))


plot(snowTab$pastSnow, (mudflatTab[snowTab$index,] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,2])




plot(tempTab[693,,1])
abline(v = arrivalTab$arrival_hist[arrivalTab$index==693])
abline(h = 0)





ggplot(tibble(snow = snowTab$currSnow, temp = apply(cbind(as.numeric(breed_polys$index), round(snowTab$currSnow)), 1, function(x) tempTab[x[1], x[2], 2])) %>% filter(snow>100),
       aes(x = snow, y = temp)) +
  geom_point() +
  geom_smooth(method='lm')

lmTab <- tibble(snow_hist = snowTab$pastSnow,
                snow_curr = snowTab$currSnow,
                temp_hist = apply(cbind(as.numeric(breed_polys$index), round(snowTab$pastSnow)), 1, function(x) tempTab[x[1], x[2], 1]),
                temp_curr = apply(cbind(as.numeric(breed_polys$index), round(snowTab$currSnow)), 1, function(x) tempTab[x[1], x[2], 2]),
                latitude  = (mudflatTab[snowTab$index,] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,2],
                start_hist = snow_hist, start_curr = snow_curr)

lmTab$start_hist[lmTab$start_hist<100] <- predict(lm(snow_hist~latitude, data = lmTab %>% filter(snow_hist>100)) ,
                                                  newdata = tibble(latitude = lmTab$latitude))[lmTab$start_hist<100]

lmTab$start_curr[lmTab$start_curr<100] <- predict(lm(snow_curr~latitude, data = lmTab %>% filter(snow_curr>100)) ,
                                                  newdata = tibble(latitude = lmTab$latitude))[lmTab$start_curr<100]

ggplot(lmTab, aes(x = latitude, y = start_hist)) +
  geom_point() +
  geom_smooth(method='lm', color = "#EEC900")+
  labs(x = "latitude" , y = "time of snow melt [day of the year]", title = "Past" ) +
  theme_bw()

ggplot(lmTab, aes(x = latitude, y = start_curr)) +
  geom_point() +
  geom_smooth(method='lm', color = "#9A32CD")+
  labs(x = "latitude" , y = "time of snow melt [day of the year]", title = "Present" ) +
  theme_bw()


ggplot() +
  geom_sf(data = eaafMap$map) +
  # geom_sf(data = breed_polys %>% mutate(snow = snowTab$currSnow), mapping = aes(fill = snow)) +
  geom_sf(data = breed_polys, mapping = aes(fill = lmTab$start_hist)) +
  geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = round(lmTab$start_hist, 0)), size = 2)


arrivalTab <- tibble(index = as.numeric(breed_polys$index), arrival_hist = lmTab$start_hist, arrival_curr = lmTab$start_curr)
hist(arrivalTab %>% filter(index!=346, index!=441) %>% pull (arrival_hist))
hist(arrivalTab %>% filter(index!=346, index!=441) %>% pull (arrival_curr))
save(arrivalTab, file = 'sdpModel/species_parameters/arrivalTab.RData')
load('sdpModel/species_parameters/arrivalTab.RData')

# ggplot()+
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = mudflatTab[arrivalTab$index,], mapping = aes(fill = arrivalTab$arrival_curr)) +
#   geom_text(data = as_tibble(st_coordinates(breed_polys %>% st_centroid())), mapping = aes(x = X, y = Y, label = round(arrivalTab$arrival_curr, 0)), size = 2)

########### RUN SIMULATION ####################################################################################################################################################

source('Simulation_setup.R')

# simulation_list <- lapply(1:length(start_end_list),
#                            function(x) start_end_list[[x]] %>% ungroup() %>% mutate(sp = names(start_end_list)[x])) %>% Reduce('rbind',.) %>%
#                           #function(x) start_end_list[[x]][sample(1:nrow(start_end_list[[x]]),10),] %>% ungroup() %>% mutate(sp = names(start_end_list)[x])) %>% Reduce('rbind',.)  %>%
#                           arrange(sample(1:nrow(.)))
# save(simulation_list, file = 'sdpModel/simulation_list.rda')

load('sdpModel/simulation_list.rda')
simulation_list <- simulation_list %>% filter(breeding_id!=346, breeding_id!=441)


## Setup parallel
library(foreach)
library(doParallel)
cl <- makeCluster(75)
registerDoParallel(cl)
invisible(clusterEvalQ(cl, {
  source('sdpModel/Simulation_setup.R')
  load('sdpModel/simulation_list.rda')
}))


sim_all <- foreach(i = 1:nrow(simulation_list), .combine = 'c') %dopar% {
  
  # if(!file.exists(paste0('simu_out/simu_', i, '.rda'))) {
  
  ### Species params
  spParms <- shorebirdScaling(spTab %>% filter(ScieName==simulation_list$sp[i]) %>% pull(LeanBodyMass))
  
  
  ### Site setup
  StartEnd_cell <- mudflatTab[c(simulation_list$winter_id[i], simulation_list$breeding_id[i]),] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates() %>%
    as_tibble() %>% setNames(c('lon', 'lat')) %>%
    st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326)
  
  size     <- (mudflatTab %>% st_centroid() %>% st_transform(4326) %>%
                 mutate(areaHist = ifelse(is.na(histArea_inner), currArea_outer, histArea_outer),
                        areaCurr = currArea_outer,
                        areaMang = mangArea_outer,
                        lake     = ifelse(is.na(lake_area), FALSE, TRUE),
                        pmax     = pmax(areaCurr, areaHist)) %>%
                 dplyr::select(pmax, areaHist, areaCurr, areaMang, lake) %>%
                 rownames_to_column(var = "index") %>% mutate(index = as.integer(index))) %>%
    filter((areaHist>0 | areaCurr>0 | areaMang>0 | lake)) %>%
    filter(st_coordinates(.)[,1] > 104 | st_coordinates(.)[,1] < -150)  %>%
    arrange(st_distance(geometry, StartEnd_cell[1,])) %>%
    bind_rows(StartEnd_cell[2,] %>% mutate(areaHist = 0, areaCurr = 0, areaMang = NA, lake = FALSE) %>% dplyr::select(names(.))) %>%
    relocate(geometry, .after = last_col()) %>%
    suppressWarnings()
  
  arrival <- as.numeric((arrivalTab %>% filter(index == simulation_list$breeding_id[i]))[,c('arrival_hist', 'arrival_curr')])
  
  ### Intake
  pmax <- apply(size %>% st_drop_geometry(), 1, function(x) pmax(x[3]+((x[5]*(1000000))/1000)*0.4, x[4]+((x[5]*(1000000))/1000)*0.4))
  hist <- apply(size %>% st_drop_geometry(), 1, function(x) x[3]+((x[5]*(1000000))/1000)*0.4)
  curr <- apply(size %>% st_drop_geometry(), 1, function(x) x[4]+((x[5]*(1000000))/1000)*0.4)
  
  intake <- tibble(index = size$index,
                   hist  = ifelse(size$lake & size$pmax==0, 0, approx(range(pmax, na.rm = T),
                                                                      c(0.25, 1), hist)$y),
                   curr  = ifelse(size$lake & size$pmax==0, 0, approx(range(pmax, na.rm = T),
                                                                      c(0.25, 1), curr)$y)) %>%
    mutate(hist = ifelse((spTab %>% filter(ScieName==simulation_list$sp[i]) %>% pull(Lake_propensity)) & hist<1 & size$lake, 0.75, hist),
           curr = ifelse((spTab %>% filter(ScieName==simulation_list$sp[i]) %>% pull(Lake_propensity)) & curr<1 & size$lake, 0.75, curr),
           intHist = hist * spParms$FDRx + spParms$EEFnc(spParms$Kesm)/spParms$X1xkJ,
           intCurr = curr * spParms$FDRx + spParms$EEFnc(spParms$Kesm)/spParms$X1xkJ)
  
  
  crds <- rbind(mudflatTab[size$index[-nrow(size)],] %>% st_centroid() %>% st_transform(4326) %>% st_coordinates(),
                StartEnd_cell[2,] %>% st_coordinates()) %>% as_tibble() %>% setNames(c("Lon", "Lat"))
  
  
  ### expend (temperature)
  daily_exp_x <- abind::abind(lapply(1:2, function(t) apply(tempTab[size$index[-nrow(size)],,t], 2, function(x)  spParms$EEFnc(x)/spParms$X1xkJ)), along = 3)
  
  
  ### simulation time
  start   <- "02-01"
  end     <- "07-30"
  doy_seq <-  as.numeric(format(seq(as.POSIXct(glue::glue("2020-{start}")), as.POSIXct(glue::glue("2020-{end}")), by = "day"), "%j"))
  
  expM  <- daily_exp_x[,doy_seq,]
  gain  <- intake %>% dplyr::select(index, intHist, intCurr)
  
  
  simuList <- lapply(c(0,1), function(s) {
    
    obj <- make_migrationSDP(init    =  list(scenario = s,
                                             minT     = min(doy_seq),
                                             maxT     = max(doy_seq),
                                             MaxX     = 100,
                                             tReward  = list(dstr = "dsnorm", mean = arrival[s+1], sd = 3, xi = 2, factor = 2),
                                             dError   = 15000),
                             species = list(B0       = 3,
                                            w        = 0.028,
                                            xc       = 10,
                                            speed    = spParms$speed,
                                            c        = spParms$c),
                             sites   = list(crds     = crds,
                                            pred     = c(1e-3, 1e-4, 1e-6, 2, 2),
                                            expend   = expM[,,s+1],
                                            gain     = gain,
                                            penalty  = c(0, 0.0015, 0)),
                             parms   = list(bearing  = TRUE, angle = 100))
    
    
    mod <- tryCatch(bwdIteration(obj), error = function(e) NULL)
    
    if(!is.null(mod)) {
      
      for(rep in 1:10) {
        simu    <- tryCatch(fwdSimulation(mod, simulation_list$n[i], start_t = 1, start_site = 1, start_x = c(30,50)), error = function(e) NULL)
        if(is.array(simu)) break
      }
      
      # x_site_plot(simu)
      # mat <- site_matrix(mudflatTab$geometry, simulation_list$breeding_id[i], simu = simu, mod = mod)
      # network_plot(mudflatTab$geometry, eaafMap$map, mat$network, species = "", linksMin = 0)
      
      
      if(!is.null(simu)) {
        index <- mod@Sites$index
        index[length(index)] <- simulation_list$breeding_id[i]
        simu[,2,] <- index[simu[,2,]+1]
        simu
      } else NULL
      
    } else NULL
  })
  
  
  save(simuList, file = paste0('sdpModel/simu_out/simu_', i, '.rda'))
  #}
  
  NULL
}

stopCluster(cl)

############ MERGE SIMULATION OUTPUT ################################################################################################################

source('sdpModel/Simulation_setup.R')


load('sdpModel/simulation_list.rda')
simulation_list <- simulation_list %>% filter(breeding_id!=346, breeding_id!=441)

mod_fls         <- tibble(file = list.files("simu_out", pattern = ".rda"))
mod_fls$id      <- as.numeric(gsub(".rda", "", gsub("simu_", "", mod_fls$file)))
mod_fls$species <- simulation_list$sp[mod_fls$id ]
mod_fls <- mod_fls %>% arrange(id)



##### Matrices (Network)
netwMat <- parallel::mclapply(unique(mod_fls$species), function(sp) {
  
  if(!file.exists(gsub(" ", "_", glue::glue("spMat_out/{sp}.rda")))) {
    out <- mod_fls %>% filter (mod_fls$species == sp)
    
    spArray <- array(0, dim = c(nrow(mudflatTab), nrow(mudflatTab), 2))
    
    for(i in 1:nrow(out)) {
      
      load(paste0("simu_out/", out$file[i]))
      
      simMat <- abind::abind(parallel::mclapply(simuList, function(a) {
        
        if(!is.null(a)){
          
          alive  <- apply(a[,5,], 1, function(x) all(x==0))
          sts0   <- a[alive,2,]
          if(sum(alive)>1) {
            
            sts    <- t(apply(sts0, 1, function(x) {
              start <- simulation_list$winter_id[out$id[i]]
              x[1:max(which(x==start))] <- start
              x}))
            
            ### Matrix
            site_index <- 1:nrow(mudflatTab)
            
            stay   <- apply(sts, 1, function(x) { 
              as.data.frame(table(x)) %>% as_tibble() %>% setNames(c('site', 'days')) %>% mutate(site = as.numeric(as.character(site))) %>%
                full_join(tibble(site = 1:length(site_index)), by = "site") %>% arrange(site) %>% mutate(days = ifelse(site==x[1] | site==x[length(x)], NA, days)) %>% pull(days)
            }) %>% apply(., 1, sum, na.rm = T)
            
            transM <- expand_grid(a = 1:length(site_index), b = 1:length(site_index)) %>% 
              left_join(apply(sts, 1, function(x) tibble(a = x[!is.na(x)][-sum(!is.na(x))], b = x[!is.na(x)][-1]) %>% filter(a!=b) %>% mutate(n = 1)) %>%
                          do.call("rbind",.) %>% group_by(a,b) %>% summarise_at("n", sum), by = join_by(a, b)) %>% filter(n>0) %>%
              mutate(a_grid = site_index[a], b_grid = site_index[b])
            
            mat <- matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
            diag(mat) <- stay
            mat[cbind(transM$a_grid, transM$b_grid)] <- transM$n
            mat
            
          } else matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
          
          
        } else matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
        
      }, mc.cores = 2), along = 3)
      
      spArray <- Reduce("+", list(spArray, simMat))
      
    }
    save(spArray, file = gsub(" ", "_", glue::glue("sdpModel/spMat_out/{sp}.rda")))
    
    NULL }
  
  
}, mc.cores = length(unique(mod_fls$species)))


##### Migration strategies
strategyTab <- parallel::mclapply(unique(mod_fls$species), function(sp) {
  
  out <- mod_fls %>% filter (mod_fls$species == sp)
  
  spTab <- parallel::mclapply(1:nrow(out), function(i) {
    
    load(paste0("sdpModel/simu_out/", out$file[i]))
    
    if(!is.null(simuList[[1]]) & !is.null(simuList[[2]])) {
      
      spList <- lapply(simuList, function(a) { 
        
        alive  <- apply(a[,5,], 1, function(x) all(x==0))
        sts0   <- a[alive,2,]
        
        if(sum(alive)>1) {
          
          sts    <- t(apply(sts0, 1, function(x) {
            start <- simulation_list$winter_id[out$id[i]]
            x[1:max(which(x==start))] <- start
            x}))
          
          
          strategy <- apply(sts, 1, function(x) {
            ID <- runif(1)
            
            mat01 <- (mudflatTab[x[!is.na(x) & !duplicated(x)],] %>% st_centroid() %>% st_distance() %>% suppressWarnings())/1000
            steps <- as.numeric(diag(as.matrix(mat01[-nrow(mat01),-1])))
            
            dur   <- min(which(x==simulation_list$breeding_id[out$id[i]])) - max(which(x==simulation_list$winter_id[out$id[i]]))
            
            mig <- ifelse(is.na(x) | (x!=simulation_list$winter_id[out$id[i]] & x!=simulation_list$breeding_id[out$id[i]]), TRUE, FALSE)
            
            if(any(!is.na(x[mig]))) {
              as.data.frame(table(x[mig])) %>% as_tibble() %>% setNames(c("Site", "Days")) %>%
                arrange(desc(Days)) %>% mutate(Prob = (Days/sum(Days))*100, Dist = NA, Dur = NA) %>% bind_rows(
                  tibble(Site = NA, Days = NA, Prob = NA, Dist = steps, Dur = NA)
                ) %>% bind_rows(tibble(Site = NA, Days = NA, Prob = NA, Dist = NA, Dur = dur)) %>% mutate(ID = ID, .before = Site) %>% mutate(Simu = out$id[i], .after = ID)
            } else {
              tibble(Site = NA, Days = NA, Prob = NA, Dist = steps, Dur = dur) %>% mutate(ID = ID, .before = Site) %>% mutate(Simu = out$id[i], .after = ID)
            }
            
          }) %>% Reduce("rbind",.)
          
        } else NULL
        
      })
      
    } else list(NULL, NULL)
    
  }, mc.cores = 40)
  
  spPast <- lapply(spTab, function(x) x[[1]]) %>% Reduce("rbind",.)
  spCurr <- lapply(spTab, function(x) x[[2]]) %>% Reduce("rbind",.)
  
  save(spPast, file = gsub(" ", "_", glue::glue("sdpModel/spMat_out/Strategy/{sp}_past.rda")))
  save(spCurr, file = gsub(" ", "_", glue::glue("sdpModel/spMat_out/Strategy/{sp}_curr.rda")))
  
}, mc.cores = 2)


################# MODEL EVALUATION #####################################################################################################################

load('sdpModel/simulation_list.rda')
load("TrackingAnalysis/TrackingData/bird_tracking_data.rda")
simulation_list <- simulation_list %>% filter(breeding_id!=346, breeding_id!=441)

source('functions.R')


#######identifying cell numbers of breeding/wintering sites from empirical data 

sci_names <- c("Calidris ferruginea", "Limosa lapponica", "Calidris tenuirostris", "Calidris canutus","Calidris alba", "Arenaria interpres")

start_end <- lapply(unique(bird_tracking_data$ID), function(id) {
  emp_sites <- bird_tracking_data %>%
    filter(ID == id, Type %in% c(0, 2)) %>%
    distinct(ID, Type, .keep_all = TRUE) %>%
    st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(crs = st_crs(mudflatTab$geometry))
  
  start <- which(st_intersects(mudflatTab$geometry, emp_sites$geometry[emp_sites$Type == 0], sparse = FALSE))
  end <- which(st_intersects(mudflatTab$geometry, emp_sites$geometry[emp_sites$Type == 2], sparse = FALSE))
  
  tibble(Species = emp_sites$Species[which(emp_sites$ID == id)[1]], ID = id, start = start, end = end) %>%
    mutate (Sci_species = sci_names[which(unique(bird_tracking_data$Species) == emp_sites$Species[which(emp_sites$ID == id)[1]])])
}) %>% do.call(rbind, .)

###### which simu datei

simulation_list_emp <- simulation_list %>% mutate (id_simu = row_number()) %>% filter(sp %in% unique(start_end$Sci_species)) %>% group_split(sp)%>%
  lapply(., function(s) {
    start_end_sp <- start_end %>% filter(Sci_species==s$sp[1])
    s %>% filter(glue::glue("{winter_id}_{breeding_id}") %in% glue::glue("{start_end_sp$start}_{start_end_sp$end}"))
  }) %>% Reduce('rbind',.)

save(simulation_list_emp, file = "sdpModel/simulation_list_emp.rda")

mod_fls         <- tibble(file = list.files("sdpModel/u_out", pattern = ".rda"))
mod_fls$id      <- as.numeric(gsub(".rda", "", gsub("simu_", "", mod_fls$file)))
mod_fls$species <- simulation_list$sp[mod_fls$id ]
mod_fls <- mod_fls %>% arrange(id) %>% filter(id %in% simulation_list_emp$id_simu)

####### Matrices for Simus for breeding and wintering sites of emperical data 

netwMat <- parallel::mclapply(unique(mod_fls$species), function(sp) {
  
  if(!file.exists(gsub(" ", "_", glue::glue("sdpModel/mod_evaluation/{sp}.rda")))) {
    out <- mod_fls %>% filter (mod_fls$species == sp)
    
    spArray <- array(0, dim = c(nrow(mudflatTab), nrow(mudflatTab), 2))
    
    for(i in 1:nrow(out)) {
      
      load(paste0("simu_out/", out$file[i]))
      
      simMat <- abind::abind(parallel::mclapply(simuList, function(a) {
        
        if(!is.null(a)){
          
          alive  <- apply(a[,5,], 1, function(x) all(x==0))
          sts0   <- a[alive,2,]
          if(sum(alive)>1) {
            
            sts    <- t(apply(sts0, 1, function(x) {
              start <- simulation_list$winter_id[out$id[i]]
              x[1:max(which(x==start))] <- start
              x}))
            
            
            ### Matrix
            site_index <- 1:nrow(mudflatTab)
            
            stay   <- apply(sts, 1, function(x) {
              as.data.frame(table(x)) %>% as_tibble() %>% setNames(c('site', 'days')) %>% mutate(site = as.numeric(as.character(site))) %>%
                full_join(tibble(site = 1:length(site_index)), by = "site") %>% arrange(site) %>% mutate(days = ifelse(site==x[1] | site==x[length(x)], NA, days)) %>% pull(days)
            }) %>% apply(., 1, sum, na.rm = T)
            
            transM <- expand_grid(a = 1:length(site_index), b = 1:length(site_index)) %>%
              left_join(apply(sts, 1, function(x) tibble(a = x[!is.na(x)][-sum(!is.na(x))], b = x[!is.na(x)][-1]) %>% filter(a!=b) %>% mutate(n = 1)) %>%
                          do.call("rbind",.) %>% group_by(a,b) %>% summarise_at("n", sum), by = join_by(a, b)) %>% filter(n>0) %>%
              mutate(a_grid = site_index[a], b_grid = site_index[b])
            
            mat <- matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
            diag(mat) <- stay
            mat[cbind(transM$a_grid, transM$b_grid)] <- transM$n
            mat
            
          } else matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
          
          
        } else matrix(0, ncol = nrow(mudflatTab), nrow = nrow(mudflatTab))
        
      }, mc.cores = 2), along = 3)
      
      spArray <- Reduce("+", list(spArray, simMat))
      
    }
    save(spArray, file = gsub(" ", "_", glue::glue("sdpModel/mod_evaluation/{sp}.rda")))
    
    NULL }
  
  
}, mc.cores = length(unique(mod_fls$species)))


#################### Plotting Networks ############################################################################################################

######### Model predictions

fls <- tibble(files = list.files('sdpModel/mod_evaluation', pattern = '.rda', full.names = T))
fls$sp <-gsub("sdpModel/mod_evaluation/(.*?)\\.rda$", "\\1", fls$files) %>% gsub("_"," ",.)


color_list <- tibble( color= c("steelblue3", "darkorange", "forestgreen", "firebrick", "darkviolet", "turquoise"),
                      sci = unique(fls$sp)) %>% left_join(spTab, by = c("sci" = "ScieName")) %>% select(color, sci, sp = Species)

plot_list_mod <- lapply(fls$files, function (file){
  load(file)
  sp <- spTab$Species[spTab$ScieName == fls$sp[fls$files == file]]
  plot1 <- network_plot(mudflatTab$geometry, eaafMap$map, spArray[,,1], "", "Model Data - Past", 0, color_list$color[color_list$sp == sp])  
  plot2 <- network_plot(mudflatTab$geometry, eaafMap$map, spArray[,,2],  "","Model Data - Present", 0, color_list$color[color_list$sp == sp])
  
  list(plot1,plot2)
})

##### Emperical Data

grid <- mudflatTab$geometry
mat_0 <- matrix(0, ncol = 1224, nrow = 1224)

spec_List <- lapply(unique(bird_tracking_data$Species)[c(2,3,4,6,1,5)], function(sp) {
  
  (bird_tracking_data%>% filter(Species == sp) %>%  st_as_sf(coords = c('Lon', 'Lat'), crs = 4326) %>%
     st_transform(crs= st_crs(grid)) %>% group_split(ID) %>% 
     lapply(., function(df) {
       first_index <- which(df$Type == 2)[1]
       tab <- df[1:first_index, ]
       tab <- tab %>% mutate(cell = apply(st_intersects(., grid, sparse = FALSE), 1, which))
       mat_out <- mat_0
       diag(mat_out)[tab$cell[tab$Type==1]] <- tab$Days[tab$Type==1]
       mat_out[cbind(tab$cell[-nrow(tab)], tab$cell[-1])] <- 1
       mat_out
     })) %>% Reduce('+', .) 
  
})

index <- unique(bird_tracking_data$Species)[c(2,3,4,6,1,5)]

plot_list_track <- lapply(1:6, function(x){
  plot <- network_plot(grid, eaafMap$map, spec_List[[x]], paste( index[x]), "Tracking Data",  0, color_list$color[color_list$sp == index[x]])
})


### combine plots 
library(cowplot)

pl1 <- plot_grid(plot_list_track[[1]],plot_list_mod[[6]][[1]],plot_list_mod[[6]][[2]], 
                 plot_list_track[[2]],plot_list_mod[[5]][[1]],plot_list_mod[[5]][[2]],
                 plot_list_track[[3]],plot_list_mod[[3]][[1]],plot_list_mod[[3]][[2]],
                 nrow=3 ,ncol =3)

pl2 <- plot_grid(plot_list_track[[4]],plot_list_mod[[1]][[1]],plot_list_mod[[1]][[2]],
                 plot_list_track[[5]],plot_list_mod[[4]][[1]],plot_list_mod[[4]][[2]],
                 plot_list_track[[6]],plot_list_mod[[2]][[1]],plot_list_mod[[2]][[2]],
                 nrow=3, ncol =3)

###########Comparing Strategies Emp/Pred Data #########################################################################################################################
library(tidyverse)
library(ggplot2)

load("TrackingAnalysis/TrackingData/step_length_df.rda")

########Tracking Data - Ratio of Time Use ##
ratio <- step_length %>% group_split(ID) %>% lapply(., function(x) {out <- x %>% mutate(Prob = (Days /(sum(Days, na.rm = T)))*100)
return(out)}) %>% Reduce("rbind",.)

color_list <- tibble( color= c("steelblue3", "darkorange", "forestgreen", "firebrick", "darkviolet", "turquoise"),
                      sp = unique(ratio$Sp)[c(6,2,3,4,5,1)])

pl_ratio_tracks <- lapply(unique((ratio %>% arrange(desc(lbm)))$Sp), function (spec){
  out <- ratio %>% filter ( Sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
    ungroup(.) %>% group_by(num_site) %>% summarize(ratio = median(Prob))
  points <- ratio %>% filter ( Sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
    ungroup(.)
  
  plot <- ggplot(data = out, aes(x= as.factor(num_site) , y=ratio))+
    geom_hline(yintercept = c(20, 40, 60, 80), linetype = "dashed", color = "gray")+ 
    geom_bar(stat = "identity", fill = color_list$color[color_list$sp == spec], color = "black")+
    geom_point(points, mapping= aes(x= as.factor(num_site), y= Prob), size = 0.7)+
    geom_line(data = points, aes(x = num_site, y = Prob, group = ID), color= "Black", alpha = 0.5)+
    labs (title = paste(spec), subtitle = "Tracking Data", x = "Sites", y = "Ratio")+
    scale_x_discrete (limits = as.factor(c(1:10))) +
    theme(  axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
    annotate("text", x = Inf, y = Inf, hjust = 1.5, vjust = 2, label = paste("N =", length(unique(points$ID))), size = 4, color = "black")+
    ylim(0, 100) + 
    theme_bw() +
    theme(panel.grid = element_blank())
  
  
  
  return(plot)
  
})

########Model Predictions - Ratio of Time Use ##

load('sdpModel/simulation_list.rda')
simulation_list <- simulation_list %>% filter(breeding_id!=346, breeding_id!=441)
load("sdpModel/simulation_list_emp.rda")
spTab <- readxl::read_xlsx('sdpModel/species_parameters/spTab2.xlsx') %>% filter(Filter==TRUE)



#### Past

fls  <- list.files('sdpModel/spMat_out/Strategy', pattern = 'past.rda', full.names = T)

sp_Past <- lapply(fls,function(file) {
  load(file)
  spCurr <- spPast %>% filter(!ID %in% ID[is.infinite(Dur)])
  return(spCurr)
}) %>% bind_rows(.)

ratio_p <- sp_Past %>% filter(!is.na(Prob)) %>% select(ID, Simu, Prob) %>% 
  filter(Simu %in% simulation_list_emp$id_simu) %>% 
  left_join(simulation_list_emp %>% select(id_simu, sp), by = c("Simu" = "id_simu"))%>%
  left_join(spTab %>% select("ScieName", "LeanBodyMass"), by = c("sp"="ScieName"))

color_list <- tibble(
  color = c("steelblue3", "darkorange", "forestgreen", "firebrick", "darkviolet", "turquoise"),
  sci = unique(ratio_p$sp)
) %>% 
  left_join(spTab, by = c("sci" = "ScieName")) %>%
  select(color, sci, sp = Species)


pl_ratio_model_past <- lapply(unique((ratio_p %>% arrange(desc(LeanBodyMass)))$sp), function (spec){
  out <- ratio_p %>% filter ( sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
    ungroup(.) %>% group_by(num_site) %>% summarize(ratio_p = median(Prob))
  
  plot <- ggplot(data = out, aes(x= as.factor(num_site) , y=ratio_p))+
    geom_hline(yintercept = c(20, 40, 60, 80), linetype = "dashed", color = "gray")+ 
    geom_bar(stat = "identity", fill = color_list$color[color_list$sci == spec], color = "black")+
    labs (title = "", subtitle = "Model Data - Past", x = "Sites", y = "Ratio")+
    scale_x_discrete (limits = as.factor(c(1:10))) +
    theme(  axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
    ylim(0, 100) + 
    theme_bw() +
    theme(panel.grid = element_blank())
  
  
  return(plot)
  
})

#### Present

fls  <- list.files('sdpModel/spMat_out/Strategy', pattern = 'curr.rda', full.names = T)

sp_Curr <- lapply(fls,function(file) {
  load(file)
  spCurr <- spCurr %>% filter(!ID %in% ID[is.infinite(Dur)])
  return(spCurr)
}) %>% bind_rows(.)

ratio_c <- sp_Curr %>% filter(!is.na(Prob)) %>% select(ID, Simu, Prob) %>% 
  filter(Simu %in% simulation_list_emp$id_simu) %>% 
  left_join(simulation_list_emp %>% select(id_simu, sp), by = c("Simu" = "id_simu"))%>%
  left_join(spTab %>% select("ScieName", "LeanBodyMass"), by = c("sp"="ScieName"))

pl_ratio_model_present <- lapply(unique((ratio_c %>% arrange(desc(LeanBodyMass)))$sp), function (spec){
  out <- ratio_c %>% filter ( sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
    ungroup(.) %>% group_by(num_site) %>% summarize(ratio_c = median(Prob))
  
  plot <- ggplot(data = out, aes(x= as.factor(num_site) , y=ratio_c))+
    geom_hline(yintercept = c(20, 40, 60, 80), linetype = "dashed", color = "gray")+ 
    geom_bar(stat = "identity", fill = color_list$color[color_list$sci == spec], color = "black")+
    labs (title = "", subtitle = "Model Data - Present", x = "Sites", y = "Ratio")+
    scale_x_discrete (limits = as.factor(c(1:10))) +
    theme(  axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
    ylim(0, 100) + 
    theme_bw() +
    theme(panel.grid = element_blank())
  
  
  return(plot)
  
})


#### combine plot 

library(cowplot)

strategies <- plot_grid(pl_ratio_tracks[[1]],pl_ratio_model_past [[1]], pl_ratio_model_present[[1]], 
                        pl_ratio_tracks[[2]],pl_ratio_model_past [[2]], pl_ratio_model_present[[2]],
                        pl_ratio_tracks[[3]],pl_ratio_model_past [[3]], pl_ratio_model_present[[3]],
                        pl_ratio_tracks[[4]],pl_ratio_model_past [[4]], pl_ratio_model_present[[4]],
                        pl_ratio_tracks[[5]],pl_ratio_model_past [[5]], pl_ratio_model_present[[5]],
                        pl_ratio_tracks[[6]],pl_ratio_model_past [[6]], pl_ratio_model_present[[6]],nrow= 6, ncol= 3)


######correlation tracks and predictions

ratioTab <- lapply(unique((ratio_c %>% arrange(desc(LeanBodyMass)))$sp), function (spec){
  
  ratioTab <- tibble (mod_past_ratio = (ratio_p %>% filter ( sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
                                          ungroup(.) %>% group_by(num_site) %>% summarize(ratio_p = median(Prob)) %>% complete(num_site = 1:10) %>% 
                                          tidyr::replace_na(list(ratio_p = 0)) %>% filter(num_site <11))[2],
                      mod_present_ratio = (ratio_c %>% filter ( sp == spec & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
                                             ungroup(.) %>% group_by(num_site) %>% summarize(ratio_c = median(Prob)) %>% complete(num_site = 1:10) %>% 
                                             tidyr::replace_na(list(ratio_c = 0)) %>% filter(num_site <11))[2],
                      
                      emp_ratio = (ratio %>% filter ( Sp == color_list$sp[color_list$sci == spec] & !is.na(Prob)) %>% group_by(ID) %>% arrange(desc(Prob)) %>% mutate (num_site = row_number()) %>% 
                                     ungroup(.) %>% group_by(num_site) %>% summarize(ratio = median(Prob))%>% complete(num_site = 1:10) %>% 
                                     tidyr::replace_na(list(ratio = 0)) %>% filter(num_site <11))[2],
                      sp = color_list$sp[color_list$sci == spec],
                      num_sites =1:10)
  
  return(ratioTab)
  
})


corr_array <- array(0,c(6,6,2))
dimnames(corr_array) <- list(unique((ratio %>% arrange(desc(lbm)))$Sp), unique((ratio %>% arrange(desc(lbm)))$Sp), c("Past", "Present"))

for ( e in 1:6){
  for (m in 1:6){
    corr <- cor(ratioTab[[e]][3], ratioTab[[m]][2])
    corr_array[e, m, 2] <- corr
  }
}


library(corrplot)
library(RColorBrewer)

corrplot(corr_array[,,2], method = "color", addCoef.col = "white", tl.col = "black", tl.srt = 45,  addgrid.col = "darkgray",
         col = colorRampPalette(brewer.pal(11, "Spectral"))(100)) 



############ Network Analysis Model Predictions ############################################################################################################################################################################################################################
library(tidyverse) # library(ggplot2) library(tidyr) library(dplyr)
library(igraph)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(scales)
load('sdpModel/site_parameter/eaafMap.rda')
load('sdpModel/site_parameter/mudflatTab.rda')


fls  <- list.files('spMat_out', pattern = '.rda', full.names = T)

all_past <- lapply(fls,function(file) {
  load(file)
  return(spArray[,,1])
}) %>% Reduce('+',.)

all_curr <- lapply(fls,function(file) {
  load(file)
  return(spArray[,,2])
}) %>% Reduce('+',.)

#### Total degree centrality

linkTab_p <- expand_grid(orig = 1:nrow(all_past), dest = 1:nrow(all_past)) %>% mutate(links=c(all_past)) %>%
  filter(orig!=dest, links>0) %>% group_split(row = 1:nrow(.)) %>%
  lapply(., function(x) tibble(orig = min(c(x$orig, x$dest)),
                               dest = max(c(x$orig, x$dest)),
                               links = x$links)) %>% Reduce('rbind',.) %>%
  group_by(orig, dest) %>% summarize(links = sum(links))


network_p <- graph_from_data_frame(d = linkTab_p, 
                                   vertices = tibble(node = sort(unique(c(linkTab_p$orig, linkTab_p$dest)))), 
                                   directed = T) 

degree_ntwrk_p <- degree(network_p, mode= "total")


linkTab_c <- expand_grid(orig = 1:nrow(all_curr), dest = 1:nrow(all_curr)) %>% mutate(links=c(all_curr)) %>%
  filter(orig!=dest, links>0) %>% group_split(row = 1:nrow(.)) %>%
  lapply(., function(x) tibble(orig = min(c(x$orig, x$dest)),
                               dest = max(c(x$orig, x$dest)),
                               links = x$links)) %>% Reduce('rbind',.) %>%
  group_by(orig, dest) %>% summarize(links = sum(links))


network_c <- graph_from_data_frame(d = linkTab_c, 
                                   vertices = tibble(node = sort(unique(c(linkTab_c$orig, linkTab_c$dest)))), 
                                   directed = T) 

degree_ntwrk_c <- degree(network_c, mode= "total")

net_p <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_p))))],
          mapping = aes(fill = degree_ntwrk_p), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title = "Degree Centrality (total)", subtitle = "Past")+
  theme_minimal()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(fill = guide_colorbar(title = "Degree 
Centrality"))

net_c <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_c))))],
          mapping = aes(fill = degree_ntwrk_c), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title=  "", subtitle = "Present") +
  theme_minimal()

library(cowplot)

legend <- get_legend(net_p)

net_plot <- plot_grid(net_p + theme(legend.position = "none"), net_c + theme(legend.position = "none"), nrow = 1) %>%
  plot_grid(.,legend, nrow =2, rel_widths = c(1, 1), rel_heights = c(1, 0.17), align = "hv", axis = "b")

ggsave("network_p_c_out.jpg",  path = "/bioing/data/PathogenTransport/Thesis/FinalPlots", plot= net_plot,width = 16, height =13.5, bg = "white", units = "cm", dpi=300)

#### Out degree centrality

degree_ntwrk_p_out <- degree(network_p, mode= "out")

degree_ntwrk_c_out <- degree(network_c, mode= "out")

net_p <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_p))))],
          mapping = aes(fill = degree_ntwrk_p_out), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title = "Degree Centrality (out)", subtitle = "Past")+
  theme_minimal()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(fill = guide_colorbar(title = "Degree 
Centrality"))

net_c <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_c))))],
          mapping = aes(fill = degree_ntwrk_c_out), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title=  "", subtitle = "Present") +
  theme_minimal()

library(cowplot)

legend <- get_legend(net_p)

net_plot <- plot_grid(net_p + theme(legend.position = "none"), net_c + theme(legend.position = "none"), nrow = 1) %>%
  plot_grid(.,legend, nrow =2, rel_widths = c(1, 1), rel_heights = c(1, 0.17), align = "hv", axis = "b")

### In degree centrality


degree_ntwrk_p_in <- degree(network_p, mode= "in")


degree_ntwrk_c_in <- degree(network_c, mode= "in")

net_p <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_p))))],
          mapping = aes(fill = degree_ntwrk_p_in), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title = "Degree Centrality (in)", subtitle = "Past")+
  theme_minimal()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(fill = guide_colorbar(title = "Degree 
Centrality"))

net_c <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey70") +
  geom_sf(data = eaafMap$map, fill = "grey70") +
  geom_sf(data = mudflatTab$geometry[as.numeric(names(unclass(V(network_c))))],
          mapping = aes(fill = degree_ntwrk_c_in), alpha = 0.5) +
  scale_fill_continuous(type = "viridis", direction = -1, name = "Value", limits = c(0, 420)) +
  labs(title=  "", subtitle = "Present") +
  theme_minimal()

library(cowplot)

legend <- get_legend(net_p)

net_plot <- plot_grid(net_p + theme(legend.position = "none"), net_c + theme(legend.position = "none"), nrow = 1) %>%
  plot_grid(.,legend, nrow =2, rel_widths = c(1, 1), rel_heights = c(1, 0.17), align = "hv", axis = "b")


####statistical tests
shapiro.test(degree_ntwrk_c)
shapiro.test(degree_ntwrk_p)
shapiro.test(degree_ntwrk_c_out)
shapiro.test(degree_ntwrk_c_in)
shapiro.test(degree_ntwrk_p_out)
shapiro.test(degree_ntwrk_p_in)

wilcox.test(degree_ntwrk_c, degree_ntwrk_p) #### p>0.05
wilcox.test(degree_ntwrk_c_in, degree_ntwrk_p_in)#### p>0.05
wilcox.test(degree_ntwrk_c_out, degree_ntwrk_p_out) ### p>0.05
wilcox.test(degree_ntwrk_c_in, degree_ntwrk_c_out) ### p = 0.0002192
wilcox.test(degree_ntwrk_p_in, degree_ntwrk_p_out) ### p = 0.0003538
