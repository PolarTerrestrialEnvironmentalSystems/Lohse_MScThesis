library(sf)
library(units)
library(tidyverse)


load("IBM/Data/report_wildbirds.rda")
load("IBM/Data/report_both_birds.rda")
load("sdpModel/site_parameter/mudflatTab.rda")
load("sdpModel/site_parameter/eaafMap.rda")


###########get location of outbreaks 

I_f <- I %>% select("Species","sero_sub_genotype_eng" ,"Outbreak_start_date","Outbreak_end_date", "Longitude","Latitude", "cases", "dead")

outbreak_sub <- I_f  %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(crs = st_crs(mudflatTab$geometry))

case_index <- lapply(1:nrow(mudflatTab), function(x) {
  sub <- mudflatTab[x, 8]
  there <- any(st_intersects(sub, outbreak_sub$geometry, sparse = FALSE))
  if (there) {
    outbreak_index <- outbreak_sub[which(st_intersects(sub, outbreak_sub$geometry, sparse = FALSE)),] %>% mutate( mud_index = x)
  }
}) %>% do.call(rbind, .)

cases <- case_index %>% group_by(mud_index) %>% summarise(cases = sum(cases)) %>% arrange(desc(cases))

ggplot() +
  geom_sf(data = eaafMap$map) +
  geom_sf(data = mudflatTab$geometry[353], alpha = 0.5)


########### prepping migration simulation 

source('sdpModel/Simulation_setup.R')

load('sdpModel/simulation_list.rda')
simulation_list <- simulation_list %>% filter(breeding_id!=346, breeding_id!=441)

mod_fls         <- tibble(file = list.files("simu_out", pattern = ".rda"))
mod_fls$id      <- as.numeric(gsub(".rda", "", gsub("simu_", "", mod_fls$file)))
mod_fls$species <- simulation_list$sp[mod_fls$id ]
mod_fls <- mod_fls %>% arrange(id)%>% filter(species != "Calidris ruficollis" )


stsTab <- parallel::mclapply(unique(mod_fls$species), function(sp) {
  
  out <- mod_fls %>% filter (mod_fls$species == sp)
  
  stsTab <- parallel::mclapply(1:nrow(out), function(i) {
    
    load(paste0("sdpModel/simu_out/", out$file[i]))
    
    if(!is.null(simuList[[1]]) & !is.null(simuList[[2]])) {
      
      stsList <- lapply(simuList, function(a) {
        
        alive  <- apply(a[,5,], 1, function(x) all(x==0))
        sts0   <- a[alive,2,]
        
        if(sum(alive)>1) {
          if(all(sts0[,179] == simulation_list$breeding_id[out$id[i]])){
            
            sts    <- t(apply(sts0, 1, function(x) {
              start <- simulation_list$winter_id[out$id[i]]
              x[1:max(which(x==start))] <- start
              x}))
          }else {
            sts    <- apply(sts0, 1, function(x) {
              if(x[179] == simulation_list$breeding_id[out$id[i]]){
                start <- simulation_list$winter_id[out$id[i]]
                x[1:max(which(x==start))] <- start
                x}}) %>% Reduce("rbind",.)
          }
          
        } else NULL
        
      })
      
    } else list(NULL, NULL)
    
    
  }, mc.cores = 30)
  
  stsPast <- lapply(stsTab, function(x) x[[1]]) %>% Reduce("rbind",.)
  stsCurr <- lapply(stsTab, function(x) x[[2]]) %>% Reduce("rbind",.)
  
  save(stsPast, file = gsub(" ", "_", glue::glue("IBM/infection_model/{sp}_sts_past.rda")))
  save(stsCurr, file = gsub(" ", "_", glue::glue("IBM/infection_model/{sp}_sts_curr.rda")))
  
}, mc.cores = 8)


####how many birds fly through cell 353?

index <- c(22,19,17,16,15,14,13,12,11,8,5,2,1)

df <- tibble(index = index,
             n= 0,
             Q1 = 0,
             Q3 = 0,
             median = 0)

for (x in 1:length(index)){
  load(fls$files[index[x]])
  i = rowSums(stsCurr == 353, na.rm = TRUE)
  df$n[x] <- sum(table(i)[2:length(table(i))])
  df$median[x] <- quantile(i[i!=0], probs = 0.5, na.rm = TRUE)
  df$Q1[x] <- quantile(i[i!=0], probs = 0.25, na.rm = TRUE)
  df$Q3[x] <- quantile(i[i!=0], probs = 0.75, na.rm = TRUE)}


######################### run infection spread model 

fls  <- tibble(files = list.files('IBM/infection_model', pattern = 'curr.rda', full.names = T))
fls$Species <- gsub("IBM/infection_model/(.*?)_sts_curr\\.rda$", "\\1", fls$files)


infection_results <- lapply(fls$Species, function(sp){
  
  if(!file.exists(paste0('IBM/infection_model/Results/', sp, '_array.rda'))) {
    
    load(fls$files[fls$Species == sp])
    
    simArray <- abind::abind(list(stsCurr, stsCurr), along = 3)
    simArray[,,2] <- 0
    
    infSite <- 353
    infProb <- 0.1
    
    medInf <- 4.61
    sdInf  <- 1.17
    
    for(t in 1:dim(simArray)[2]) {
      
      ind <- which(simArray[,t,1] ==infSite & simArray[,t,2]==0 & runif(dim(simArray)[1]) < infProb)
      
      if(length(ind)>0) {
        
        newInf <- parallel::mclapply(ind, function(x) {
          out <- simArray[x,,2]
          out[t:min((t+round(rnorm(1, medInf, sdInf), 0)),dim(simArray)[2])] <- 1
          out
        }, mc.cores = parallel::detectCores()) %>% Reduce("rbind",.)
        
        simArray[ind,,2] <- newInf
      }
      
    }
    
    save(simArray, file = glue::glue("infection_model/Results3/{sp}_array.rda"))
    
    load(glue::glue("infection_model/Results2/{sp}_array.rda"))
    
    cl <- parallel::makeCluster(parallel::detectCores()-2)
    invisible(parallel::clusterEvalQ(cl, {
      library(dplyr)
      library(sf)
    }))
    
    parallel::clusterExport(cl, c("simArray", "mudflatTab"))
    
    load(glue::glue("IBM/infection_model/Results2/{sp}_array.rda"))
    
    if(any(simArray[,,2]!=0)) {
      
      tracks <- lapply(which(apply(simArray, 1, function(u) any(u==1))), function(x) {
        
        tab <- tibble(sites = simArray[x,,1], inf = simArray[x,,2]) %>% mutate(time = 1:nrow(.)) %>% filter(inf>0)
        if(any(diff(tab$time)>1))  tab <- tab[1:min(which(diff(tab$time)>1)),]
        
        (mudflatTab %>% st_centroid())[tab$sites[!is.na(tab$sites)],] %>% mutate(id = x, time = tab$time[!is.na(tab$sites)]) %>%
            dplyr::select(id, time) %>% filter(!duplicated(tab$sites[!is.na(tab$sites)])) %>% suppressWarnings()
        
      })
      parallel::stopCluster(cl)
      
      tracks_filter <- sapply(tracks, function(x) nrow(x)>1)
      
      track <- tracks[tracks_filter] %>% do.call("rbind",.)
      
      linestr <- tracks[tracks_filter] %>% do.call("rbind",.) %>% group_by(id) %>%
        summarise(do_union = FALSE) %>%
        st_cast("LINESTRING")
      
      save(track,file = glue::glue("IBM/infection_model/Results2/{sp}_tracks.rda") )
      save(linestr,file = glue::glue("IBM/infection_model/Results2/{sp}_linestr.rda") )
    } else NULL
    
  }else NULL
  
})


######### create df from model output
spTab <- readxl::read_xlsx('sdpModel/species_parameters/spTab2.xlsx') %>% filter(Filter==TRUE)

# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data =mudflatTab$geometry[infection_df$InfectedPositions[infection_df$InfectedPositions != 353]], mapping = aes(fill=infection_df$n[infection_df$InfectedPositions != 353]), alpha = 0.5) +
#   geom_sf(data = mudflatTab$geometry[353],fill = "red",color = "red", alpha = 0.5)


fls  <- tibble(files = list.files('IBM/infection_model/Results2', pattern = 'tracks.rda', full.names = T))
fls$sp <- gsub("IBM/infection_model/Results2/(.*?)_tracks\\.rda$", "\\1", fls$files)
fls$sp <- gsub ("_", " ", fls$sp)

all_tracks<- lapply(fls$files,function(file) {
  load(file)
  tracks <- track %>% mutate(lbm = spTab$LeanBodyMass[spTab$ScieName == fls$sp[fls$files == file]],
                             sp = fls$sp[fls$files == file],
                             mud_id = sapply(track$geometry, function (x) { mud <- which(st_intersects(x,mudflatTab$geometry, sparse = FALSE))})) 
  return(tracks)
}) %>% Reduce("rbind",.)

### calculate directions
library(units)

track_bear <- all_tracks %>%
  st_transform(crs = st_crs("+proj=longlat +datum=WGS84")) %>% 
  group_split(id,lbm) %>% lapply(., function(x) { 
    x <- x %>% mutate(bearing = NA, dist= NA)
    for (i in 1:nrow(x)) {
      if(i>1) {
        out <-Reduce("rbind", list(c(x$geometry[1], x$geometry[i]))) 
        x[i,]$bearing <- lwgeom::st_geod_azimuth(out) %>% units::set_units(., "degrees")
        x[i,8] <- drop_units(st_distance(x$geometry[i], x$geometry[1]))/1000
        x}
    }
    return(x)
  }) %>% Reduce("rbind", .)

#### how many breeding sites with infected individuals


breed_infect <- length(unique(simulation_list$breeding_id[simulation_list$breeding_id %in% (all_tracks %>% filter(mud_id != 353))$mud_id]))

total_infect <- length(unique((all_tracks %>% filter(mud_id != 353))$mud_id))


########### distance to bm 

plot_data <- track_bear %>% group_by(id) %>% filter(row_number() == n()) %>% select(lbm,dist, sp) %>% mutate(sp_eng = spTab$Species[match(sp,spTab$ScieName)]) %>% arrange(lbm)

mod <- glm(dist ~ lbm, data = plot_data, family = quasipoisson)
summary(mod)


distance_spread <- ggplot(data = plot_data, aes(x = lbm, y = dist)) +
  geom_point(aes(color = sp_eng), position = position_dodge(width = 2), alpha = 0.6) +
  geom_smooth(method = "glm", formula = y~x,method.args = list(family = quasipoisson), se = TRUE, color = "black") +
  scale_color_viridis_d(option = "plasma", limits = unique(plot_data$sp_eng)) +
  labs(x = "Lean Body Mass in g", y = "Distance in km", color = "Species") +
  ggtitle("Max. Distances from Point of Infection") +
  theme_bw()


plot_data[which.min(plot_data$dist), ]
plot_data[which.max(plot_data$dist), ]

summary(plot_data$dist)

hist_data <- hist(plot_data$dist, plot = FALSE)

histogram <- hist(plot_data$dist[plot_data$sp_eng == "Bar-tailed Godwit"], main = "Distribution of Distances Travled While Infected for Bar-tailed Godwit", 
                  xlab ="Distances in km")

####### circular plot
library(circular)
angles  <- circular((track_bear$bearing), units = "degrees")

anglVec <- tibble(bearing = as.numeric(angles)) %>% filter(!is.na(bearing)) %>%
  mutate(bearing = ifelse(bearing>0, bearing, bearing+360)) %>%
  mutate(bin = cut(bearing, seq(0,360,20), lables = T))%>% group_by(bin) %>%
  summarise(freq = n()) %>%
  mutate(x_bin = sapply(strsplit(gsub("]", "", gsub("[(]","", as.character(bin))), ","), 
                        function(x) mean(as.numeric(c(x[[2]], x[[1]])))))

spread_direc <- ggplot(anglVec, aes(x = x_bin, y = freq)) +
  geom_bar(stat = 'identity', col = "black", fill= "#CD5555") +
  coord_polar()+
  scale_x_continuous(breaks = c(0,90,180,360), labels = c("N","E","S","W"))+
  labs(title= "Directions of Infection Spread")+
  theme_bw()

###### map of spread
fls  <- tibble(files = list.files('IBM/infection_model/Results2', pattern = 'linestr.rda', full.names = T))
fls$sp <- gsub("IBM/infection_model/Results2/(.*?)_linestr\\.rda$", "\\1", fls$files)
fls$sp <- gsub ("_", " ", fls$sp)

all_linestr <- lapply(fls$files,function(file) {
  load(file)
  linestr <- linestr %>% mutate(lbm = spTab$LeanBodyMass[spTab$ScieName == fls$sp[fls$files == file]])
  return(linestr)
}) %>% Reduce("rbind", . )


spread_map <- ggplot() +
  geom_sf(data = eaafMap$map) +
  geom_sf(data = all_linestr, alpha = 0.2) +
  geom_sf(data = track_bear %>% st_transform(., st_crs(all_tracks)) %>% filter(!is.na(bearing))%>% 
            group_by(geometry) %>% summarize(count = n()), color = "#CD5555", 
          aes(size = count), alpha = 0.8) +
  geom_sf(data = mudflatTab$geometry[353], color = "black", fill = "#EEDC82") +
  labs(size = "Number 
of Birds") +
  theme_bw() +
  theme(legend.position = c(0.89,0.12), legend.box = "horizontal" ,
        legend.background = element_rect(fill = alpha("white", 0.5)),
        text = element_text(size = 10))    # Position und Box der Legende ändern


