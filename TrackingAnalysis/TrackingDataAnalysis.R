#### Script for Tracking Data Analysis #######



### MERGE TRACKS ################################################################################################################################
library(tidyverse)
library(sf)

path <- "TrackingAnalysis/TrackingData" 
load ("TrackingAnalysis/TrackingData/breeedTab.RData")

breedTab <- breedTab  %>%
  select(species, id, arr_breed, breed_lon, breed_lat, -geometry) %>%
  mutate(Departure = arr_breed + days(2), Type = 2) %>% mutate (Days = as.numeric(difftime(Departure, arr_breed, units = 'days')))%>%
  rename(Species = species, Arrival = arr_breed, Lon = breed_lon, Lat = breed_lat, ID = id)


###### Ruddy turnstone 

fls <- c(list.files(glue::glue("{path}/Turnstone"), pattern = "movementSummary", recursive = T, full.names = T))

tracks <- lapply(fls, function(x) {
  id <- gsub("_Grouped_movementSummary.csv", "", sapply(strsplit(x, "/"), function(y) y[length(y)]))
  read.csv(x) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(StartTime), Departure = as.POSIXct(EndTime)) %>% ### transfer into dates
    mutate(Species = "Ruddy Turnstone", ID = id) %>% rename(Lon = Lon.50., Lat= Lat.50.) %>% mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days)))  %>% filter(Days>0 | is.na(Days), ID != "WA",ID != "WT", ID != "YD", ID!= "ZSC", ID != "ZSJ", ID != "ZV", ID != "ZY") %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) %>% mutate(Days = ifelse(Type==0, NA, Days))
}) %>% do.call("rbind", .)

save(tracks, file = glue::glue("{path}/Turnstone_Tracks.rda"))

###### Sanderling 

fls <- list.files(glue::glue("{path}/Sanderling"), pattern = "MovementSummary", recursive = T, full.names = T)

tracks <- lapply(fls, function(x) {
  id <- gsub("_Grouped_MovementSummary.csv", "", sapply(strsplit(x, "/"), function(y) y[length(y)]))
  read.csv(x) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(Arrival), Departure = as.POSIXct(Departure),
           Days = as.numeric(difftime(Departure, Arrival, units = 'days'))) %>% ### transfer into dates
    mutate(Species = "Sanderling", ID = id) %>% rename(Type = type) %>% mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days)))  %>% filter(Days>0 | is.na(Days), ID != "2030") %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) 
}) %>% do.call("rbind", .)

save(tracks, file = glue::glue("{path}/Sanderling_Tracks.rda"))

###### Curlew Sandpiper 

fls <- list.files(glue::glue("{path}/CurlewSandpiper"), pattern = "Grouped_movementSummary", recursive = T, full.names = T)

tracks <- lapply(fls, function(x) {
  id <- gsub("_Grouped_movementSummary.csv", "", sapply(strsplit(x, "/"), function(y) y[length(y)]))
  read.csv(x) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(StartTime), Departure = as.POSIXct(EndTime),
           Days = as.numeric(difftime(EndTime, StartTime, units = 'days'))) %>% ### transfer into dates
    mutate(Species = "Curlew Sandpiper", ID = id) %>% rename(Lon = Lon.50., Lat= Lat.50.) %>%mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days))) %>% filter(is.na(Days) | Days > 0) %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) 
}) %>% do.call("rbind", .)

save(tracks, file = glue::glue("{path}/CurlewSandpiper_Tracks.rda"))

###### Godwit 

fls <- c(list.files(glue::glue("{path}/Godwit"), pattern = "Grouped_MovementSummary", recursive = T, full.names = T),
         list.files(glue::glue("{path}/Godwit"), pattern = "Grouped_movementSummary", recursive = T, full.names = T))

tracks <- lapply(fls, function(x) {
  id  <- sapply(strsplit(sapply(strsplit(x, "/"), function(y) y[length(y)]), "_"), function(y) y[1])
  tmp <- read.csv(x) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(Arrival.50.), Departure = as.POSIXct(Departure.50.),
           Days = as.numeric(difftime(Departure.50., Arrival.50., units = 'days'))) %>% ### transfer into dates
    mutate(Species = "Bar-tailed Godwit", ID = id) %>% rename(Lon = Lon.50., Lat= Lat.50.) %>%mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days))) %>% filter(is.na(Days) | Days > 0, ID!= "H373") %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) 
}) %>% do.call("rbind", .)

save(tracks, file = glue::glue("{path}/Godwit_Tracks.rda"))


###### Great knot 

fls <- list.files(glue::glue("{path}/GreatKnot"), pattern = "MovementSummary", recursive = T, full.names = T)

tracks <- lapply(fls, function(x) {
  id <- gsub("_Grouped_MovementSummary.csv", "", sapply(strsplit(x, "/"), function(y) y[length(y)]))
  read.csv(x) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(Arrival), Departure = as.POSIXct(Departure),
           Days = as.numeric(difftime(Departure, Arrival, units = 'days'))) %>% ### transfer into dates
    mutate(Species = "Great Knot", ID = id) %>% rename(Type = type) %>% mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days))) %>% filter(is.na(Days) | Days > 0) %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) 
}) %>% do.call("rbind", .)

save(tracks, file = glue::glue("{path}/GreatKnot_Tracks.rda"))


###### Red knot 

fls <- list.files(glue::glue("{path}/RedKnot"), pattern = "MovementSummary", recursive = T, full.names = T)

tracks <- lapply(fls, function(x) {
  id <- gsub("_MovementSummary.csv", "", sapply(strsplit(x, "/"), function(y) y[length(y)]))
  read.csv(x, check.names = TRUE) %>% as_tibble() %>%
    mutate(Arrival = as.POSIXct(Arrival.50., format = "%Y-%m-%d %H:%M:%S"), Departure = as.POSIXct(Departure.50., format = "%Y-%m-%d %H:%M:%S"),
           Days = as.numeric(difftime(Departure.50., Arrival.50., units = 'days'))) %>% ### transfer into dates
    mutate(Species = "Red Knot", ID = id) %>% rename(Lon =Lon.50., Lat=Lat.50.) %>% mutate(Days = ifelse(Type==0, NA,ifelse (Type==2,NA,Days))) %>% filter(is.na(Days) | Days > 0) %>%
    dplyr::select(Species, ID, Arrival, Departure, Days, Lon, Lat, Type) 
}) %>% do.call("rbind", .)

common_ids <- intersect(tracks$ID, breedTab$ID)

breedTabRK <- st_drop_geometry(breedTab) %>% filter(Species == "RedKnot") %>% mutate (Species = "Red Knot")

# adding artificial breeding sites to missing ids from breedTab
for (id in common_ids) {
  breed_row <- breedTabRK[breedTabRK$ID == id, ]
  if (!any(tracks$Type == 2 & tracks$ID == id)) {
    type1_idx <- max(which(tracks$Type == 1 & tracks$ID == id))
    type3_idx <- min(which(tracks$Type == 3 & tracks$ID == id))
    tracks <- rbind(tracks[1:type1_idx, ], breed_row, tracks[(type1_idx + 1):type3_idx, ], tracks[(type3_idx + 1):nrow(tracks), ])
  }
}

save(tracks , file = glue::glue("{path}/RedKnot _Tracks.rda"))


## merge everything

track_fls <- list.files(path, pattern = "Tracks.rda", full.names = T) 
bird_tracking_data <- lapply(track_fls, function(x) {load(x); tracks})  %>% do.call("rbind", .) %>% filter(ID != 'BQ916', ID != "D148", ID!= "W552", ID!= "W539") 
#filter doubeling IDs
bird_tracking_data <-bird_tracking_data[!grepl("\\(.*\\)", bird_tracking_data$ID), ]

save(bird_tracking_data, file = glue::glue("{path}/bird_tracking_data.rda"))


###### ANALYSIS OF TRACKING DATA #############################################################################################################
library(dplyr)
library(geosphere)
library(tidyr)

path <- "TrackingAnalysis/TrackingData"  

load(glue::glue("{path}/bird_tracking_data.rda"))
spTab <- readxl::read_xlsx('sdpModel/species_parameters/spTab2.xlsx') %>% filter(Filter==TRUE)


### Preparation/Calculation of Parameters #####

### Step length
step_length <- bird_tracking_data %>%
  group_split(ID) %>%
  lapply(., function(df) {
    out <- df %>% mutate(dist = suppressWarnings(distHaversine(cbind(Lon, Lat), cbind(lag(Lon), lag(Lat)))/1000),
                         lbm = spTab$LeanBodyMass[spTab$Species == .$Species[1]],
                         sci_sp = spTab$ScieName[spTab$Species == .$Species[1]])
    do.call('rbind', lapply(1:2, function(x) {
      if (x == 1) {
        if (any(out$Type == 2)) {
          tibble(ID = out$ID[1], Sp = out$Species[1], Sci_sp =out$sci_sp[1], lbm = out$lbm[1], Dir = x, Step = out$dist[min(which(out$Type == 1)):min(which(out$Type == 2))], Days = out$Days[min(which(out$Type == 1)):min(which(out$Type == 2))])
        } else NULL
      } else {
        if (any(out$Type == 2) & out$Type[nrow(out)] == 0) {
          sm <- out[(max(which(out$Type == 2)) + 1):nrow(out), ]
          tibble(ID = sm$ID[1], Sp = sm$Species[1], Sci_sp =sm$sci_sp[1], lbm = sm$lbm[1], Dir = x, Step = sm$dist[1:min(which(sm$Type == 0))], Days = sm$Days[1:min(which(sm$Type == 0))])
        } else NULL
      }
    }))
  }) %>%
  do.call('rbind', .)


### migration duration
mig_dur <- bird_tracking_data %>%
  group_split(ID) %>%
  lapply(., function(df) {
    out <- df %>% mutate(lbm = spTab$LeanBodyMass[spTab$Species == .$Species[1]],
                         sci_sp = spTab$ScieName[spTab$Species == .$Species[1]])
    do.call('rbind', lapply(1:2, function(x) {
      if(x==1) {
        if(any(out$Type==2)) {
          tibble(ID = out$ID[1], Sp = out$Species[1], Dir = x, Sci_sp = out$sci_sp[1], lbm = out$lbm[1],
                 Dur = as.numeric(difftime(out$Arrival[min(which(out$Type==2))], out$Departure[min(which(out$Type==1))-1], units = 'days')))
        } else NULL
      } else {
        if(any(out$Type==2) & out$Type[nrow(out)]==0) {
          sm <- out[(max(which(out$Type==2))):nrow(out),]
          tibble(ID = out$ID[1], Sp = out$Species[1], Dir = x, Sci_sp = sm$sci_sp[1], lbm = sm$lbm[1],
                 Dur = as.numeric(difftime(sm$Arrival[min(which(sm$Type==0))], sm$Departure[1], units = 'days')))
        } else NULL
      }
    }))
  }) %>%
  do.call('rbind', .)


mig_dist <- aggregate(Step ~ ID+Sp+lbm+Dir, data = step_length, FUN = sum)
colnames(mig_dist) <- c("ID","Sp", "lean_body_mass","Dir", "dist")
#hist(mig_dist$dist,main = "Distribution of migration distance", xlab = "distance in km", ylab = "frequency")

step_length <- merge(step_length, mig_dist, by = c("ID", "Sp", "Dir")) %>% select (c("ID","Sp","Sci_sp","lbm","Step","Days","dist", "Dir"))

mig_dur <- merge(mig_dur, mig_dist, by= c("ID","Sp", "Dir")) %>% select (c("ID","Sp","Sci_sp","lbm","dist","Dur", "Dir"))


save(step_length, file = glue::glue("{path}/step_length_df.rda"))
save(mig_dur, file = glue::glue("{path}/mig_dur_df.rda"))

load(glue::glue("{path}/step_length_df.rda"))
load(glue::glue("{path}/mig_dur_df.rda"))

####### PLOTTING STRATEGIES ################################################################################################################

library(ggplot2)
library(dplyr)


##### Step length ######################################

data <- step_length 

median <- data %>%
  group_by(lbm, Dir) %>%
  summarise(Q1 = quantile(Step, probs = 0.25, na.rm = TRUE),
            Q3 = quantile(Step, probs = 0.75, na.rm = TRUE),
            Step = quantile(Step, probs = 0.5, na.rm = TRUE))

p1 <- ggplot(data, aes(x = lbm, y = Step)) +
  geom_point(
    aes(color = as.factor(Dir)),
    position = position_jitterdodge(jitter.width = 1, dodge.width = 6),
    alpha = 0.4, show.legend = TRUE, size = 1
  ) +
  geom_errorbar(
    data = median,
    aes(x = lbm, ymin = Q1, ymax = Q3, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    width = 0.5, size = 0.7, show.legend = TRUE
  ) +
  geom_point(
    data = median,
    aes(x = lbm, y = Step, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    size = 4, shape = 18, stroke = 0
  ) +
  labs(x = "Lean Body Mass in g", y = "Step Length in km", title = "Step Length", subtitle = "Tracking Data") +
  scale_color_manual(values = c("steelblue3", "sandybrown"), name= "Direction", labels = c("Northward", "Southward")) + 
  scale_y_continuous(limits = c(0, 12000)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(size = 0.25, linetype = 1),
    # panel.grid.major.x = element_line(size = 0.75, linetype = 2),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
## Above 90th  percentile of Step Length###################################

data <- step_length %>%
  group_by(ID,Dir) %>%
  filter(Step>quantile(Step, probs = 0.9))


median <- data %>%
  group_by(lbm, Dir) %>%
  summarise(Q1 = quantile(Step, probs = 0.25, na.rm = TRUE),
            Q3 = quantile(Step, probs = 0.75, na.rm = TRUE),
            Step = quantile(Step, probs = 0.5, na.rm = TRUE))

p3 <- ggplot(data, aes(x = lbm, y = Step)) +
  geom_point(
    aes(color = as.factor(Dir)),
    position = position_jitterdodge(jitter.width = 1, dodge.width = 6),
    alpha = 0.4, show.legend = TRUE, size = 1
  ) +
  geom_errorbar(
    data = median,
    aes(x = lbm, ymin = Q1, ymax = Q3, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    width = 0.5, size = 0.7, show.legend = TRUE
  ) +
  geom_point(
    data = median,
    aes(x = lbm, y = Step, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    size = 4, shape = 18, stroke = 0
  ) +
  labs(x = "Lean Body Mass in g", y = "Step Length in km", title = "Upper 10% of Step Lengths", subtitle ="Tracking Data") +
  scale_color_manual(values = c("steelblue3", "sandybrown"), name= "Direction", labels = c("Northward", "Southward")) + 
  scale_y_continuous(limits = c(0, 12000)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(size = 0.25, linetype = 1),
    # panel.grid.major.x = element_line(size = 0.75, linetype = 2),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12) 
  )

###Duration##################################
data <- mig_dur 

median <- data %>%
  group_by(lbm, Dir) %>%
  summarise(Q1 = quantile(Dur, probs = 0.25, na.rm = TRUE),
            Q3 = quantile(Dur, probs = 0.75, na.rm = TRUE),
            Dur = quantile(Dur, probs = 0.5, na.rm = TRUE))

p5 <- ggplot(data, aes(x = lbm, y = Dur)) +
  geom_point(
    aes(color = as.factor(Dir)),
    position = position_jitterdodge(jitter.width = 1, dodge.width = 6),
    alpha = 0.4, show.legend = TRUE, size = 1
  ) +
  geom_errorbar(
    data = median,
    aes(x = lbm, ymin = Q1, ymax = Q3, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    width = 0.5, size = 0.7, show.legend = TRUE
  ) +
  geom_point(
    data = median,
    aes(x = lbm, y = Dur, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    size = 4, shape = 18, stroke = 0
  ) +
  labs(x = "Lean Body Mass in g", y = "Time in days", title = "Migration Duration", subtitle ="Tracking Data") +
  scale_color_manual(values = c("steelblue3", "sandybrown"), name= "Direction", labels = c("Northward", "Southward")) + 
  scale_y_continuous(limits = c(0, 170)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(size = 0.25, linetype = 1),
    # panel.grid.major.x = element_line(size = 0.75, linetype = 2),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )


##Stopover Time ####################################################

data <- step_length 

median <- data %>%
  group_by(lbm, Dir)  %>%
  summarise(Q1 = quantile(Days, probs = 0.25, na.rm = TRUE),
            Q3 = quantile(Days, probs = 0.75, na.rm = TRUE),
            Days = quantile(Days, probs = 0.5, na.rm = TRUE))

p7 <- ggplot(data, aes(x = lbm, y = Days)) +
  geom_point(
    aes(color = as.factor(Dir)),
    position = position_jitterdodge(jitter.width = 1, dodge.width = 6),
    alpha = 0.4, show.legend = TRUE, size = 1
  ) +
  geom_errorbar(
    data = median,
    aes(x = lbm, ymin = Q1, ymax = Q3, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    width = 0.5, size = 0.7, show.legend = TRUE
  ) +
  geom_point(
    data = median,
    aes(x = lbm, y = Days, color = as.factor(Dir)),
    position = position_dodge(width = 6),
    size = 4, shape = 18, stroke = 0
  ) +
  labs(x = "Lean Body Mass in g", y = "Time in days", title = "Stopover Time", subtitle = "Tracking Data") +
  scale_color_manual(values = c("steelblue3", "sandybrown"), name= "Direction", labels = c("Northward", "Southward")) + 
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(size = 0.25, linetype = 1),
    # panel.grid.major.x = element_line(size = 0.75, linetype = 2),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
####### STATISTICAL ANALYSIS###############################################################################################################
path <- "TrackingAnalysis/TrackingData"  

load(glue::glue("{path}/bird_tracking_data.rda"))
spTab <- readxl::read_xlsx('sdpModel/species_parameters/spTab2.xlsx') %>% filter(Filter==TRUE)

library(tidyverse)
library(ggplot2)

load(glue::glue("{path}/step_length_df.rda"))
load(glue::glue("{path}/mig_dur_df.rda"))

##step length

data <- data.frame(step_length) %>%
  setNames(c("ID", "SP", "Scie_SP", "BM", "Step", "Days", "Dist", "Dir")) %>%
  as_tibble() %>%
  left_join(mig_dur %>% distinct(ID, Dir, Dur), by = c("ID", "Dir")) %>%
  mutate(logBM = log(BM))

bodySize <- sort(unique(data$BM))

model <- list(nw = lm(Step ~ logBM*Dist, data = data %>% filter(Dir==1)),
              sw = lm(Step ~ logBM*Dist, data = data %>% filter(Dir==2)))

summary(model[[1]]) ##model[[2]]
plot(model[[1]])##model[[2]]

pred <- lapply(1:2, function(x) {
  newdat  = data.frame(logBM = log(seq(15, 250, length=100)), Dist = rep(median(data$Dist), 100))
  cbind(newdat, Dir = x, sig = summary(model[[x]])$coefficients[2,4], predict(model[[x]], newdata = newdat, interval="confidence")) %>% as_tibble()
}) %>% Reduce('rbind',.)

pred$rounded_sig <- ifelse(pred$sig < 0.009, sprintf("%.3e", as.numeric(pred$sig)), sprintf("%.6f", as.numeric(pred$sig)))
pred <- filter(pred, pred$Dir == 1)

p2 <- ggplot(pred, aes(x = exp(logBM), y = fit)) +
  geom_line(mapping = aes(linetype = as.factor(rounded_sig), color = as.factor(Dir)), size = 1) +
  geom_ribbon(data = pred, mapping = aes(ymin = lwr, max = upr, fill = as.factor(Dir)), alpha=0.2, colour= NA) +
  scale_x_log10(breaks = bodySize[bodySize != 84]) +
  scale_color_manual(values=c("steelblue3"), guide = guide_legend(title = "Direction"), labels = c("Northward"))+
  scale_fill_manual(values=c("steelblue3"),guide = guide_legend(title = "Direction"),labels = c("Northward"))+
  labs(x= "Lean Body Mass in g", y= "Step Length in km", subtitle = "Linear Model Predictions", title ="") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.25,linetype = 1),
        panel.grid.major.x = element_line(size = 0.75, linetype = 2),
        axis.text.x = element_text(size = 12),
        axis.title = element_text (size =12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size =12))+
  coord_cartesian(ylim = c(0, max(pred$upr) + 0.8))





###max step_length

data <- data.frame(step_length) %>%
  setNames(c("ID", "SP", "Scie_SP", "BM", "Step", "Days", "Dist", "Dir")) %>%
  as_tibble() %>%
  left_join(mig_dur %>% distinct(ID, Dir, Dur), by = c("ID", "Dir")) %>%
  mutate(logBM = log(BM))

data <- data %>%
  group_by(ID, Dir) %>%
  filter(Step>quantile(Step, probs = 0.9))



bodySize <- sort(unique(data$BM))

model <- list(nw = lm(Step ~ logBM*Dist, data = data %>% filter(Dir==1)),
              sw = lm(Step ~ logBM*Dist, data = data %>% filter(Dir==2)))


summary(model[[1]]) ##model[[2]]
plot(model[[1]])##model[[2]]

pred <- lapply(1:2, function(x) {
  newdat  = data.frame(logBM = log(seq(15, 250, length=100)), Dist = rep(median(data$Dist), 100))
  cbind(newdat, Dir = x, sig = summary(model[[x]])$coefficients[2,4], predict(model[[x]], newdata = newdat, interval="confidence")) %>% as_tibble()
}) %>% Reduce('rbind',.)

pred$rounded_sig <- ifelse(pred$sig < 0.009, sprintf("%.3e", as.numeric(pred$sig)), sprintf("%.6f", as.numeric(pred$sig)))
pred <- filter(pred, pred$Dir == 1)

p4 <- ggplot(pred, aes(x = exp(logBM), y = fit)) +
  geom_line(mapping = aes(linetype = as.factor(rounded_sig), color = as.factor(Dir)), size = 1) +
  geom_ribbon(data = pred, mapping = aes(ymin = lwr, max = upr, fill = as.factor(Dir)), alpha=0.2, colour= NA) +
  scale_x_log10(breaks = bodySize[bodySize != 84]) +
  scale_color_manual(values=c("steelblue3"), guide = guide_legend(title = "Direction"), labels = c("Northward"))+
  scale_fill_manual(values=c("steelblue3"),guide = guide_legend(title = "Direction"),labels = c("Northward"))+
  labs(x= "Lean Body Mass in g", y= "Step Length in km", subtitle = "Linear Model Predictions", title ="") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.25,linetype = 1),
        panel.grid.major.x = element_line(size = 0.75,linetype = 2),
        axis.text.x = element_text(size = 12),
        axis.title = element_text (size =12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size =12))+
  coord_cartesian(ylim = c(0, max(pred$upr) + 0.8))

###migration duration

data <- data.frame(step_length) %>%
  setNames(c("ID", "SP", "Scie_SP", "BM", "Step", "Days", "Dist", "Dir")) %>%
  as_tibble() %>%
  left_join(mig_dur %>% distinct(ID, Dir, Dur), by = c("ID", "Dir")) %>%
  mutate(logBM = log(BM))


model <- list(nw = lm(Dur ~ logBM*Dist, data = data %>% filter(Dir==1)),
              sw = lm(Dur ~ logBM*Dist, data = data %>% filter(Dir==2)))
summary(model[[1]])

pred <- lapply(1:2, function(x) {
  newdat  = data.frame(logBM = log(seq(15, 250, length=100)), Dist = rep(median(data$Dist), 100))
  cbind(newdat, Dir = x, sig = summary(model[[x]])$coefficients[2,4], predict(model[[x]], newdata = newdat, interval="confidence")) %>% as_tibble()
}) %>% Reduce('rbind',.)

pred$rounded_sig <- ifelse(pred$sig < 0.009, sprintf("%.3e", as.numeric(pred$sig)), sprintf("%.6f", as.numeric(pred$sig)))
pred <- filter(pred, pred$Dir == 2)

p6 <- ggplot(pred, aes(x = exp(logBM), y = fit)) +
  geom_line(mapping = aes(linetype = as.factor(rounded_sig), color = as.factor(Dir)), size = 1) +
  geom_ribbon(data = pred, mapping = aes(ymin = lwr, max = upr, fill = as.factor(Dir)), alpha=0.2, colour= NA) +
  scale_x_log10(breaks = bodySize[bodySize != 84]) +
  scale_color_manual(values=c("sandybrown"), guide = guide_legend(title = "Direction"), labels = c("Southward"))+
  scale_fill_manual(values=c("sandybrown"),guide = guide_legend(title = "Direction"),labels = c("Southward"))+
  labs(x= "Lean Body Mass in g", y= "Time in days", subtitle = "Linear Model Predictions", title ="") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.25,linetype = 1),
        panel.grid.major.x = element_line(size = 0.75,linetype = 2),
        axis.text.x = element_text(size = 12),
        axis.title = element_text (size =12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size =12))+
  coord_cartesian(ylim = c(0, max(pred$upr) + 0.8))

###Stop over time

data <- data.frame(step_length) %>%
  setNames(c("ID", "SP", "Scie_SP", "BM", "Step", "Days", "Dist", "Dir")) %>%
  as_tibble() %>%
  left_join(mig_dur %>% distinct(ID, Dir, Dur), by = c("ID", "Dir")) %>%
  mutate(logBM = log(BM))

model <- list(nw = lm(Days ~ logBM*scale(Dist), data = data %>% filter(Dir==1)),
              sw = lm(Days ~ logBM*scale(Dist), data = data %>% filter(Dir==2)))
summary(model[[1]]) ##model[[2]]
plot(model[[1]])##model[[2]]

pred <- lapply(1:2, function(x) {
  newdat  = data.frame(logBM = log(seq(15, 250, length=100)), Dist =rep(median(data$Dist), 100))
  cbind(newdat, Dir = x, sig = summary(model[[x]])$coefficients[2,4], predict(model[[x]], newdata = newdat, interval="confidence")) %>% as_tibble()
}) %>% Reduce('rbind',.)

pred$rounded_sig <- sprintf("%.3e", as.numeric(pred$sig))

p8<- ggplot(pred, aes(x = exp(logBM), y = fit, group = Dir, color = Dir)) +
  geom_line(mapping = aes(linetype = as.factor(rounded_sig), color = as.factor(Dir)), size = 1) +
  geom_ribbon(data = pred, mapping = aes(ymin = lwr, max = upr, fill = as.factor(Dir)), alpha=0.2, colour= NA) +
  scale_x_log10(breaks = bodySize[bodySize != 84]) +
  scale_color_manual(values=c("steelblue3", "sandybrown"), guide = guide_legend(title = "Direction"), labels = c("Northward", "Southward"))+
  scale_fill_manual(values=c("steelblue3", "sandybrown"),guide = guide_legend(title = "Direction"),labels = c("Northward", "Southward"))+
  scale_linetype_manual(values = c("solid", "solid"), guide = guide_legend(title = "Significance")) +
  labs(x= "Lean Body Mass in g", y= "Time in days", subtitle = "Linear Model Predictions", title ="") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.25,linetype = 1),
        panel.grid.major.x = element_line(size = 0.75,linetype = 2),
        axis.text.x = element_text( size = 12),
        axis.title = element_text (size =12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size =12))+
  coord_cartesian(ylim = c(0, max(pred$upr) + 0.8))

########### merge Plot 

library(cowplot)
library(grid)
library(gridExtra)

mylegend<-get_legend((p1+ theme(legend.direction = "horizontal")))

strategy_plot <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p2 + theme(legend.position="none"),
                                          p3 + theme(legend.position="none"),
                                          p4 + theme(legend.position="none"),
                                          p5 + theme(legend.position="none"),
                                          p6 + theme(legend.position="none"),
                                          p7 + theme(legend.position="none"),
                                          p8 + theme(legend.position="none"),nrow=4),
                              mylegend, nrow=2,heights=c(10, 1), top=  textGrob("Effects of Body Size on Migration Factors",
                                                                                gp = gpar(fontsize = 15, fontface = "bold")))

print(strategy_plot)

