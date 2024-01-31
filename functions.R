library(bbmle)
library(scales)

## Scaling functions

shorebirdScaling     <- function(lbm) {
  
  out <- list(lbm   = lbm,
              mbm   = predict(lm(mbm~lbm, data = tibble(lbm=c(15, 300), mbm = c(23.05, 643.4))),
                              newdata = tibble(lbm = lbm)))
  
  ### Fuel deposition rate
  out$X1gX  <- (out$mbm-out$lbm)/100
  out$X1xkJ <- out$X1gX*32.77
  out$FDR   <- 2.37*(out$lbm/1000)^-0.27
  out$FDRg  <- (out$FDR*(out$lbm))/100
  out$FDRx  <- out$FDRg/out$X1gX
  out$daysToRefuel <- (out$mbm-out$lbm)/out$FDRg
  
  ### Metabolic rate
  out$DER    <- 912*(lbm/1000)^0.704
  out$BMR    <- 5.06*(lbm/1000)^0.729
  out$cond   <- 0.0067*(lbm)^0.4384
  out$Kesm   <- 41-(out$BMR/(1.2*out$cond))
  out$W      <- 57.3*(lbm/1000)^0.81
  
  ### Flight parameters
  x1  <- c(250, 144, 110, 20)
  dst <- c(14000, 7500, 6500, 3000) 
  fit1 <- lm(log(dst)~log(x1))
  pr   <- exp(predict(fit1))
  out$FlightCap <- exp(predict(fit1, newdata = data.frame(x1 = lbm))) + 1000
  
  # Energy exp ~ temp
  t  <- seq(-35, 45, length = 100)
  te <- 41 - (out$BMR/(1.2*out$cond))
  b = out$BMR - (-out$cond)*te
  
  tm <- data.frame(tm = seq(-35, te, length = 100), W = -out$cond*seq(-25, te, length = 100) + b)
  tm <- rbind(tm, data.frame(tm = seq(te, 45, length = 100), W = rep(out$BMR, 100)))
  
  out$EEFnc  <- suppressWarnings(approxfun(x = tm[,1], y = 2*tm[,2]*86.4, rule = 3))
  
  out$speed  <- 1440
  out$c      <- out$FlightCap / (1-(1+99/100)^-0.5)
  
  return(out)
}

#### Plotting
x_site_plot <- function(simu) {
  
  require(ggplot2)
  require(zoo)
  require(patchwork)
  
  alive <- apply(simu[,5,], 1, function(x) all(x==0))
  
  x_traj <- t(apply(simu[alive,3,], 1, function(x) zoo::na.approx(x, rule = 2)))
  x_traj_tab <- tibble(id = 1:sum(alive)) %>% bind_cols(x_traj %>% as_tibble()) %>%
    pivot_longer(cols = -id, values_to = "x") %>% mutate(time = as.numeric(gsub("V", "", name))) %>% dplyr::select(id, time, x)
  
  pl1 <- ggplot(x_traj_tab, aes(x = time, y = x, group = id, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    scale_color_manual(values = rainbow(sum(alive))) +
    theme_light()
  
  site_traj <- simu[alive,2,] + 1
  site_traj_tab <- tibble(id = 1:sum(alive)) %>% bind_cols(site_traj %>% as_tibble()) %>%
    pivot_longer(cols = -id, values_to = "site") %>% mutate(time = as.numeric(gsub("V", "", name))) %>% 
    dplyr::select(id, time, site) %>% filter(!is.na(site))
  
  pl2 <- ggplot(site_traj_tab, aes(x = time, y = site, group = id, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    scale_color_manual(values = rainbow(sum(alive))) +
    theme_light()
  
  pl1 / pl2
  
}

site_matrix <- function(grid, breeding_index, simu, mod) {
  
  # grid <- mudflatTab$geometry
  # breeding_index <- simulation_list$breeding_id[i]
  
  require(ggplot2)

  alive  <- apply(simu[,5,], 1, function(x) all(x==0))
  sts0   <- simu[alive,2,] + 1  
  sts    <- t(apply(sts0, 1, function(x) {
    x[1:max(which(x==1))] <- 1
    x
  }))
  
  strategy <- t(apply(sts, 1, function(x) {
      x[x==1 | x==max(x, na.rm = T)] <- NA
      t <- as.data.frame(table(x[!is.na(x)])) %>% as_tibble()
      out <- sort(t$Freq, decreasing = T)
      c(out[1:25], rep(NA, 25-length(out[1:25])))
    }))
  
  steps_ind <- t(apply(sts, 1, function(sts_ind) {
    step  <- which(abs(diff(sts_ind[!is.na(sts_ind)]))>0)
    dists <- as.numeric(mod@Sites$dist[cbind(sts_ind[!is.na(sts_ind)][step], sts_ind[!is.na(sts_ind)][step+1])])
    c(dists[1:25], rep(NA, 25-length(dists[1:25])))
  }))
  
  mig_dur <- t(apply(sts, 1, function(x) c(min(which(x==max(x, na.rm = T))) - max(which(x==1)), rep(NA, 24))))
  
  site_index <- mod@Sites$index
  site_index[length(site_index)] <- breeding_index
  
  stay   <- apply(sts, 1, function(x) {
    as.data.frame(table(x)) %>% as_tibble() %>% setNames(c('site', 'days')) %>% mutate(site = as.numeric(as.character(site))) %>%
      full_join(tibble(site = 1:length(site_index)), by = "site") %>% arrange(site) %>% mutate(days = ifelse(site==1 | site==max(site), NA, days)) %>% pull(days)
  }) %>% apply(., 1, sum, na.rm = T)
  
  
  
  
  transM <- expand_grid(a = 1:length(site_index), b = 1:length(site_index)) %>% 
    left_join(apply(sts, 1, function(x) tibble(a = x[!is.na(x)][-sum(!is.na(x))], b = x[!is.na(x)][-1]) %>% filter(a!=b) %>% mutate(n = 1)) %>%
                do.call("rbind",.) %>% group_by(a,b) %>% summarise_at("n", sum), by = join_by(a, b)) %>% filter(n>0) %>%
    mutate(a_grid = site_index[a], b_grid = site_index[b])
    
    
  mat <- matrix(0, ncol = length(grid), nrow = length(grid))
  diag(mat)[site_index] <- stay
  
  mat[cbind(transM$a_grid, transM$b_grid)] <- transM$n
    
  list(network = mat, strategy = abind(list(strategy, steps_ind, mig_dur), along = 3))
}

network_plot <- function(grid, map, trans, title, subtitle, linksMin, col) {
  
  crds <- grid %>% st_centroid() %>% st_coordinates()
  
  lines<- trans
  diag(lines) <- 0  
  links <- tibble(a = rep(1:nrow(trans), each =  nrow(trans)), b = rep(1:nrow(trans), times = nrow(trans)), n = c(trans)) %>% 
  filter(n>0, a!=b, n>linksMin) %>% mutate(lon_1 = crds[a,1], lat_1 = crds[a,2], lon_2 = crds[b,1], lat_2 = crds[b,2])

  
  ggplot() +
    geom_sf(data = eaafMap$grid, color = 'grey90') +
    geom_sf(data = eaafMap$map, fill = 'grey90') +
    geom_curve(
      data = links,
      aes(
        x = lon_1, y = lat_1, xend = lon_2, yend = lat_2,
        linewidth = as.numeric(n), alpha = as.numeric(n)
      ),
      color = col, curvature = 0.1, lineend = "round"
    ) +
    geom_point(
      data = crds[diag(trans)>0,] %>% as_tibble(),
      aes(x = X, y = Y, size = diag(trans)[diag(trans)>0]),
      color = "black",shape = 21, fill = "grey30", stroke = 0.5
    ) +
    labs(x = "", y = "", title = paste(title), subtitle = paste(subtitle)) +
    scale_linewidth_continuous(range = c(0.2, 3), name = "number of individuals", labels = comma_format(big.mark =" "), guide = FALSE)+
    scale_alpha_continuous(range = c(0.4, 3), guide =FALSE)+
    scale_size_continuous(name = "Stopover time", labels = comma_format(big.mark =" "), guide = FALSE) +
    theme_minimal()
        
  
    
}


gaussMLE <- function(day, size, prob, thresh) {
  
  tab <- data.frame(day = day, size = size, p = prob)
  
  gauss.curve <- function(parms, intv = 1) {
    t <- seq(1, 366, intv)
    parms <- as.list(parms)
    fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
    fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
    c(fit1, fit2[-1])
  }
  
  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), 1)
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\r")
    -sum(dbinom(x = round(tab[,3]*tab[,2],0), size = rep(100, length(fit)), prob = fit[day], log=TRUE), na.rm=T)
  }
  
  mle <- suppressWarnings(mle2(gauss.loglik, method="L-BFGS-B",
                               start=list(a1 = 225, a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
                               lower=list(a1 = 50,  a2 = 5,   a3 = 0.5,a4=5,  a5 = 0.5),
                               upper=list(a1 = 225, a2 =  Inf,  a3 =  Inf, a4 =  Inf, a5 =  Inf),
  ))
  
  t <- seq(1, 366, 1)
  fit <- gauss.curve(coef(mle), intv = 1)
  
  start <- t[min(which(fit<thresh))]
  
  list(result = data.frame(week = t, fit = fit), start = start, end = end)
}


# serial_network_plot <- function(mod, simu) {
#   
#   require(ggplot2)
#   
#   alive  <- apply(simu[,5,], 1, function(x) all(x==0))
#   sts0   <- simu[alive,2,] + 1  
#   sts    <- t(apply(sts0, 1, function(x) {
#     x[1:max(which(x==1))] <- 1
#     x
#   }))
#   
#   stay   <- apply(sts, 1, function(x) {
#     as.data.frame(table(x)) %>% as_tibble() %>% setNames(c('site', 'days')) %>% mutate(site = as.numeric(as.character(site))) %>%
#       full_join(tibble(site = 1:nrow(mod@Sites$crds)), by = "site") %>% arrange(site) %>% mutate(days = ifelse(site==1 | site==max(site), NA, days)) %>% pull(days)
#   }) %>% apply(., 1, sum, na.rm = T)
#   
#   transM <- expand_grid(a = 1:nrow(mod@Sites$crds), b = 1:nrow(mod@Sites$crds)) %>% 
#     left_join(apply(sts, 1, function(x) tibble(a = x[!is.na(x)][-sum(!is.na(x))], b = x[!is.na(x)][-1]) %>% filter(a!=b) %>% mutate(n = 1)) %>%
#                 do.call("rbind",.) %>% group_by(a,b) %>% summarise_at("n", sum), by = join_by(a, b)) %>% filter(n>0) %>%
#     mutate(lon1 = mod@Sites$crd$Lon[a], lat1 = mod@Sites$crd$Lat[a], lon2 = mod@Sites$crd$Lon[b], lat2 = mod@Sites$crd$Lat[b])
#   
#   # ggplot() +
#   #   geom_curve(data = transM, mapping = aes(x = lon1, y = lat1, xend = lon2, yend = lat2, linewidth = n), 
#   #              alpha = 0.7, lineend='round') +
#   #   geom_point(data = crds %>% mutate(days = stay), mapping = aes(x = Lon, y = Lat, size = days), 
#   #              shape = 16, color = "orange") +
#   #   theme_light()
#   
#   invisible(list(mod@Sites$crds %>% mutate(days = stay), transM))
#   
# }
  

