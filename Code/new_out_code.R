########
######
library(MASS, lib="myR")
library(RColorBrewer, lib="myR")
library(reshape2, lib="myR")
library(ggstance, lib="myR")
library(ggplot2, lib="myR")
library(ggformula, lib="myR")
library(mosaicData, lib="myR")
library(mosaic, lib="myR")
library(latticeExtra, lib="myR")
library(gridExtra, lib="myR")
library(gridBase, lib="myR")
library(grid, lib="myR")
library(mgcv, lib="myR")
library(bbmle, lib="myR")
load('yes_changing_100_sens_act.RData')                 # ^ 0 Y

###########################
###########################
###### POP LEVEL OUTS ######
########
cumHumanInfections_yes_sens = matrix(0, Nreps, T + 1)
humanInfections_yes_sens = matrix(0, Nreps, T + 1)
for(ii in 1 : Nreps)
{
  cumHumanInfections_yes_sens[ii, ] = data_yes_changing_sens_act[1, ii][[1]]
  humanInfections_yes_sens[ii, ] = data_yes_changing_sens_act[2, ii][[1]]
}
humanInfections_yes_sens <- humanInfections_yes_sens/cumHumanInfections_yes_sens[1,1]
cumHumanInfections_yes_sens <- cumHumanInfections_yes_sens/cumHumanInfections_yes_sens[1,1]

e <- melt(as.data.frame(humanInfections_yes_sens))
e$time <- as.numeric(rep(1:201, each=Nreps))
e$rep <- as.numeric(rep(1:Nreps, 201))
e <- e[,c(4,3,2)]
colnames(e) <- c("rep", "time", "inf_prev")
no_inf <- e$rep[which(e$time==30 & e$inf_prev < 0.003)]
e_new<- e[-c(which(e$rep %in% no_inf)),]

e_cum <- melt(as.data.frame(cumHumanInfections_yes_sens))
e_cum$time <- as.numeric(rep(1:201, each=Nreps))
e_cum$rep <- as.numeric(rep(1:Nreps, 201))
e_cum <- e_cum[,c(4,3,2)]
colnames(e_cum) <- c("rep", "time", "inf_prev")
no_inf <- e_cum$rep[which(e_cum$time==30 & e_cum$inf_prev < 0.003)]
e_cum_new<- e_cum[-c(which(e_cum$rep %in% no_inf)),]

e_levels <- sort(unique(e_new$rep))
e_prev_max <- numeric(length=length(e_levels))
e_prev_max_date <- numeric(length=length(e_levels))
e_end_epi <- numeric(length=length(e_levels))
e_prev_start_5 <- numeric(length=length(e_levels))
e_prev_end_5 <- numeric(length=length(e_levels))
e_total_inf <- numeric(length=length(e_levels))
e_10_inf <- numeric(length=length(e_levels))
e_30_inf <- numeric(length=length(e_levels))
e_50_inf <- numeric(length=length(e_levels))
e_65_inf <- numeric(length=length(e_levels))
for(xx in e_levels){ 
  e_prev_max[which(e_levels==xx)]<- max(e_new$inf_prev[which(e_new$rep == xx)])
  temp_time <- which.max(e_new$inf_prev[which(e_new$rep == xx)])
  e_prev_max_date[which(e_levels==xx)]<- e_new$time[which(e_new$rep == xx)][temp_time]
  e_end_epi[which(e_levels==xx)] <- which.max(cumsum(e_new$inf_prev[which(e_new$rep == xx)])) + 1
  start <- which(e_new$inf_prev[which(e_new$rep == xx)] > 0.05)[1]
  e_prev_start_5[which(e_levels==xx)] <- start
  e_prev_end_5[which(e_levels==xx)] <- which(e_new$inf_prev[which(e_new$rep == xx)] < 0.05 & e_new$time[which(e_new$rep == xx)] > start)[1]
  e_total_inf[which(e_levels==xx)] <- e_cum_new$inf_prev[which(e_cum_new$rep == xx & e_cum_new$time == 201)]
  e_10_inf[which(e_levels==xx)] <- which(e_cum_new$inf_prev[which(e_cum_new$rep == xx)] >= 0.10)[2]
  e_30_inf[which(e_levels==xx)] <- which(e_cum_new$inf_prev[which(e_cum_new$rep == xx)] >= 0.30)[2]
  e_50_inf[which(e_levels==xx)] <- which(e_cum_new$inf_prev[which(e_cum_new$rep == xx)] >= 0.50)[2]
  e_65_inf[which(e_levels==xx)] <- which(e_cum_new$inf_prev[which(e_cum_new$rep == xx)] >= 0.65)[2]
}

prev_max <- rbind(favstats(e_prev_max))
prev_max_date <- rbind(favstats(e_prev_max_date))
end_epi <- rbind(favstats(e_end_epi))
prev_start_5 <- rbind(favstats(e_prev_start_5))
prev_end_5 <- rbind(favstats(e_prev_end_5))
total_inf <- rbind(favstats(e_total_inf))
all_10_inf <- rbind(favstats(e_10_inf))
all_30_inf <- rbind(favstats(e_30_inf))
all_50_inf <- rbind(favstats(e_50_inf))
all_65_inf <- rbind(favstats(e_65_inf))

#######
inputs_sens <- list(e_new=e_new,e_cum_new=e_cum_new)

outputs_sens <- list(prev_max=prev_max, prev_max_date=prev_max_date, end_epi=end_epi, prev_start_5=prev_start_5, prev_end_5=prev_end_5, 
                     total_inf=total_inf, all_10_inf=all_10_inf, all_30_inf=all_30_inf, all_50_inf=all_50_inf, all_65_inf=all_65_inf)

save(inputs_sens, outputs_sens, file = "metric_inputs_outputs_sens.RData")
save(outputs_sens, file = "metric_outputs_sens.RData")

save(e_new, e_cum_new, file="metric_inputs_sens.RData")
save(prev_max, prev_max_date, end_epi, prev_start_5, prev_end_5, total_inf, all_10_inf, all_30_inf, all_50_inf, all_65_inf, file="metric_outputs2_sens.RData")
###### END POP LEVEL OUTS ######
###########################
###########################
#### GET DATA SETUP 100 ####

reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

total_bites_per_rep_yes_changing_house_100 <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max+1, ncol=Nh_temp)
rel_bites_per_rep_yes_changing_house_100 <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max, ncol=Nh_temp)
perc_rel_bites_per_rep_yes_changing_house_100 <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  norm_bites <- colSums(data_yes_changing_sens_act[3,xx][[1]][[(rho_max+1)]])
  total_bites_per_rep_yes_changing_house_100[[x]][6,] <- norm_bites
  for(y in 1:(rho_max)){
    total_day_y <- colSums(data_yes_changing_sens_act[3,xx][[1]][[y]]) ## xx is replicate; y is the day of infectiousness
    total_bites_per_rep_yes_changing_house_100[[x]][y,] <- total_day_y
    rel_bites_per_rep_yes_changing_house_100[[x]][y,] <- total_day_y - norm_bites
    perc_rel_bites_per_rep_yes_changing_house_100[[x]][y,] <- ((total_day_y - norm_bites)/norm_bites)
  }
}

## to get distribution/frequency plot with  expect # mosq contacts, 
## for each day of infectiousness (row), get a list of values from all possible replicates (ncol=Nreps_yes_changing_house_100*Nh_temp)
total_bites_per_day_abs_yes_changing_house_100 <- matrix(NA, nrow=(rho_max + 1), ncol=(Nreps_yes_changing_house_100*Nh_temp))
rel_bites_per_day_yes_changing_house_100 <- matrix(NA, nrow=(rho_max), ncol=(Nreps_yes_changing_house_100*Nh_temp))
perc_rel_bites_per_day_yes_changing_house_100 <- matrix(NA, nrow=(rho_max), ncol=(Nreps_yes_changing_house_100*Nh_temp))
for(y in 1:(rho_max+1)){ # for each row of matrix
  for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
    index <- c((Nh_temp*(x-1)+1) : (Nh_temp*x))
    total_bites_per_day_abs_yes_changing_house_100[y,index] <- total_bites_per_rep_yes_changing_house_100[[x]][y,]
    if(y %in% c(1:rho_max)){
      rel_bites_per_day_yes_changing_house_100[y,index] <- rel_bites_per_rep_yes_changing_house_100[[x]][y,]
      perc_rel_bites_per_day_yes_changing_house_100[y,index] <- perc_rel_bites_per_rep_yes_changing_house_100[[x]][y,]
    }
  }
}

total_bites_per_day_abs2_yes_changing_house_100 <- as.data.frame(t(total_bites_per_day_abs_yes_changing_house_100))
total_bites_per_day_abs3_yes_changing_house_100 <- melt(total_bites_per_day_abs2_yes_changing_house_100)
total_bites_per_day_abs3_yes_changing_house_100$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
total_bites_per_day_abs3_yes_changing_house_100$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(total_bites_per_day_abs3_yes_changing_house_100) <- c("day_inf", "mosq_contacts","rep","person")
total_bites_per_day_abs3_yes_changing_house_100$day_inf <- as.character(total_bites_per_day_abs3_yes_changing_house_100$day_inf)
total_bites_per_day_abs3_yes_changing_house_100$day_inf <- unlist(strsplit(total_bites_per_day_abs3_yes_changing_house_100$day_inf,"V"))[(seq(0,length(total_bites_per_day_abs3_yes_changing_house_100$day_inf)*2,2))]
total_bites_per_day_abs3_yes_changing_house_100$day_inf <- factor(total_bites_per_day_abs3_yes_changing_house_100$day_inf, levels=c(1:6))
total_bites_per_day_abs4_yes_changing_house_100 <- droplevels(total_bites_per_day_abs3_yes_changing_house_100)
total_bites_per_day_abs4_yes_changing_house_100$mosq_contacts <- round(total_bites_per_day_abs4_yes_changing_house_100$mosq_contacts, digits=2)
total_bites_per_day_abs4_yes_changing_house_100$day_inf <- factor(total_bites_per_day_abs4_yes_changing_house_100$day_inf, levels=c(6,1:5))
total_bites_per_day_abs4_yes_changing_house_100 <- total_bites_per_day_abs4_yes_changing_house_100[,c(3,4,1,2)]


rel_bites_per_day_2_yes_changing_house_100 <- as.data.frame(t(rel_bites_per_day_yes_changing_house_100))
rel_bites_per_day_3_yes_changing_house_100 <- melt(rel_bites_per_day_2_yes_changing_house_100)
rel_bites_per_day_3_yes_changing_house_100$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
rel_bites_per_day_3_yes_changing_house_100$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(rel_bites_per_day_3_yes_changing_house_100) <- c("day_inf", "rel_mosq_contacts","rep","person")
rel_bites_per_day_3_yes_changing_house_100$day_inf <- as.character(rel_bites_per_day_3_yes_changing_house_100$day_inf)
rel_bites_per_day_3_yes_changing_house_100$day_inf <- unlist(strsplit(rel_bites_per_day_3_yes_changing_house_100$day_inf,"V"))[(seq(0,length(rel_bites_per_day_3_yes_changing_house_100$day_inf)*2,2))]
rel_bites_per_day_3_yes_changing_house_100$day_inf <- factor(rel_bites_per_day_3_yes_changing_house_100$day_inf, levels=c(1:5))
rel_bites_per_day_4_yes_changing_house_100 <- droplevels(rel_bites_per_day_3_yes_changing_house_100)
rel_bites_per_day_4_yes_changing_house_100$rel_mosq_contacts <- round(rel_bites_per_day_4_yes_changing_house_100$rel_mosq_contacts, digits=2)
rel_bites_per_day_4_yes_changing_house_100$day_inf <- factor(rel_bites_per_day_4_yes_changing_house_100$day_inf, levels=c(1:5))
rel_bites_per_day_4_yes_changing_house_100 <- rel_bites_per_day_4_yes_changing_house_100[,c(3,4,1,2)]

perc_rel_bites_per_day_2_yes_changing_house_100 <- as.data.frame(t(perc_rel_bites_per_day_yes_changing_house_100))
perc_rel_bites_per_day_3_yes_changing_house_100 <- melt(perc_rel_bites_per_day_2_yes_changing_house_100)
perc_rel_bites_per_day_3_yes_changing_house_100$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
perc_rel_bites_per_day_3_yes_changing_house_100$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(perc_rel_bites_per_day_3_yes_changing_house_100) <- c("day_inf", "perc_rel_mosq_contacts","rep","person")
perc_rel_bites_per_day_3_yes_changing_house_100$day_inf <- as.character(perc_rel_bites_per_day_3_yes_changing_house_100$day_inf)
perc_rel_bites_per_day_3_yes_changing_house_100$day_inf <- unlist(strsplit(perc_rel_bites_per_day_3_yes_changing_house_100$day_inf,"V"))[(seq(0,length(perc_rel_bites_per_day_3_yes_changing_house_100$day_inf)*2,2))]
perc_rel_bites_per_day_3_yes_changing_house_100$day_inf <- factor(perc_rel_bites_per_day_3_yes_changing_house_100$day_inf, levels=c(1:5))
perc_rel_bites_per_day_4_yes_changing_house_100 <- droplevels(perc_rel_bites_per_day_3_yes_changing_house_100)
perc_rel_bites_per_day_4_yes_changing_house_100$perc_rel_mosq_contacts <- round(perc_rel_bites_per_day_4_yes_changing_house_100$perc_rel_mosq_contacts, digits=2)
perc_rel_bites_per_day_4_yes_changing_house_100$day_inf <- factor(perc_rel_bites_per_day_4_yes_changing_house_100$day_inf, levels=c(1:5))
perc_rel_bites_per_day_4_yes_changing_house_100 <- perc_rel_bites_per_day_4_yes_changing_house_100[,c(3,4,1,2)]

###### ####### ######
no_inf_by_rep_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  no_inf <- which(is.na(total_bites_per_rep_yes_changing_house_100[[x]][1,]))
  extra <- Nh_temp - length(no_inf)
  no_inf_by_rep_100[x,] <- c( which(is.na(total_bites_per_rep_yes_changing_house_100[[x]][1,])) , rep(NA,extra))
}

total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes <- numeric(length=length(total_bites_per_day_abs4_yes_changing_house_100$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                                     (total_bites_per_day_abs4_yes_changing_house_100$person %in% inf_no_yes))] <- 0
  total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                                     !(total_bites_per_day_abs4_yes_changing_house_100$person %in% inf_no_yes))] <- 1
}
total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes <- as.factor(total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes)


rel_bites_per_day_4_yes_changing_house_100$inf_no_yes <- numeric(length=length(rel_bites_per_day_4_yes_changing_house_100$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  rel_bites_per_day_4_yes_changing_house_100$inf_no_yes[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                (rel_bites_per_day_4_yes_changing_house_100$person %in% inf_no_yes))] <- 0
  rel_bites_per_day_4_yes_changing_house_100$inf_no_yes[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                !(rel_bites_per_day_4_yes_changing_house_100$person %in% inf_no_yes))] <- 1
}
rel_bites_per_day_4_yes_changing_house_100$inf_no_yes <- as.factor(rel_bites_per_day_4_yes_changing_house_100$inf_no_yes)

perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes <- numeric(length=length(perc_rel_bites_per_day_4_yes_changing_house_100$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                     (perc_rel_bites_per_day_4_yes_changing_house_100$person %in% inf_no_yes))] <- 0
  perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                     !(perc_rel_bites_per_day_4_yes_changing_house_100$person %in% inf_no_yes))] <- 1
}
perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes <- as.factor(perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes)


total_bites_per_day_abs4_yes_changing_house_100b <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes==1),])
rel_bites_per_day_4_yes_changing_house_100b <- droplevels(rel_bites_per_day_4_yes_changing_house_100[which(rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1),])
perc_rel_bites_per_day_4_yes_changing_house_100b <- droplevels(perc_rel_bites_per_day_4_yes_changing_house_100[which(perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1),])



rel_bites_per_day_4_yes_changing_house_100$top10 <- numeric(length=length(rel_bites_per_day_4_yes_changing_house_100$rep))
total_bites_per_day_abs4_yes_changing_house_100$top10 <- numeric(length=length(total_bites_per_day_abs4_yes_changing_house_100$rep))
rel_bites_per_day_4_yes_changing_house_100$top10_1_3 <- numeric(length=length(rel_bites_per_day_4_yes_changing_house_100$rep))
total_bites_per_day_abs4_yes_changing_house_100$top10_1_3 <- numeric(length=length(total_bites_per_day_abs4_yes_changing_house_100$rep))
rel_bites_per_day_4_yes_changing_house_100$top10_4_6 <- numeric(length=length(rel_bites_per_day_4_yes_changing_house_100$rep))
total_bites_per_day_abs4_yes_changing_house_100$top10_4_6 <- numeric(length=length(total_bites_per_day_abs4_yes_changing_house_100$rep))
rel_bites_per_day_4_yes_changing_house_100$top2 <- numeric(length=length(rel_bites_per_day_4_yes_changing_house_100$rep))
total_bites_per_day_abs4_yes_changing_house_100$top2 <- numeric(length=length(total_bites_per_day_abs4_yes_changing_house_100$rep))

perc_rel_bites_per_day_4_yes_changing_house_100$top10 <- numeric(length=length(perc_rel_bites_per_day_4_yes_changing_house_100$rep))
perc_rel_bites_per_day_4_yes_changing_house_100$top10_1_3 <- numeric(length=length(perc_rel_bites_per_day_4_yes_changing_house_100$rep))
perc_rel_bites_per_day_4_yes_changing_house_100$top10_4_6 <- numeric(length=length(perc_rel_bites_per_day_4_yes_changing_house_100$rep))
perc_rel_bites_per_day_4_yes_changing_house_100$top2 <- numeric(length=length(perc_rel_bites_per_day_4_yes_changing_house_100$rep))

total_bites_100_norm <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$day_inf ==6),])
total_bites_100_1_3 <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$day_inf ==2),])
total_bites_100_4_6 <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$day_inf ==3),])
for(xx in 1:Nreps_yes_changing_house_100){
  norm_rep1 <- total_bites_100_norm[which(total_bites_100_norm$rep == xx),]
  norm_rep <- total_bites_100_norm$mosq_contacts[which(total_bites_100_norm$rep == xx)]
  top10length <- round(length(norm_rep)*0.10)
  index <- which(norm_rep %in% sort(norm_rep, decreasing=TRUE)[1:top10length])
  index2 <- norm_rep1$person[index]
  total_bites_per_day_abs4_yes_changing_house_100$top10[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                                total_bites_per_day_abs4_yes_changing_house_100$person %in% index2)] <- 1
  rel_bites_per_day_4_yes_changing_house_100$top10[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                           rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  perc_rel_bites_per_day_4_yes_changing_house_100$top10[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                perc_rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  
  top2length <- round(length(norm_rep)*0.025)
  index <- which(norm_rep %in% sort(norm_rep, decreasing=TRUE)[1:top2length])
  index2 <- norm_rep1$person[index]
  total_bites_per_day_abs4_yes_changing_house_100$top2[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                               total_bites_per_day_abs4_yes_changing_house_100$person %in% index2)] <- 1
  rel_bites_per_day_4_yes_changing_house_100$top2[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                          rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  perc_rel_bites_per_day_4_yes_changing_house_100$top2[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                               perc_rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  
  
  day_1_3_rep1 <- total_bites_100_1_3[which(total_bites_100_1_3$rep == xx),]
  day_1_3_rep <- total_bites_100_1_3$mosq_contacts[which(total_bites_100_1_3$rep == xx)]
  top10length <- round(length(day_1_3_rep)*0.10)
  index <- which(day_1_3_rep %in% sort(day_1_3_rep, decreasing=TRUE)[1:top10length])
  index2 <- day_1_3_rep1$person[index]
  total_bites_per_day_abs4_yes_changing_house_100$top10_1_3[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                                    total_bites_per_day_abs4_yes_changing_house_100$person %in% index2)] <- 1
  rel_bites_per_day_4_yes_changing_house_100$top10_1_3[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                               rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  perc_rel_bites_per_day_4_yes_changing_house_100$top10_1_3[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                    perc_rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  
  
  day_4_6_rep1 <- total_bites_100_4_6[which(total_bites_100_4_6$rep == xx),]
  day_4_6_rep <- total_bites_100_4_6$mosq_contacts[which(total_bites_100_4_6$rep == xx)]
  top10length <- round(length(day_4_6_rep)*0.10)
  index <- which(day_4_6_rep %in% sort(day_4_6_rep, decreasing=TRUE)[1:top10length])
  index2 <- day_4_6_rep1$person[index]
  total_bites_per_day_abs4_yes_changing_house_100$top10_4_6[which(total_bites_per_day_abs4_yes_changing_house_100$rep == xx & 
                                                                    total_bites_per_day_abs4_yes_changing_house_100$person %in% index2)] <- 1
  rel_bites_per_day_4_yes_changing_house_100$top10_4_6[which(rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                               rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
  perc_rel_bites_per_day_4_yes_changing_house_100$top10_4_6[which(perc_rel_bites_per_day_4_yes_changing_house_100$rep == xx & 
                                                                    perc_rel_bites_per_day_4_yes_changing_house_100$person %in% index2)] <- 1
}

###################################################################
##### HOUSE MOSQUITO COUNTS ####
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

house_mosq_per_rep_yes_changing_house_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
house_mosq_per_rep_when_1_3_yes_changing_house_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  house_mosq <- data_yes_changing_sens_act[9,xx][[1]][[(rho_max+1)]][4,]
  house_mosq_1_3 <- data_yes_changing_sens_act[9,xx][[1]][[2]][4,]
  house_mosq_per_rep_yes_changing_house_100[x,] <- house_mosq
  house_mosq_per_rep_when_1_3_yes_changing_house_100[x,] <- house_mosq_1_3
}

house_mosq_per_rep_yes_changing_house_100b <- as.data.frame(t(house_mosq_per_rep_yes_changing_house_100))
house_mosq_per_rep_when_1_3_yes_changing_house_100b <- as.data.frame(t(house_mosq_per_rep_when_1_3_yes_changing_house_100))
house_mosq_per_rep_yes_changing_house_100c <- melt(house_mosq_per_rep_yes_changing_house_100b)
house_mosq_per_rep_when_1_3_yes_changing_house_100c <- melt(house_mosq_per_rep_when_1_3_yes_changing_house_100b)
house_mosq_per_rep_yes_changing_house_100c$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
house_mosq_per_rep_when_1_3_yes_changing_house_100c$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
house_mosq_per_rep_yes_changing_house_100c$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
house_mosq_per_rep_when_1_3_yes_changing_house_100c$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(house_mosq_per_rep_yes_changing_house_100c) <- c("variable", "house_mosq","rep","person")
colnames(house_mosq_per_rep_when_1_3_yes_changing_house_100c) <- c("variable", "house_mosq","rep","person")
house_mosq_per_rep_yes_changing_house_100d <- house_mosq_per_rep_yes_changing_house_100c[,c(3,4,2)]
house_mosq_per_rep_when_1_3_yes_changing_house_100d <- house_mosq_per_rep_when_1_3_yes_changing_house_100c[,c(3,4,2)]


house_mosq_yes_changing_house_100 <- merge(house_mosq_per_rep_yes_changing_house_100d, house_mosq_per_rep_when_1_3_yes_changing_house_100d, by = c("rep","person"))
colnames(house_mosq_yes_changing_house_100) <- c("rep","person", "norm_mosq", "mosq_1_3")

total_bites_per_day_abs4_yes_changing_house_100b <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes==1),])
rel_bites_per_day_4_yes_changing_house_100b <- droplevels(rel_bites_per_day_4_yes_changing_house_100[which(rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1),])
perc_rel_bites_per_day_4_yes_changing_house_100b <- droplevels(perc_rel_bites_per_day_4_yes_changing_house_100[which(perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1),])

total_bites_per_day_abs4_yes_changing_house_100c <- droplevels(total_bites_per_day_abs4_yes_changing_house_100[which(total_bites_per_day_abs4_yes_changing_house_100$inf_no_yes==1 & total_bites_per_day_abs4_yes_changing_house_100$top10 == 1),])
rel_bites_per_day_4_yes_changing_house_100c <- droplevels(rel_bites_per_day_4_yes_changing_house_100[which(rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1 & rel_bites_per_day_4_yes_changing_house_100$top10==1),])
perc_rel_bites_per_day_4_yes_changing_house_100c <- droplevels(perc_rel_bites_per_day_4_yes_changing_house_100[which(perc_rel_bites_per_day_4_yes_changing_house_100$inf_no_yes==1 & perc_rel_bites_per_day_4_yes_changing_house_100$top10==1),])


total_mosq_bites_house_100 <- merge(total_bites_per_day_abs4_yes_changing_house_100b, house_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
total_mosq_bites_house_100 <- total_mosq_bites_house_100[order(total_mosq_bites_house_100[,1],total_mosq_bites_house_100[,2],total_mosq_bites_house_100[,3]),]

rel_mosq_bites_house_100 <- merge(rel_bites_per_day_4_yes_changing_house_100b, house_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
rel_mosq_bites_house_100 <- rel_mosq_bites_house_100[order(rel_mosq_bites_house_100[,1],rel_mosq_bites_house_100[,2],rel_mosq_bites_house_100[,3]),]

perc_rel_mosq_bites_house_100 <- merge(perc_rel_bites_per_day_4_yes_changing_house_100b, house_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
perc_rel_mosq_bites_house_100 <- perc_rel_mosq_bites_house_100[order(perc_rel_mosq_bites_house_100[,1],perc_rel_mosq_bites_house_100[,2],perc_rel_mosq_bites_house_100[,3]),]

total_mosq_bites_house_100$top10 <- factor(total_mosq_bites_house_100$top10, levels=c(0,1))
rel_mosq_bites_house_100$top10 <- factor(rel_mosq_bites_house_100$top10, levels=c(0,1))
total_mosq_bites_house_100$top10_1_3 <- factor(total_mosq_bites_house_100$top10_1_3, levels=c(0,1))
rel_mosq_bites_house_100$top10_1_3 <- factor(rel_mosq_bites_house_100$top10_1_3, levels=c(0,1))
total_mosq_bites_house_100$top10_4_6 <- factor(total_mosq_bites_house_100$top10_4_6, levels=c(0,1))
rel_mosq_bites_house_100$top10_4_6 <- factor(rel_mosq_bites_house_100$top10_4_6, levels=c(0,1))
total_mosq_bites_house_100$top2 <- factor(total_mosq_bites_house_100$top2, levels=c(0,1))
rel_mosq_bites_house_100$top2 <- factor(rel_mosq_bites_house_100$top2, levels=c(0,1))

perc_rel_mosq_bites_house_100$top10 <- factor(perc_rel_mosq_bites_house_100$top10, levels=c(0,1))
perc_rel_mosq_bites_house_100$top10_1_3 <- factor(perc_rel_mosq_bites_house_100$top10_1_3, levels=c(0,1))
perc_rel_mosq_bites_house_100$top10_4_6 <- factor(perc_rel_mosq_bites_house_100$top10_4_6, levels=c(0,1))
perc_rel_mosq_bites_house_100$top2 <- factor(perc_rel_mosq_bites_house_100$top2, levels=c(0,1))

rel_mosq_bites_house_100_1_3 <- droplevels(rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$day_inf == 2),])

rel_mosq_bites_house_100_1_3b <- rel_mosq_bites_house_100_1_3
rel_mosq_bites_house_100_1_3b$norm_mosq <- as.factor(as.numeric(rel_mosq_bites_house_100_1_3b$norm_mosq))

perc_rel_mosq_bites_house_100_1_3 <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf == 2),])
perc_rel_mosq_bites_house_100_1_3b <- perc_rel_mosq_bites_house_100_1_3
perc_rel_mosq_bites_house_100_1_3b$norm_mosq <- as.factor(as.numeric(perc_rel_mosq_bites_house_100_1_3b$norm_mosq))

######## MOSQ BASED ON HOME #########
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

house_yes_changing_house_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  house <- data_yes_changing_sens_act[10,xx][[1]]
  house_yes_changing_house_100[x,] <- house
}

house_yes_changing_house_100b <- as.data.frame(t(house_yes_changing_house_100))
house_yes_changing_house_100c <- melt(house_yes_changing_house_100b)
house_yes_changing_house_100c$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
house_yes_changing_house_100c$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(house_yes_changing_house_100c) <- c("variable", "house_index","rep","person")
house_yes_changing_house_100d <- house_yes_changing_house_100c[,c(3,4,2)]

house_index_mosq_yes_changing_house_100 <- merge(house_mosq_yes_changing_house_100, house_yes_changing_house_100d, by = c("rep","person"))


total_mosq_bites_house_100 <- merge(total_bites_per_day_abs4_yes_changing_house_100, house_index_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
total_mosq_bites_house_100 <- total_mosq_bites_house_100[order(total_mosq_bites_house_100[,1],total_mosq_bites_house_100[,2],total_mosq_bites_house_100[,3]),]

rel_mosq_bites_house_100 <- merge(rel_bites_per_day_4_yes_changing_house_100, house_index_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
rel_mosq_bites_house_100 <- rel_mosq_bites_house_100[order(rel_mosq_bites_house_100[,1],rel_mosq_bites_house_100[,2],rel_mosq_bites_house_100[,3]),]

perc_rel_mosq_bites_house_100 <- merge(perc_rel_bites_per_day_4_yes_changing_house_100, house_index_mosq_yes_changing_house_100, by=c("rep","person"), all.x = TRUE)
perc_rel_mosq_bites_house_100 <- perc_rel_mosq_bites_house_100[order(perc_rel_mosq_bites_house_100[,1],perc_rel_mosq_bites_house_100[,2],perc_rel_mosq_bites_house_100[,3]),]


total_mosq_bites_house_100$top10 <- factor(total_mosq_bites_house_100$top10, levels=c(0,1))
rel_mosq_bites_house_100$top10 <- factor(rel_mosq_bites_house_100$top10, levels=c(0,1))
total_mosq_bites_house_100$top10_1_3 <- factor(total_mosq_bites_house_100$top10_1_3, levels=c(0,1))
rel_mosq_bites_house_100$top10_1_3 <- factor(rel_mosq_bites_house_100$top10_1_3, levels=c(0,1))
total_mosq_bites_house_100$top10_4_6 <- factor(total_mosq_bites_house_100$top10_4_6, levels=c(0,1))
rel_mosq_bites_house_100$top10_4_6 <- factor(rel_mosq_bites_house_100$top10_4_6, levels=c(0,1))
total_mosq_bites_house_100$top2 <- factor(total_mosq_bites_house_100$top2, levels=c(0,1))
rel_mosq_bites_house_100$top2 <- factor(rel_mosq_bites_house_100$top2, levels=c(0,1))

total_mosq_bites_house_100_1_3 <- droplevels(total_mosq_bites_house_100[which(total_mosq_bites_house_100$day_inf == 2),])
rel_mosq_bites_house_100_1_3 <- droplevels(rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$day_inf == 2),])

perc_rel_mosq_bites_house_100$top10 <- factor(perc_rel_mosq_bites_house_100$top10, levels=c(0,1))
perc_rel_mosq_bites_house_100$top10_1_3 <- factor(perc_rel_mosq_bites_house_100$top10_1_3, levels=c(0,1))
perc_rel_mosq_bites_house_100$top10_4_6 <- factor(perc_rel_mosq_bites_house_100$top10_4_6, levels=c(0,1))
perc_rel_mosq_bites_house_100$top2 <- factor(perc_rel_mosq_bites_house_100$top2, levels=c(0,1))
perc_rel_mosq_bites_house_100_1_3 <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf == 2),])

total_mosq_bites_house_100_1_3b <- droplevels(total_mosq_bites_house_100_1_3[which(!(is.na(total_mosq_bites_house_100_1_3$mosq_contacts))),])


rel_mosq_bites_house_100_4_6 <- droplevels(rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$day_inf == 3),])
rel_mosq_bites_house_100_7_9 <- droplevels(rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$day_inf == 4),])
rel_mosq_bites_house_100_10_12 <- droplevels(rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$day_inf == 5),])
perc_rel_mosq_bites_house_100_4_6 <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf == 3),])
perc_rel_mosq_bites_house_100_7_9 <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf == 4),])
perc_rel_mosq_bites_house_100_10_12 <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf == 5),])


total_mosq_bites_house_100_1_3b <- total_mosq_bites_house_100_1_3[which(total_mosq_bites_house_100_1_3$top10==1 & total_mosq_bites_house_100_1_3$top10_1_3==1),]


#########
favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$top10==1 & total_mosq_bites_house_100$day_inf == 6)])
favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$day_inf == 6)])
favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$top10==1 & total_mosq_bites_house_100$day_inf == 2)])
favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$day_inf == 2)])

total_mosq_stats <- rbind(favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$top10==1 & total_mosq_bites_house_100$day_inf == 6)]),
                          favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$day_inf == 6)]),
                          favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$top10==1 & total_mosq_bites_house_100$day_inf == 2)]),
                          favstats(total_mosq_bites_house_100$mosq_contacts[which(total_mosq_bites_house_100$day_inf == 2)]))
row.names(total_mosq_stats) <- c("healthy_top", "healthy", "1_3_top", "1_3")

perc_rel_mosq_bites_house_100c <- droplevels(perc_rel_mosq_bites_house_100[which(perc_rel_mosq_bites_house_100$day_inf %in% c("2","3","4")),])
perc_rel_mosq_bites_house_100c$perc_rel_mosq_contacts <- perc_rel_mosq_bites_house_100c$perc_rel_mosq_contacts*100

percents <- perc_rel_mosq_bites_house_100c$perc_rel_mosq_contacts[which(perc_rel_mosq_bites_house_100c$top10==1 & perc_rel_mosq_bites_house_100c$day_inf == 2)]
percents2 <- perc_rel_mosq_bites_house_100c$perc_rel_mosq_contacts[which(perc_rel_mosq_bites_house_100c$day_inf == 2)]
length(which(percents > 0))/length(which(!(is.na(percents))))
length(which(percents == 0))/length(which(!(is.na(percents))))
length(which(percents < 0))/length(which(!(is.na(percents))))

length(which(percents2 > 0))/length(which(!(is.na(percents2))))
length(which(percents2 == 0))/length(which(!(is.na(percents2))))
length(which(percents2 < 0))/length(which(!(is.na(percents2))))
length(which(percents == -100))/length(which(!(is.na(percents))))
length(which(percents < (-90)))/length(which(!(is.na(percents))))


percent_stats <- rbind(c(">0_top", "==0_top", "<0_top",">0", "==0", "<0", "==-100", "==-90" ), c(length(which(percents > 0))/length(which(!(is.na(percents)))),
                                                                                                 length(which(percents == 0))/length(which(!(is.na(percents)))),
                                                                                                 length(which(percents < 0))/length(which(!(is.na(percents)))),
                                                                                                 length(which(percents2 > 0))/length(which(!(is.na(percents2)))),
                                                                                                 length(which(percents2 == 0))/length(which(!(is.na(percents2)))),
                                                                                                 length(which(percents2 < 0))/length(which(!(is.na(percents2)))),
                                                                                                 length(which(percents == -100))/length(which(!(is.na(percents)))),
                                                                                                 length(which(percents < (-90)))/length(which(!(is.na(percents))))))

write.csv(total_mosq_stats,  "total_mosq_stats.csv")
write.csv(percent_stats, "percent_stats.csv")
##### ACTIVITY SPACE MOSQUITO COUNTS ####
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

mosq_per_rep_yes_changing_house_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  for(j in 1:Nh_temp){
    house_list <- which(HM2[j,] > 0)
    house_ind <- numeric(length=length(house_list))
    for(i in 1:(length(house_ind))){
      house_ind[i] <- which(data_yes_changing_sens_act[10,xx][[1]] %in% house_list[i])[1]
    }
    mosq_counts <- data_yes_changing_sens_act[9,xx][[1]][[(rho_max+1)]][4,house_ind]
    mosq_per_rep_yes_changing_house_100[x,j] <- sum(mosq_counts)
  }
}


mosq_per_rep_yes_changing_house_100b <- as.data.frame(t(mosq_per_rep_yes_changing_house_100))
mosq_per_rep_yes_changing_house_100c <- melt(mosq_per_rep_yes_changing_house_100b)
mosq_per_rep_yes_changing_house_100c$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
mosq_per_rep_yes_changing_house_100c$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(mosq_per_rep_yes_changing_house_100c) <- c("variable", "total_mosq","rep","person")
mosq_per_rep_yes_changing_house_100d <- mosq_per_rep_yes_changing_house_100c[,c(3,4,2)]

###################################################################
##### BITE PERCENT AT HOME ####
bite_percent <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max+1, ncol=Nh_temp)
rel_bite_percent <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  house <- data_yes_changing_sens_act[10,xx][[1]]
  house_mat <- as.matrix(cbind(1:Nh,house))
  norm <- t(data_yes_changing_sens_act[3,xx][[1]][[(rho_max+1)]])
  house_mosq <- norm[house_mat]
  total_mosq <- rowSums(norm)
  fraction_home_bites <- house_mosq/total_mosq
  bite_percent[[x]][6,] <- fraction_home_bites
  for(y in 1:(rho_max)){
    house <- data_yes_changing_sens_act[10,xx][[1]]
    house_mat <- as.matrix(cbind(1:Nh,house))
    norm <- t(data_yes_changing_sens_act[3,xx][[1]][[y]])
    house_mosq <- norm[house_mat]
    total_mosq <- rowSums(norm)
    fraction_home_bites_day <- house_mosq/total_mosq
    bite_percent[[x]][y,] <- fraction_home_bites_day
    
    rel_bite_percent[[x]][y,] <- fraction_home_bites_day - fraction_home_bites
  }
}


bite_percent_per_day <- matrix(NA, nrow=(rho_max + 1), ncol=(Nreps_yes_changing_house_100*Nh_temp))
rel_bite_percent_per_day <- matrix(NA, nrow=(rho_max), ncol=(Nreps_yes_changing_house_100*Nh_temp))
for(y in 1:(rho_max+1)){ # for each row of matrix
  for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
    index <- c((Nh_temp*(x-1)+1) : (Nh_temp*x))
    bite_percent_per_day[y,index] <- bite_percent[[x]][y,]
    if(y %in% c(1:rho_max)){
      rel_bite_percent_per_day[y,index] <- rel_bite_percent[[x]][y,]
    }
  }
}

bite_percent_per_day2 <- as.data.frame(t(bite_percent_per_day))
bite_percent_per_day3 <- melt(bite_percent_per_day2)
bite_percent_per_day3$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
bite_percent_per_day3$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(bite_percent_per_day3) <- c("day_inf", "mosq_contacts","rep","person")
bite_percent_per_day3$day_inf <- as.character(bite_percent_per_day3$day_inf)
bite_percent_per_day3$day_inf <- unlist(strsplit(bite_percent_per_day3$day_inf,"V"))[(seq(0,length(bite_percent_per_day3$day_inf)*2,2))]
bite_percent_per_day3$day_inf <- factor(bite_percent_per_day3$day_inf, levels=c(1:6))
bite_percent_per_day4 <- droplevels(bite_percent_per_day3)
bite_percent_per_day4$mosq_contacts <- round(bite_percent_per_day4$mosq_contacts, digits=4)
bite_percent_per_day4$day_inf <- factor(bite_percent_per_day4$day_inf, levels=c(6,1:5))
bite_percent_per_day4 <- bite_percent_per_day4[,c(3,4,1,2)]


rel_bite_percent_per_day_2 <- as.data.frame(t(rel_bite_percent_per_day))
rel_bite_percent_per_day_3 <- melt(rel_bite_percent_per_day_2)
rel_bite_percent_per_day_3$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
rel_bite_percent_per_day_3$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(rel_bite_percent_per_day_3) <- c("day_inf", "rel_mosq_contacts","rep","person")
rel_bite_percent_per_day_3$day_inf <- as.character(rel_bite_percent_per_day_3$day_inf)
rel_bite_percent_per_day_3$day_inf <- unlist(strsplit(rel_bite_percent_per_day_3$day_inf,"V"))[(seq(0,length(rel_bite_percent_per_day_3$day_inf)*2,2))]
rel_bite_percent_per_day_3$day_inf <- factor(rel_bite_percent_per_day_3$day_inf, levels=c(1:5))
rel_bite_percent_per_day_4 <- droplevels(rel_bite_percent_per_day_3)
rel_bite_percent_per_day_4$rel_mosq_contacts <- round(rel_bite_percent_per_day_4$rel_mosq_contacts, digits=4)
rel_bite_percent_per_day_4$day_inf <- factor(rel_bite_percent_per_day_4$day_inf, levels=c(1:5))
rel_bite_percent_per_day_4 <- rel_bite_percent_per_day_4[,c(3,4,1,2)]
rel_bite_percent_per_day_4$percent_change_mosq_contacts <- rel_bite_percent_per_day_4$rel_mosq_contacts*100


no_inf_by_rep_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  no_inf <- which(is.na(bite_percent[[x]][1,]))
  extra <- Nh_temp - length(no_inf)
  no_inf_by_rep_100[x,] <- c(which(is.na(bite_percent[[x]][1,])) , rep(NA,extra))
}

bite_percent_per_day4$inf_no_yes <- numeric(length=length(bite_percent_per_day4$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  bite_percent_per_day4$inf_no_yes[which(bite_percent_per_day4$rep == xx & 
                                           (bite_percent_per_day4$person %in% inf_no_yes))] <- 0
  bite_percent_per_day4$inf_no_yes[which(bite_percent_per_day4$rep == xx & 
                                           !(bite_percent_per_day4$person %in% inf_no_yes))] <- 1
}
bite_percent_per_day4$inf_no_yes <- as.factor(bite_percent_per_day4$inf_no_yes)


rel_bite_percent_per_day_4$inf_no_yes <- numeric(length=length(rel_bite_percent_per_day_4$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  rel_bite_percent_per_day_4$inf_no_yes[which(rel_bite_percent_per_day_4$rep == xx & 
                                                (rel_bite_percent_per_day_4$person %in% inf_no_yes))] <- 0
  rel_bite_percent_per_day_4$inf_no_yes[which(rel_bite_percent_per_day_4$rep == xx & 
                                                !(rel_bite_percent_per_day_4$person %in% inf_no_yes))] <- 1
}
rel_bite_percent_per_day_4$inf_no_yes <- as.factor(rel_bite_percent_per_day_4$inf_no_yes)

bite_percent_per_day4b <- droplevels(bite_percent_per_day4[which(bite_percent_per_day4$inf_no_yes==1),])
rel_bite_percent_per_day_4b <- droplevels(rel_bite_percent_per_day_4[which(rel_bite_percent_per_day_4$inf_no_yes==1),])

bite_percent_per_day4c <- bite_percent_per_day4b[-c(which(is.na(bite_percent_per_day4b$mosq_contacts))),]
bite_percent_per_day4c$mosq_contacts <- bite_percent_per_day4c$mosq_contacts*100

##### BITE PERCENT GRAPHS ####
summary(bite_percent_per_day4c$mosq_contacts[which(bite_percent_per_day4c$day_inf == 6)])
temp <- bite_percent_per_day4c$mosq_contacts[which(bite_percent_per_day4c$day_inf == 6)]
(length(which(temp == 100)))/length(temp)


rel_bite_percent_per_day_4c <- rel_bite_percent_per_day_4b[-c(which(is.na(rel_bite_percent_per_day_4b$percent_change_mosq_contacts))),]
#######
top10_rep_person <- rel_bites_per_day_4_yes_changing_house_100[,c(1,2,6)]
top10_rep_person <- top10_rep_person[!duplicated(top10_rep_person[,c(1,2)]),]

bite_percent_per_day4c_top10 <- merge(bite_percent_per_day4c,top10_rep_person, by=c("rep","person"), all.x = TRUE)
rel_bite_percent_per_day_4c_top10 <- merge(rel_bite_percent_per_day_4c,top10_rep_person, by=c("rep","person"), all.x = TRUE)

summary(bite_percent_per_day4c_top10$mosq_contacts[which(bite_percent_per_day4c_top10$day_inf == 6 & bite_percent_per_day4c_top10$top10 == 1)])
temp <- bite_percent_per_day4c_top10$mosq_contacts[which(bite_percent_per_day4c_top10$day_inf == 6 & bite_percent_per_day4c_top10$top10 == 1)]
(length(which(temp == 100)))/length(temp) 


#######################
###########################
##### PRIMARY #######
#######

reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

R_norm_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_norm_home_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_norm_rest_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  R_norm_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[6,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_norm_home_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[13,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_norm_rest_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[12,xx][[1]]) ## xx is replicate; y is the day of infectiousness
}

R_norm_yes_changing_100_melt <- melt(R_norm_yes_changing_100, na.rm=TRUE)
colnames(R_norm_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_yes_changing_100_melt$rep <- factor(as.character(R_norm_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_home_yes_changing_100_melt <- melt(R_norm_home_yes_changing_100, na.rm=TRUE)
colnames(R_norm_home_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_home_yes_changing_100_melt$rep <- factor(as.character(R_norm_home_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_rest_yes_changing_100_melt <- melt(R_norm_rest_yes_changing_100, na.rm=TRUE)
colnames(R_norm_rest_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_rest_yes_changing_100_melt$rep <- factor(as.character(R_norm_rest_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_yes_changing_100_melt$type <- as.character("all")
R_norm_home_yes_changing_100_melt$type <- as.character("home")
R_norm_rest_yes_changing_100_melt$type <- as.character("rest")

R_norm_yes_changing_melt <- rbind(R_norm_yes_changing_100_melt, R_norm_home_yes_changing_100_melt, R_norm_rest_yes_changing_100_melt)
R_norm_yes_changing_melt$type <- factor(R_norm_yes_changing_melt$type, levels=c( "home", "rest", "all"))

#######
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

R_move_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_move_home_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_move_rest_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  R_move_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[7,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_move_home_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[14,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_move_rest_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[15,xx][[1]]) ## xx is replicate; y is the day of infectiousness
}

R_move_yes_changing_100_melt <- melt(R_move_yes_changing_100, na.rm=TRUE)
colnames(R_move_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_yes_changing_100_melt$rep <- factor(as.character(R_move_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_home_yes_changing_100_melt <- melt(R_move_home_yes_changing_100, na.rm=TRUE)
colnames(R_move_home_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_home_yes_changing_100_melt$rep <- factor(as.character(R_move_home_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_rest_yes_changing_100_melt <- melt(R_move_rest_yes_changing_100, na.rm=TRUE)
colnames(R_move_rest_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_rest_yes_changing_100_melt$rep <- factor(as.character(R_move_rest_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_yes_changing_100_melt$type <- as.character("all")
R_move_home_yes_changing_100_melt$type <- as.character("home")
R_move_rest_yes_changing_100_melt$type <- as.character("rest")

R_move_yes_changing_melt <- rbind(R_move_yes_changing_100_melt, R_move_home_yes_changing_100_melt, R_move_rest_yes_changing_100_melt)
R_move_yes_changing_melt$type <- factor(R_move_yes_changing_melt$type, levels=c( "home", "rest", "all"))

###
inf_top10_rep_person <- total_bites_per_day_abs4_yes_changing_house_100[,c(1,2,5,6)]
inf_top10_rep_person <- inf_top10_rep_person[!duplicated(inf_top10_rep_person[,c(1,2)]),]

R_norm_yes_changing_melt10 <- merge(R_norm_yes_changing_melt,inf_top10_rep_person, by=c("rep","person"), all.x = TRUE)
R_move_yes_changing_melt10 <- merge(R_move_yes_changing_melt,inf_top10_rep_person, by=c("rep","person"), all.x = TRUE)


R_norm_yes_changing_melt10$type2 <- as.character("norm")
R_move_yes_changing_melt10$type2 <- as.character("move")
colnames(R_norm_yes_changing_melt10) <- c("rep","person","R","type","inf_no_yes", "top10", "type2")
colnames(R_move_yes_changing_melt10) <- c("rep","person","R","type", "inf_no_yes", "top10", "type2")


R_move_yes_changing_melt2 <- rbind(R_norm_yes_changing_melt10, R_move_yes_changing_melt10)
R_move_yes_changing_melt2$type2 <- factor(R_move_yes_changing_melt2$type2, levels=c("norm","move"))



R_move_yes_changing_melt2c <- droplevels(R_move_yes_changing_melt2[-c(which(R_move_yes_changing_melt2$type == "all")),])

first_code_R_graph <- ggplot(data=R_move_yes_changing_melt2c, aes(x=R, ..count.., color=type, linetype=type2)) +
  geom_density(size=0.75) +
  theme_classic() +
  coord_cartesian(xlim=c(0,16)) +
  scale_x_continuous("Expected R value", breaks=seq(0,16,2)) +
  scale_y_continuous("Counts", breaks=seq(0,300000,25000)) +
  scale_color_discrete("Where Primary Infectious\nIndividual was Biten", labels=c("Home","Other Houses")) +
  scale_linetype_discrete("How R Value\nCalculated", labels=c("Pre-Epidemic","With Movement\nChange")) +
  #labs(title="Effect of Movement Change on Expected R Value") +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = c(0.80)), 
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14, lineheight=0.80),
        legend.title.align=0.5,
        legend.position = c(0.8,0.6),
        legend.key.size = unit(1.5,"line"),
        legend.text = element_text(size=12))


tiff("first_code_R_graph.tiff", units="in", width=6, height=6, res=300)
first_code_R_graph
dev.off()
######

R_move_yes_changing_melt3 <- droplevels(R_move_yes_changing_melt2[-c(which(R_move_yes_changing_melt2$inf_no_yes == 0)),])

norm_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="norm")]
move_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="move")]
norm_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "home" & R_move_yes_changing_melt3$type2=="norm")]
move_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "home" & R_move_yes_changing_melt3$type2=="move")]
norm_rest <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "rest" & R_move_yes_changing_melt3$type2=="norm")]
move_rest <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "rest" & R_move_yes_changing_melt3$type2=="move")]

primary_R_stats <- rbind(favstats(norm_all), favstats(move_all), 
                         favstats(norm_home), favstats(move_home), 
                         favstats(norm_rest), favstats(move_rest))
row.names(primary_R_stats) <- c("norm_all","move_all",
                                "norm_home", "move_home",
                                "norm_rest", "move_rest")

write.csv(primary_R_stats, "primary_R_stats.csv")
##########
R_move_yes_changing_melt3$top10 <- factor(R_move_yes_changing_melt3$top10, levels=c("0","1"))
R_move_yes_changing_melt3a <- R_move_yes_changing_melt3[order(R_move_yes_changing_melt3[,1], R_move_yes_changing_melt3[,2], R_move_yes_changing_melt3[,4]),]
R_move_yes_changing_melt4 <- R_move_yes_changing_melt3a[which(R_move_yes_changing_melt3a$type2=="norm"),]
R_move <- R_move_yes_changing_melt3a$R[which(R_move_yes_changing_melt3a$type2 == "move")]
R_move_yes_changing_melt4$R_move <- R_move
R_move_yes_changing_melt4$R_change <- R_move_yes_changing_melt4$R_move - R_move_yes_changing_melt4$R
R_move_yes_changing_melt4$R_change_percent <- ((R_move_yes_changing_melt4$R_move - R_move_yes_changing_melt4$R)/R_move_yes_changing_melt4$R)*100

#######

all <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all")]
home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home")]
rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest")]
low_all <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)]
high_all <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)]
low_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==0)]
high_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==1)]
low_rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==0)]
high_rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==1)]

####
all_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all")]))*100
home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="home")]))*100
rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="rest")]))*100
low_all_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)]))*100
high_all_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)]))*100
low_home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==0)]))*100
high_home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="home"& R_move_yes_changing_melt4$top10==1)]))*100
low_rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==0)]))*100
high_rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="rest"& R_move_yes_changing_melt4$top10==1)]))*100

primary_R_change_stats <- rbind(favstats(all), favstats(home), favstats(rest),
                                favstats(low_all), favstats(low_home), favstats(low_rest),
                                favstats(high_all), favstats(high_home), favstats(high_rest),
                                favstats(all_perc), favstats(home_perc), favstats(rest_perc),
                                favstats(low_all_perc), favstats(low_home_perc), favstats(low_rest_perc),
                                favstats(high_all_perc), favstats(high_home_perc), favstats(high_rest_perc))
row.names(primary_R_change_stats) <- c("all","home","rest",
                                       "low_all","low_home","low_rest",
                                       "high_all","high_home","high_rest",
                                       "all_perc","home_perc","rest_perc",
                                       "low_all_perc","low_home_perc","low_rest_perc",
                                       "high_all_perc","high_home_perc","high_rest_perc")

write.csv(primary_R_change_stats, "primary_R_change_stats.csv")
#################
###################################################################
##################################################################
##### SECOND CODE #####
###################################################################
###################################################################
##########
total_mosq_bites_house_100a <- total_mosq_bites_house_100

total_mosq_bites_house_100a$top20 <- numeric(length=length(total_mosq_bites_house_100a$rep))
total_mosq_bites_house_100a$top20_sick <- numeric(length=length(total_mosq_bites_house_100a$rep))
total_bites_100_norm <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 6),])
total_bites_100_norm_sick <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 6 & total_mosq_bites_house_100a$inf_no_yes == 1),])

for(xx in 1:Nreps_yes_changing_house_100){
  norm_rep1 <- total_bites_100_norm[which(total_bites_100_norm$rep == xx),]
  norm_rep <- total_bites_100_norm$mosq_contacts[which(total_bites_100_norm$rep == xx)]
  top20length <- round(length(norm_rep)*0.20)
  index <- which(norm_rep %in% sort(norm_rep, decreasing=TRUE)[1:top20length])
  index2 <- norm_rep1$person[index]
  total_mosq_bites_house_100a$top20[which(total_mosq_bites_house_100a$rep == xx & 
                                            total_mosq_bites_house_100a$person %in% index2)] <- 1
  
  norm_rep1_sick <- total_bites_100_norm_sick[which(total_bites_100_norm_sick$rep == xx),]
  norm_rep_sick <- total_bites_100_norm_sick$mosq_contacts[which(total_bites_100_norm_sick$rep == xx)]
  top20length_sick <- round(length(norm_rep_sick)*0.20)
  index_sick <- which(norm_rep_sick %in% sort(norm_rep_sick, decreasing=TRUE)[1:top20length_sick])
  index2_sick <- norm_rep1_sick$person[index_sick]
  total_mosq_bites_house_100a$top20_sick[which(total_mosq_bites_house_100a$rep == xx & 
                                                 total_mosq_bites_house_100a$person %in% index2_sick)] <- 1
}
total_mosq_bites_house_100a$top20 <- as.factor(as.character(total_mosq_bites_house_100a$top20))
total_mosq_bites_house_100a$top20_sick <- as.factor(as.character(total_mosq_bites_house_100a$top20_sick))

total_mosq_bites_house_100b <- droplevels(total_mosq_bites_house_100a[-c(which(total_mosq_bites_house_100a$day_inf == 6)),])

#####
a <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1),]
a2 <- droplevels(a[-c(which(is.na(a$mosq_contacts))),])
out_pre<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==1)])$out
out_1_3<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==2)])$out
out_4_6<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==3)])$out
out_7_9<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==4)])$out
out_10_12<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==5)])$out

total_short <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1),]
outs <- c(which(total_short$day_inf == 1 & total_short$mosq_contacts %in% out_pre),
          which(total_short$day_inf == 2 & total_short$mosq_contacts %in% out_1_3),
          which(total_short$day_inf == 3 & total_short$mosq_contacts %in% out_4_6),
          which(total_short$day_inf == 4 & total_short$mosq_contacts %in% out_7_9),
          which(total_short$day_inf == 5 & total_short$mosq_contacts %in% out_10_12))
total_short2 <- droplevels(total_short[-outs,])

bites_graph <- ggplot(data=total_short2, aes(day_inf, mosq_contacts)) +
  geom_violin() +
  scale_x_discrete("Day of Symptoms",  labels=c("Presymptomatic", "Days 1-3","Days 4-6","Days 7-9","Days 10-12"))  +
  scale_y_continuous("Number of Expected Mosquito Contacts", breaks=seq(0,10,2))

tiff("bites_graph.tiff", units="in", width=8, height=6, res=300)
bites_graph
dev.off()

perc_rel_mosq_bites_house_100a <- perc_rel_mosq_bites_house_100
rel_mosq_bites_house_100$perc_rel_mosq_contacts <- (perc_rel_mosq_bites_house_100a$perc_rel_mosq_contacts)*100

a <- rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$inf_no_yes == 1),]
a2 <- droplevels(a[-c(which(is.na(a$rel_mosq_contacts))),])
out_pre<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==1)])$out
out_1_3<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==2)])$out
out_4_6<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==3)])$out
out_7_9<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==4)])$out
out_10_12<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==5)])$out

rel_short <- rel_mosq_bites_house_100[which(rel_mosq_bites_house_100$inf_no_yes == 1),]
outs <- c(which(rel_short$day_inf == 1 & rel_short$rel_mosq_contacts %in% out_pre),
          which(rel_short$day_inf == 2 & rel_short$rel_mosq_contacts %in% out_1_3),
          which(rel_short$day_inf == 3 & rel_short$rel_mosq_contacts %in% out_4_6),
          which(rel_short$day_inf == 4 & rel_short$rel_mosq_contacts %in% out_7_9),
          which(rel_short$day_inf == 5 & rel_short$rel_mosq_contacts %in% out_10_12))
rel_short2 <- droplevels(rel_short[-outs,])
##

#### REMOVE OUTLIERS FOR VIOLIN PLOT #####
a <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1 & total_mosq_bites_house_100b$top20 == 0),]
a2 <- droplevels(a[-c(which(is.na(a$mosq_contacts))),])
out_pre<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==1)])$out
out_1_3<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==2)])$out
out_4_6<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==3)])$out
out_7_9<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==4)])$out
out_10_12<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==5)])$out

mosq_80_with_outs <- favstats(a2$mosq_contacts ~ a2$day_inf + a2$top20)

a <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1 & total_mosq_bites_house_100b$top20 == 1),]
a2 <- droplevels(a[-c(which(is.na(a$mosq_contacts))),])
out20_pre<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==1)])$out
out20_1_3<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==2)])$out
out20_4_6<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==3)])$out
out20_7_9<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==4)])$out
out20_10_12<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==5)])$out

mosq_20_with_outs <- favstats(a2$mosq_contacts ~ a2$day_inf + a2$top20)

total_short <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1),]
outs <- c(which(total_short$day_inf == 1 & total_short$top20 == 0 & total_short$mosq_contacts %in% out_pre),
          which(total_short$day_inf == 2 & total_short$top20 == 0 & total_short$mosq_contacts %in% out_1_3),
          which(total_short$day_inf == 3 & total_short$top20 == 0 & total_short$mosq_contacts %in% out_4_6),
          which(total_short$day_inf == 4 & total_short$top20 == 0 & total_short$mosq_contacts %in% out_7_9),
          which(total_short$day_inf == 5 & total_short$top20 == 0 & total_short$mosq_contacts %in% out_10_12),
          which(total_short$day_inf == 1 & total_short$top20 == 1 & total_short$mosq_contacts %in% out20_pre),
          which(total_short$day_inf == 2 & total_short$top20 == 1 & total_short$mosq_contacts %in% out20_1_3),
          which(total_short$day_inf == 3 & total_short$top20 == 1 & total_short$mosq_contacts %in% out20_4_6),
          which(total_short$day_inf == 4 & total_short$top20 == 1 & total_short$mosq_contacts %in% out20_7_9),
          which(total_short$day_inf == 5 & total_short$top20 == 1 & total_short$mosq_contacts %in% out20_10_12))
total_short2 <- droplevels(total_short[-outs,])

mosq_20_80_no_outs <- favstats(total_short2$mosq_contacts ~ total_short2$day_inf + total_short2$top20)

bites_20_graph <- ggplot(data=total_short2, aes(day_inf, mosq_contacts)) +
  geom_violin() +
  scale_x_discrete("Day of Symptoms",  labels=c("Presymptomatic", "Days 1-3","Days 4-6","Days 7-9","Days 10-12"))  +
  scale_y_continuous("Number of Expected Mosquito Contacts", breaks=seq(0,20,5)) +
  facet_grid(. ~ top20, labeller = labeller(top20 = c('0' = "Bottom 80% bites pre-exposure",
                                                      '1' = "Top 20% bites pre-exposure"))) +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = 0.85), 
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12))

tiff("bites_20_graph.tiff", units="in", width=12, height=6, res=300)
bites_20_graph
dev.off()

write.csv(mosq_80_with_outs, "mosq_80_with_outs.csv")
write.csv(mosq_20_with_outs, "mosq_20_with_outs.csv")
write.csv(mosq_20_80_no_outs, "mosq_20_80_no_outs.csv")


#####
rel_mosq_bites_house_100a <- rel_mosq_bites_house_100
rel_mosq_bites_house_100a$top20 <- total_mosq_bites_house_100b$top20

a <- rel_mosq_bites_house_100a[which(rel_mosq_bites_house_100a$inf_no_yes == 1 & rel_mosq_bites_house_100a$top20 == 0),]
a2 <- droplevels(a[-c(which(is.na(a$rel_mosq_contacts))),])
out_pre<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==1)])$out
out_1_3<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==2)])$out
out_4_6<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==3)])$out
out_7_9<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==4)])$out
out_10_12<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==5)])$out

a <- rel_mosq_bites_house_100a[which(rel_mosq_bites_house_100a$inf_no_yes == 1 & rel_mosq_bites_house_100a$top20 == 1),]
a2 <- droplevels(a[-c(which(is.na(a$rel_mosq_contacts))),])
out20_pre<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==1)])$out
out20_1_3<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==2)])$out
out20_4_6<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==3)])$out
out20_7_9<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==4)])$out
out20_10_12<- boxplot.stats(a2$rel_mosq_contacts[which(a2$day_inf==5)])$out

rel_short <- rel_mosq_bites_house_100a[which(rel_mosq_bites_house_100a$inf_no_yes == 1),]
outs <- c(which(rel_short$day_inf == 1 & rel_short$top20 == 0 & rel_short$rel_mosq_contacts %in% out_pre),
          which(rel_short$day_inf == 2 & rel_short$top20 == 0 & rel_short$rel_mosq_contacts %in% out_1_3),
          which(rel_short$day_inf == 3 & rel_short$top20 == 0 & rel_short$rel_mosq_contacts %in% out_4_6),
          which(rel_short$day_inf == 4 & rel_short$top20 == 0 & rel_short$rel_mosq_contacts %in% out_7_9),
          which(rel_short$day_inf == 5 & rel_short$top20 == 0 & rel_short$rel_mosq_contacts %in% out_10_12),
          which(rel_short$day_inf == 1 & rel_short$top20 == 1 & rel_short$rel_mosq_contacts %in% out20_pre),
          which(rel_short$day_inf == 2 & rel_short$top20 == 1 & rel_short$rel_mosq_contacts %in% out20_1_3),
          which(rel_short$day_inf == 3 & rel_short$top20 == 1 & rel_short$rel_mosq_contacts %in% out20_4_6),
          which(rel_short$day_inf == 4 & rel_short$top20 == 1 & rel_short$rel_mosq_contacts %in% out20_7_9),
          which(rel_short$day_inf == 5 & rel_short$top20 == 1 & rel_short$rel_mosq_contacts %in% out20_10_12))
rel_short2 <- droplevels(rel_short[-outs,])

rel_mosq_20_80 <- favstats(rel_short2$rel_mosq_contacts ~ rel_short2$day_inf + rel_short2$top20)
perc_rel_mosq_20_80 <- favstats(rel_short2$perc_rel_mosq_contacts ~ rel_short2$day_inf + rel_short2$top20)

rel_bite_graph <- ggplot(data=rel_short2, aes(day_inf, rel_mosq_contacts)) +
  geom_violin(scale="width") +
  scale_x_discrete("Day of Symptoms",  labels=c("Presymptomatic", "Days 1-3","Days 4-6","Days 7-9","Days 10-12"))  +
  scale_y_continuous("Change in Number of Expected Mosquito Contacts", breaks=seq(-12,12,2)) +
  facet_grid(. ~ top20, labeller = labeller(top20 = c('0' = "Bottom 80% bites pre-exposure",
                                                      '1' = "Top 20% bites pre-exposure"))) +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = 0.85), 
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12))

tiff("rel_bite_graph.tiff", units="in", width=10, height=6, res=300)
rel_bite_graph
dev.off()

write.csv(rel_mosq_20_80, "rel_mosq_20_80.csv")
write.csv(perc_rel_mosq_20_80, "perc_rel_mosq_20_80.csv")



perc_rel_bite_graph <- ggplot(data=rel_short2, aes(day_inf, perc_rel_mosq_contacts)) +
  geom_violin(scale="width") +
  scale_x_discrete("Day of Symptoms",  labels=c("Presymptomatic", "Days 1-3","Days 4-6","Days 7-9","Days 10-12"))  +
  scale_y_continuous("Percent Change in Expected Mosquito Contacts", breaks=seq(-100,300,50)) +
  facet_grid(. ~ top20, labeller = labeller(top20 = c('0' = "Bottom 80% bites pre-exposure",
                                                      '1' = "Top 20% bites pre-exposure"))) +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = 0.85), 
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12))

tiff("perc_rel_bite_graph.tiff", units="in", width=12, height=6, res=300)
perc_rel_bite_graph
dev.off()


#
#### REMOVE OUTLIERS FOR VIOLIN PLOT top20_sick #####
a <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1 & total_mosq_bites_house_100b$top20_sick == 0),]
a2 <- droplevels(a[-c(which(is.na(a$mosq_contacts))),])
out_pre<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==1)])$out
out_1_3<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==2)])$out
out_4_6<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==3)])$out
out_7_9<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==4)])$out
out_10_12<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==5)])$out

a <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1 & total_mosq_bites_house_100b$top20_sick == 1),]
a2 <- droplevels(a[-c(which(is.na(a$mosq_contacts))),])
out20_pre<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==1)])$out
out20_1_3<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==2)])$out
out20_4_6<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==3)])$out
out20_7_9<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==4)])$out
out20_10_12<- boxplot.stats(a2$mosq_contacts[which(a2$day_inf==5)])$out

total_short <- total_mosq_bites_house_100b[which(total_mosq_bites_house_100b$inf_no_yes == 1),]
outs <- c(which(total_short$day_inf == 1 & total_short$top20_sick == 0 & total_short$mosq_contacts %in% out_pre),
          which(total_short$day_inf == 2 & total_short$top20_sick == 0 & total_short$mosq_contacts %in% out_1_3),
          which(total_short$day_inf == 3 & total_short$top20_sick == 0 & total_short$mosq_contacts %in% out_4_6),
          which(total_short$day_inf == 4 & total_short$top20_sick == 0 & total_short$mosq_contacts %in% out_7_9),
          which(total_short$day_inf == 5 & total_short$top20_sick == 0 & total_short$mosq_contacts %in% out_10_12),
          which(total_short$day_inf == 1 & total_short$top20_sick == 1 & total_short$mosq_contacts %in% out20_pre),
          which(total_short$day_inf == 2 & total_short$top20_sick == 1 & total_short$mosq_contacts %in% out20_1_3),
          which(total_short$day_inf == 3 & total_short$top20_sick == 1 & total_short$mosq_contacts %in% out20_4_6),
          which(total_short$day_inf == 4 & total_short$top20_sick == 1 & total_short$mosq_contacts %in% out20_7_9),
          which(total_short$day_inf == 5 & total_short$top20_sick == 1 & total_short$mosq_contacts %in% out20_10_12))
total_short2 <- droplevels(total_short[-outs,])


##### BITE PERCENT AT HOME ####
home_bites_count <- lapply(1:Nreps_yes_changing_house_100, matrix, data=NA, nrow=rho_max+1, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  house <- data_yes_changing_sens_act[10,xx][[1]]
  house_mat <- as.matrix(cbind(1:Nh,house))
  norm <- t(data_yes_changing_sens_act[3,xx][[1]][[(rho_max+1)]])
  house_mosq <- norm[house_mat]
  home_bites <- house_mosq
  home_bites_count[[x]][6,] <- home_bites
  for(y in 1:(rho_max)){
    house <- data_yes_changing_sens_act[10,xx][[1]]
    house_mat <- as.matrix(cbind(1:Nh,house))
    norm <- t(data_yes_changing_sens_act[3,xx][[1]][[y]])
    house_mosq <- norm[house_mat]
    home_bites_day <- house_mosq
    home_bites_count[[x]][y,] <- home_bites_day
  }
}
home_bites_count_per_day <- matrix(NA, nrow=(rho_max + 1), ncol=(Nreps_yes_changing_house_100*Nh_temp))
for(y in 1:(rho_max+1)){ # for each row of matrix
  for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
    index <- c((Nh_temp*(x-1)+1) : (Nh_temp*x))
    home_bites_count_per_day[y,index] <- home_bites_count[[x]][y,]
  }
}

home_bites_count_per_day2 <- as.data.frame(t(home_bites_count_per_day))
home_bites_count_per_day3 <- melt(home_bites_count_per_day2)
home_bites_count_per_day3$rep <- as.numeric(rep(1:Nreps_yes_changing_house_100, each=Nh_temp))
home_bites_count_per_day3$person <- as.numeric(rep(1:Nh_temp, Nreps_yes_changing_house_100))
colnames(home_bites_count_per_day3) <- c("day_inf", "mosq_contacts","rep","person")
home_bites_count_per_day3$day_inf <- as.character(home_bites_count_per_day3$day_inf)
home_bites_count_per_day3$day_inf <- unlist(strsplit(home_bites_count_per_day3$day_inf,"V"))[(seq(0,length(home_bites_count_per_day3$day_inf)*2,2))]
home_bites_count_per_day3$day_inf <- factor(home_bites_count_per_day3$day_inf, levels=c(1:6))
home_bites_count_per_day4 <- droplevels(home_bites_count_per_day3)
home_bites_count_per_day4$mosq_contacts <- round(home_bites_count_per_day4$mosq_contacts, digits=4)
home_bites_count_per_day4$day_inf <- factor(home_bites_count_per_day4$day_inf, levels=c(6,1:5))
home_bites_count_per_day4 <- home_bites_count_per_day4[,c(3,4,1,2)]


no_inf_by_rep_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  no_inf <- which(is.na(home_bites_count[[x]][1,]))
  extra <- Nh_temp - length(no_inf)
  no_inf_by_rep_100[x,] <- c(which(is.na(home_bites_count[[x]][1,])) , rep(NA,extra))
}

home_bites_count_per_day4$inf_no_yes <- numeric(length=length(home_bites_count_per_day4$rep))
for(xx in 1:Nreps_yes_changing_house_100){
  inf_no_yes <- no_inf_by_rep_100[xx,which(!(is.na(no_inf_by_rep_100[xx,])))]
  home_bites_count_per_day4$inf_no_yes[which(home_bites_count_per_day4$rep == xx & 
                                               (home_bites_count_per_day4$person %in% inf_no_yes))] <- 0
  home_bites_count_per_day4$inf_no_yes[which(home_bites_count_per_day4$rep == xx & 
                                               !(home_bites_count_per_day4$person %in% inf_no_yes))] <- 1
}
home_bites_count_per_day4$inf_no_yes <- as.factor(home_bites_count_per_day4$inf_no_yes)


#####
home_bites_count_per_day4 <- droplevels(home_bites_count_per_day4[order(home_bites_count_per_day4$rep, 
                                                                        home_bites_count_per_day4$person, 
                                                                        home_bites_count_per_day4$day_inf),])

home_bites_count_per_day4$percentile <- numeric(length=length(home_bites_count_per_day4$rep))
home_bites_count_per_day4$percentile_sick_presymp <- character(length=length(home_bites_count_per_day4$rep))
home_bites_count_per_day4$percentile_sick_days <- character(length=length(home_bites_count_per_day4$rep))
home_bites_count_per_day4$mosq_presymp <- numeric(length=length(home_bites_count_per_day4$rep))

#### percentile_sick done based on mosquito count at time point 1 (presymp <- -> similar to when infectious bite)
home_bites_count_per_day4_norm <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 6),])
home_bites_count_per_day4_norm_sick_presymp <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 1 & 
                                                                                            !(is.na(home_bites_count_per_day4$mosq_contacts))),])
home_bites_count_per_day4_norm_sick_1_3 <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 2 & 
                                                                                        !(is.na(home_bites_count_per_day4$mosq_contacts))),])
home_bites_count_per_day4_norm_sick_4_6 <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 3 &
                                                                                        !(is.na(home_bites_count_per_day4$mosq_contacts))),])
home_bites_count_per_day4_norm_sick_7_9 <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 4 &
                                                                                        !(is.na(home_bites_count_per_day4$mosq_contacts))),])
home_bites_count_per_day4_norm_sick_10_12 <- droplevels(home_bites_count_per_day4[which(home_bites_count_per_day4$day_inf == 5 & 
                                                                                          !(is.na(home_bites_count_per_day4$mosq_contacts))),])

for(xx in 1:Nreps_yes_changing_house_100){
  norm_rep <- home_bites_count_per_day4$mosq_contacts[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$day_inf == 6)]
  home_bites_count_per_day4$percentile[which(home_bites_count_per_day4$rep == xx)] <- rep(ntile(norm_rep, n=100), each=6)
  
  norm_rep_sick_presymp <- home_bites_count_per_day4_norm_sick_presymp$mosq_contacts[which(home_bites_count_per_day4_norm_sick_presymp$rep == xx)]
  home_bites_count_per_day4$percentile_sick_presymp[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$inf_no_yes==1)] <- 
    rep(ntile(norm_rep_sick_presymp, n=100), each=6)
  
  norm_rep_sick_presymp <- home_bites_count_per_day4_norm_sick_presymp$mosq_contacts[which(home_bites_count_per_day4_norm_sick_presymp$rep == xx)]
  home_bites_count_per_day4$mosq_presymp[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$inf_no_yes==1)] <- 
    rep(norm_rep_sick_presymp, each=6)
  
  
  norm_rep_sick_1_3 <- home_bites_count_per_day4_norm_sick_1_3$mosq_contacts[which(home_bites_count_per_day4_norm_sick_1_3$rep == xx)]
  home_bites_count_per_day4$percentile_sick_days[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$day_inf == 2 & 
                                                         !(is.na(home_bites_count_per_day4$mosq_contacts)))]<- ntile(norm_rep_sick_1_3, n=100)
  
  norm_rep_sick_4_6 <- home_bites_count_per_day4_norm_sick_4_6$mosq_contacts[which(home_bites_count_per_day4_norm_sick_4_6$rep == xx)]
  home_bites_count_per_day4$percentile_sick_days[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$day_inf == 3 & 
                                                         !(is.na(home_bites_count_per_day4$mosq_contacts)))]<- ntile(norm_rep_sick_4_6, n=100)
  
  norm_rep_sick_7_9 <- home_bites_count_per_day4_norm_sick_7_9$mosq_contacts[which(home_bites_count_per_day4_norm_sick_7_9$rep == xx)]
  home_bites_count_per_day4$percentile_sick_days[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$day_inf == 4 & 
                                                         !(is.na(home_bites_count_per_day4$mosq_contacts)))]<- ntile(norm_rep_sick_7_9, n=100)
  
  norm_rep_sick_10_12 <- home_bites_count_per_day4_norm_sick_10_12$mosq_contacts[which(home_bites_count_per_day4_norm_sick_10_12$rep == xx)]
  home_bites_count_per_day4$percentile_sick_days[which(home_bites_count_per_day4$rep == xx & home_bites_count_per_day4$day_inf == 5 & 
                                                         !(is.na(home_bites_count_per_day4$mosq_contacts)))]<- ntile(norm_rep_sick_10_12, n=100)
}

home_bites_count_per_day4$percentile_sick_presymp[which(home_bites_count_per_day4$percentile_sick_presymp == "")] <- NA
home_bites_count_per_day4$percentile_sick_days[which(home_bites_count_per_day4$percentile_sick_days == "")] <- NA
home_bites_count_per_day4$percentile_sick_presymp <- as.numeric(as.character(home_bites_count_per_day4$percentile_sick_presymp))
home_bites_count_per_day4$percentile_sick_days <- as.numeric(as.character(home_bites_count_per_day4$percentile_sick_days))

home_bites_count_per_day4b <- droplevels(home_bites_count_per_day4[c(which(home_bites_count_per_day4$day_inf == 6)),])

##########

total_mosq_bites_house_100a$percentile <- numeric(length=length(total_mosq_bites_house_100a$rep))
total_mosq_bites_house_100a$percentile_sick_presymp <- character(length=length(total_mosq_bites_house_100a$rep))
total_mosq_bites_house_100a$percentile_sick_days <- character(length=length(total_mosq_bites_house_100a$rep))
total_mosq_bites_house_100a$mosq_presymp <- numeric(length=length(total_mosq_bites_house_100a$rep))

#### percentile_sick done based on mosquito count at time point 1 (presymp <- -> similar to when infectious bite)
total_bites_100_norm <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 6),])
total_bites_100_norm_sick_presymp <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 1 & !(is.na(total_mosq_bites_house_100a$mosq_contacts))),])
total_bites_100_norm_sick_1_3 <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 2 & !(is.na(total_mosq_bites_house_100a$mosq_contacts))),])
total_bites_100_norm_sick_4_6 <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 3 & !(is.na(total_mosq_bites_house_100a$mosq_contacts))),])
total_bites_100_norm_sick_7_9 <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 4 & !(is.na(total_mosq_bites_house_100a$mosq_contacts))),])
total_bites_100_norm_sick_10_12 <- droplevels(total_mosq_bites_house_100a[which(total_mosq_bites_house_100a$day_inf == 5 & !(is.na(total_mosq_bites_house_100a$mosq_contacts))),])

for(xx in 1:Nreps_yes_changing_house_100){
  norm_rep <- total_mosq_bites_house_100a$mosq_contacts[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$day_inf == 6)]
  total_mosq_bites_house_100a$percentile[which(total_mosq_bites_house_100a$rep == xx)] <- rep(ntile(norm_rep, n=100), each=6)
  
  norm_rep_sick_presymp <- total_bites_100_norm_sick_presymp$mosq_contacts[which(total_bites_100_norm_sick_presymp$rep == xx)]
  total_mosq_bites_house_100a$percentile_sick_presymp[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$inf_no_yes==1)] <- rep(ntile(norm_rep_sick_presymp, n=100), each=6)
  
  norm_rep_sick_presymp <- total_bites_100_norm_sick_presymp$mosq_contacts[which(total_bites_100_norm_sick_presymp$rep == xx)]
  total_mosq_bites_house_100a$mosq_presymp[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$inf_no_yes==1)] <- rep(norm_rep_sick_presymp, each=6)
  
  
  norm_rep_sick_1_3 <- total_bites_100_norm_sick_1_3$mosq_contacts[which(total_bites_100_norm_sick_1_3$rep == xx)]
  total_mosq_bites_house_100a$percentile_sick_days[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$day_inf == 2 & !(is.na(total_mosq_bites_house_100a$mosq_contacts)))]<- ntile(norm_rep_sick_1_3, n=100)
  
  norm_rep_sick_4_6 <- total_bites_100_norm_sick_4_6$mosq_contacts[which(total_bites_100_norm_sick_4_6$rep == xx)]
  total_mosq_bites_house_100a$percentile_sick_days[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$day_inf == 3 & !(is.na(total_mosq_bites_house_100a$mosq_contacts)))]<- ntile(norm_rep_sick_4_6, n=100)
  
  norm_rep_sick_7_9 <- total_bites_100_norm_sick_7_9$mosq_contacts[which(total_bites_100_norm_sick_7_9$rep == xx)]
  total_mosq_bites_house_100a$percentile_sick_days[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$day_inf == 4 & !(is.na(total_mosq_bites_house_100a$mosq_contacts)))]<- ntile(norm_rep_sick_7_9, n=100)
  
  norm_rep_sick_10_12 <- total_bites_100_norm_sick_10_12$mosq_contacts[which(total_bites_100_norm_sick_10_12$rep == xx)]
  total_mosq_bites_house_100a$percentile_sick_days[which(total_mosq_bites_house_100a$rep == xx & total_mosq_bites_house_100a$day_inf == 5 & !(is.na(total_mosq_bites_house_100a$mosq_contacts)))]<- ntile(norm_rep_sick_10_12, n=100)
}

total_mosq_bites_house_100a$percentile_sick_presymp[which(total_mosq_bites_house_100a$percentile_sick_presymp == "")] <- NA
total_mosq_bites_house_100a$percentile_sick_days[which(total_mosq_bites_house_100a$percentile_sick_days == "")] <- NA
total_mosq_bites_house_100a$percentile_sick_presymp <- as.numeric(as.character(total_mosq_bites_house_100a$percentile_sick_presymp))
total_mosq_bites_house_100a$percentile_sick_days <- as.numeric(as.character(total_mosq_bites_house_100a$percentile_sick_days))

total_mosq_bites_house_100b <- droplevels(total_mosq_bites_house_100a[-c(which(total_mosq_bites_house_100a$day_inf == 6)),])

#########
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1

attractiveness <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
attractiveness_percent <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  att <- data_yes_changing_sens_act[8,xx][[1]]
  attractiveness[x,] <- att
  attractiveness_percent[x,] <- percent_rank(att)
}

attractiveness2 <- melt(attractiveness)
attractiveness2 <- attractiveness2[order(attractiveness2$Var1, attractiveness2$Var2),] 
attractiveness <- attractiveness2$value
attractiveness_percent2 <- melt(attractiveness_percent)
attractiveness_percent2 <- attractiveness_percent2[order(attractiveness_percent2$Var1, attractiveness_percent2$Var2),] 
attractiveness_percent <- attractiveness_percent2$value
rel_mosq_bites_house_100_1_3_new <- cbind(rel_mosq_bites_house_100_1_3, attractiveness)
rel_mosq_bites_house_100_1_3_new2 <- cbind(rel_mosq_bites_house_100_1_3_new, attractiveness_percent)
top20 <- total_mosq_bites_house_100a$top20[which(total_mosq_bites_house_100a$day_inf == 1)]
top20_sick <- total_mosq_bites_house_100a$top20_sick[which(total_mosq_bites_house_100a$day_inf == 1)]
rel_mosq_bites_house_100_1_3_new3 <- cbind(rel_mosq_bites_house_100_1_3_new2, top20)
rel_mosq_bites_house_100_1_3_new4 <- cbind(rel_mosq_bites_house_100_1_3_new3, top20_sick)

perc_rel_mosq_bites_house_100_1_3_new <- cbind(perc_rel_mosq_bites_house_100_1_3, attractiveness)
perc_rel_mosq_bites_house_100_1_3_new2 <- cbind(perc_rel_mosq_bites_house_100_1_3_new, attractiveness_percent)
top20 <- total_mosq_bites_house_100a$top20[which(total_mosq_bites_house_100a$day_inf == 1)]
top20_sick <- total_mosq_bites_house_100a$top20_sick[which(total_mosq_bites_house_100a$day_inf == 1)]
perc_rel_mosq_bites_house_100_1_3_new3 <- cbind(perc_rel_mosq_bites_house_100_1_3_new2, top20)
perc_rel_mosq_bites_house_100_1_3_new4 <- cbind(perc_rel_mosq_bites_house_100_1_3_new3, top20_sick)



percent_home <- (bite_percent_per_day4$mosq_contacts[which(bite_percent_per_day4$day_inf == 6)]) * 100
rel_mosq_bites_house_100_1_3_new5 <- cbind(rel_mosq_bites_house_100_1_3_new4, percent_home)
perc_rel_mosq_bites_house_100_1_3_new5 <- cbind(perc_rel_mosq_bites_house_100_1_3_new4, percent_home)


#####
rel_mosq_bites_house_100_1_3_new6 <- rel_mosq_bites_house_100_1_3_new5[-c(which(is.na(rel_mosq_bites_house_100_1_3_new5$rel_mosq_contacts))),]
perc_rel_mosq_bites_house_100_1_3_new6 <- perc_rel_mosq_bites_house_100_1_3_new5[-c(which(is.na(perc_rel_mosq_bites_house_100_1_3_new5$perc_rel_mosq_contacts))),]

home_count_bites_norm <- home_bites_count_per_day4b$mosq_contacts
rel_mosq_bites_house_100_1_3_new_new <- cbind(rel_mosq_bites_house_100_1_3_new4, home_count_bites_norm)
rel_mosq_bites_house_100_1_3_new_new2 <- rel_mosq_bites_house_100_1_3_new_new[-c(which(is.na(rel_mosq_bites_house_100_1_3_new_new$rel_mosq_contacts))),]

rel_mosq_bites_house_100_1_3_update <- cbind(rel_mosq_bites_house_100_1_3_new4, home_count_bites_norm, percent_home)
rel_mosq_bites_house_100_1_3_update2 <- rel_mosq_bites_house_100_1_3_update[-c(which(is.na(rel_mosq_bites_house_100_1_3_update$rel_mosq_contacts))),]

perc_rel_mosq_bites_house_100_1_3_new_new <- cbind(perc_rel_mosq_bites_house_100_1_3_new4, home_count_bites_norm)
perc_rel_mosq_bites_house_100_1_3_new_new2 <- perc_rel_mosq_bites_house_100_1_3_new_new[-c(which(is.na(perc_rel_mosq_bites_house_100_1_3_new_new$perc_rel_mosq_contacts))),]

perc_rel_mosq_bites_house_100_1_3_update <- cbind(perc_rel_mosq_bites_house_100_1_3_new4, home_count_bites_norm, percent_home)
perc_rel_mosq_bites_house_100_1_3_update2 <- perc_rel_mosq_bites_house_100_1_3_update[-c(which(is.na(perc_rel_mosq_bites_house_100_1_3_update$perc_rel_mosq_contacts))),]

#### R #######
#### REMOVE OUTLIERS FOR VIOLIN PLOT #####
a <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "norm"),])
out_home<- boxplot.stats(a$R[which(a$type == "home")])$out
out_rest<- boxplot.stats(a$R[which(a$type == "rest")])$out
out_all<- boxplot.stats(a$R[which(a$type == "all")])$out

a <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "move"),])
out_home_move<- boxplot.stats(a$R[which(a$type == "home")])$out
out_rest_move<- boxplot.stats(a$R[which(a$type == "rest")])$out
out_all_move<- boxplot.stats(a$R[which(a$type == "all")])$out

R_short <- R_move_yes_changing_melt2
outs <- c(which(R_short$type2 == "norm" & R_short$type == "home" & R_short$R %in% out_home),
          which(R_short$type2 == "norm" & R_short$type == "rest" & R_short$R %in% out_rest),
          which(R_short$type2 == "norm" & R_short$type == "all" & R_short$R %in% out_all),
          which(R_short$type2 == "move" & R_short$type == "home" & R_short$R %in% out_home_move),
          which(R_short$type2 == "move" & R_short$type == "rest" & R_short$R %in% out_rest_move),
          which(R_short$type2 == "move" & R_short$type == "all" & R_short$R %in% out_all_move))
R_short2 <- droplevels(R_short[-outs,])


primary_R_graph <- ggplot(data=droplevels(R_short2[which(R_short2$type != "all"),]), aes(type2, R, fill=type)) +
  geom_boxplot(outlier.color = NA) + 
  scale_x_discrete("Did R Value Calculation Included Movement Change", labels= c("No", "Yes"))  +
  scale_y_continuous("R value") +
  scale_fill_discrete("Where Primary Infectious\nIndividual was Biten",labels=c("Home","Other Houses")) +
  coord_cartesian(ylim=c(0,20)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = c(0.70)),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.key.size = unit(2,"lines"))

favstats(R_short2$R[which(R_short2$type2=="norm")]~R_short2$type[which(R_short2$type2=="norm")])
favstats(R_short2$R[which(R_short2$type2=="move")]~R_short2$type[which(R_short2$type2=="move")])


norm_all <- R_short2$R[which(R_short2$type == "all" & R_short2$type2=="norm")]
move_all <- R_short2$R[which(R_short2$type == "all" & R_short2$type2=="move")]
norm_home <- R_short2$R[which(R_short2$type == "home" & R_short2$type2=="norm")]
move_home <- R_short2$R[which(R_short2$type == "home" & R_short2$type2=="move")]
norm_rest <- R_short2$R[which(R_short2$type == "rest" & R_short2$type2=="norm")]
move_rest <- R_short2$R[which(R_short2$type == "rest" & R_short2$type2=="move")]

primary_R_stats_short <- rbind(favstats(norm_all), favstats(move_all), 
                         favstats(norm_home), favstats(move_home), 
                         favstats(norm_rest), favstats(move_rest))
row.names(primary_R_stats_short) <- c("norm_all","move_all",
                                "norm_home", "move_home",
                                "norm_rest", "move_rest")

write.csv(primary_R_stats_short, "primary_R_stats_short.csv")


##### Overall R Boxplot #####
overall_R_graph <- ggplot(data=droplevels(R_short2[which(R_short2$type == "all"),]), aes(type2, R)) +
  geom_boxplot(outlier.color = NA) + 
  scale_x_discrete("Did R Value Calculation Included Movement Change", labels= c("No", "Yes"))  +
  scale_y_continuous("R value") +
  coord_cartesian(ylim=c(0,20)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = c(0.70)),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.key.size = unit(2,"lines"))

tiff("overall_R_graph.tiff", units="in", width=6, height=6, res=300)
overall_R_graph
dev.off()


R_temp <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type == "all" & R_move_yes_changing_melt2$inf_no_yes==1),])
favstats(R_temp$R[which(R_temp$type2=="norm")]) # norm
favstats(R_temp$R[which(R_temp$type2=="move")]) 

favstats(R_short2$R[which(R_short2$type == "all" & R_short2$type2 == "norm")])
favstats(R_short2$R[which(R_short2$type == "all" & R_short2$type2 == "move")])


R_temp2 <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$inf_no_yes==1),])
favstats(R_temp2$R[which(R_temp2$type2=="norm")]~R_temp2$type[which(R_temp2$type2=="norm")]) # norm
favstats(R_temp2$R[which(R_temp2$type2=="move")]~R_temp2$type[which(R_temp2$type2=="move")]) 

favstats(R_short2$R[which(R_short2$type2 == "norm")] ~ R_short2$type[which(R_short2$type2 == "norm")])
favstats(R_short2$R[which(R_short2$type2 == "move")] ~ R_short2$type[which(R_short2$type2 == "move")])

###
R_move_yes_changing_melt3 <- droplevels(R_move_yes_changing_melt2[-c(which(R_move_yes_changing_melt2$inf_no_yes == 0)),])
R_move_yes_changing_melt3a <- R_move_yes_changing_melt3[order(R_move_yes_changing_melt3[,1], R_move_yes_changing_melt3[,2], R_move_yes_changing_melt3[,4]),]
R_move_yes_changing_melt4 <- R_move_yes_changing_melt3a[which(R_move_yes_changing_melt3a$type2=="norm"),]
R_move <- R_move_yes_changing_melt3a$R[which(R_move_yes_changing_melt3a$type2 == "move")]
R_move_yes_changing_melt4$R_move <- R_move
R_move_yes_changing_melt4$R_change <- R_move_yes_changing_melt4$R_move - R_move_yes_changing_melt4$R

R_move_yes_changing_melt2_new <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2=="move"),])
R_move_yes_changing_melt2_new <- droplevels(R_move_yes_changing_melt2_new[order(R_move_yes_changing_melt2_new$rep, R_move_yes_changing_melt2_new$person, R_move_yes_changing_melt2_new$type),])


R_move_yes_changing_melt4a <- droplevels(R_move_yes_changing_melt4[which(R_move_yes_changing_melt4$type=="all"),])
favstats(R_move_yes_changing_melt4a$R) # norm
favstats(R_move_yes_changing_melt4a$R_change)
favstats(((R_move_yes_changing_melt4a$R_change)/(R_move_yes_changing_melt4a$R))*100) # percent change

favstats(R_move_yes_changing_melt4$R~R_move_yes_changing_melt4$type) # norm
favstats(R_move_yes_changing_melt4$R_change~R_move_yes_changing_melt4$type)
R_move_yes_changing_melt4$R_change_percent <- (((R_move_yes_changing_melt4$R_change)/(R_move_yes_changing_melt4$R))*100)
favstats(R_move_yes_changing_melt4$R_change_percent~R_move_yes_changing_melt4$type) # percent change

###
a <- droplevels(R_move_yes_changing_melt4)
out_home<- boxplot.stats(a$R[which(a$type == "home")])$out
out_rest<- boxplot.stats(a$R[which(a$type == "rest")])$out
out_all<- boxplot.stats(a$R[which(a$type == "all")])$out

R_short_new <- R_move_yes_changing_melt4
outs <- c(which(R_short_new$type2 == "norm" & R_short_new$type == "home" & R_short_new$R %in% out_home),
          which(R_short_new$type2 == "norm" & R_short_new$type == "rest" & R_short_new$R %in% out_rest),
          which(R_short_new$type2 == "norm" & R_short_new$type == "all" & R_short_new$R %in% out_all))
R_short_new2 <- droplevels(R_short_new[-outs,])


R_short_new2a <- droplevels(R_short_new2[which(R_short_new2$type=="all"),])

prim_R_stats2 <- favstats(R_short_new2$R~R_short_new2$type) # norm
prim_R_rel_stats2 <- favstats(R_short_new2$R_change~R_short_new2$type)
R_short_new2$R_change_percent <- (((R_short_new2$R_change)/(R_short_new2$R))*100)
prim_R_perc_rel_stats2 <- favstats(R_short_new2$R_change_percent~R_short_new2$type) # percent change

write.csv(prim_R_stats2, "prim_R_stats2.csv")
write.csv(prim_R_rel_stats2, "prim_R_rel_stats2.csv")
write.csv(prim_R_perc_rel_stats2, "prim_R_perc_rel_stats2.csv")


#####
change_R_graph<- ggplot(data=R_short_new2, aes(type,R_change)) +
  geom_hline(yintercept = 0, col="red") +
  geom_boxplot() +
  scale_x_discrete("Where Primary Bites Occured", labels= c("Home", "Other Houses", "All"))  +
  scale_y_continuous("Change in R value", breaks=seq(-20,30,10)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = c(0.70)),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.key.size = unit(2,"lines"))

tiff("change_R_graph.tiff", units="in", width=6, height=6, res=300)
change_R_graph
dev.off()

###########################################
######## MODELS ###########
###########################################
#### MOVE ####
R_move_yes_changing_melt_new <- R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "move"),]
R_move_yes_changing_melt_new$person <- as.numeric(as.character(R_move_yes_changing_melt_new$person))
R_move_yes_changing_melt_new2 <- droplevels(R_move_yes_changing_melt_new[order(R_move_yes_changing_melt_new$rep, R_move_yes_changing_melt_new$person),])
R_move_yes_changing_melt_new_home <- droplevels(R_move_yes_changing_melt_new2[which(R_move_yes_changing_melt_new2$type=="home"),])
R_move_yes_changing_melt_new_all <- droplevels(R_move_yes_changing_melt_new2[which(R_move_yes_changing_melt_new2$type=="all"),])

H_deg <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  H_deg1 <- data_yes_changing_sens_act[5,x][[1]][,1]
  H_deg[x,] <- H_deg1
}
H_deg2 <- melt(H_deg)
H_deg2$rep <- H_deg2$Var1
H_deg2$person <- H_deg2$Var2
H_deg2$deg <- H_deg2$value
H_deg3 <- H_deg2[,c(4:6)]
H_deg3a <- H_deg3[order(H_deg3$rep, H_deg3$person),]

place_deg <- H_deg3a$deg
place_deg2 <- place_deg[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_mosq_num <- total_mosq_bites_house_100_1_3$norm_mosq[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_count_bites_norm <- home_bites_count_per_day4b$mosq_contacts[which(home_bites_count_per_day4b$inf_no_yes==1)]
percent_home <- (bite_percent_per_day4$mosq_contacts[which(bite_percent_per_day4$day_inf == 6 & bite_percent_per_day4$inf_no_yes==1)]) * 100
attractiveness_percent<- rel_mosq_bites_house_100_1_3_new5$attractiveness_percent[which(rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1)]
R_change <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all")]
R_change_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="home")]
R_change_perc <- R_move_yes_changing_melt4$R_change_percent[which(R_move_yes_changing_melt4$type=="all")]
R_change_home_perc <- R_move_yes_changing_melt4$R_change_percent[which(R_move_yes_changing_melt4$type=="home")]

temp <- cbind(R_move_yes_changing_melt_new_all$rep, R_move_yes_changing_melt_new_all$person)
temp2 <- merge(temp, mosq_per_rep_yes_changing_house_100d, by.x=c("V1","V2"), by.y=c("rep","person"), all.x=TRUE)
temp2 <- temp2[order(temp2$V1, temp2$V2),]
total_bites <- temp2$total_mosq

R_move_yes_changing_melt_new_all2 <- cbind(R_move_yes_changing_melt_new_all, percent_home, home_count_bites_norm, attractiveness_percent, place_deg2, home_mosq_num, R_change, R_change_perc, total_bites)
R_move_yes_changing_melt_new_home2 <- cbind(R_move_yes_changing_melt_new_home, percent_home, home_count_bites_norm, attractiveness_percent, place_deg2, home_mosq_num, R_change_home, R_change_home_perc, total_bites)

R_move_yes_changing_melt_new_all2_short <- droplevels(R_move_yes_changing_melt_new_all2[which(R_move_yes_changing_melt_new_all2$home_mosq_num < 55),]) 
# 


qb <- gam(R~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num,k=9), data=R_move_yes_changing_melt_new_all2_short, method="REML")
qb0<- gam(R~s(percent_home), data=R_move_yes_changing_melt_new_all2_short, method="REML")
qb1<- gam(R~s(attractiveness_percent), data=R_move_yes_changing_melt_new_all2_short, method="REML")
qb2<- gam(R~s(home_mosq_num, k=9), data=R_move_yes_changing_melt_new_all2_short, method="REML")

qb6aa<- gam(R~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num, k=9) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_all2_short, method="REML")

qb5baa<- gam(R~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num, k=9)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_all2_short, method="REML")
qb5baa1<- gam(R~ti(home_mosq_num,percent_home)+
                s(home_mosq_num, k=9)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
qb5baa2<- gam(R~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num, k=9)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
# 
deviance <- c(summary(qb0)$dev.expl, summary(qb2)$dev.expl, summary(qb1)$dev.expl,
              summary(qb)$dev.expl, summary(qb5baa)$dev.expl, summary(qb5baa2)$dev.expl, 
              summary(qb5baa1)$dev.expl,summary(qb6aa)$dev.expl)
R.sq <- c(summary(qb0)$r.sq, summary(qb2)$r.sq, summary(qb1)$r.sq,
          summary(qb)$r.sq, summary(qb5baa)$r.sq, summary(qb5baa2)$r.sq, 
          summary(qb5baa1)$r.sq, summary(qb6aa)$r.sq)
form <- c(as.character(summary(qb0)$formula)[3],as.character(summary(qb2)$formula)[3],as.character(summary(qb1)$formula)[3],
          as.character(summary(qb)$formula)[3],as.character(summary(qb5baa)$formula)[3],as.character(summary(qb5baa2)$formula)[3],
          as.character(summary(qb5baa1)$formula)[3],as.character(summary(qb6aa)$formula)[3])
temp_tab <- cbind(AIC(qb0,qb2,qb1,qb,qb5baa,qb5baa2,qb5baa1,qb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_move_AICc <- AICctab(qb0,qb2,qb1,qb,qb5baa,qb5baa2,qb5baa1,qb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_all2_short,"R_move_model_data.csv")
write.csv(R_move_AICc, "R_move_model_AICc.csv")
write.csv(temp_tab, "R_move_model.csv")
###### NORM #######
R_move_yes_changing_melt_new_norm1 <- R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "norm"),]
R_move_yes_changing_melt_new_norm <- R_move_yes_changing_melt_new_norm1[which(R_move_yes_changing_melt_new_norm1$inf_no_yes == 1),]
R_move_yes_changing_melt_new_norm$person <- as.numeric(as.character(R_move_yes_changing_melt_new_norm$person))
R_move_yes_changing_melt_new_norm2 <- droplevels(R_move_yes_changing_melt_new_norm[order(R_move_yes_changing_melt_new_norm$rep, R_move_yes_changing_melt_new_norm$person),])
R_move_yes_changing_melt_new_norm_home <- droplevels(R_move_yes_changing_melt_new_norm2[which(R_move_yes_changing_melt_new_norm2$type=="home"),])
R_move_yes_changing_melt_new_norm_all <- droplevels(R_move_yes_changing_melt_new_norm2[which(R_move_yes_changing_melt_new_norm2$type=="all"),])

H_deg <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  H_deg1 <- data_yes_changing_sens_act[5,x][[1]][,1]
  H_deg[x,] <- H_deg1
}
H_deg2 <- melt(H_deg)
H_deg2$rep <- H_deg2$Var1
H_deg2$person <- H_deg2$Var2
H_deg2$deg <- H_deg2$value
H_deg3 <- H_deg2[,c(4:6)]
H_deg3a <- H_deg3[order(H_deg3$rep, H_deg3$person),]

place_deg <- H_deg3a$deg
place_deg2 <- place_deg[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_mosq_num <- total_mosq_bites_house_100_1_3$norm_mosq[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_count_bites_norm <- home_bites_count_per_day4b$mosq_contacts[which(home_bites_count_per_day4b$inf_no_yes==1)]
percent_home <- (bite_percent_per_day4$mosq_contacts[which(bite_percent_per_day4$day_inf == 6 & bite_percent_per_day4$inf_no_yes==1)]) * 100
attractiveness_percent<- rel_mosq_bites_house_100_1_3_new5$attractiveness_percent[which(rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1)]

temp <- cbind(R_move_yes_changing_melt_new_norm_all$rep, R_move_yes_changing_melt_new_norm_all$person)
temp2 <- merge(temp, mosq_per_rep_yes_changing_house_100d, by.x=c("V1","V2"), by.y=c("rep","person"), all.x=TRUE)
temp2 <- temp2[order(temp2$V1, temp2$V2),]
total_bites <- temp2$total_mosq

R_move_yes_changing_melt_new_norm_all2 <- cbind(R_move_yes_changing_melt_new_norm_all, percent_home, home_count_bites_norm, attractiveness_percent, place_deg2, home_mosq_num,total_bites)
R_move_yes_changing_melt_new_norm_home2 <- cbind(R_move_yes_changing_melt_new_norm_home, percent_home, home_count_bites_norm, attractiveness_percent, place_deg2, home_mosq_num, total_bites)

R_move_yes_changing_melt_new_norm_all2_short <- droplevels(R_move_yes_changing_melt_new_norm_all2[which(R_move_yes_changing_melt_new_norm_all2$home_mosq_num < 55),]) 
R_move_yes_changing_melt_new_norm_all2_short_total <- droplevels(R_move_yes_changing_melt_new_norm_all2[which(R_move_yes_changing_melt_new_norm_all2$total_bites < 165),]) 


pb <- gam(R~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num,k=9)+s(total_bites), data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb0<- gam(R~s(percent_home), data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb1<- gam(R~s(attractiveness_percent), data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb2<- gam(R~s(home_mosq_num, k=9), data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb3<- gam(R~s(total_bites), data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")

pb6aa<- gam(R~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num, k=9) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")

pb5baa<- gam(R~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num, k=9)+s(attractiveness_percent)+
               s(percent_home)+s(total_bites), 
             data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb5baa1<- gam(R~ti(home_mosq_num,percent_home)+
                s(home_mosq_num, k=9)+s(attractiveness_percent)+
                s(percent_home)+s(total_bites), 
              data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
pb5baa2<- gam(R~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num, k=9)+s(attractiveness_percent)+
                s(percent_home)+s(total_bites), 
              data=R_move_yes_changing_melt_new_norm_all2_short_total, method="REML")
# 
deviance <- c(summary(pb0)$dev.expl, summary(pb2)$dev.expl, summary(pb3)$dev.expl, summary(pb1)$dev.expl,
              summary(pb)$dev.expl, summary(pb5baa)$dev.expl, summary(pb5baa2)$dev.expl, 
              summary(pb5baa1)$dev.expl,summary(pb6aa)$dev.expl)
R.sq <- c(summary(pb0)$r.sq, summary(pb2)$r.sq, summary(pb3)$r.sq, summary(pb1)$r.sq,
          summary(pb)$r.sq, summary(pb5baa)$r.sq, summary(pb5baa2)$r.sq, 
          summary(pb5baa1)$r.sq, summary(pb6aa)$r.sq)
form <- c(as.character(summary(pb0)$formula)[3],as.character(summary(pb2)$formula)[3],as.character(summary(pb3)$formula)[3], as.character(summary(pb1)$formula)[3],
          as.character(summary(pb)$formula)[3],as.character(summary(pb5baa)$formula)[3],as.character(summary(pb5baa2)$formula)[3],
          as.character(summary(pb5baa1)$formula)[3],as.character(summary(pb6aa)$formula)[3])
temp_tab <- cbind(AIC(pb0,pb2,pb3,pb1,pb,pb5baa,pb5baa2,pb5baa1,pb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_norm_AICc <- AICctab(pb0,pb2,pb3,pb1,pb,pb5baa,pb5baa2,pb5baa1,pb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_norm_all2_short_total,"R_norm_model_data.csv")
write.csv(R_norm_AICc, "R_norm_model_AICc.csv")
write.csv(temp_tab, "R_norm_model.csv")

#### MOVE HOME ####
R_move_yes_changing_melt_new_home2$rep <- factor(R_move_yes_changing_melt_new_home2$rep, levels=c(1:Nreps_yes_changing_house_100))
R_move_yes_changing_melt_new_home2_short <- droplevels(R_move_yes_changing_melt_new_home2[which(R_move_yes_changing_melt_new_home2$home_mosq_num < 55),]) 

vb <- gam(R~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num, k=20), data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb0<- gam(R~s(percent_home), data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb1<- gam(R~s(attractiveness_percent, k=15), data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb2<- gam(R~s(home_mosq_num), data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb6aa<- gam(R~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb5baa<- gam(R~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb5baa1<- gam(R~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")
vb5baa2<- gam(R~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")

deviance <- c(summary(vb0)$dev.expl, summary(vb2)$dev.expl, summary(vb1)$dev.expl,
              summary(vb)$dev.expl, summary(vb5baa)$dev.expl, summary(vb5baa2)$dev.expl, 
              summary(vb5baa1)$dev.expl,summary(vb6aa)$dev.expl)
R.sq <- c(summary(vb0)$r.sq, summary(vb2)$r.sq, summary(vb1)$r.sq,
          summary(vb)$r.sq, summary(vb5baa)$r.sq, summary(vb5baa2)$r.sq, 
          summary(vb5baa1)$r.sq, summary(vb6aa)$r.sq)
form <- c(as.character(summary(vb0)$formula)[3],as.character(summary(vb2)$formula)[3], as.character(summary(vb1)$formula)[3],
          as.character(summary(vb)$formula)[3],as.character(summary(vb5baa)$formula)[3],as.character(summary(vb5baa2)$formula)[3],
          as.character(summary(vb5baa1)$formula)[3],as.character(summary(vb6aa)$formula)[3])
temp_tab <- cbind(AIC(vb0,vb2,vb1,vb,vb5baa,vb5baa2,vb5baa1,vb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_move_home_AICc <- AICctab(vb0,vb2,vb1,vb,vb5baa,vb5baa2,vb5baa1,vb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_home2_short,"R_move_home_model_data.csv")
write.csv(R_move_home_AICc, "R_move_home_model_AICc.csv")
write.csv(temp_tab, "R_move_home_model.csv")

##### R CHANGE #####
#### ALL ####
R_move_yes_changing_melt_new_all2_short <- droplevels(R_move_yes_changing_melt_new_all2[which(R_move_yes_changing_melt_new_all2$home_mosq_num < 55),]) 
R_move_yes_changing_melt_new_all2_short2 <- droplevels(R_move_yes_changing_melt_new_all2_short[which(R_move_yes_changing_melt_new_all2_short$percent_home >9),]) 


gb <- gam(R_change~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb0<- gam(R_change~s(percent_home), data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb1<- gam(R_change~s(attractiveness_percent), data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb2<- gam(R_change~s(home_mosq_num), data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb5baa<- gam(R_change~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb5baa1<- gam(R_change~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb5baa2<- gam(R_change~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
gb6aa<- gam(R_change~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_all2_short, method="REML")

deviance <- c(summary(gb0)$dev.expl, summary(gb2)$dev.expl, summary(gb1)$dev.expl,
              summary(gb)$dev.expl, summary(gb5baa)$dev.expl, summary(gb5baa2)$dev.expl, 
              summary(gb5baa1)$dev.expl,summary(gb6aa)$dev.expl)
R.sq <- c(summary(gb0)$r.sq, summary(gb2)$r.sq,summary(gb1)$r.sq,
          summary(gb)$r.sq, summary(gb5baa)$r.sq, summary(gb5baa2)$r.sq, 
          summary(gb5baa1)$r.sq, summary(gb6aa)$r.sq)
form <- c(as.character(summary(gb0)$formula)[3],as.character(summary(gb2)$formula)[3], as.character(summary(gb1)$formula)[3],
          as.character(summary(gb)$formula)[3],as.character(summary(gb5baa)$formula)[3],as.character(summary(gb5baa2)$formula)[3],
          as.character(summary(gb5baa1)$formula)[3],as.character(summary(gb6aa)$formula)[3])
temp_tab <- cbind(AIC(gb0,gb2,gb1,gb,gb5baa,gb5baa2,gb5baa1,gb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_change_AICc <- AICctab(gb0,gb2,gb1,gb,gb5baa,gb5baa2,gb5baa1,gb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_all2_short,"R_change_model_data.csv")
write.csv(R_change_AICc, "R_change_model_AICc.csv")
write.csv(temp_tab, "R_change_model.csv")


#
#### HOME ####
R_move_yes_changing_melt_new_home2_short <- droplevels(R_move_yes_changing_melt_new_home2[which(R_move_yes_changing_melt_new_home2$home_mosq_num < 55),]) 

fb <- gam(R_change_home~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb0<- gam(R_change_home~s(percent_home), data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb1<- gam(R_change_home~s(attractiveness_percent), data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb2<- gam(R_change_home~s(home_mosq_num), data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb5baa<- gam(R_change_home~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb5baa1<- gam(R_change_home~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb5baa2<- gam(R_change_home~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")
fb6aa<- gam(R_change_home~ti(home_mosq_num,attractiveness_percent,percent_home) + 
              ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) +
              ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_home2_short, method="REML")

deviance <- c(summary(fb0)$dev.expl, summary(fb2)$dev.expl,  summary(fb1)$dev.expl,
              summary(fb)$dev.expl, summary(fb5baa)$dev.expl, summary(fb5baa2)$dev.expl, 
              summary(fb5baa1)$dev.expl,summary(fb6aa)$dev.expl)
R.sq <- c(summary(fb0)$r.sq, summary(fb2)$r.sq,  summary(fb1)$r.sq,
          summary(fb)$r.sq, summary(fb5baa)$r.sq, summary(fb5baa2)$r.sq, 
          summary(fb5baa1)$r.sq, summary(fb6aa)$r.sq)
form <- c(as.character(summary(fb0)$formula)[3],as.character(summary(fb2)$formula)[3], as.character(summary(fb1)$formula)[3],
          as.character(summary(fb)$formula)[3],as.character(summary(fb5baa)$formula)[3],as.character(summary(fb5baa2)$formula)[3],
          as.character(summary(fb5baa1)$formula)[3],as.character(summary(fb6aa)$formula)[3])
temp_tab <- cbind(AIC(fb0,fb2,fb1,fb,fb5baa,fb5baa2,fb5baa1,fb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_change_home_AICc <- AICctab(fb0,fb2,fb1,fb,fb5baa,fb5baa2,fb5baa1,fb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_home2_short,"R_change_home_model_data.csv")
write.csv(R_change_home_AICc, "R_change_home_model_AICc.csv")
write.csv(temp_tab, "R_change_home_model.csv")


#
##### R PERCENT CHANGE #####
#### ALL ####
R_move_yes_changing_melt_new_all2_short <- droplevels(R_move_yes_changing_melt_new_all2[which(R_move_yes_changing_melt_new_all2$home_mosq_num < 55),]) 

hb <- gam(R_change_perc~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb0<- gam(R_change_perc~s(percent_home), data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb1<- gam(R_change_perc~s(attractiveness_percent), data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb2<- gam(R_change_perc~s(home_mosq_num), data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb5baa<- gam(R_change_perc~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb5baa1<- gam(R_change_perc~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb5baa2<- gam(R_change_perc~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_all2_short, method="REML")
hb6aa<- gam(R_change_perc~ti(home_mosq_num,attractiveness_percent,percent_home) + 
              ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) +
              ti(attractiveness_percent, percent_home) +
              s(home_mosq_num)+s(attractiveness_percent)+
              s(percent_home), 
            data=R_move_yes_changing_melt_new_all2_short, method="REML")

deviance <- c(summary(hb0)$dev.expl, summary(hb2)$dev.expl, summary(hb1)$dev.expl,
              summary(hb)$dev.expl, summary(hb5baa)$dev.expl, summary(hb5baa2)$dev.expl, 
              summary(hb5baa1)$dev.expl,summary(hb6aa)$dev.expl)
R.sq <- c(summary(hb0)$r.sq, summary(hb2)$r.sq, summary(hb1)$r.sq,
          summary(hb)$r.sq, summary(hb5baa)$r.sq, summary(hb5baa2)$r.sq, 
          summary(hb5baa1)$r.sq, summary(hb6aa)$r.sq)
form <- c(as.character(summary(hb0)$formula)[3],as.character(summary(hb2)$formula)[3],as.character(summary(hb1)$formula)[3],
          as.character(summary(hb)$formula)[3],as.character(summary(hb5baa)$formula)[3],as.character(summary(hb5baa2)$formula)[3],
          as.character(summary(hb5baa1)$formula)[3],as.character(summary(hb6aa)$formula)[3])
temp_tab <- cbind(AIC(hb0,hb2,hb1,hb,hb5baa,hb5baa2,hb5baa1,hb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_change_perc_AICc <- AICctab(hb0,hb2,hb1,hb,hb5baa,hb5baa2,hb5baa1,hb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_all2_short,"R_change_perc_model_data.csv")
write.csv(R_change_perc_AICc, "R_change_perc_model_AICc.csv")
write.csv(temp_tab, "R_change_perc_model.csv")

#
#### HOME ####
R_move_yes_changing_melt_new_home2_short <- droplevels(R_move_yes_changing_melt_new_home2[which(R_move_yes_changing_melt_new_home2$home_mosq_num < 55),]) 

nnb1 <- gam(R_change_home_perc~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num)+s(place_deg2), 
            data=R_move_yes_changing_melt_new_home2_short, method="REML")

nb <- gam(R_change_home_perc~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb0<- gam(R_change_home_perc~s(percent_home), data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb1<- gam(R_change_home_perc~s(attractiveness_percent), data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb2<- gam(R_change_home_perc~s(home_mosq_num), data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb5baa<- gam(R_change_home_perc~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb5baa1<- gam(R_change_home_perc~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")
nb5baa2<- gam(R_change_home_perc~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=R_move_yes_changing_melt_new_home2_short, method="REML")

nb6aa<- gam(R_change_home_perc~ti(home_mosq_num,attractiveness_percent,percent_home) + 
              ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) +
              ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=R_move_yes_changing_melt_new_home2_short, method="REML")

deviance <- c(summary(nb0)$dev.expl, summary(nb2)$dev.expl, summary(nb1)$dev.expl,
              summary(nb)$dev.expl, summary(nb5baa)$dev.expl, summary(nb5baa2)$dev.expl, 
              summary(nb5baa1)$dev.expl,summary(nb6aa)$dev.expl)
R.sq <- c(summary(nb0)$r.sq, summary(nb2)$r.sq, summary(nb1)$r.sq,
          summary(nb)$r.sq, summary(nb5baa)$r.sq, summary(nb5baa2)$r.sq, 
          summary(nb5baa1)$r.sq, summary(nb6aa)$r.sq)
form <- c(as.character(summary(nb0)$formula)[3],as.character(summary(nb2)$formula)[3],as.character(summary(nb1)$formula)[3],
          as.character(summary(nb)$formula)[3],as.character(summary(nb5baa)$formula)[3],as.character(summary(nb5baa2)$formula)[3],
          as.character(summary(nb5baa1)$formula)[3],as.character(summary(nb6aa)$formula)[3])
temp_tab <- cbind(AIC(nb0,nb2,nb1,nb,nb5baa,nb5baa2,nb5baa1,nb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
R_change_home_perc_AICc <- AICctab(nb0,nb2,nb1,nb,nb5baa,nb5baa2,nb5baa1,nb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(R_move_yes_changing_melt_new_home2_short,"R_change_home_perc_model_data.csv")
write.csv(R_change_home_perc_AICc, "R_change_home_perc_model_AICc.csv")
write.csv(temp_tab, "R_change_home_perc_model.csv")

#
##############################
###### REL BITES ######
H_deg <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  H_deg1 <- data_yes_changing_sens_act[5,x][[1]][,1]
  H_deg[x,] <- H_deg1
}
H_deg2 <- melt(H_deg)
H_deg2$rep <- H_deg2$Var1
H_deg2$person <- H_deg2$Var2
H_deg2$deg <- H_deg2$value
H_deg3 <- H_deg2[,c(4:6)]
H_deg3a <- H_deg3[order(H_deg3$rep, H_deg3$person),]
place_deg <- H_deg3a$deg
place_deg2 <- place_deg[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_mosq_num <- total_mosq_bites_house_100_1_3$norm_mosq[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_count_bites_norm <- home_bites_count_per_day4b$mosq_contacts[which(home_bites_count_per_day4b$inf_no_yes==1)]
percent_home <- (bite_percent_per_day4$mosq_contacts[which(bite_percent_per_day4$day_inf == 6 & bite_percent_per_day4$inf_no_yes==1)]) * 100
attractiveness_percent<- rel_mosq_bites_house_100_1_3_new5$attractiveness_percent[which(rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1)]

rel_mosq_bites_house_100_1_3_new5_new <- droplevels(rel_mosq_bites_house_100_1_3_new5[which(rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1),])
rel_bites_new <- cbind(rel_mosq_bites_house_100_1_3_new5_new, percent_home, home_count_bites_norm, attractiveness_percent, place_deg2, home_mosq_num)

rel_bites_new_short <- droplevels(rel_bites_new[which(rel_bites_new$home_mosq_num < 55),]) 
# 

jb <- gam(rel_mosq_contacts~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=rel_bites_new_short, method="REML")
jb0<- gam(rel_mosq_contacts~s(percent_home), data=rel_bites_new_short, method="REML")
jb1<- gam(rel_mosq_contacts~s(attractiveness_percent), data=rel_bites_new_short, method="REML")
jb2<- gam(rel_mosq_contacts~s(home_mosq_num), data=rel_bites_new_short, method="REML")
jb5baa<- gam(rel_mosq_contacts~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=rel_bites_new_short, method="REML")
jb5baa1<- gam(rel_mosq_contacts~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=rel_bites_new_short, method="REML")
jb5baa2<- gam(rel_mosq_contacts~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=rel_bites_new_short, method="REML")
jb6aa<- gam(rel_mosq_contacts~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=rel_bites_new_short, method="REML")

deviance <- c(summary(jb0)$dev.expl, summary(jb2)$dev.expl, summary(jb1)$dev.expl,
              summary(jb)$dev.expl, summary(jb5baa)$dev.expl, summary(jb5baa2)$dev.expl, 
              summary(jb5baa1)$dev.expl,summary(jb6aa)$dev.expl)
R.sq <- c(summary(jb0)$r.sq, summary(jb2)$r.sq, summary(jb1)$r.sq,
          summary(jb)$r.sq, summary(jb5baa)$r.sq, summary(jb5baa2)$r.sq, 
          summary(jb5baa1)$r.sq, summary(jb6aa)$r.sq)
form <- c(as.character(summary(jb0)$formula)[3],as.character(summary(jb2)$formula)[3],as.character(summary(jb1)$formula)[3],
          as.character(summary(jb)$formula)[3],as.character(summary(jb5baa)$formula)[3],as.character(summary(jb5baa2)$formula)[3],
          as.character(summary(jb5baa1)$formula)[3],as.character(summary(jb6aa)$formula)[3])
temp_tab <- cbind(AIC(jb0,jb2,jb1,jb,jb5baa,jb5baa2,jb5baa1,jb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
rel_bites_AICc <- AICctab(jb0,jb2,jb1,jb,jb5baa,jb5baa2,jb5baa1,jb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(rel_bites_new_short,"rel_bites_model_data.csv")
write.csv(rel_bites_AICc, "rel_bites_model_AICc.csv")
write.csv(temp_tab, "rel_bites_model.csv")
##############################
###### PERC REL BITES ######
H_deg <- matrix(NA, nrow=Nreps_yes_changing_house_100 , ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){ ## get values from all reps
  H_deg1 <- data_yes_changing_sens_act[5,x][[1]][,1]
  H_deg[x,] <- H_deg1
}
H_deg2 <- melt(H_deg)
H_deg2$rep <- H_deg2$Var1
H_deg2$person <- H_deg2$Var2
H_deg2$deg <- H_deg2$value
H_deg3 <- H_deg2[,c(4:6)]
H_deg3a <- H_deg3[order(H_deg3$rep, H_deg3$person),]
place_deg <- H_deg3a$deg
place_deg2 <- place_deg[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_mosq_num <- total_mosq_bites_house_100_1_3$norm_mosq[which(total_mosq_bites_house_100_1_3$inf_no_yes=="1")]
home_count_bites_norm <- home_bites_count_per_day4b$mosq_contacts[which(home_bites_count_per_day4b$inf_no_yes==1)]
percent_home <- (bite_percent_per_day4$mosq_contacts[which(bite_percent_per_day4$day_inf == 6 & bite_percent_per_day4$inf_no_yes==1)]) * 100
attractiveness_percent<- rel_mosq_bites_house_100_1_3_new5$attractiveness_percent[which(rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1)]

perc_rel_mosq_bites_house_100_1_3_new5_new <- droplevels(perc_rel_mosq_bites_house_100_1_3_new5[which(perc_rel_mosq_bites_house_100_1_3_new5$inf_no_yes==1),])
perc_rel_bites_new <- cbind(perc_rel_mosq_bites_house_100_1_3_new5_new, home_count_bites_norm, place_deg2, home_mosq_num)

perc_rel_bites_new_short <- droplevels(perc_rel_bites_new[which(perc_rel_bites_new$home_mosq_num < 55),]) 
perc_rel_bites_new_short$perc_rel_mosq_contacts <- (perc_rel_bites_new_short$perc_rel_mosq_contacts)*100


kb <- gam(perc_rel_mosq_contacts~s(attractiveness_percent)+s(percent_home)+s(home_mosq_num), data=perc_rel_bites_new_short, method="REML")
kb0<- gam(perc_rel_mosq_contacts~s(percent_home), data=perc_rel_bites_new_short, method="REML")
kb1<- gam(perc_rel_mosq_contacts~s(attractiveness_percent), data=perc_rel_bites_new_short, method="REML")
kb2<- gam(perc_rel_mosq_contacts~s(home_mosq_num), data=perc_rel_bites_new_short, method="REML")
kb5baa<- gam(perc_rel_mosq_contacts~ti(home_mosq_num,attractiveness_percent)+
               s(home_mosq_num)+s(attractiveness_percent)+
               s(percent_home), 
             data=perc_rel_bites_new_short, method="REML")
kb5baa1<- gam(perc_rel_mosq_contacts~ti(home_mosq_num,percent_home)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=perc_rel_bites_new_short, method="REML")
kb5baa2<- gam(perc_rel_mosq_contacts~ti(percent_home,attractiveness_percent)+
                s(home_mosq_num)+s(attractiveness_percent)+
                s(percent_home), 
              data=perc_rel_bites_new_short, method="REML")
kb6aa<- gam(perc_rel_mosq_contacts~ti(home_mosq_num,attractiveness_percent,percent_home) + ti(home_mosq_num, attractiveness_percent) +
              ti(home_mosq_num, percent_home) + ti(attractiveness_percent, percent_home) +
              s(home_mosq_num) + s(attractiveness_percent) + s(percent_home),
            data=perc_rel_bites_new_short, method="REML")

deviance <- c(summary(kb0)$dev.expl, summary(kb2)$dev.expl, summary(kb1)$dev.expl,
              summary(kb)$dev.expl, summary(kb5baa)$dev.expl, summary(kb5baa2)$dev.expl, 
              summary(kb5baa1)$dev.expl,summary(kb6aa)$dev.expl)
R.sq <- c(summary(kb0)$r.sq, summary(kb2)$r.sq, summary(kb1)$r.sq,
          summary(kb)$r.sq, summary(kb5baa)$r.sq, summary(kb5baa2)$r.sq, 
          summary(kb5baa1)$r.sq, summary(kb6aa)$r.sq)
form <- c(as.character(summary(kb0)$formula)[3],as.character(summary(kb2)$formula)[3],as.character(summary(kb1)$formula)[3],
          as.character(summary(kb)$formula)[3],as.character(summary(kb5baa)$formula)[3],as.character(summary(kb5baa2)$formula)[3],
          as.character(summary(kb5baa1)$formula)[3],as.character(summary(kb6aa)$formula)[3])
temp_tab <- cbind(AIC(kb0,kb2,kb1,kb,kb5baa,kb5baa2,kb5baa1,kb6aa), R.sq,(deviance*100), form) 
temp_tab2 <- temp_tab[order(temp_tab$R.sq),]
perc_rel_bites_AICc <- AICctab(kb0,kb2,kb1,kb,kb5baa,kb5baa2,kb5baa1,kb6aa, weights=TRUE, base=TRUE, sort=FALSE)

write.csv(perc_rel_bites_new_short,"perc_rel_bites_model_data.csv")
write.csv(perc_rel_bites_AICc, "perc_rel_bites_model_AICc.csv")
write.csv(temp_tab, "perc_rel_bites_model.csv")
#
###################
###################################################################
###################################################################
##### END SECOND CODE #####
###################################################################
##################################################################
##### SECONDARY #######
#######
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1


R_norm_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_norm_secondary_home_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_norm_secondary_rest_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  R_norm_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[6,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_norm_secondary_home_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[17,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_norm_secondary_rest_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[18,xx][[1]]) ## xx is replicate; y is the day of infectiousness
}

R_norm_yes_changing_100_melt <- melt(R_norm_yes_changing_100, na.rm=TRUE)
colnames(R_norm_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_yes_changing_100_melt$rep <- factor(as.character(R_norm_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_secondary_home_yes_changing_100_melt <- melt(R_norm_secondary_home_yes_changing_100, na.rm=TRUE)
colnames(R_norm_secondary_home_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_secondary_home_yes_changing_100_melt$rep <- factor(as.character(R_norm_secondary_home_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_secondary_rest_yes_changing_100_melt <- melt(R_norm_secondary_rest_yes_changing_100, na.rm=TRUE)
colnames(R_norm_secondary_rest_yes_changing_100_melt) <- c("rep", "person","R_norm")
R_norm_secondary_rest_yes_changing_100_melt$rep <- factor(as.character(R_norm_secondary_rest_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_norm_yes_changing_100_melt$type <- as.character("all")
R_norm_secondary_home_yes_changing_100_melt$type <- as.character("secondary_home")
R_norm_secondary_rest_yes_changing_100_melt$type <- as.character("secondary_rest")

R_norm_yes_changing_melt <- rbind(R_norm_yes_changing_100_melt, R_norm_secondary_home_yes_changing_100_melt, R_norm_secondary_rest_yes_changing_100_melt)
R_norm_yes_changing_melt$type <- factor(R_norm_yes_changing_melt$type, levels=c( "secondary_home", "secondary_rest", "all"))

#######
reps <- numeric(length=Nreps)
for(x in 1:Nreps){
  if(data_yes_changing_sens_act[1,x][[1]][T] > 4){
    reps[x] <- x
  }
}
reps2 <- reps[which(reps!=0)]
Nreps_yes_changing_house_100 <- length(reps2)
Nh_temp <- (data_yes_changing_sens_act[1,1][[1]][1]) + 1


R_move_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_move_secondary_home_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
R_move_secondary_rest_yes_changing_100 <- matrix(NA, nrow=Nreps_yes_changing_house_100, ncol=Nh_temp)
for(x in 1:Nreps_yes_changing_house_100){
  xx <- reps2[x]
  R_move_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[7,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_move_secondary_home_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[19,xx][[1]]) ## xx is replicate; y is the day of infectiousness
  R_move_secondary_rest_yes_changing_100[x,] <- rowSums(data_yes_changing_sens_act[20,xx][[1]]) ## xx is replicate; y is the day of infectiousness
}

R_move_yes_changing_100_melt <- melt(R_move_yes_changing_100, na.rm=TRUE)
colnames(R_move_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_yes_changing_100_melt$rep <- factor(as.character(R_move_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_secondary_home_yes_changing_100_melt <- melt(R_move_secondary_home_yes_changing_100, na.rm=TRUE)
colnames(R_move_secondary_home_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_secondary_home_yes_changing_100_melt$rep <- factor(as.character(R_move_secondary_home_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_secondary_rest_yes_changing_100_melt <- melt(R_move_secondary_rest_yes_changing_100, na.rm=TRUE)
colnames(R_move_secondary_rest_yes_changing_100_melt) <- c("rep", "person","R_move")
R_move_secondary_rest_yes_changing_100_melt$rep <- factor(as.character(R_move_secondary_rest_yes_changing_100_melt$rep), levels=c(1:Nreps_yes_changing_house_100))

R_move_yes_changing_100_melt$type <- as.character("all")
R_move_secondary_home_yes_changing_100_melt$type <- as.character("secondary_home")
R_move_secondary_rest_yes_changing_100_melt$type <- as.character("secondary_rest")

R_move_yes_changing_melt <- rbind(R_move_yes_changing_100_melt, R_move_secondary_home_yes_changing_100_melt, R_move_secondary_rest_yes_changing_100_melt)
R_move_yes_changing_melt$type <- factor(R_move_yes_changing_melt$type, levels=c( "secondary_home", "secondary_rest", "all"))


inf_top10_rep_person <- total_bites_per_day_abs4_yes_changing_house_100[,c(1,2,5,6)]
inf_top10_rep_person <- inf_top10_rep_person[!duplicated(inf_top10_rep_person[,c(1,2)]),]

R_norm_yes_changing_melt10 <- merge(R_norm_yes_changing_melt,inf_top10_rep_person, by=c("rep","person"), all.x = TRUE)
R_move_yes_changing_melt10 <- merge(R_move_yes_changing_melt,inf_top10_rep_person, by=c("rep","person"), all.x = TRUE)

R_norm_yes_changing_melt10$type2 <- as.character("norm")
R_move_yes_changing_melt10$type2 <- as.character("move")
colnames(R_norm_yes_changing_melt10) <- c("rep","person","R","type","inf_no_yes", "top10", "type2")
colnames(R_move_yes_changing_melt10) <- c("rep","person","R","type", "inf_no_yes", "top10", "type2")

R_move_yes_changing_melt2 <- rbind(R_norm_yes_changing_melt10, R_move_yes_changing_melt10)
R_move_yes_changing_melt2$type2 <- factor(R_move_yes_changing_melt2$type2, levels=c("norm","move"))

R_move_yes_changing_melt2c <- droplevels(R_move_yes_changing_melt2[-c(which(R_move_yes_changing_melt2$type == "all")),])
# 
R_move_yes_changing_melt3 <- droplevels(R_move_yes_changing_melt2[-c(which(R_move_yes_changing_melt2$inf_no_yes == 0)),])



favstats(R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="norm")])
favstats(R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="move")])
# 
norm_secondary_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="norm")]
move_secondary_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="move")]
# 
norm_secondary_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_home" & R_move_yes_changing_melt3$type2=="norm")]
move_secondary_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_home" & R_move_yes_changing_melt3$type2=="move")]

favstats(R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_rest" & R_move_yes_changing_melt3$type2=="norm")])
favstats(R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_rest" & R_move_yes_changing_melt3$type2=="move")])


norm_secondary_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="norm")]
move_secondary_all <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "all" & R_move_yes_changing_melt3$type2=="move")]
norm_secondary_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_home" & R_move_yes_changing_melt3$type2=="norm")]
move_secondary_home <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_home" & R_move_yes_changing_melt3$type2=="move")]
norm_secondary_rest <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_rest" & R_move_yes_changing_melt3$type2=="norm")]
move_secondary_rest <- R_move_yes_changing_melt3$R[which(R_move_yes_changing_melt3$type == "secondary_rest" & R_move_yes_changing_melt3$type2=="move")]

secondary_R_stats <- rbind(favstats(norm_secondary_all), favstats(move_secondary_all), 
                           favstats(norm_secondary_home), favstats(move_secondary_home), 
                           favstats(norm_secondary_rest), favstats(move_secondary_rest))
row.names(secondary_R_stats) <- c("norm_secondary_all","move_secondary_all",
                                  "norm_secondary_home", "move_secondary_home",
                                  "norm_secondary_rest", "move_secondary_rest")

write.csv(secondary_R_stats, "secondary_R_stats.csv")
#####
R_move_yes_changing_melt3$top10 <- factor(R_move_yes_changing_melt3$top10, levels=c("0","1"))
R_move_yes_changing_melt3a <- R_move_yes_changing_melt3[order(R_move_yes_changing_melt3[,1], R_move_yes_changing_melt3[,2], R_move_yes_changing_melt3[,4]),]
R_move_yes_changing_melt4 <- R_move_yes_changing_melt3a[which(R_move_yes_changing_melt3a$type2=="norm"),]
R_move <- R_move_yes_changing_melt3a$R[which(R_move_yes_changing_melt3a$type2 == "move")]
R_move_yes_changing_melt4$R_move <- R_move
R_move_yes_changing_melt4$R_change <- R_move_yes_changing_melt4$R_move - R_move_yes_changing_melt4$R
# 

#######
low_secondary_all <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)]
high_secondary_all <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)]
low_secondary_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==0)]
high_secondary_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==1)]
low_secondary_rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==0)]
high_secondary_rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==1)]
secondary_all  <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all")]
secondary_home <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home")]
secondary_rest <- R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest")]


low_secondary_all_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==0)]))*100
high_secondary_all_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all"& R_move_yes_changing_melt4$top10==1)]))*100
low_secondary_home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==0)]))*100
high_secondary_home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_home"& R_move_yes_changing_melt4$top10==1)]))*100
low_secondary_rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==0)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==0)]))*100
high_secondary_rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==1)])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_rest"& R_move_yes_changing_melt4$top10==1)]))*100
secondary_all_perc  <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="all")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="all")]))*100
secondary_home_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_home")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_home")]))*100
secondary_rest_perc <- ((R_move_yes_changing_melt4$R_change[which(R_move_yes_changing_melt4$type=="secondary_rest")])/(R_move_yes_changing_melt4$R[which(R_move_yes_changing_melt4$type=="secondary_rest")]))*100

secondary_R_change_stats <- rbind(favstats(secondary_all), favstats(secondary_home), favstats(secondary_rest),
                                  favstats(low_secondary_all), favstats(low_secondary_home), favstats(low_secondary_rest),
                                  favstats(high_secondary_all), favstats(high_secondary_home), favstats(high_secondary_rest),
                                  favstats(secondary_all_perc), favstats(secondary_home_perc), favstats(secondary_rest_perc),
                                  favstats(low_secondary_all_perc), favstats(low_secondary_home_perc), favstats(low_secondary_rest_perc),
                                  favstats(high_secondary_all_perc), favstats(high_secondary_home_perc), favstats(high_secondary_rest_perc))
row.names(secondary_R_change_stats) <- c("secondary_all","secondary_home","secondary_rest",
                                         "low_secondary_all","low_secondary_home","low_secondary_rest",
                                         "high_secondary_all","high_secondary_home","high_secondary_rest",
                                         "secondary_all_perc","secondary_home_perc","secondary_rest_perc",
                                         "low_secondary_all_perc","low_secondary_home_perc","low_secondary_rest_perc",
                                         "high_secondary_all_perc","high_secondary_home_perc","high_secondary_rest_perc")


write.csv(secondary_R_change_stats, "secondary_R_change_stats.csv")
###################################################################
##################################################################
##### SECOND CODE #####
###################################################################
###################################################################
#### secondary -- REMOVE OUTLIERS FOR VIOLIN PLOT #####
a <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "norm"),])
out_secondary_home<- boxplot.stats(a$R[which(a$type == "secondary_home")])$out
out_secondary_rest<- boxplot.stats(a$R[which(a$type == "secondary_rest")])$out
out_all<- boxplot.stats(a$R[which(a$type == "all")])$out

a <- droplevels(R_move_yes_changing_melt2[which(R_move_yes_changing_melt2$type2 == "move"),])
out_secondary_home_move<- boxplot.stats(a$R[which(a$type == "secondary_home")])$out
out_secondary_rest_move<- boxplot.stats(a$R[which(a$type == "secondary_rest")])$out
out_all_move<- boxplot.stats(a$R[which(a$type == "all")])$out

R_short <- R_move_yes_changing_melt2
outs <- c(which(R_short$type2 == "norm" & R_short$type == "secondary_home" & R_short$R %in% out_secondary_home),
          which(R_short$type2 == "norm" & R_short$type == "secondary_rest" & R_short$R %in% out_secondary_rest),
          which(R_short$type2 == "norm" & R_short$type == "all" & R_short$R %in% out_all),
          which(R_short$type2 == "move" & R_short$type == "secondary_home" & R_short$R %in% out_secondary_home_move),
          which(R_short$type2 == "move" & R_short$type == "secondary_rest" & R_short$R %in% out_secondary_rest_move),
          which(R_short$type2 == "move" & R_short$type == "all" & R_short$R %in% out_all_move))
R_short2 <- droplevels(R_short[-outs,])

favstats(R_short$R[which(R_short$type2 == "norm")]~R_short$type[which(R_short$type2=="norm")])
favstats(R_short2$R[which(R_short2$type2 == "norm")]~R_short2$type[which(R_short2$type2=="norm")])

favstats(R_short$R[which(R_short$type2 == "move")]~R_short$type[which(R_short$type2=="move")])
favstats(R_short2$R[which(R_short2$type2 == "move")]~R_short2$type[which(R_short2$type2=="move")])

secondary_R_graph <- ggplot(data=droplevels(R_short2[which(R_short2$type != "all"),]), aes(type2, R, fill=type)) +
  geom_boxplot(outlier.color = NA) + 
  scale_x_discrete("Did R Value Calculation Included Movement Change", labels= c("No", "Yes"))  +
  scale_y_continuous("R value") +
  scale_fill_discrete("Where Secondary\nInfection Occured",labels=c("Home of Primary\nInfected Individual","Other Houses")) +
  coord_cartesian(ylim=c(0,20)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size=16, lineheight = c(0.70)),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.key.size = unit(2,"lines"))


favstats(R_short2$R[which(R_short2$type2 == "norm")]~R_short2$type[which(R_short2$type2=="norm")])
favstats(R_short2$R[which(R_short2$type2 == "move")]~R_short2$type[which(R_short2$type2=="move")])

norm_all <- R_short2$R[which(R_short2$type == "all" & R_short2$type2=="norm")]
move_all <- R_short2$R[which(R_short2$type == "all" & R_short2$type2=="move")]
norm_home <- R_short2$R[which(R_short2$type == "secondary_home" & R_short2$type2=="norm")]
move_home <- R_short2$R[which(R_short2$type == "secondary_home" & R_short2$type2=="move")]
norm_rest <- R_short2$R[which(R_short2$type == "secondary_rest" & R_short2$type2=="norm")]
move_rest <- R_short2$R[which(R_short2$type == "secondary_rest" & R_short2$type2=="move")]

secondary_R_stats_short <- rbind(favstats(norm_all), favstats(move_all), 
                               favstats(norm_home), favstats(move_home), 
                               favstats(norm_rest), favstats(move_rest))
row.names(secondary_R_stats_short) <- c("norm_all","move_all",
                                      "norm_home", "move_home",
                                      "norm_rest", "move_rest")

write.csv(secondary_R_stats_short, "secondary_R_stats_short.csv")

##### GRAPH BOTH #######
primary_R_graph2 <- primary_R_graph +   
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=20, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=18),
        legend.key.size = unit(2,"lines"))

secondary_R_graph2 <- secondary_R_graph +   
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=20, lineheight=0.75),
        legend.title.align=0.5,
        legend.text = element_text(size=18),
        legend.key.size = unit(2,"lines"))

vps <- baseViewports()
pushViewport(vps$figure)
vp1 <-plotViewport(c(1,3,1,1)) ## create new vp with margins, you play with this values

tiff("R_primary_secondary_graphs.tiff", units="in", width=16, height=8, res=300)
grid.arrange(primary_R_graph2, secondary_R_graph2, ncol=2, widths=c(1.05,1), vp=vp1)
dev.off()

#####
R_short_new <- R_move_yes_changing_melt4
outs <- c(which(R_short_new$type2 == "norm" & R_short_new$type == "secondary_home" & R_short_new$R %in% out_home),
          which(R_short_new$type2 == "norm" & R_short_new$type == "secondary_rest" & R_short_new$R %in% out_rest),
          which(R_short_new$type2 == "norm" & R_short_new$type == "all" & R_short_new$R %in% out_all))
R_short_new2 <- droplevels(R_short_new[-outs,])


R_short_new2a <- droplevels(R_short_new2[which(R_short_new2$type=="all"),])

sec_R_stats2 <- favstats(R_short_new2$R~R_short_new2$type) # norm
sec_R_rel_stats2 <- favstats(R_short_new2$R_change~R_short_new2$type)
R_short_new2$R_change_percent <- (((R_short_new2$R_change)/(R_short_new2$R))*100)
q <- R_short_new2[-which(R_short_new2$R_change_percent==Inf),]
sec_R_perc_rel_stats2 <- favstats(R_short_new2$R_change_percent~R_short_new2$type) # percent change

write.csv(sec_R_stats2, "sec_R_stats2.csv")
write.csv(sec_R_rel_stats2, "sec_R_rel_stats2.csv")
