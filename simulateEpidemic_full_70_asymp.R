simulateEpidemic_full_70_asymp = function(hh, T = 200)
{
  
  asymp_symp <- character(length=Nh) ## 
  asymp_symp <- factor(asymp_symp, levels=c(0,1)) ## 
  asymp <- character(length=Nh) ## 
  asymp <- factor(asymp, levels=c(0,1)) ## 
  
  
  #### CONFIGURE SN/HM NETOWRK ####
  houses <- rep(c(1:m), House_counts)
  counts <- integer(length=Nh)
  i=1
  counter <- House_counts[i]
  counts[1:(counter)] <- c(1:counter)
  for(i in 2:m){
    counter <- House_counts[i]
    index <- sum(House_counts[1:(i-1)])
    counts[(index+1):(index+counter)] <- c(1:counter)
  }
  names <- character(length=Nh)
  for(i in 1:Nh){
    names[i] <- paste(houses[i],"_", counts[i])
  }
  degs <- rpois(Nh,2.8)
  degs <- rpois(Nh,2.8)
  more_home_body_count <- (round(.15*Nh)) - (length(which(degs == 0)))
  home_body_index <- sample(1:Nh, more_home_body_count, replace=FALSE)
  degs[home_body_index] <- 0
  if(sum(degs) %% 2 != 0){degs[1] <- degs[1] + 1}
  ## RANDOMLY GENERATING NETWORK ##
  stubs <- rep(names[1:Nh], degs)
  stub_length <- length(stubs)
  i = 0
  while( i <= (stub_length-2)/2){
    if(i == 0){
      stub_perm <- sample(stubs)
    }
    temp_a <- 2*i + 1
    temp_b <- 2*i + 2
    pair <- stub_perm[c(temp_a,temp_b)]
    temp_homes <- as.integer(gsub(" _.*", pair, replacement=""))
    i <- i+1
    if(pair[1] == pair[2]){
      i = 0
      next
    }
    if(temp_homes[1] == temp_homes[2]){
      i = 0
      next
    }
  }
  SN <- matrix(0,nrow = Nh, ncol = Nh, dimnames = list(c(names[1:Nh]), c(names[1:Nh])))
  for(i in 0:((stub_length-2)/2)){
    temp_a <- 2*i + 1
    temp_b <- 2*i + 2
    pair <- stub_perm[c(temp_a,temp_b)]
    index1 <- which(names == pair[1])
    index2 <- which(names == pair[2])
    SN[index1, index2] <- 1
    SN[index2, index1] <- 1
  }
  SN1 <- graph_from_adjacency_matrix(SN, mode="undirected")
  ## HOME NETWORK ##
  name_homes <- as.integer(gsub(" _.*", names, replacement=""))
  Homes <- matrix(0,nrow = Nh, ncol = m, dimnames = list(c(names[1:Nh]), c(1:m)))
  for(i in 1:m){
    ## for each home
    rows <- which(name_homes == i)
    Homes[rows,i] <- 1
  }
  ## HUMAN MOVEMENT TO HOUSES --> based on SN and Homes of friends ##
  HM <- SN %*% Homes
  SNOld <- SN
  HMOld <- HM
  ## CHANGE HM to have rowSums(HM) == 1
  home_body_old <- which(degs == 0)
  for(i in 1:length(home_body_old)){
    name <- home_body_old[i]
    house_col <- name_homes[name]
    HM[name, house_col] <- 1
  }
  HM2 <- HM
  not_home_body <- which(degs != 0)
  percent_time_home <- numeric(length=Nh)
  percent_time_home[home_body_old] <- 1
  percent_time_home[not_home_body] <- 0.5
  for(i in 1:length(not_home_body)){
    name <- not_home_body[i]
    house_col <- name_homes[name]
    house_visit_cols <- as.numeric(which(HM[name,] > 0))
    house_visit_amount <- as.numeric(HM[name, house_visit_cols])
    sum_amount <- sum(house_visit_amount)
    time_left <- 1- percent_time_home[name]
    fract <- time_left/sum_amount
    percent_visit <- house_visit_amount*fract
    HM2[name, house_col] <- percent_time_home[name]
    HM2[name, house_visit_cols] <- percent_visit
  }
  H <- as.matrix(HM2, nrow=Nh, ncol=Nf)
  
  HpoorHostMixing <- as.matrix(HM2, nrow=Nh, ncol=Nf)
  UpoorHostMixingUnevenBiting = makeUunevenBiting(HpoorHostMixing)
  H = HpoorHostMixing
  U = UpoorHostMixingUnevenBiting[[2]]
  U_omega <- UpoorHostMixingUnevenBiting[[1]]
  
  ##### ####
  source('variables.R')
  ####
  
  SN_out_deg <- matrix(nrow=Nh, ncol=(rho_max+tau+1))
  ##############################################
  #### Get SN metrics for HEALTHY TIMEPOINT ####
  SN_net <- graph_from_adjacency_matrix(SN, mode="directed", add.rownames = TRUE)
  tt=1
  SN_out_deg[,1] <- degree(SN_net, mode="out") 
  #############################
  H_deg <- matrix(nrow=Nh, ncol=(rho_max+tau+1))
  ##############################################
  #### Get H metrics for HEALTHY TIMEPOINT ####
  H_net <- graph_from_incidence_matrix(H, weighted=TRUE)
  people <- which(V(H_net)$type == FALSE)
  H_deg[,1] <- degree(H_net)[people]
  #############################
  SN_out_deg[hh,2] <- degree(SN_net, mode="out")[hh]
  H_deg[hh,2] <- degree(H_net)[hh]
  ##############################################
  
  
  ####
  houses_index <- which(V(H_net)$type == TRUE)
  H_deg_houses <- numeric(length=Nh)
  
  
  
  
  ## create list with B matrices for each time point during illness
  B_mats <- vector("list", rho_max+1)
  for(j in 1:(rho_max+1)){
    B_mats[[j]] <- matrix(nrow=Nf,ncol=Nh)
  }
  
  ## time point "41" is the normal B values before epidemic (keep at 41 so days of infectiousness 1-40 can have matrices 1-40)
  B_norm <- diag(as.vector(M), Nf, Nf) %*% t(U)
  B_mats[[6]] <- B_norm
  ####
  
  ## create list with B matrices for each time point during illness
  B_rest_mats <- vector("list", rho_max+1)
  for(j in 1:(rho_max+1)){
    B_rest_mats[[j]] <- matrix(nrow=Nf,ncol=Nh)
  }
  ## create list with B matrices for each time point during illness
  B_home_mats <- vector("list", rho_max+1)
  for(j in 1:(rho_max+1)){
    B_home_mats[[j]] <- matrix(nrow=Nf,ncol=Nh)
  }
  
  
  home_U <- matrix(0,nrow=Nh, ncol=Nf)
  rest_U <- U
  for(i in 1:Nh){
    home_U[i,name_homes[i]] <- U[i,name_homes[i]]
    rest_U[i,name_homes[i]] <- 0
  }
  ## time point "41" is the normal B values before epidemic (keep at 41 so days of infectiousness 1-40 can have matrices 1-40)
  B_rest_norm = diag(as.vector(M), Nf, Nf) %*% t(rest_U)
  B_rest_mats[[6]] <- B_rest_norm
  ## time point "41" is the normal B values before epidemic (keep at 41 so days of infectiousness 1-40 can have matrices 1-40)
  B_home_norm = diag(as.vector(M), Nf, Nf) %*% t(home_U)
  B_home_mats[[6]] <- B_home_norm
  
  
  mosq_house_mats <- vector("list", rho_max+1)
  for(j in 1:(rho_max+1)){
    mosq_house_mats[[j]] <- matrix(nrow=4,ncol=Nh)
  }
  ##
  mosq_house_norm <- matrix(nrow=4,ncol=Nh)
  house_index_temp <- name_homes[1:Nh]
  mosq_house_norm[1,] <- c(Smf[house_index_temp])
  mosq_house_norm[2,] <- c(Emf[house_index_temp])
  mosq_house_norm[3,] <- c(Imf[house_index_temp])
  mosq_house_norm[4,] <- c(Smf[house_index_temp] + Emf[house_index_temp] + Imf[house_index_temp])
  mosq_house_mats[[6]] <- mosq_house_norm
  ####
  
  
  
  
  
  
  ## create R_nonlinear_norm matrix
  ## This version has the c_1 values for each time period of infectiousness (5)
  R_nonlinear_norm_symp = 1 - exp(-b * (0.4+0.7+0.4+0.1+0.01) * t(B_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_symp) = 0
  R_nonlinear_norm_asymp = 1 - exp(-b * (0.35+0.6+0.35+0.1+0.01) * t(B_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_asymp) = 0
  
  R_nonlinear_norm_rest_symp = 1 - exp(-b * (0.4+0.7+0.4+0.1+0.01) * t(B_rest_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_rest_symp) = 0
  R_nonlinear_norm_home_symp = 1 - exp(-b * (0.4+0.7+0.4+0.1+0.01) * t(B_home_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_home_symp) = 0
  
  R_nonlinear_norm_rest_asymp = 1 - exp(-b * (0.35+0.6+0.35+0.1+0.01) * t(B_rest_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_rest_asymp) = 0
  R_nonlinear_norm_home_asymp = 1 - exp(-b * (0.35+0.6+0.35+0.1+0.01) * t(B_home_norm) %*% Q %*% t(U))
  diag(R_nonlinear_norm_home_asymp) = 0
  
  
  
  #############
  V_norm_secondary_home_symp <- t(B_norm) %*% Q %*% t(U)
  V_norm_secondary_rest_symp <- t(B_norm) %*% Q %*% t(U)
  V_norm_secondary_home_asymp <- t(B_norm) %*% Q %*% t(U)
  V_norm_secondary_rest_asymp <- t(B_norm) %*% Q %*% t(U)
  for(ii in 1:Nh){
    V_norm_secondary_home_symp[ii,which(name_homes != name_homes[ii])] <- 0
    V_norm_secondary_rest_symp[ii,which(name_homes == name_homes[ii])] <- 0
    
    V_norm_secondary_home_asymp[ii,which(name_homes != name_homes[ii])] <- 0
    V_norm_secondary_rest_asymp[ii,which(name_homes == name_homes[ii])] <- 0
  }
  
  R_nonlinear_norm_secondary_home_symp = 1 - exp(-b * (0.4+0.7+0.4+0.1+0.01) * V_norm_secondary_home_symp)
  diag(R_nonlinear_norm_secondary_home_symp) = 0
  R_nonlinear_norm_secondary_rest_symp = 1 - exp(-b * (0.4+0.7+0.4+0.1+0.01) * V_norm_secondary_rest_symp)
  diag(R_nonlinear_norm_secondary_rest_symp) = 0
  R_nonlinear_norm_secondary_home_asymp = 1 - exp(-b * (0.35+0.6+0.35+0.1+0.01) * V_norm_secondary_home_asymp)
  diag(R_nonlinear_norm_secondary_home_asymp) = 0
  R_nonlinear_norm_secondary_rest_asymp = 1 - exp(-b * (0.35+0.6+0.35+0.1+0.01) * V_norm_secondary_rest_asymp)
  diag(R_nonlinear_norm_secondary_rest_asymp) = 0
  
  
  
  
  ## set up R_nonlinear_movement
  R_nonlinear_movement <- matrix(nrow=Nh, ncol=Nh)
  
  R_nonlinear_movement_rest <- matrix(nrow=Nh, ncol=Nh)
  R_nonlinear_movement_home <- matrix(nrow=Nh, ncol=Nh)
  
  R_nonlinear_movement_secondary_rest <- matrix(nrow=Nh, ncol=Nh)
  R_nonlinear_movement_secondary_home <- matrix(nrow=Nh, ncol=Nh)
  
  ## matrix with Nh rows and 5 columns --> V_k_l values 
  ## for each person at each inf. level to calculate R_movement
  V_k_l_vecs <- vector("list", rho_max)
  for(j in 1:(rho_max)){
    V_k_l_vecs[[j]] <- matrix(nrow=Nh,ncol=Nh)
  }
  
  V_move_rest <- vector("list", rho_max)
  for(j in 1:(rho_max)){
    V_move_rest[[j]] <- matrix(nrow=Nh,ncol=Nh)
  }
  V_move_home <- vector("list", rho_max)
  for(j in 1:(rho_max)){
    V_move_home[[j]] <- matrix(nrow=Nh,ncol=Nh)
  }
  
  V_move_secondary_home <- vector("list", rho_max)
  for(j in 1:(rho_max)){
    V_move_secondary_home[[j]] <- matrix(nrow=Nh,ncol=Nh)
  }
  V_move_secondary_rest <- vector("list", rho_max)
  for(j in 1:(rho_max)){
    V_move_secondary_rest[[j]] <- matrix(nrow=Nh,ncol=Nh)
  }
  
  
  
  
  Sh[hh] = 0
  Ih[hh, 1] = 1
  EIR_cum = rep(0, T + 1)
  EIR_cum[1] = Nh - 1
  Ih_sum = rep(0, T + 1)
  Ih_sum[1] = 1
  
  asymp_symp[hh] <- rbinom(1,1,0.3) ## prob=0.5; prob=0.1 
  asymp[which(asymp_symp == 0)] <- 1
  asymp[which(asymp_symp == 1)] <- 0
  if(asymp_symp[hh] == 0){asymp_index <- hh}
  if(asymp_symp[hh] == 1){asymp_index <- 0}
  
  
  
  ####
  B_temp <- diag(as.vector(M), Nf, Nf) %*% t(U)
  B_mats[[1]][,hh] <- B_temp[,hh]
  
  
  ####
  V_k_l_vecs[[1]][hh,] <- (t(B_temp) %*% Q %*% t(U))[hh,]
  
  ###########
  home_U <- matrix(0,nrow=Nh, ncol=Nf)
  rest_U <- U
  for(i in 1:Nh){
    home_U[i,name_homes[i]] <- U[i,name_homes[i]]
    rest_U[i,name_homes[i]] <- 0
  }
  B_rest_temp = diag(as.vector(M), Nf, Nf) %*% t(rest_U)
  B_home_temp = diag(as.vector(M), Nf, Nf) %*% t(home_U)
  B_rest_mats[[1]][,hh] <- B_rest_temp[,hh]
  B_home_mats[[1]][,hh] <- B_home_temp[,hh]
  
  V_move_rest[[1]][hh,] = (t(B_rest_temp) %*% Q %*% t(U))[hh,]
  V_move_home[[1]][hh,] = (t(B_home_temp) %*% Q %*% t(U))[hh,]
  
  #############
  V_move_secondary_home[[1]][hh,] <- (t(B_temp) %*% Q %*% t(U))[hh,]
  V_move_secondary_rest[[1]][hh,] <- (t(B_temp) %*% Q %*% t(U))[hh,]
  V_move_secondary_home[[1]][hh,which(name_homes != name_homes[hh])] <- 0
  V_move_secondary_rest[[1]][hh,which(name_homes == name_homes[hh])] <- 0
  
  
  
  
  mosq_house_mats[[1]][1,hh] <- c(Smf[name_homes[hh]])
  mosq_house_mats[[1]][2,hh] <- c(Emf[name_homes[hh]])
  mosq_house_mats[[1]][3,hh] <- c(Imf[name_homes[hh]])
  mosq_house_mats[[1]][4,hh] <- c(Smf[name_homes[hh]] + Emf[name_homes[hh]] + Imf[name_homes[hh]])
  
  time_bite <- numeric(length=Nh)
  time_bite[hh] <- 1
  
  
  # loop over time
  for(tt in 2 : (T + 1)){
    ############################
    #### Get network-level SN metrics for each timepoint ####
    SN_net <- graph_from_adjacency_matrix(SN, mode="directed", add.rownames = TRUE)
    H_net <- graph_from_incidence_matrix(H, weighted=TRUE)
    people <- which(V(H_net)$type == FALSE)
    ###################
    #### Get node-level SN metrics for DURING INFECTIOUSNESS (columns 3-7) ####
    inf_1 <- which(Ih[,1] > 0)
    inf_2 <- which(Ih[,2] > 0)
    inf_3 <- which(Ih[,3] > 0)
    inf_4 <- which(Ih[,4] > 0)
    inf_5 <- which(Ih[,5] > 0)
    if(length(inf_1) > 0){
      SN_out_deg[inf_1,3] <- degree(SN_net, mode="out")[inf_1]
      H_deg[inf_1,3] <- degree(H_net)[inf_1]
    }
    if(length(inf_2) > 0){
      SN_out_deg[inf_2,4] <- degree(SN_net, mode="out")[inf_2]
      H_deg[inf_2,4] <- degree(H_net)[inf_2]
    }
    if(length(inf_3) > 0){
      SN_out_deg[inf_3,5] <- degree(SN_net, mode="out")[inf_3]
      H_deg[inf_3,5] <- degree(H_net)[inf_3]
    }
    if(length(inf_4) > 0){
      SN_out_deg[inf_4,6] <- degree(SN_net, mode="out")[inf_4]
      H_deg[inf_4,6] <- degree(H_net)[inf_4]
    }
    if(length(inf_5) > 0){
      SN_out_deg[inf_5,7] <- degree(SN_net, mode="out")[inf_5]
      H_deg[inf_5,7] <- degree(H_net)[inf_5]
    }
    ############################
    
    # mosquito movement from blood feeding habitats to larval habitats
    Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L)
    for(ii in 1 : sigma)
      Eml[, ii] = matrix(rmultinomMatrix(Emf[, ii], L), nrow = Nl, ncol = 1)
    Iml = rmultinomMatrix(Imf, L)
    
    # egg laying and advancement through larval stages
    M_larvae = cbind(rpois(Nl, v * (Sml + rowSums(Eml) + Iml)), M_larvae[, c(1 : xi)[-xi]])
    M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi)
    
    # mosquito movement from larval habitats to blood feeding habitats
    Smf = rmultinomMatrix(Sml, F)
    for(ii in 1 : (sigma - 1))
      Emf[, ii + 1] = matrix(rmultinomMatrix(Eml[, ii], F), nrow = Nf, ncol = 1)
    Imf = rmultinomMatrix(Iml, F) + rmultinomMatrix(Eml[, sigma], F)
    
    ####
    # transmission from infectious hosts to susceptible mosquitoes
    c_1 = rep(0,Nh) # transmission from host to mosquito for each person (length = Nh)
    
    ##### ASYMP #####
    c_1[which((Ih[,1] == 1) & asymp_symp == 0)] = 0.35 ## presymptomatic period
    c_1[which((Ih[,2] == 1) & asymp_symp == 0)] = 0.6
    c_1[which((Ih[,3] == 1) & asymp_symp == 0)] = 0.35
    c_1[which((Ih[,4] == 1) & asymp_symp == 0)] = 0.1
    c_1[which((Ih[,5] == 1) & asymp_symp == 0)] = 0.01
    
    ##### SYMP #####
    c_1[which((Ih[,1] == 1) & asymp_symp == 1)] = 0.4 ## presymptomatic period
    c_1[which((Ih[,2] == 1) & asymp_symp == 1)] = 0.7
    c_1[which((Ih[,3] == 1) & asymp_symp == 1)] = 0.4
    c_1[which((Ih[,4] == 1) & asymp_symp == 1)] = 0.1
    c_1[which((Ih[,5] == 1) & asymp_symp == 1)] = 0.01
    
    if(sum(Ih) > 1)
      Emf[,1] = sapply(1:Nf, function(ii) sum(rbinom(length(which(rowSums(Ih)>0)),rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)],c_1[which(rowSums(Ih) > 0)])))
    if(sum(Ih) == 1)
      Emf[,1] = sapply(1:Nf, function(ii) sum(rbinom(length(which(rowSums(Ih)>0)),rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)],c_1[which(rowSums(Ih) > 0)])))
    ####
    
    if(sum(Ih) == 0)
      Emf[, 1] = 0
    Smf = Smf - Emf[, 1]
    
    # infectious host recovery
    IhOld = Ih
    RhOld = Rh
    for(ii in seq(rho_max - 1, 1, -1))
      Ih[, ii + 1] = rbinom(Nh, Ih[, ii], 1 - rho_fail[ii])
    Rh = Rh + rowSums(IhOld) - rowSums(Ih[, 2 : rho_max])
    rm(IhOld)
    
    # host progression through the pathogen incubation period
    Ih[, 1] = Eh[, tau]
    # for(ii in seq(tau - 1, 1, -1))
    #   Eh[, ii + 1] = Eh[, ii]
    
    # transmission from infectious mosquitoes to susceptible hosts
    secBites_h = rowSums(sapply(1 : Nf, function(ii) rmultinom(1, Imf[ii], U[, ii])))
    Eh[, 1] = rbinom(Nh, Sh, 1 - (1 - b) ^ secBites_h)
    Sh = Sh - Eh[, 1]
    
    
    ############################
    #### Get node-level SN metrics for TIME WHEN BIT (column 2) ####
    SN_net <- graph_from_adjacency_matrix(SN, mode="directed", add.rownames = TRUE)
    H_net <- graph_from_incidence_matrix(H, weighted=TRUE)
    new_bite <- which(Eh[,1] > 0)
    if(length(new_bite) > 0){
      SN_out_deg[new_bite,2] <- degree(SN_net, mode="out")[new_bite]
      H_deg[new_bite,2] <- degree(H_net)[new_bite]
      
      time_bite[new_bite] <- tt
    }
    ############################
    
    EIR_cum[tt] = Nh - sum(Sh)
    Ih_sum[tt] = length(which(rowSums(Ih) > 0))
    
    
    ###### ###### ###### ASYMP vs SYMP ###### ###### ###### 
    ## prob = prob of symptomatic (1)
    ## Kyle 2008 --> Asymptomatic rate 50-90%
    ## Morisson 2010 --> Inapparent rate in Iquitos: mean=.9206, low=.5, high=.9756
    asymp_symp[which(Ih[,1] == 1)] <- rbinom(length(which(Ih[,1] == 1)),1,0.3) ## prob=0.3; prob=0.1 
    asymp[which(asymp_symp == 0)] <- 1
    asymp_index <- c(asymp_index,as.numeric(which(asymp == 1)))
    
    
    ###### ###### ###### MOVEMENT CHANGE ###### ###### ######
    total_houses <- numeric(length=Nh)
    num_houses_replace_3 <- numeric(length=Nh)
    num_houses_replace_4 <- numeric(length=Nh)
    
    for(ii in 1:Nh){
      total_houses[ii] <- length(as.numeric(which(HMOld[ii,] >= 1)))
      num_houses_replace_3[ii] <- round(total_houses[ii]/3)
      num_houses_replace_4[ii] <- round((2*total_houses[ii])/3)
    }
    num_houses_add_back_4 <- num_houses_replace_4 - num_houses_replace_3
    
    for(ii in 1:Nh){
      if(ii %in% home_body_old){}
      else{
        if((Ih[ii,2] == 1) & (asymp_symp[ii] == 1)){ # Ih[ii,1] is presymptomatic period
          percent_time_home[ii] <- 1
          SN[ii,] <- 0
          HM[ii,] <- 0
        }
        if((Ih[ii,3] == 1) & (asymp_symp[ii] == 1)){
          if(num_houses_replace_3[ii] == 0){ ## still going 0 places
            percent_time_home[ii] <- 1
            SN[ii,] <- 0
            HM[ii,] <- 0
          }
          else{
            percent_time_home[ii] <- 0.80
            house_col <- name_homes[ii]
            house_visit_cols_old <- as.numeric(which(HMOld[ii,] >= 1))
            house_visit_amount_old <- as.numeric(HMOld[ii,house_visit_cols_old])
            prob_replace <- (house_visit_amount_old/sum(house_visit_amount_old))
            replace_house <- sample(house_visit_cols_old,num_houses_replace_3[ii],prob=prob_replace)
            people_visit_cols_old <- as.numeric(which(SNOld[ii,] == 1))
            people_visit_houses_old <- name_homes[people_visit_cols_old]
            people_replace <- which(people_visit_houses_old %in% replace_house)
            replace_SN <- people_visit_cols_old[people_replace]
            SN[ii,replace_SN] <- SNOld[ii,replace_SN]
            HM[ii,] <- SN[ii,] %*% Homes
            house_visit_cols_new <- as.numeric(which(HM[ii,] >= 1))
            house_visit_amount_new <- as.numeric(HM[ii,house_visit_cols_new])
          }
        }
        if((Ih[ii,4] == 1) & (asymp_symp[ii] == 1)){
          if(num_houses_add_back_4[ii] == 0){ ## still going 0 places
            percent_time_home[ii] <- 0.70
          }
          else{
            percent_time_home[ii] <- 0.70    
            house_col <- name_homes[ii]
            house_visit_cols_old <- as.numeric(which(HMOld[ii,] >= 1))
            house_visit_amount_old <- as.numeric(HMOld[ii,house_visit_cols_old])
            house_visit_cols <- as.numeric(which(HM[ii,] >= 1))
            house_visit_amount <- as.numeric(HM[ii,house_visit_cols])
            remove_3 <- house_visit_cols_old[which(!(house_visit_cols_old %in% house_visit_cols))] ## houses removed in time step 3
            ## of these that were removed at step 3, the number we want to replace in movement matrix is num_houses_add_back_4
            ## probability of replacing in movement matrix depends on number of house members, i.e. house_visit_amount_old
            houses_removed_3 <- house_visit_cols_old[which(house_visit_cols_old %in% remove_3)]
            houses_removed_3_amount <- house_visit_amount_old[which(house_visit_cols_old %in% remove_3)]
            prob_replace <- houses_removed_3_amount/sum(houses_removed_3_amount) 
            ## not doing "1-value" because want probability of replacing, not removing
            ## more likely to replace those with more people living there
            if(num_houses_add_back_4[ii] == 1){
              replace_house <- houses_removed_3
            }
            else{
              replace_house <- sample(houses_removed_3,num_houses_add_back_4[ii],prob=prob_replace)
            }
            people_visit_cols_old <- as.numeric(which(SNOld[ii,] == 1))
            people_visit_houses_old <- name_homes[people_visit_cols_old]
            people_replace <- which(people_visit_houses_old %in% replace_house)
            replace_SN <- people_visit_cols_old[people_replace]
            SN[ii,replace_SN] <- SNOld[ii,replace_SN]
            HM[ii,] <- SN[ii,] %*% Homes
            house_visit_cols_new <- as.numeric(which(HM[ii,] >= 1))
            house_visit_amount_new <- as.numeric(HM[ii,house_visit_cols_new])
          }
        }
        if((Ih[ii,5] == 1) & (asymp_symp[ii] == 1)){
          percent_time_home[ii] <- 0.5
          SN[ii,] <- SNOld[ii,]
          HM[ii,] <- HMOld[ii,]
          house_visit_cols_new <- as.numeric(which(HM[ii,] >= 1))
          house_visit_amount_new <- as.numeric(HM[ii,house_visit_cols_new])
          ##
          asymp_symp[ii] <- NA
          ##
        }
        if(Rh[ii] > RhOld[ii]){
          percent_time_home[ii] <- 0.5
          SN[ii,] <- SNOld[ii,]
          HM[ii,] <- HMOld[ii,]
          house_visit_cols_new <- as.numeric(which(HM[ii,] >= 1))
          house_visit_amount_new <- as.numeric(HM[ii,house_visit_cols_new])
          ##
          asymp_symp[ii] <- NA
          ##
        }
      }
    }
    home_body_new <- which(rowSums(SN) == 0)
    for(i in 1:length(home_body_new)){
      name <- home_body_new[i]
      house_col <- name_homes[name]
      HM[name, house_col] <- 1
    }
    HM2 <- HM
    
    not_home_body_new <- which(rowSums(SN) != 0)
    for(ii in 1:(length(not_home_body_new))){
      name <- not_home_body_new[ii]
      house_col <- name_homes[name]
      HM2[name, house_col] <- percent_time_home[name]
      house_visit_cols <- as.numeric(which(HM[name,] > 0))
      house_visit_amount <- as.numeric(HM[name, house_visit_cols])
      sum_amount <- sum(house_visit_amount)
      time_left <- 1-percent_time_home[name]
      fract <- time_left/sum_amount
      percent_visit <- house_visit_amount*fract
      HM2[name, house_visit_cols] <- percent_visit
    }
    
    H <- as.matrix(HM2, nrow=Nh, ncol=Nf)
    U <- (diag(U_omega) %*% H %*% diag(1 / colSums(diag(U_omega) %*% H)))
    
    
    
    #### Calculate B matrix for current time point ####
    ## M calculated in "runInfectionsOverTime" before epidemic simulation run
    B_temp <- diag(as.vector(M), Nf, Nf) %*% t(U)
    
    ## for 1,2, ..., rho_max --> --> which people in this infection stage
    ## put their B_temp row into the respective B_I# matrix
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        ## assign B matrix values for those participants ##
        B_mats[[ii]][,ind_temp] <- B_temp[,ind_temp]
      }
    }
    
    
    #### Calculate B matrix for current time point ####
    home_U <- matrix(0,nrow=Nh, ncol=Nf)
    rest_U <- U
    for(i in 1:Nh){
      home_U[i,name_homes[i]] <- U[i,name_homes[i]]
      rest_U[i,name_homes[i]] <- 0
    }
    B_rest_temp = diag(as.vector(M), Nf, Nf) %*% t(rest_U)
    B_home_temp = diag(as.vector(M), Nf, Nf) %*% t(home_U)
    
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        ## assign B matrix values for those participants ##
        B_rest_mats[[ii]][,ind_temp] <- B_rest_temp[,ind_temp]
      }
    }
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        ## assign B matrix values for those participants ##
        B_home_mats[[ii]][,ind_temp] <- B_home_temp[,ind_temp]
      }
    }
    ####
    
    #### Calculate B matrix for current time point ####
    ## M calculated in "runInfectionsOverTime" before epidemic simulation run
    mosq_house_temp <- matrix(nrow=4,ncol=Nh)
    house_index_temp <- name_homes[1:Nh]
    mosq_house_temp[1,] <- c(Smf[house_index_temp])
    mosq_house_temp[2,] <- c(Emf[house_index_temp])
    mosq_house_temp[3,] <- c(Imf[house_index_temp])
    mosq_house_temp[4,] <- c(Smf[house_index_temp] + Emf[house_index_temp] + Imf[house_index_temp])
    ## for 1,2, ..., rho_max --> --> which people in this infection stage
    ## put their B_temp row into the respective B_I# matrix
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        ## assign B matrix values for those participants ##
        mosq_house_mats[[ii]][,ind_temp] <- mosq_house_temp[,ind_temp]
      }
    }
    
    #####
    houses_index <- which(V(H_net)$type == TRUE)
    ind_temp <- which(Ih[,2] == 1)
    ind_homes <- name_homes[ind_temp]
    house_temp <- houses_index[ind_homes]
    if(length(ind_temp) > 0){
      house_deg <- degree(H_net)[house_temp]
      H_deg_houses[ind_temp] <- house_deg
    }
    
    
    
    #### For anyone (hh) in infectiousness stages 1--> 5, caluclate V_k_l_vecs
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        V_k_l_vecs[[ii]][ind_temp,] <- (t(B_temp) %*% Q %*% t(U))[ind_temp,]
      }
    }
    ###########
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        V_move_rest[[ii]][ind_temp,] <- (t(B_rest_temp) %*% Q %*% t(U))[ind_temp,]
      }
    }
    for(ii in 1:rho_max){
      ind_temp <- which(Ih[,ii] == 1)
      if(length(ind_temp) > 0){
        V_move_home[[ii]][ind_temp,] <- (t(B_home_temp) %*% Q %*% t(U))[ind_temp,]
      }
    }
    ####
  }
  
  
  asymp_index2 <- sort(unique(asymp_index))
  if(asymp_index2[1] == 0){asymp_index2 <- asymp_index2[-1]}
  
  
  
  ## Make R_nonlinear_norm, where the people who were asymptomatic have ASYMP infectiousness
  ## and the people who were symptomatic have SYMP infectiousness
  ## set as symp. values
  R_nonlinear_norm <- R_nonlinear_norm_symp 
  ## if person was asymp --> change to asymp. values matrix row
  R_nonlinear_norm[asymp_index2,] <-   R_nonlinear_norm_asymp[asymp_index2,]
  diag(R_nonlinear_norm) <- 0
  
  R_nonlinear_norm_home <- R_nonlinear_norm_home_symp 
  ## if person was asymp --> change to asymp. values matrix row
  R_nonlinear_norm_home[asymp_index2,] <-   R_nonlinear_norm_home_asymp[asymp_index2,]
  diag(R_nonlinear_norm_home) <- 0
  
  R_nonlinear_norm_rest <- R_nonlinear_norm_rest_symp 
  ## if person was asymp --> change to asymp. values matrix row
  R_nonlinear_norm_rest[asymp_index2,] <-   R_nonlinear_norm_rest_asymp[asymp_index2,]
  diag(R_nonlinear_norm_rest) <- 0
  
  
  R_nonlinear_norm_secondary_home <- R_nonlinear_norm_secondary_home_symp
  ## if person was asymp --> change to asymp. values matrix row
  R_nonlinear_norm_secondary_home[asymp_index2,] <-   R_nonlinear_norm_secondary_home_asymp[asymp_index2,]
  diag(R_nonlinear_norm_secondary_home) <- 0
  
  R_nonlinear_norm_secondary_rest <- R_nonlinear_norm_secondary_rest_symp
  ## if person was asymp --> change to asymp. values matrix row
  R_nonlinear_norm_secondary_rest[asymp_index2,] <-   R_nonlinear_norm_secondary_rest_asymp[asymp_index2,]
  diag(R_nonlinear_norm_secondary_rest) <- 0
  
  
  
  
  
  
  ## set up R_nonlinear_movement
  for(ii in 1:Nh){
    ## SYMPTOMATIC
    if(length(which(asymp_index2 == ii)) == 0){
      if(sum(is.na(V_k_l_vecs[[5]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_k_l_vecs[[1]][ii,]) + (0.7 * V_k_l_vecs[[2]][ii,]) + (0.4 * V_k_l_vecs[[3]][ii,]) +
                        (0.1 * V_k_l_vecs[[4]][ii,]) + (0.01 * V_k_l_vecs[[5]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[4]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_k_l_vecs[[1]][ii,]) + (0.7 * V_k_l_vecs[[2]][ii,]) + (0.4 * V_k_l_vecs[[3]][ii,]) +
                        (0.1 * V_k_l_vecs[[4]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[3]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_k_l_vecs[[1]][ii,]) + (0.7 * V_k_l_vecs[[2]][ii,]) + (0.4 * V_k_l_vecs[[3]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[2]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_k_l_vecs[[1]][ii,]) + (0.7 * V_k_l_vecs[[2]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[1]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_k_l_vecs[[1]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
    }
    ## ASYMPTOMATIC
    else if(length(which(asymp_index2 == ii)) == 1){
      if(sum(is.na(V_k_l_vecs[[5]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_k_l_vecs[[1]][ii,]) + (0.6 * V_k_l_vecs[[2]][ii,]) + (0.35 * V_k_l_vecs[[3]][ii,]) +
                        (0.1 * V_k_l_vecs[[4]][ii,]) + (0.01 * V_k_l_vecs[[5]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[4]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_k_l_vecs[[1]][ii,]) + (0.6 * V_k_l_vecs[[2]][ii,]) + (0.35 * V_k_l_vecs[[3]][ii,]) +
                        (0.1 * V_k_l_vecs[[4]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[3]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_k_l_vecs[[1]][ii,]) + (0.6 * V_k_l_vecs[[2]][ii,]) + (0.35 * V_k_l_vecs[[3]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[2]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_k_l_vecs[[1]][ii,]) + (0.6 * V_k_l_vecs[[2]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_k_l_vecs[[1]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_k_l_vecs[[1]][ii,]))
        R_nonlinear_movement[ii,] = 1 - exp(-b * rate_temp)
      }
    }
  }
  
  
  
  ## set up R_nonlinear_movement_home
  for(ii in 1:Nh){
    ## SYMPTOMATIC
    if(length(which(asymp_index2 == ii)) == 0){
      if(sum(is.na(V_move_home[[5]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_home[[1]][ii,]) + (0.7 * V_move_home[[2]][ii,]) + (0.4 * V_move_home[[3]][ii,]) +
                        (0.1 * V_move_home[[4]][ii,]) + (0.01 * V_move_home[[5]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[4]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_home[[1]][ii,]) + (0.7 * V_move_home[[2]][ii,]) + (0.4 * V_move_home[[3]][ii,]) +
                        (0.1 * V_move_home[[4]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[3]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_home[[1]][ii,]) + (0.7 * V_move_home[[2]][ii,]) + (0.4 * V_move_home[[3]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[2]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_home[[1]][ii,]) + (0.7 * V_move_home[[2]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[1]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_home[[1]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
    }
    ## ASYMPTOMATIC
    else if(length(which(asymp_index2 == ii)) == 1){
      if(sum(is.na(V_move_home[[5]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_home[[1]][ii,]) + (0.6 * V_move_home[[2]][ii,]) + (0.35 * V_move_home[[3]][ii,]) +
                        (0.1 * V_move_home[[4]][ii,]) + (0.01 * V_move_home[[5]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[4]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_home[[1]][ii,]) + (0.6 * V_move_home[[2]][ii,]) + (0.35 * V_move_home[[3]][ii,]) +
                        (0.1 * V_move_home[[4]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[3]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_home[[1]][ii,]) + (0.6 * V_move_home[[2]][ii,]) + (0.35 * V_move_home[[3]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[2]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_home[[1]][ii,]) + (0.6 * V_move_home[[2]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_home[[1]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_home[[1]][ii,]))
        R_nonlinear_movement_home[ii,] = 1 - exp(-b * rate_temp)
      }
    }
  }
  
  
  ## set up R_nonlinear_movement_rest
  for(ii in 1:Nh){
    ## SYMPTOMATIC
    if(length(which(asymp_index2 == ii)) == 0){
      if(sum(is.na(V_move_rest[[5]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_rest[[1]][ii,]) + (0.7 * V_move_rest[[2]][ii,]) + (0.4 * V_move_rest[[3]][ii,]) +
                        (0.1 * V_move_rest[[4]][ii,]) + (0.01 * V_move_rest[[5]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[4]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_rest[[1]][ii,]) + (0.7 * V_move_rest[[2]][ii,]) + (0.4 * V_move_rest[[3]][ii,]) +
                        (0.1 * V_move_rest[[4]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[3]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_rest[[1]][ii,]) + (0.7 * V_move_rest[[2]][ii,]) + (0.4 * V_move_rest[[3]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[2]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_rest[[1]][ii,]) + (0.7 * V_move_rest[[2]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[1]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_rest[[1]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
    }
    ## ASYMPTOMATIC
    else if(length(which(asymp_index2 == ii)) == 1){
      if(sum(is.na(V_move_rest[[5]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_rest[[1]][ii,]) + (0.6 * V_move_rest[[2]][ii,]) + (0.35 * V_move_rest[[3]][ii,]) +
                        (0.1 * V_move_rest[[4]][ii,]) + (0.01 * V_move_rest[[5]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[4]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_rest[[1]][ii,]) + (0.6 * V_move_rest[[2]][ii,]) + (0.35 * V_move_rest[[3]][ii,]) +
                        (0.1 * V_move_rest[[4]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[3]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_rest[[1]][ii,]) + (0.6 * V_move_rest[[2]][ii,]) + (0.35 * V_move_rest[[3]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[2]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_rest[[1]][ii,]) + (0.6 * V_move_rest[[2]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_rest[[1]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_rest[[1]][ii,]))
        R_nonlinear_movement_rest[ii,] = 1 - exp(-b * rate_temp)
      }
    }
  }
  
  
  ###########
  ## set up R_nonlinear_movement_secondary_home
  for(ii in 1:Nh){
    ## SYMPTOMATIC
    if(length(which(asymp_index2 == ii)) == 0){
      if(sum(is.na(V_move_secondary_home[[5]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_home[[1]][ii,]) + (0.7 * V_move_secondary_home[[2]][ii,]) + (0.4 * V_move_secondary_home[[3]][ii,]) +
                        (0.1 * V_move_secondary_home[[4]][ii,]) + (0.01 * V_move_secondary_home[[5]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[4]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_home[[1]][ii,]) + (0.7 * V_move_secondary_home[[2]][ii,]) + (0.4 * V_move_secondary_home[[3]][ii,]) +
                        (0.1 * V_move_secondary_home[[4]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[3]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_home[[1]][ii,]) + (0.7 * V_move_secondary_home[[2]][ii,]) + (0.4 * V_move_secondary_home[[3]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[2]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_home[[1]][ii,]) + (0.7 * V_move_secondary_home[[2]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[1]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_home[[1]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
    }
    
    ## ASYMPTOMATIC
    else if(length(which(asymp_index2 == ii)) == 1){
      if(sum(is.na(V_move_secondary_home[[5]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_home[[1]][ii,]) + (0.6 * V_move_secondary_home[[2]][ii,]) + (0.35 * V_move_secondary_home[[3]][ii,]) +
                        (0.1 * V_move_secondary_home[[4]][ii,]) + (0.01 * V_move_secondary_home[[5]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[4]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_home[[1]][ii,]) + (0.6 * V_move_secondary_home[[2]][ii,]) + (0.35 * V_move_secondary_home[[3]][ii,]) +
                        (0.1 * V_move_secondary_home[[4]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[3]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_home[[1]][ii,]) + (0.6 * V_move_secondary_home[[2]][ii,]) + (0.35 * V_move_secondary_home[[3]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[2]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_home[[1]][ii,]) + (0.6 * V_move_secondary_home[[2]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_home[[1]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_home[[1]][ii,]))
        R_nonlinear_movement_secondary_home[ii,] = 1 - exp(-b * rate_temp)
      }
    }
  }
  
  
  
  ## set up R_nonlinear_movement_secondary_rest
  for(ii in 1:Nh){
    ## SYMPTOMATIC
    if(length(which(asymp_index2 == ii)) == 0){
      if(sum(is.na(V_move_secondary_rest[[5]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_rest[[1]][ii,]) + (0.7 * V_move_secondary_rest[[2]][ii,]) + (0.4 * V_move_secondary_rest[[3]][ii,]) +
                        (0.1 * V_move_secondary_rest[[4]][ii,]) + (0.01 * V_move_secondary_rest[[5]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[4]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_rest[[1]][ii,]) + (0.7 * V_move_secondary_rest[[2]][ii,]) + (0.4 * V_move_secondary_rest[[3]][ii,]) +
                        (0.1 * V_move_secondary_rest[[4]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[3]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_rest[[1]][ii,]) + (0.7 * V_move_secondary_rest[[2]][ii,]) + (0.4 * V_move_secondary_rest[[3]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[2]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_rest[[1]][ii,]) + (0.7 * V_move_secondary_rest[[2]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[1]][ii,])) < Nh){
        rate_temp <- ((0.4 * V_move_secondary_rest[[1]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
    }
    ## ASYMPTOMATIC
    else if(length(which(asymp_index2 == ii)) == 1){
      if(sum(is.na(V_move_secondary_rest[[5]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_rest[[1]][ii,]) + (0.6 * V_move_secondary_rest[[2]][ii,]) + (0.35 * V_move_secondary_rest[[3]][ii,]) +
                        (0.1 * V_move_secondary_rest[[4]][ii,]) + (0.01 * V_move_secondary_rest[[5]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[4]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_rest[[1]][ii,]) + (0.6 * V_move_secondary_rest[[2]][ii,]) + (0.35 * V_move_secondary_rest[[3]][ii,]) +
                        (0.1 * V_move_secondary_rest[[4]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[3]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_rest[[1]][ii,]) + (0.6 * V_move_secondary_rest[[2]][ii,]) + (0.35 * V_move_secondary_rest[[3]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[2]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_rest[[1]][ii,]) + (0.6 * V_move_secondary_rest[[2]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
      else if(sum(is.na(V_move_secondary_rest[[1]][ii,])) < Nh){
        rate_temp <- ((0.35 * V_move_secondary_rest[[1]][ii,]))
        R_nonlinear_movement_secondary_rest[ii,] = 1 - exp(-b * rate_temp)
      }
    }
  }
  
  
  #######################
  
  
  diag(R_nonlinear_movement_rest) = 0
  diag(R_nonlinear_movement_home) = 0
  diag(R_nonlinear_movement_secondary_rest) = 0
  diag(R_nonlinear_movement_secondary_home) = 0
  diag(R_nonlinear_movement) = 0
  
  return(list(EIR_cum, Ih_sum, B_mats, asymp_index2, SN_out_deg, 
              H_deg, R_nonlinear_norm_symp, R_nonlinear_norm, R_nonlinear_movement, U_omega,
              mosq_house_mats, name_homes, H_deg_houses, 
              R_nonlinear_norm_rest_symp, R_nonlinear_norm_home_symp,
              R_nonlinear_norm_rest, R_nonlinear_norm_home,
              R_nonlinear_movement_home, R_nonlinear_movement_rest, time_bite,
              R_nonlinear_norm_secondary_home_symp, R_nonlinear_norm_secondary_rest_symp,
              R_nonlinear_norm_secondary_home, R_nonlinear_norm_secondary_rest,
              R_nonlinear_movement_secondary_home, R_nonlinear_movement_secondary_rest))
}

