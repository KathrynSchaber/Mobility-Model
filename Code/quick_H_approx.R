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
  ## CHANGE HM to have rowSums(HM) == 1
  home_body <- which(degs == 0)
  for(i in 1:length(home_body)){
    name <- home_body[i]
    house_col <- name_homes[name]
    HM[name, house_col] <- 1
  }
  HM2 <- HM
  not_home_body <- which(degs != 0)
  percent_time_home <- numeric(length=600)
  percent_time_home[home_body] <- 1
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
  UpoorHostMixingRossMacdonald = makeUrossMacdonald(HpoorHostMixing)
  U <- UpoorHostMixingRossMacdonald
