
distanceMatrix = function(x, y)
{
  distances = t(sapply(1 : length(x), function(ii) sqrt((x[ii] - x) ^ 2 + (y[ii] - y) ^ 2)))
  
  for(ii in 1 : length(x))
  {
    distances[ii, ii] = 0
  } # for diagonal element of the matrix
  
  return(distances)
} # end function distanceMatrix definition



estimateLarvalEquilib = function()
{
  MM = array(0, dim = c(Nl, xi, T))
  Smf = rep(1, Nf)
  
  for (tt in 1 : T)
  {
    # mosquito movement and reproduction
    Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L) # uninfected 
    Smf = rmultinomMatrix(Sml, F) # uninfected recruitment, mortality, and movement
    M_larvae = cbind(rpois(Nl, v * Sml), M_larvae[, c(1 : xi)[-xi]]) # new eggs and advancement through larval stages
    M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi) # larval mortality
    
    # store variables over time
    MM[,, tt] = M_larvae
  } # end for loop over time
  
  return(matrix(sapply(1 : xi, function(ii) sapply(1 : Nl, function(jj) mean(MM[jj, ii, floor(T / 2) : T]))), Nl, xi))
} # end function estimateMosquitoEquilibLarval definition



estimateMosquitoEquilibFeeding = function()
{
  SSmf = matrix(0, Nf, T)
  Smf = rep(1, Nf)
  
  for (tt in 1 : T)
  {
    # mosquito movement and reproduction
    Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L) # uninfected 
    Smf = rmultinomMatrix(Sml, F) # uninfected recruitment, mortality, and movement
    M_larvae = cbind(rpois(Nl, v * Sml), M_larvae[, c(1 : xi)[-xi]]) # new eggs and advancement through larval stages
    M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi) # larval mortality
    
    # store variables over time
    SSmf[, tt] = Smf
  } # end for loop over time
  
  return(sapply(1 : Nf, function(ii) mean(SSmf[ii, floor(T / 2) : T])))
} # end function estimateMosquitoEquilibFeeding definition



estimateMosquitoEquilibLarval = function()
{
  SSml = matrix(0, Nl, T)
  Smf = rep(1, Nf)
  
  for (tt in 1 : T)
  {
    # mosquito movement and reproduction
    Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L) # uninfected 
    Smf = rmultinomMatrix(Sml, F) # uninfected recruitment, mortality, and movement
    M_larvae = cbind(rpois(Nl, v * Sml), M_larvae[, c(1 : xi)[-xi]]) # new eggs and advancement through larval stages
    M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi) # larval mortality
    
    # store variables over time
    SSml[, tt] = Sml
  } # end for loop over time
  
  return(sapply(1 : Nl, function(ii) mean(SSml[ii, floor(T / 2) : T])))
} # end function estimateMosquitoEquilibLarval definition




makeHpoorHostMixing = function(Nf, Nh, Hh, percentTimeAtHome = .5, distances = distanceMatrix(XFcoord, YFcoord))
{ 
  # Define each host to spend percentTimeAtHome at their home
  Hindex = rep(1 : Nf, Hh) 
  H = matrix (data = 0, nrow = Nh, ncol = Nf)
  H[cbind(1 : Nh, Hindex)] = percentTimeAtHome
  
  # Determine how many places besides home a host spends time
  placesPerHost = ceiling(rweibull(Nh, 2, 1.5))
  
  # Define probabilities that inhabitants of i spend some amount of time at j according to the distance between i and j
  probByDistance = dnorm(distances, 0, .1 * 2 * sqrt(1 / pi))
  
  # Randomly pick places where each host spends time other than home
  for(hh in 1 : Nh)
  {
    H[hh, sample((1 : Nf)[-Hindex[hh]], placesPerHost[hh], prob = probByDistance[Hindex[hh], -Hindex[hh]])] = (1 - percentTimeAtHome) / placesPerHost[hh]
  } # end for
  
  return(H)
} # end function makeHW definition



makeHrossMacdonald = function(Nf, Nh, Hh)
{ 
  Hindex = rep(1 : Nf, Hh) # a list of each person's home 
  H = matrix(data = 1 / Nf, nrow = Nh, ncol = Nf) # humans spend an equal amount of time at each feeding station
  
  return(H)
} # end function makeHrossMacdonald definition



makeUrossMacdonald = function(H)
{
  omega = normalize(rep(1 / Nh, Nh))
  
  return(diag(omega) %*% H %*% diag(1 / colSums(diag(omega) %*% H)))
} # end function makeUrossMacdonald definition



makeUunevenBiting = function(H)
{
  omega = normalize(rexp(Nh, rate = bitingRate))
  
  return(list(omega, diag(omega) %*% H %*% diag(1 / colSums(diag(omega) %*% H))))
} # end function makeUunevenBiting definition


####### I changed "dnorm(dd, 0, .025 * 2 * sqrt(1 / pi))" to "dnorm(dd, 0, .01 * 2 * sqrt(1 / pi))"
####### 2*sqrt(1/pi)) is the diameter of the circle with x/y coords
####### dnorm gives density for quantiles, dd, with mean 0 and sd= ____ (0.025*diameter or 0.01*diameter)
####### since there are 100 sites, put sd as 1/100 with diameter
mosquitoMovement = function(x, y, X, Y, surv, mixing = 'poor')
# Generates a movement matrix M between the points located at
# coordinates (x, y) and those at (X, Y). Movement occurs between
# a point (x, y) and its N nearest neigbhors in the X-Y set.
# A proportion (1 - surv) perish.
{
  di = function(i)
  # Returns a movement vector from (x[i], y[i]) to all (X, Y)
  {
    m = rep(0, length(X))
    dd = sqrt((x[i] - X) ^ 2 + (y[i] - Y) ^ 2)
    if(mixing == 'poor')
    {
      m = surv * normalize(dnorm(dd, 0, .01 * 2 * sqrt(1 / pi)))
    } # end if
    if(mixing == 'well')
    {
      m = rep(surv / length(dd), length(dd))
    } # end if
    
    return(m)
  } # end function di definition
  
  return(t(sapply(1 : length(x), di)))
} # end function mosquitoMovement definition



normalize = function(x)
{
  return(x / sum(x))
} # end function normalize definition



points.clustered = function(n, meanParents = 10, clusteredness = .25)
{
  meanDist = clusteredness / sqrt(meanParents)
  meanChildren = Nf / meanParents
  
  ps = rMatClust(meanParents, meanDist, meanChildren, win = disc(radius = sqrt(1 / pi)))
  while(ps$n != n)
  {
    ps = rMatClust(meanParents, meanDist, meanChildren, win = disc(radius = sqrt(1 / pi)))
  } # end while
  
  return(ps)
} # end function points.clustered definition



points.overdispersed = function(n)
{
  inhibitionFactor = 1
  
  ps = rSSI(inhibitionFactor / sqrt(n), n, win = disc(radius = sqrt(1 / pi)))
  while(ps$n != n)
  {
    inhibitionFactor = inhibitionFactor - .01
    ps = rSSI(inhibitionFactor / sqrt(n), n, win = disc(radius = sqrt(1 / pi)))
  } # end while
  
  return(ps)
} # end function points.overdispersed definition



points.poisson = function(n)
{
  ps = rpoispp(n, win = disc(radius = sqrt(1 / pi)))
  while(ps$n != n)
  {
    ps = rpoispp(n, win = disc(radius = sqrt(1 / pi)))
  } # end while
  
  return(ps)
} # end function points.poisson definition



rmultinomMatrix = function(X, M)
# Redistributes a vector X into an output vector of the same size
# plus one, where the extra element accounts for the fact that the
# probabilities in M may not sum to one.  Redistribution occurs
# according to probability mass functions defined on the rows of M.
#
# The dimension of each object is
#   [X] = 1 x m
# 	[M] = m x n
# 	[P] = m x n + 1
# 	[A] = m x n + 1
# 	[output] = 1 x n + 1
{
  lengthIn = nrow(M)
  lengthOut = ncol(M)
  indexIn = 1 : lengthIn
  probability_remainder = 1 - rowSums(M)
  P = cbind(M, probability_remainder) # add column for remainder so total probability == 1
  
  # Distributes the number contained in X[i] into lengthOut + 1 bins and
  # outputs a vector of length lengthOut.
  rmultinomVector = function(i)
  {
    if (X[i] == 0) {
      return(rep(0, lengthOut))
    } # end if
    else
    {
      return(rmultinom(1, size = X[i], prob = P[i,])[1 : lengthOut])
    } # end else
  } # end function rm definition
  
  A = sapply(indexIn, rmultinomVector)
  return(rowSums(A))
} # end function rmultinomMatrix definition

