lambda_equilib = estimateLarvalEquilib()[, xi]

nSum = 100
Q = matrix(data = 0, nrow = Nf, ncol = Nf)
for (qq in sigma : nSum)
{
  Q = Q + mtx.exp((L %*% F), qq)
} # summing over time since infection in mosquitoes

nSum = 100
M = rep(0, Nf)
for (bb in 0 : (sigma + nSum))
{
  M = M + F %*% mtx.exp((L %*% F), bb)
} # summing over adult mosquito ages
M = lambda_equilib %*% M

B = diag(as.vector(M), Nf, Nf) %*% t(U)

V = t(B) %*% Q %*% t(U)

R_nonlinear = 1 - exp(-b * c * t(B) %*% Q %*% t(U) * rho_mean)
diag(R_nonlinear) = 0

R0_nonlinear = lambda(R_nonlinear)
