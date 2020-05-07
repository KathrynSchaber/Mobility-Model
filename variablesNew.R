M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvae = M_larvaeNew
MM = array(0, dim = c(Nl, xi, T))

Sml = SmlNew
Smf = SmfNew
Eml = matrix(0, nrow = Nl, ncol = sigma + 1)
Emf = matrix(0, nrow = Nf, ncol = sigma + 1)
Iml = rep(0, Nl)
Imf = rep(0, Nf)

Sh = rep(1, Nh)
Eh = matrix(0, nrow = Nh, ncol = tau)
Ih = matrix(0, nrow = Nh, ncol = rho_max + 1)
Rh = rep(0, Nh)

