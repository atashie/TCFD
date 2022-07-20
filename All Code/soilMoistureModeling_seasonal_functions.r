
infiltration_f = function(PPT = NA, crn = 50, Smax = 100, Ia_scalar = 0.2){
	Ia = Ia_scalar * Smax
	infil = PPT - ifelse(PPT > Ia, (PPT - Ia)^2 / (PPT - Ia + Smax), 0)
	return(infil)
	}
	
PET_f = function(	)	{
	
	}

	
SM_routine_f = function(infiltration = NA, PET = NA, 
	rcScl = 0.1,	# water retention curve scalar
	rcExp = 1.3,	# water retention curve exponent
	PETexp = 2, 	# exponent of PET decay
	Zr = 1000,	# root depth
	n = 0.5,	# soil porosity
	smhp = 0.00,	# soil moisture at hydroscopic point
	smwp = 0.10,	# soil moisture at wilting point
	smfc = 0.25,	# soil moisture field capacity
	s0 = .5)	# initial soil moisture 
	{
	nt = length(infiltration)
	SM_store = c(s0, rep(NA, nt))
	smhp_stor = Zr * smhp
	smwp_stor = Zr * smwp
	smfc_stor = Zr * smfc
	max_stor = Zr * n
	
	for(i in 1:length(infiltration))	{
		thisSMstor = SM_store[i]
		AET = ifelse(thisSMstor > smhp_stor, 
			min(thisSMstor - smhp_stor, PET[i] * (thisSMstor / max_stor) ^ PETexp),
			0)
		thisSMstor = thisSMstor - AET
		
		deepDrainage = ifelse(thisSMstor > smfc_stor,
			min(thisSMstor - smfc_stor, rcScl * thisSMstor * (thisSMstor / smfc_stor)^rcExp),
			0)
		
		SM_store[i+1] = min(max_stor, thisSMstor - min(c(thisSMstor, deepDrainage)) + infiltration[i]) 
		}
	return(SM_store / n)
	}
	
	
SM_routine_f(infiltration = infil, PET = PET)
	

rain = 90 * RainPoisson(ndays = 365, lambda = 0.05, alpha = 0.60)








sm_out_1 = swb_f(R = rain,
	Rstar = 3,	# max amount that the cannopy intercepts
	Emax = seq(4,5,length.out=365),	# max ET rate
	Ew = 0.5,	# min ET rate
	Ks = 2000,	# Ksat
	b = 4.38,	# exponent of water retention curve
	Zr = 400,	# root depth
	n = 0.5,	# soil porosity
	sh = 0.01,	# soil moisture at hydroscopic point
	sw = 0.10,	# soil moisture at wilting point
	sstar = 0.25,	# soil moisture field capacity
	s0 = 0.10,	# initial soil moisture (else S0 = sh)
	nsteps = 1,	# number of steps/division for the numerical colution #???
	gr = TRUE)	# logical to show graphics of results
