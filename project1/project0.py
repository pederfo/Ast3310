from numpy import *

# Mass fractions of atomic species

X = 0.7 # Hydrogen
Y3 = 10**(-10) # Helium 3
Y = 0.29 # Helium 4
Z = 0.01 # Other metals
Z_Li = 10**(-7) # Lithium 7
Z_Be = 10**(-7) # Beryllium 7

# Other constants
ro = 1.62*10**5
mu = 1.660538921*10**(-27) # Unit mass
T = 1.5*10**7 # Core temperature

# Densities of atomic species
n_p = ro*X/(1*mu)
n_He = ro*Y/(4*mu)
n_e = (ro/2*mu)*(1 + X)

n_d = ro*X/(2*mu)
n_3He = ro*Y3/(3*mu)
n_Li = ro*Z_Li/(7*mu)
n_Be = ro*Z_Be/(7*mu)
n_8Be = ro*Z/(8*mu)
n_12C = ro*Z/(12*mu)
n_13C = ro*Z/(13*mu)
n_14N = ro*Z/(14*mu)
n_15N = ro*Z/(15*mu)

## Avogadros konstant

Na = 6.02214129*10**23 # 1/mol

def reaction(lamda, n1, n2):
	r = lamda*n1*n2
	return r

def supertrooper(ro, T):
	## Units of temperature
	T9 = T*1E-9
	T_mongo = T9/(1+4.95*10**(-2)*T9)
	T_mongo2 = T9/(1+0.759*T9)

	## Reaction rates
	#2x H
	lam_pp = (4.01E-15*T9**(-2./3)*exp(-3.380/(T9**(1./3)))*(1 + 0.123*T9**(1./3) + 1.09*T9**(2./3) + 0.938*T9))/Na*1E-6

	r_pp = lam_pp*n_p**2
	Q1 = 0.15E6*1.602176565E-19 # Free energy [Joule]
	Q2 = 5.49E6*1.602176565E-19 # (From Deuterium reaction) 

	#2x 3He
	lam_33 = (6.04E10*T9**(-2./3)*exp(-12.276/(T9**(1./3)))*(1 + 0.034*T9**(1./3) - 0.522*T9**(2./3) - 0.124*T9 + 0.353*T9**(4./3) + 0.213*T9**(-5./3)))/Na*1E-6

	r_33 = lam_33*n_3He**2
	Q3 = 12.86E6*1.602176565E-19 # Free energy [Joule]

	#3He 7Be
	lam_34 = (5.61E6*T_mongo**(5./6)*T9**(-3./2)*exp(-12.826*T_mongo**(-1./3)))/Na*1E-6
	
	r_34 = lam_34*n_3He*n_He
	Q4 = 1.59E6*1.602176565E-19 # Free energy [Joule]

	#7Be 7Li
	#lam_e7 = (1.34E-10*T9**(-1./2)*(1 - 0.537*T9**(1./3) + 3.86*T9**(2./3) + 0.0027*T9**(-1.)*exp(2.515E-3*T9**(-1.))))/Na*1E-6
	lam_e7 = 1.34*10**(-10)*T9**(-0.5)*(1.0-0.537*T9**(1/3.)+3.86*T9**(2/3.)+0.0027*T9**(-1.)*exp(2.515*10**(-3)*T9**(-1.)))/Na*1E-6
	r_7e = lam_e7*n_Be*n_e
	Q7e = 0.05E6*1.602176565E-19 # Free energy [Joule]

	#7Li
	#lamd_17 = (1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) -4.830E8*T_mongo2**(5./6)*T9**(-2./3)*exp(-8.472**T_mongo2**(-1./3)) + 1.06E10*T9**(-2./3)*exp(-30.442*T9**(-1)))/Na*1E-6
	lamd_17 = (1.096*10**9*T9**(-2./3)*exp(-8.472*T9**(-1./3))-4.830*10**8*T_mongo2**(5./6)*T9**(-3./2)*exp(-8.472*T_mongo2**(-1./3))+1.06*10**(10)*T9**(-3./2)*exp(-30.442*T9**(-1.)))/Na*1E-6
	r_71 = lamd_17*n_Li*n_p
	Q71 = 17.35E6*1.602176565E-19 # Free energy [Joule]	

	#7Be 8Be
	lam_17 = (3.11E5*T9**(-2./3)*exp(-10.262*T9**(-1./3)) + 2.53E3*T9**(-3./2)*exp(-7.306*T9**(-1)))/Na*1E-6

	#14N 15O
	lam_p14 = (4.90E7*T9**(-2./3)*exp(-15.228*T9**(-1./3) - 0.092*T9**2)*(1 + 0.027*T9**(1./3) - 0.778*T9**(2./3) - 0.149*T9 + 0.261*T9**(4./3) + 0.127*T9**(5./3)) + 2.37E3*T9**(2./3)*exp(-3.011*T9**(-1)) + 2.19E4*exp(-12.53*T9**(-1)))/Na*1E-6

	if r_7e > r_34:
		r_7e = r_34

	if r_71 > r_7e:
		r_17 = r_7e

	### Energy produced

	Epp = r_pp*(Q1 + Q2) #PP
	E33 = r_33*Q3 #PPI
	E34 = r_34*Q4	#PPII
	E7e = r_7e*Q7e  #PPII
	E71 = r_71*Q71  #PPII

	return Epp, E33, E34, E7e, E71


[Epp, E33, E34, E7e, E71] = supertrooper(ro, T)

print Epp, E33, E34, E7e, E71
