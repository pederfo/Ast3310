from numpy import *

ro = 1.62E5
T = 1.57E7
mu = 1.6605E-27
Na = 6.02214129E23

X = 0.7 # Hydrogen
Y3 = 1E-10 # Helium 3
Y = 0.29 # Helium 4
Z = 0.01 # Other metals
Z_Li = 1E-7 # Lithium 7
Z_Be = 1E-7 # Beryllium 7


def lamdaNa(T):
	T9 = T*1E-9		
	T9_star1 = T9/(1+4.95E-2*T9)	
	T9_star2 = T9/(1+.759*T9)

	#2xH
	l_1   = 4.01E-15*T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1 +\
                     .123*T9**(1./3) + 1.09*T9**(2./3) + .938*T9)

	#2x3He
	l_2	= 6.04E10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 +\
                     .034*T9**(1./3) - .522*T9**(2./3)-.124*T9 + \
                     .353*T9**(4./3)+ .213*T9**(-5./3))

	#3He Be
	l_3  = 5.61E6*T9_star1**(5./6)*T9**(-3./2)*exp(-12.826*\
                     T9_star1**(-1./3))

	#Be Li
	l_4   = 1.34E-10*T9**(-1./2)*(1 - .537*T9**(1./3)  + \
                     3.86*T9**(2./3) +  .0027*T9**-1*exp(2.515E-3*T9**-1))

	#Li H
	l_5	= 1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - \
                     4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*\
                     T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)\
                     *exp(-30.442*T9**-1)

	return l_1, l_2, l_3, l_4, l_5

def efunk(ro, T):
	# Densities of atomic species
	n_p = ro*X/(1.*mu)
	n_He = ro*Y/(4.*mu)
	n_e = (ro/2.*mu)*(1. + X)
	n_d = ro*X/(2.*mu)
	n_3He = ro*Y3/(3.*mu)
	n_Li = ro*Z_Li/(7.*mu)
	n_Be = ro*Z_Be/(7.*mu)

	
	[l_1, l_2, l_3, l_4, l_5] = lamdaNa(T)

	#complete lambda functions
	lm_1 = l_1/Na*1E-6
	lm_2 = l_2/Na*1E-6
	lm_3 = l_3/Na*1E-6
	lm_4 = l_4/Na*1E-6
	lm_5 = l_5/Na*1E-6


	#reaktion rates
	r_pp = n_p**2/(ro*2)*lm_1
	r_pd = r_pp
	r_33 = n_3He**2  /(ro*2)*lm_2
	r_34 = n_3He*n_He/ro*lm_3
	r_e7 = n_Li*n_Be/ro*lm_4
	r_71 = n_Li*n_p/ro*lm_5

	#Check and correct for reactions that overstep their fuel
	if 2.*r_33 + r_34 > r_pd:
		r_33 = 2./3*r_pd
		r_34 = 1./3*r_pd
	if r_e7 > r_34:
		r_7e = r_34
	if r_71 > r_e7:
		r_71 = r_7e

	
	# free energy
	Q1 = 0.15E6*1.602176565E-19 
	Q2 = 5.49E6*1.602176565E-19 
	Q3 = 12.86E6*1.602176565E-19
	Q4 = 1.59E6*1.602176565E-19
	Q7e = 0.05E6*1.602176565E-19
	Q71 = 17.35E6*1.602176565E-19

	# full energy
	e_1 = r_pp*(Q1+Q2)
	e_2 = r_33*Q3
	e_3 = r_34*Q4
	e_4 = r_e7*Q7e
	e_5 = r_71*Q71
	
	
	E = e_1+e_2+e_3+e_4+e_5
	return e_1, e_2, e_3, e_4, e_5


[e_1, e_2, e_3, e_4, e_5] = efunk(ro, T)
print e_1, e_2, e_3, e_4, e_5
