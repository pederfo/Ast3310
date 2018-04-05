from numpy import *
ro = 1.62E5
T = 1.57E7
m_u = 1.6605E-27
N_A = 6.02214129E23

X = 0.7 # Hydrogen
Y3 = 1E-10 # Helium 3
Y4 = 0.29 # Helium 4
Z = 0.01 # Other metals
Z_73Li = 1E-7 # Lithium 7
Z_74Be = 1E-7 # Beryllium 7


#Computes number of atoms of the different elements.
def atomnumbers(rho):
	n_e   = rho*(1+X) /(2.*m_u)    	#Number of electrons
	n_d   = 0			#Number of deuterium atoms
	n_p   = rho*X	  /m_u     	#Number of Hydrogen atoms
	n_3   = rho*Y3	  /(3.*m_u)     #Number of Helium 3 atoms
	n_4   = rho*Y4	  /(4.*m_u)     #Number of Helium 4 atoms
	n_Li  = rho*Z_73Li/(7.*m_u) 	#Numver of Beryllium 7 atoms
	n_Be  = rho*Z_74Be/(7.*m_u) 	#Number of Lithium 7 atoms
	atomnumbers = [n_e,n_d,n_p,n_3,n_4,n_Li,n_Be]
	return atomnumbers

#Reactions
"""In this part the reaction rates of the processes in the PPI and PPII \
chain are calculated. I neglect the other processes since I choose to \
look at a star with low temperature at it's core, and in that case those \
processes are very slow compared to PPI and PPII."""

def reactions(T):
	T9 		= T*1E-9		#Convert T to form used in reactions.
	T9_star1 	= T9/(1+4.95E-2*T9)	#Other form used.
	T9_star2 	= T9/(1+.759*T9)	#Yet another form used.

	H_D   	= 	4.01E-15*T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1 +\
                     .123*T9**(1./3) + 1.09*T9**(2./3) + .938*T9)

	He3_p	= 	6.04E10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 +\
                     .034*T9**(1./3) - .522*T9**(2./3)-.124*T9 + \
                     .353*T9**(4./3)+ .213*T9**(-5./3))

	He3_Be  = 	5.61E6*T9_star1**(5./6)*T9**(-3./2)*exp(-12.826*\
                     T9_star1**(-1./3))

	Be_Li   = 	1.34E-10*T9**(-1./2)*(1 - .537*T9**(1./3)  + \
                     3.86*T9**(2./3) +  .0027*T9**-1*exp(2.515E-3*T9**-1))

	Li_p	= 	1.096E9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - \
                     4.830E8*T9_star2**(5./6)*T9**(-3./2)*exp(-8.472*\
                     T9_star2**(-1./3)) + 1.06E10*T9**(-3./2)\
                     *exp(-30.442*T9**-1)
	"""
	"These two reactions are not used in this assignment, but may be useful \
     at a later point, so I've kept them here."
     
	Berylliym_Boron = 	3.11E5*T9**(-2./3)*exp(-10.262*T9**(-1./3))+2.53E3\
                         *T9**(-3./2)*exp(-7.306*T9**-1)

	Nitrogen_Oxygen = 	4.9E7*T9**(-2./3)*exp(-15.228*T9**(-1./3)-\
                         .092*T9**2)*(1+.027*T9**(1./3)-.778*T9**(2./3)\
                         -.149*T9+.261*T9**(4./3)+.127*T9**(5./3)) + \
                         2.37E3*T9**(2./3)*exp(-3.011*T9**-1) + \
                         2.19E4*exp(-12.53*T9**-1)
	"""
	rrs = [H_D,He3_p,He3_Be,Be_Li,Li_p]#,Berylliym_Boron,Nitrogen_Oxygen
	return rrs

#Time evolution
"""
This segment computes the reactionrates for the different processes in the sun.
There is no time evolution here 
"""
def Reactionrates(T,rho):
	n = atomnumbers(rho)
	rrs = reactions(T)
	#Reaction rates
	"""All of the following rates are calculated by the equation
	r_ik = n_i*n_k/(rho*(1+delta_ik)*lambda_ik)
	where i,k define the elements and delta_ik =1 if i=k else 0.
	"""
	lambda_pp = rrs[0]/N_A*1E-6		#[m^3/s]
	#lambda_pd = 1				#This one is actually unknown so it's just \
						#set to 1 since the reaction it is used in 
						#includes the number of deuterium atoms and \
						#that number is set to 0 earlier.
	lambda_33 = rrs[1]/N_A*1E-6		#[m^3/s]
	lambda_34 = rrs[2]/N_A*1E-6		#[m^3/s]
	lambda_7e = rrs[3]/N_A*1E-6		#[m^3/s]
	lambda_71prime  = rrs[4]/N_A*1E-6	#[m^3/s]
	#lambda_71  = Berylliym_Boron/N_A*1E-6	#[m^3/s]#Not used in this assigment
	#lambda_p14 = Nitrogen_Oxygen/N_A*1E-6	#[m^3/s]#Not used in this assigment

	r_pp 	  = n[2]**2/(rho*2)*lambda_pp	#[kg-1*s-1]
	r_pd  	  = r_pp			#[kg-1*s-1] Assume this reaction happens \
						#instantly such that it happens at the same \
						#rate as the elements it needs become available.
						#Thus it must be the same as the reaction that \
						#creates those elements (in this case Deuterium)
	r_33 	  = n[3]**2  /(rho*2)*lambda_33	#[kg-1*s-1]
	r_34 	  = n[3]*n[4]/rho*lambda_34	#[kg-1*s-1]
	r_7e 	  = n[5]*n[6]/rho*lambda_7e	#[kg-1*s-1]
	r_71prime = n[5]*n[2]/rho*lambda_71prime#[kg-1*s-1]
	#r_71  	 = n_Li*n_p	/rho*lambda_71	#Not used in this assigment
	#r_p14 	 = n_p*n_4	/rho*lambda_p14 #Not used in this assigment
		
	#This part makes sure that reaction rates which rely on other rates \
	#don't use elements that are not yet present.
	if 2.*r_33 + r_34 > r_pd:
		r_33 = 2./3*r_pd
		r_34 = 1./3*r_pd
	if r_7e > r_34:
		r_7e = r_34
	if r_71prime > r_7e:
		r_71prime = r_7e
	reactionrate = r_pp,r_pd,r_33,r_34,r_7e,r_71prime

	return reactionrate

#Solving the equation for the change in Luminosity

def Energy(T,rho):
	
	Q1 = 0.15E6*1.602176565E-19 
	Q2 = 5.49E6*1.602176565E-19 
	Q3 = 12.86E6*1.602176565E-19
	Q4 = 1.59E6*1.602176565E-19
	Q7e = 0.05E6*1.602176565E-19
	Q71 = 17.35E6*1.602176565E-19

	reactionrate = Reactionrates(T,rho)
	e_1  = reactionrate[0]*(Q1+Q2)
	e_2  = reactionrate[2]*Q3
	e_3  = reactionrate[3]*Q4
	e_4  = reactionrate[4]*Q7e
	e_5  = reactionrate[5]*Q71
	e    = e_1+e_2+e_3+e_4+e_5
	return e_1, e_2, e_3, e_4, e_5

[e_1, e_2, e_3, e_4, e_5] = Energy(T,ro)
print e_1, e_2, e_3, e_4, e_5
