


dr_dm = 1./(4*pi*r[i]**2*ro)
dP_dm = - G*M[i]/(4*pi*r[i]**4)
dL_dm = E
dT_dm = -3.*k_SI*L[i]/(256.*pi**2*sb*r[i]**4*T[i]**3)



## Fra variablesteplength.pdf



f1 = dr_dm
dm1 = abs(r[i]*0.01/f1)
dm_list.append(dm1)

f2 = dP_dm
dm2 = abs(P[i]*0.01/f2)
dm_list.append(dm2)

f3 = dL_dm
dm3 = abs(L[i]*0.01/f3)
dm_list.append(dm3)

f4 = dT_dm
dm4 = abs(T[i]*0.01/f4)
dm_list.append(dm4)

dm = -min(dm_list)

