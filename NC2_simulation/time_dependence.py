import numpy as np
import photochem_turnover as hv
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import product_formation as pf
from scipy import interpolate as interp
import time


fIN = open('fit_params.txt','r').readlines()
params = [float(line) for line in fIN];
[kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0] = params

Dconc = 50.0e-3
#Aconc = 0e-6

kET2 = kET
kBET2 = kBET

ratevals = [['D',Dconc],['kET',kET],['kHT',kHT],['kAET',kAET],['kBET',kBET]]
secvals = [['kET',kET2],['kBET',kBET2]]

ttot = 55 # min
pN2 = 1.00
pH2 = 0.00



start = time.time()

print('Starting Integration...\n')

[time_array,t2,vcat_array,vox_array,H2_array,NH3_array,D_array,A_array,D2_array,A2_array,population]=pf.product_formation(khp2,khp3,khp4,kre,koa,ttot*60,readrates=ratevals,secondary=secvals,paramfile='photochem_params.txt',N2=pN2,H2atm=pH2,kto=k0)

end = time.time(); elapsed=(end-start)/60
print('Integration complete: time elapsed = %.3f min' % elapsed)
print('\nPlotting kinetics over %.2f min' % ttot)


## Print results
population=population.transpose()

fOUT = open('product_formation_results.txt','w')
print('time(s)\tH2(nmol)\tNH3(nmol)\t[D](mM)\t[A](mM)\t[D2](mM)\t[A2](mM)\tInact\tMox\tE0\tE1\tE2\tE3\tE44H\tE42n2h\tE5\tE6\tE7\tE8',file=fOUT)
for i in range(len(time_array)):
	print('%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f' % (time_array[i],H2_array[i],NH3_array[i],D_array[i],A_array[i],D2_array[i],A2_array[i],population[i,0],population[i,1],population[i,2],population[i,3],population[i,4],population[i,5],population[i,6],population[i,7],population[i,8],population[i,9],population[i,10],population[i,11]),file=fOUT)

fOUT.close()


fOUT = open('charge_separation_rates.txt','w')
for i in range(len(t2)):
	print('%.5f\t%.5f\t%.5f' % (t2[i],vcat_array[i],vox_array[i]),file=fOUT)
fOUT.close()




print('NH3 = %.3f nmol' % NH3_array[-1])
print('H2 = %.3f nmol' % H2_array[-1])


# READ EXPT DATA
fIN = open('expt_data_NC2.txt','r').readlines()
data = [[] for i in range(len(fIN))]
i = 0
for line in fIN:
	values = line.split()
	data[i] = [float(values[0]),float(values[1]),float(values[2])]
	i+=1
data = np.array(data).transpose()

#RUN INTERPOLATION
cs = interp.CubicSpline(time_array/60,H2_array)
H2_model = cs(data[0,:])

cs = interp.CubicSpline(time_array/60,NH3_array)
NH3_model = cs(data[0,:])

diff1 = data[1,:] - H2_model; diff2 = data[2,:] - NH3_model
diffs = np.array(list(diff1)+list(diff2))

rmsd = np.sqrt(np.mean(diffs**2))
print('RMSD = %.3f nmol' % rmsd)


fig=plt.figure()
ax=fig.gca()
ax.scatter(data[0,:],data[1,:],marker='o',s=105,color='g',edgecolors='k',linewidths=2)
ax.scatter(data[0,:],data[2,:],marker='o',s=105,color='c',edgecolors='k',linewidths=2)
ax.plot(time_array/60,NH3_array,linewidth=3,c='c')
ax.plot(time_array/60,H2_array,linewidth=3,c='g')

if ttot >= 10:
	plt.show()
