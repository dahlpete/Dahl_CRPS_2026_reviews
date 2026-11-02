import numpy as np
import photochem_turnover_test3 as hv
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import product_formation6 as pf
from scipy import interpolate as interp
import time


fIN = open('params_temp0.txt','r').readlines()
params = [float(line) for line in fIN];
[kET,kHT2,kAET,kBET,khp2,khp3,khp4,kre,koa,k0] = params


# Run time dependence

Dconc = 50.0e-3
#Aconc = 0e-6

kET2 = kET
kBET2 = kBET

ratevals = [['D',Dconc],['kET',kET],['kHT2',kHT2],['kAET',kAET],['kBET',kBET]]
secvals = [['kET',kET2],['kBET',kBET2]]

ttot =  55.5# min
pN2 = 1.00
pH2 = 0.00


start = time.time()

print('Starting Integration...\n')

[time_array,t2,vcat_array,vox_array,H2_array,NH3_array,D_array,A_array,D2_array,A2_array,population]=pf.product_formation(khp2,khp3,khp4,kre,koa,ttot*60,readrates=ratevals,secondary=secvals,paramfile='photochem_params_test2.txt',N2=pN2,H2atm=pH2,kto=k0)

end = time.time(); elapsed=(end-start)/60
print('Integration complete: time elapsed = %.3f min' % elapsed)

time_array = time_array / 60; # convert to minutes


# extract the timepoints that align with experiment
expt_times = np.array([2,5,8,12,18,25,32,38,45,55]);
#expt_times = np.array([2,5,8])
ixarray = np.zeros(len(expt_times))
i = 0
for et in expt_times:
	tarray = np.abs(time_array - et)
	minval = np.min(tarray); tarray=list(tarray)
	ixval = tarray.index(minval)
	ixarray[i] = ixval
	i += 1
ixarray = ixarray.astype(int)

nh3_vals = [NH3_array[i] for i in ixarray]
h2_vals = [H2_array[i] for i in ixarray]

# print time dependent values
fOUT = open('nh3_values.txt','w')
for i in range(len(nh3_vals)):
	print('%.5E' % nh3_vals[i],file=fOUT)
fOUT.close()

fOUT = open('h2_values.txt','w')
for i in range(len(h2_vals)):
	print('%.5E' % h2_vals[i],file=fOUT)
fOUT.close()


# run HQ dependence
hq_vals = np.array([0,1,2,5,7,10,20,30,40,50,65,80,100])

evals = np.zeros(len(hq_vals))
nh3vals = np.zeros(len(hq_vals))
h2vals = np.zeros(len(hq_vals))
for i in range(len(evals)):
	Dconc = hq_vals[i] * 1e-3

	ratevals = [['D',Dconc],['kET',kET],['kHT2',kHT2],['kAET',kAET],['kBET',kBET]]
	secvals = [['kET',kET2],['kBET',kBET2]]

	ttot = 7 # min
	pN2 = 1.00
	pH2 = 0.00

	[time_array,t2,vcat_array,vox_array,H2_array,NH3_array,D_array,A_array,D2_array,A2_array,population]=pf.product_formation(khp2,khp3,khp4,kre,koa,ttot*60,readrates=ratevals,secondary=secvals,paramfile='photochem_params_test2.txt',N2=pN2,H2atm=pH2,kto=k0)

	evals[i] = 3*NH3_array[-1] + 2*H2_array[-1]; nh3vals[i] = NH3_array[-1]; h2vals[i] = H2_array[-1]

fOUT = open('hq_dependence.txt','w')
for i in range(len(evals)):
	print('%.5E\t%.5E' % (nh3vals[i],h2vals[i]),file=fOUT)
fOUT.close()
	

