import numpy as np
import photochem_turnover as hv
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy import interpolate as interp
import product_formation as pf

fIN = open('expt_data_NC2.txt','r').readlines()
data = [[] for i in range(len(fIN))]
i = 0
for line in fIN:
        values = line.split()
        data[i] = [float(values[0]),float(values[1]),float(values[2])]
        i+=1
data = np.array(data).transpose()

hqIN = open('hq_dependence.txt','r').readlines()
hq_data = [[] for i in range(len(hqIN))]
i = 0
for line in fIN:
        values = line.split()
        hq_data[i] = [float(values[0]),float(values[1])]
        i+=1
hq_data = np.array(hq_data).transpose()

ttot = 55 # min

def rmsd(x):
        global kET, kET2, kHT, kAET, kBET, kBET2, khp2, kre, khp3, khp4, koa, k0
        [kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0] = list(x);

        Dconc = 50e-3
        ratevals = [['D',Dconc],['kET',kET],['kHT',kHT],['kAET',kAET],['kBET',kBET]]
        secvals = [['kET',kET2],['kBET',kBET2]]

        # run first using the parameters for RT product assays
        [time_array,t2,vcat_array,vox_array,H2_array,NH3_array,D_array,A_array,D2_array,A2_array,population]=pf.product_formation(khp2,khp3,khp4,kre,koa,ttot*60,readrates=ratevals,secondary=secvals,paramfile='photochem_params.txt',N2=1.0,H2atm=0.0,kto=k0)

        # interpolate H2 and NH3 data
        cs = interp.CubicSpline(time_array/60,H2_array)
        H2_model = cs(data[0,:])

        cs = interp.CubicSpline(time_array/60,NH3_array)
        NH3_model = cs(data[0,:])

        diff1 = data[1,:] - H2_model; diff2 = data[2,:] - NH3_model

        # run HQ concentration dependence
        hq_vals = hq_data[0,:]
        #hq_vals = np.array([0,1,2,5,7,10,20,30,40,50,65,80,100])
        
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
        
                [time_array,t2,vcat_array,vox_array,H2_array,NH3_array,D_array,A_array,D2_array,A2_array,population]=pf.product_formation(khp2,khp3,khp4,kre,koa,ttot*60,readrates=ratevals,secondary=secvals,paramfile='photochem_params.txt',N2=pN2,H2atm=pH2,kto=k0)
        
                evals[i] = 3*NH3_array[-1] + 2*H2_array[-1]; nh3vals[i] = NH3_array[-1]; h2vals[i] = H2_array[-1]

        diff3 = hq_data[1,:] - e_vals


        diffs = np.array(list(diff1)+list(diff2)+list(diff3))

        rmsd = np.sqrt(np.mean(diffs**2))

        fOUT = open('params_temp.txt','w')
        for l in [kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0]:
                print('%.8e' % l,file=fOUT)
        fOUT.close()

        print('%.8e' % rmsd,file=fOUT2)
        print(rmsd)
        return rmsd



param_read = False
if param_read:
        fIN = open('params_temp.txt','r').readlines()
        params = [float(line) for line in fIN];
        [kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0] = params
else:
        [kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0] = [5.0e+04,3.5e+06,4.5e+04,8.0e+05,2.0e-01,2.0e-01,2.0e-01,2.0e-01,4.0e+00,2.5e-02]

kET2 = kET
kBET2 = kBET

fOUT2 = open('rmsd.txt','w')

x0 = [kET,kHT,kAET,kBET,khp2,khp3,khp4,kre,koa,k0]
bnds = [(1e3,1e8),(1e3,1e8),(1e3,1e8),(1e3,1e8),(0.0001,2.0),(0.0001,2.0),(0.0001,2.0),(0.0001,2.0),(0.0001,10),(0.0001,2.0)]


mthd = 'nelder-mead'
res = minimize(rmsd,x0,method=mthd,bounds=bnds,options={'adaptive':True,'maxiter':5000,'disp':True})

results = res.x;

print(results)
fOUT = open('params.txt','w')
for m in range(len(results)):
        print('%.4e' % (results[m]),file=fOUT)
fOUT.close()


