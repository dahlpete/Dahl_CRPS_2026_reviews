import numpy as np

# construct matrix for time dependence

n_timepoints = 10
n_hqs = 13
n_ks = 10

J1 = np.zeros([n_timepoints,n_ks])
J2 = np.zeros([n_timepoints,n_ks])
J3 = np.zeros([n_hqs,n_ks])
J4 = np.zeros([n_hqs,n_ks])

# go column by column

param_sets = np.zeros([n_ks,n_ks+1])

for i in range(n_ks):
    kval = i+1
    file1 = 'r%s/params_temp0.txt' % kval
    fIN = open(file1,'r').readlines()
    values1 = np.array([float(line) for line in fIN])

    file2 = 'r0/params_temp0.txt'
    fIN2 = open(file2,'r').readlines()
    values2 = np.array([float(line) for line in fIN2])

    diffvals = values1 - values2

    dk = diffvals[i]

    if i == 0:
        param_sets[:,i] = values2
    param_sets[:,kval] = values1

    # nh3 data
    file1 = 'r%s/nh3_values.txt' % kval
    fIN = open(file1,'r').readlines()
    nh3_vals1 = np.array([float(line) for line in fIN])

    file2 = 'r0/nh3_values.txt'
    fIN2 = open(file2,'r').readlines()
    nh3_vals2 = np.array([float(line) for line in fIN2])

    dC = nh3_vals1 - nh3_vals2
    dCdk = dC / dk
    J1[:,i] = dCdk

    # h2 data
    file1 = 'r%s/h2_values.txt' % kval
    fIN = open(file1,'r').readlines()
    h2_vals1 = np.array([float(line) for line in fIN])

    file2 = 'r0/h2_values.txt'
    fIN2 = open(file2,'r').readlines()
    h2_vals2 = np.array([float(line) for line in fIN2])

    dC = h2_vals1 - h2_vals2
    dCdk = dC / dk
    J2[:,i] = dCdk

    # HQ data
    file1 = 'r%s/hq_dependence.txt' % kval
    fIN = open(file1,'r').readlines()
    hq_vals1 = np.array([np.array(line.split()).astype(float) for line in fIN])
    #print(hq_vals1)
    nh3_vals1 = hq_vals1[:,0]; h2_vals1 = hq_vals1[:,1]

    file2 = 'r0/hq_dependence.txt'
    fIN2 = open(file2,'r').readlines()
    hq_vals2 = np.array([np.array(line.split()).astype(float) for line in fIN2])
    nh3_vals2 = hq_vals2[:,0]; h2_vals2 = hq_vals2[:,1]

    dC = nh3_vals1 - nh3_vals2
    dCdk = dC / dk
    J3[:,i] = dCdk

    dC = h2_vals1 - h2_vals2
    dCdk = dC / dk
    J4[:,i] = dCdk

#print(J1)

# construct weights
# NH3
expt_var = np.array([0.00000239454728,1.40501498382368,0.345772627232,0.32427917084672,3.430928894,0.009808066,3.144797256,0.81388028,0.001618714,2.211414699])
inv_var = 1/expt_var

weight_matrix1 = np.zeros([len(inv_var),len(inv_var)])
for i in range(len(inv_var)):
    weight_matrix1[i,i] = inv_var[i]

# H2
expt_var = np.array([0.244340603,9.531239429,2.67587352,0.12661075,0.012698333,51.09500819,87.55198214,0.723748119,113.8953377,0.006408889])
inv_var = 1/expt_var

weight_matrix2 = np.zeros([len(inv_var),len(inv_var)])
for i in range(len(inv_var)):
    weight_matrix2[i,i] = inv_var[i]

# HQ - nh3
expt_var = np.array([0.263218657,0.043362769,0.008197121,0.002950963,0.308357469,1.486879766,0.914535075,0.102603956,2.34215698,0.224215829,0.208385477,0.005749308,1.08778536])
inv_var = 1/expt_var

weight_matrix3 = np.zeros([len(inv_var),len(inv_var)])
for i in range(len(inv_var)):
    weight_matrix3[i,i] = inv_var[i]

# HQ - h2
expt_var = np.array([0.269633802,0.044504725,0.662845509,1.562316316,0.055275301,0.176521567,0.008041149,1.85988389,0.148966154,0.436837679,0.90308516,1.445999403,8.927577751])
inv_var = 1/expt_var

weight_matrix4 = np.zeros([len(inv_var),len(inv_var)])
for i in range(len(inv_var)):
    weight_matrix4[i,i] = inv_var[i]


# compute FID
#F = np.linalg.matmul(np.linalg.matmul(np.transpose(J1),weight_matrix),J1)
F1 = J1.T @ weight_matrix1 @ J1
F2 = J2.T @ weight_matrix2 @ J2
F3 = J3.T @ weight_matrix3 @ J3
F4 = J4.T @ weight_matrix4 @ J4

Ftot = F1+F2+F3+F4



# print the fisher information matrix
np.savetxt('FIM.txt',Ftot)

np.savetxt('param_sets.txt',param_sets)



C_matrix = np.linalg.inv(Ftot)
np.savetxt('covar.txt',C_matrix)

stdev_vals = np.zeros(n_ks)
for i in range(n_ks):
    stdev_vals[i] = np.sqrt(C_matrix[i,i])


D = np.zeros([n_ks,n_ks])
for i in range(n_ks):
    D[i,i] = stdev_vals[i]
Dinv = np.linalg.inv(D)

R = Dinv @ C_matrix @ Dinv

np.savetxt('correlation.txt',R)

print(stdev_vals)
