This directory contains the files necessary to simulate the kinetics of ammonia and H2 formation with the NC1:MoFe protein biohybrids described in Figure 2. To simulate the kinetics of ammonia and H2 accumulation, as presented in Figure 2b, run

python time_dependence.py

This will simulate the reaction for 55 minutes, and will plot the simulated products along with the experimentally determined values.

To simulate the HQ concentration dependence, the file time_dependence.py will need to be modified. The variable ttot, which specifies the total simulation time, should be set to 7 minutes, and the variable Dconc, which specifies the HQ concentration in units of M, should be updated to the desired concentration. The simulation is run the same way as specified above. Because the total time is less than 10 minutes, the code will not produce a graph of the kinetics. If a graph is desired, this can be altered in the last line of the time_dependence.py script.

Similarly, to simulate the pN2 dependence the time_dependence.py script should be modified. The variable ttot should be set to 7 minutes, the variable Dconc should be set to 50 x 10-3 M, and the variable pN2 should be set to the desired N2 partial pressure in units of atm. 

To simulate the H2 inhibition, the variable ttot should be set to 15 min, the variable Dconc should be set to 50 x 10-3 M, the variable pN2 should be set to 0.5 atm, and the variable pH2 should be set to the desired H2 partial pressure in units of atm.
