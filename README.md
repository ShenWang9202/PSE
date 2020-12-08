# PSE
The code is for Probabilistic State Estimation in WDN
# Main file
MonteCarloNew.m is the main file.
# Intro
This work is to explore how the uncertainty from demand, measurement noise, and pipe roughness coefficients propagate to the uncertainty of system states, i.e., heads and flows.

AnalycalDistribution_FOSM.m is using the FOSM method and is good for any network.

"over" in the filename considers an over-determined scenario
This file is only for the 3-node network and would give an error when trying the other networks.

The keyword "Pipe"  in the filename considers the uncertainty from pipe roughness coefficients.

All AnalycalDistribution_MC based files consider demand uncertainty.

AnalycalDistribution_MC_over_pipe.m is considered an over-determined scenario and pipe uncertainty demand uncertainty at the same time

AnalycalDistribution_Kxx.m means we are using the second version of our derived formula, which is much easier to understand and the no need to construct Vech(Kxx).

# Note
Note that for PES1 Network and BAK network, their unit system is based on LPS. So when you Linearize Pipes or Pumps, Headloss_pipe_R should be 10.66 instead of 4.727; Please be careful to convert unit, otherwise, the covariance result is wrong even no error or warnings pop out.

As for how to obtain (i) the probability density function of the linear combinations of independent random variables with different (known or unknown) distributions and (ii) its corresponding confidence interval numerically, please see PSE/tests/PDF_Sum_Random_Var.m file.
