# RhoActinRD (Rho and actin reaction-diffusion)

This repository contains the Matlab code for the publications
1. "Actin network heterogeneity tunes activator-inhibitor dynamics at the cell cortex," by O. Maxian, A. Dinner, and E. Munro, [PNAS](https://www.pnas.org/doi/abs/10.1073/pnas.2520485122), 2025. 
2. "Learning activator-inhibitor dynamics at the cell cortex with neural likelihood ratio estimation," by O. Maxian and A. Dinner. 

The main code is divided into two folders, corresponding to the "forward model" (RhoAndActin), 
and the neural-network based inference tools (SimBasedInf). 
For the forward model, the first file to run 
can be found in [MakeMoviesBestFit](https://github.com/omaxian/RhoActinRD/blob/master/RhoAndActin/MakeMoviesBestFit.m). 
This will generate Fig. 2 in the PNAS paper (all of the simulations best fit to the experimental data), and more generally is the interface for the discrete model coupling actin filaments to continuum Rho dynamics. Codes for continuum model (Fig. 5 in the text) can be found [here](https://github.com/omaxian/RhoActinRD/blob/master/RhoAndActin/ContinuumModels/RhoAndActinTauPDEs.m)

The main driver file for the hybrid model is [RhoAndActinBasalNuc](https://github.com/omaxian/RhoActinRD/blob/master/RhoAndActin/RhoAndActinBasalNuc.m). 

For simulation-based inference, the main interface is in [InferFromData](https://github.com/omaxian/RhoActinRD/blob/master/SimBasedInf/InferFromData.m), which gives you the likelihood functions assuming the classifier is already trained. 