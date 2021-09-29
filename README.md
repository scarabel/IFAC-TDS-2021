# IFAC-TDS-2021
 Material related to the talk "Bifurcation Analysis of Delay Equations Using Software for ODEs" in session TB3T1 of the IFAC Workshop on Time Delay Systems 2021 (September 30th, 2021).

 The repository contains:
 - PDF introductory slides for the presentation  "Bifurcation Analysis of Delay Equations Using Software for ODEs"
 - MATLAB codes for the analysis of the logistic Delay Differential Equation ("PS_logistic.m" and "MC_logistic.m") and for the Renewal Equation for cannibalism ("PS_renewal" and "MC_renewal"); the files named "PS_" contain the definition of the right-hand side of the pseudospectral discretization of each equation, while the files named "MC_" contain the command line instructions for the numerical bifurcation using the software package MatCont for MATLAB.

The package MatCont can be downloaded from https://sourceforge.net/projects/matcont/

The code "init.m" in the MatCont folder should be run before using MatCont, and the files "PS_" should be saved in an active folder in the MATLAB path. 

## Main references

- Breda, D., Diekmann, O., Gyllenberg, M., Scarabel, F., Vermiglio, R. (2016) Pseudospectral discretization of nonlinear delay equations: new prospects for numerical bifurcation analysis, SIAM Journal on Applied Dynamical Systems, 15(1), 1–23. https://doi.org/10.1137/15M1040931 
- Breda, D., Diekmann, O., Liessi, D., Scarabel, F. (2016) Numerical bifurcation analysis of a class of nonlinear renewal equations, Electronic Journal of Qualitative Theory of Differential Equations, 65, 1–24. https://doi.org/10.14232/ejqtde.2016.1.65
- Scarabel, F., Diekmann, O., Vermiglio, R. (2021) Numerical bifurcation analysis of renewal equations via pseudospectral approximation, Journal of Computational and Applied Mathematics, 397:113611. https://doi.org/10.1016/j.cam.2021.113611


