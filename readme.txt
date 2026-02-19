Programs and datasets for the paper "Time-Varying Risk Premium in Large Cross-Sectional Equity Datasets",
by P. Gagliardini, E. Ossola and O. Scaillet.

(1) Monthly data used for Section 4 from Jul 64 to Dec 09:
Wspace_25FF.mat    25 Fama-French portfolios 
Wspace_44Indu.mat  44 Indu. portfolios 
Wspace_CRSPCMST_ret.mat Individual stocks matching CRSP and Compustat 
Wspace_Fact.mat  factors: mkt,smb,hml,mom 
RiskFree.mat  free-risk rate: monthly 30-day T-bill 
Wspace_Instrum.mat  instruments: default spread, term spread

(2) MATLAB codes for the empirical exercises in Section 4.
GOS_main.m is the main code that allows you to perform the exercises in Section 4.
You need to insert the INPUT (dataset, number of factors, model) at the beginning of the code.
Each MATLAB code is annotated and contains a step-by-step description of the successive computations. 
We use EstTimeInvariantMod.m to estimate time-invariant models and get Tables 1 and 3.
We use EstTimeVariantMod.m to estimate time-varying models and get Tables 2 and 3, and Figures 1-2, and Figures in Appendices 9-11. We also get results for the time-variation tests shown in Section 4.3.

