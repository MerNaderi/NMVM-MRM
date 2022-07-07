# NMVM-MRM
Robust mixture regression modeling based on the normal mean-variance mixture distributions

Please copy the files to the "current working directory" of the R package.

The getwd() function shall determine an absolute pathname of the "current working directory". 


./Functions

	contains 

        (1) GHST-mix.r: EM-based estimating program for the GHST-MRM;

        (2) VG-mix.r: EM-based estimating program for the VG-MRM;

        (3) NIG-mix.r: EM-based estimating program for the NIG-MRM;

        (4) SL-mix.r: EM-based estimating program for the SL-MRM;

        (5) NMVBS-mix.r: EM-based estimating program for the NMVBS-MRM;

        (6) NMVL-mix.r: EM-based estimating program for the NMVL-MRM;

        (7) normal-mix.r: EM-based estimating program for the normal-MRM;

        (8) Additional.r: function required for random data generation and EM implementation;

        (9) SMSN-fiting.r: EM-based estimating program for the SMSN-MRMs (including ST-MRM, SSL-MRM and SCN-MRM);


./Code

       contains

        (1) Simulation 1 - NMVBS.r: main script for conducting simulation 1 for the NMVBS-MRM;

        (2) Simulation 1 - SL.r: main script for conducting simulation 1 for the SL-MRM; (Results do not reported in the main text);

        (3) Simulation 2-- Parallel.r: main script for conducting simulation 2 for the parallel configuration;

        (4) Simulation 3-- ESN---Parallel.r: main script for conducting simulation 3 for the parallel configuration with 50 noise added;

        (5) Simulation 4.r: main script for conducting simulation 4;

        (6) Fig 8.r: main script for reproducing the heat-map compairsion plots;

        (7) Tone.r: code for fitting the special cases to the Tone perception dataset;

        (8) Tone-sensivity.r: code for creating Figure 2 which shows the LD measures for the Tone perception data fitted by four selected NMVM-MRMs;

        (9) Fig 7.r: main script for reproducing the Figure 7;

Notes: 
    
      1. The "sn", "ghyp", "aricode", "fpc", "ClusterR", "VGAM", and "GIGrvg" R packages are required;

      2. The "Tone" data set is avalible in the R package "fpc";

      3. The authors do not have permission to distribute "Pinus nigra" data. Request should be submitted to the correcponding author of Garcia-Escudero et al. (2010);
