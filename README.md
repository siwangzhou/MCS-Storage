# MCS-Storage
## Test datasets
"T_normal.mat" and "T_normal600.mat" are ocean temperature datasets. The size of "T_normal.mat" is 16*16*300, and that of "T_normal600.mat" is 16*16*600. "PM_normal.mat" and "PM_normal600.mat" are PM2.5 concentration datasets, the sizes of which are 16*16*300 and 16*16*600, respectively.

## Figures
Seven .m files, using which the figures shown in Figs. 4-7 are plotted.

## Source codes
“mcs_storage_period.m” and “mcs_storage_period_withoutPDR.m” are used to do the recovery performance comparison for 12 sensing periods with and without employing PDR algorithm, respectively. 
“mcs_storage_proposed.m” and “mcs_storage_proposed_withoutPDR.m” are used to compute the relative square errors at various decoding ratios with and without employing PDR algorithm, respectively. 
“DAMPblock.m” and “denoiseblock.m” are used to implement our PDR algorithm.
“DAMP.m”, “denoise.m”, and the folder “D-AMP_toolbox-master” refers to the reference [16], and are used to obtain the experimental results without the technique of PDR.

## Trajectory
“trajectory_3d.m” is used to create trajectory matrixes of the participants. “Phi_3d_300up_30000.mat” is the created trajectory matrix with more than 300 steps. The other two .mat files are the trajectory matrixes corresponding to 400 and 500 steps, respectively.
