This project provides the code corresponding to the reconstruction method 
presented in the paper 

Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
      Expectation-Maximization based approach to 3D reconstruction from single-
      waveform multispectral Lidar data. IEEE Transactions on Computational 
      Imaging, 6, pp.1033-1043.

A demo is available in 'main_synthetic_multi.m'.
This demo uses pregenerated data. However, two codes are also provided to generate
synthetic data:

   - 'gener_data4git.m' is the code used to generate the ground truth depth and 
	reflectivity profiles

   - 'initdata.m' generates histograms of photon counts from the ground truth 
	generated in 'gener_data4git.m'.

other dataset can be generated using these two functions. If you are interested
 by other ground truths you should consider modifying 'gener_data4git.m'.
If you want to use other impulse response functions, or to generate histograms
with other signal-to-backgorund ratio or average signal count you should modify 
'initdata.m'.

The main algorithm is separated into three functions:

   - 'Estim_W.m' allows for the estimation of the mixture weights

   - 'Estim_T.m' allows for the estimation of the depth profile

   - 'Estim_R.m' allows for the estimation of the spectral reflectivity

Note that in order to be able to run the algorithm, you should first compile 
the mex file 'MyHess_multi.c' using the command:

mex MyHess_multi.c

The BM3D folder contains the BM3D package for Poisson imaging that can be 
downloaded in
http://www.cs.tut.fi/~foi/invansc/
Note that this files are necessary to run 'Estim_R.m'.


