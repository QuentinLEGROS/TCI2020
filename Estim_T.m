function T_out=Estim_T(Y,P,Fc,w_out,ind0,Nrow,Ncol,c,step_Nx,Titer)


% Main algorithm: estimate the depth profile
% 
% INPUT:
% Y         : Histograms of photon count
% P         : Depth Posterior distribution
% Fc        : Impulse response functions
% w_out     : Mixture weight estimates
% ind0      : Non-empty time bins
% Nrow      : Number of rows
% Ncol      : Number of columns
% c         : Depth TV prior weight
% Titer     : Number of MC samples
% 
% OUTPUT:
% T_out     : Depth estimates
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

[~,T0]=max(P,[],2); % hot start depth for sampling
T0 = T0.*step_Nx;
C=compute_plik(Y,Fc,w_out,ind0); % Likelihood
P_out=compute_P_VEM_4n_rest(C,reshape(T0,Nrow,Ncol),c,Titer); % Gibbs sampling

% MAP estimation
[~,T_out]=max(P_out,[],2);

