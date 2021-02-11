function R_out=Estim_R(Y,Fc,W_out,Nrow,Ncol,Ns)

% Main algorithm: estimate the reflectivity profile
% 
% INPUT:
% Y         : Histograms of photon count
% Fc        : Impulse response functions
% W_out     : Estimated mixture weight 
% Nrow      : Number of rows
% Ncol      : Number of columns
% Ns        : Number of spectral component
% alp1      : weight relative to the time exposure : from reflevtivity to
%             intensity
%
% OUTPUT:
% R_out    : Reflectivity estimate
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

alp1 = sum(Fc(:,1,1)); % from intensity to reflevtivity
addpath('BM3D');

flat = @(x)x(:); 
flatI = @(x)reshape(x,[Nrow,Ncol]);

F_su=squeeze(sum(Fc(:,1,:),1)); % Integral of the IRFs
EY_bar = iterVSTpoisson(flatI(sum(Y,2))); % denoising E[y]

R_out = zeros(size(W_out));
for ll = 1:Ns
    R_out(:,ll) = (W_out(:,ll).*EY_bar(:))./(F_su(ll));
%     R_out(:,ll) = flat(iterVSTpoisson(flatI(R_out(:,ll))));
end
