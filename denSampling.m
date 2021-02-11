function Y = denSampling(X,P)

% Compute depth samples from P
% 
% INPUT:
% X        : Density support
% P        : Probability density function
%
% OUTPUT:
% Y        : samples
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

[N,D]=size(P);
cum_P = cumsum(P,2);
z = rand(N,1);
A=(z*ones(1,D)<cum_P);
ind=D-sum(A,2)+1;
Y=X(ind);
Y=Y(:);
