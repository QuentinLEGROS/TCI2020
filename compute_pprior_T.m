function pprior=compute_pprior_T(T,Nz,V,c)


% Compute the log of the depth TV prior
% 
% INPUT:
% T        : Current depth estimate
% Nz       : Number of admissible position
% V        : indices of pixel neighbors
% c        : TV prior weight
%
% OUTPUT:
% pprior   : Log TV prior
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.


[N,NN]=size(V);
Xm=200;
T1=ones(N,1)*(1:Nz); %N x Nz
pprior=zeros(N,Nz);
for t=1:NN
    X=abs(T1-T(V(:,t))*ones(1,Nz));
    X(isnan(X))=0;
    X(X>Xm)=Xm;
    pprior=pprior + X;
end

pprior=-2*c*pprior;