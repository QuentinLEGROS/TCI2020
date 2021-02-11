function P_out=compute_P_VEM_4n_rest(plik0,T,c,Titer)

% Estimate the posterior distribution using samples
% 
% INPUT:
% plik0     : Log likelihood
% T         : Current depth estimate
% c         : TV prior weight
% Titer     : Number of MC samples
%
% OUTPUT:
% P         : Posterior distribution
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.


Nx=size(plik0,2);
[Nrow,Ncol]=size(T);
Z=NaN*ones(Nrow+2,Ncol+2);
P1=zeros(Nrow+2,Ncol+2,Nx);
P1(2:end-1,2:end-1,:)=reshape(plik0,Nrow,Ncol,Nx);
P1=reshape(P1,(Nrow+2)*(Ncol+2),Nx);
T1=NaN*ones(Nrow+2,Ncol+2);
T1(2:end-1,2:end-1)=T;
Z(2:2:end-1,2:2:end-1)=0;
Z(3:2:end-1,3:2:end-1)=0;
Z(2:2:end-1,3:2:end-1)=1;
Z(3:2:end-1,2:2:end-1)=1;

P_out = zeros(Nx,Nrow*Ncol);
Pindvec = (0:Nx:Nrow*Ncol*Nx - 1)';

flat = @(x)x(:);


T2=T1;
ind0=cell(2,1);
for t=1:Titer
    for r=randperm(2)
        ind0{r}=find(Z==r-1);
        V=[ind0{r}-1 ind0{r}+1 ind0{r}+Nrow+2  ind0{r}-(Nrow+2)];
        plik=P1(ind0{r},:);
        pprior=compute_pprior_T(T2(:),Nx,V,c);
        ppost=plik+pprior;
        ppost=ppost-max(ppost,[],2)*ones(1,Nx);
        post=exp(ppost);
        ppost=post./(sum(post,2)*ones(1,Nx));
        T2(ind0{r})=denSampling(1:Nx,ppost);
    end
T3 = flat(T2(2:end-1,2:end-1)) + Pindvec;
P_out(T3) = P_out(T3) + 1;

end
P_out = transpose(reshape(P_out,Nx,Nrow*Ncol));

