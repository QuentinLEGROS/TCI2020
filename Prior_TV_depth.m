function [P,T_out]=Prior_TV_depth(plik0,T,c)

% Compute the posterior distribution and sample the depth
% 
% INPUT:
% plik0     : Log likelihood
% T         : Current depth estimate
% c         : TV prior weight
%
% OUTPUT:
% P         : Posterior distribution
% T_out     : New depth estimates
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

%% Initialization
Nx=size(plik0,2);
[Nrow,Ncol]=size(T);
% Create temporary images
Z=NaN*ones(Nrow+2,Ncol+2); % for pixel indices
P1=zeros(Nrow+2,Ncol+2,Nx); % Log likelihood
P1(2:end-1,2:end-1,:)=reshape(plik0,Nrow,Ncol,Nx);
P1=reshape(P1,(Nrow+2)*(Ncol+2),Nx);
T1=NaN*ones(Nrow+2,Ncol+2); % Current depth
T1(2:end-1,2:end-1)=T;
Z(2:2:end-1,2:2:end-1)=0; % Odd pixels
Z(3:2:end-1,3:2:end-1)=0; % Odd pixels
Z(2:2:end-1,3:2:end-1)=1; % Even pixels
Z(3:2:end-1,2:2:end-1)=1; % Even pixels
C=zeros((Nrow+2)*(Ncol+2),Nx); % Log likelihood to estimate

%% Sampling
%T2=T1;
ind0=cell(2,1);
for t=1:2 % Sample twice
    for r=randperm(2) %Iterate on even or odd pixels
        ind0{r}=find(Z==r-1); % Select even or odd pixels
        V=[ind0{r}-1 ind0{r}+1 ind0{r}+Nrow+2  ind0{r}-(Nrow+2)]; % Neighborhood of each selected pixel
        plik=P1(ind0{r},:); % log likelihood associated to selected pixels
        pprior=compute_pprior_T(T1(:),Nx,V,c); % Compute log TV prior
        ppost=plik+pprior; % Compute log posterior
        ppost=ppost-max(ppost,[],2)*ones(1,Nx);
        post=exp(ppost); %Compute posterior
        ppost=post./(sum(post,2)*ones(1,Nx)); % Normalization
        T1(ind0{r})=denSampling(1:Nx,ppost); % Sample from posterior
    end
end

%% Compute posterior from samples
T_out=T1;
for r=randperm(2)
    ind0{r}=find(Z==r-1);
    V=[ind0{r}-1 ind0{r}+1 ind0{r}+Nrow+2  ind0{r}-(Nrow+2)];
    plik=P1(ind0{r},:);
    pprior=compute_pprior_T(T1(:),Nx,V,c);
    C((ind0{r}),:)=plik+pprior;
end

C=reshape(C,Nrow+2,Ncol+2,Nx);
C=C(2:end-1,2:end-1,:);
C=reshape(C,Nrow*Ncol,Nx);
C=C-max(C,[],2)*ones(1,Nx);
P=exp(C); % Compute posterior
P=P./(sum(P,2)*ones(1,Nx)); % Normalization
T_out=T_out(2:end-1,2:end-1);
T_out=T_out(:);