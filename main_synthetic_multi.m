clear all                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     % clear all
close all
clc
warning('off','all')


%% Initialize Parameters
c=0.05; % Range regularization parameter
step_Nx = 10;
si=1; % Neighborhood for clustering -> patches of size (2*si+1)x(2*si+1)
nb_label = 5; % Number of cluster
Titer = 200; % : Number of MC samples for MC depth sampling

%% Initialize data
data=load('exampleCirc.mat'); 
Y=data.Y;
Fc=data.Fc;
Nrow=data.Nrow;
Ncol=data.Ncol;
A0=data.A0;
T0=data.T0;
 
%% Prior choice
Chpw{1} = 'W-Dirichlet';
Chpw{2} = 'C-Dirichlet';

met_w=Chpw{1};

%% EM algorithm - estimate W
disp('Estimation : W')
[W_out,P,ind0]=Estim_W(Y,Fc,Nrow,Ncol,c,met_w,si,step_Nx,nb_label);

%% Sampling T from the posterior
disp('Estimation : T')
T_out=Estim_T(Y,P,Fc,W_out,ind0,Nrow,Ncol,c,step_Nx,Titer);

%% Compute R from W - BM3D algorithm
disp('Estimation : R')
R_out=Estim_R(Y,Fc,W_out,Nrow,Ncol,size(W_out,2));

figure(1)
subplot(2,1,1)
imagesc(T0)
title('true T')
subplot(2,1,2)
imagesc(reshape(T_out,Nrow,Ncol))
title('estim T')

figure(2)
subplot(3,2,1)
imagesc(squeeze(A0(:,:,1)))
title('Current r1')
subplot(3,2,3)
imagesc(squeeze(A0(:,:,2)))
title('Current r2')
subplot(3,2,5)
imagesc(squeeze(A0(:,:,3)))
title('Current r3')
subplot(3,2,2)
imagesc(reshape(R_out(:,1),Nrow,Ncol))
title('Estim r1')
subplot(3,2,4)
imagesc(reshape(R_out(:,2),Nrow,Ncol))
title('Estim r2')
subplot(3,2,6)
imagesc(reshape(R_out(:,3),Nrow,Ncol))
title('Estim r3')
