function [wc,alphaloc]=Mstep_CDir_cube(Y,Ypos,P,Fc,wc,alphaloc,ind0,WLab,nb_label)

% Compute the M-step using weak-gamma prior
% 
% INPUT:
% Y           : Histograms of photon count
% Ypos        : Photon counts in non-empty time bins
% P           : Posterior distribution
% Fc          : Impulse response functions
% wc          : current mixture weights
% alphaloc    : C-Dirichlet prior hyperparameter
% ind0        : Set of non-empty time bins
% WLab        : clusters Label
% nb_label    : Number of cluster for patch k-means
% 
% OUTPUT:
% wc          : Mixture weights estimate
% alphaloc    : Updated C-Dirichlet prior hyperparameter
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

[N,~] = size(wc);
WLab = WLab(:);
parfor i = 1:N
    wc(i,:)=update_Mstep(Y(i,:),Ypos(i),P(i,:),Fc,wc(i,:),alphaloc(WLab(i),:),ind0(i));
end

for i = 1:nb_label
    % restriction aux elements de l'image du label i
    alphaloc(i,:) = update_Dir_parameters(wc((WLab==i),:),alphaloc(i,:));
end
