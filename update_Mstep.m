function wc=update_Mstep(Y,Ypos,P,Fc,wc,alpha,ind0)


% Compute the M-step using weak-Dirichlet prior
% 
% INPUT:
% Y           : Histograms of photon count
% Ypos        : photon counts in non-empty time bins
% Fc          : Impulse response functions
% wc          : current mixture weights
% alpha       : Dirichlet prior hyperparameter
% ind0        : set of non-empty time bins
% P           : Posterior distribution
% 
% OUTPUT:
% wc         : estimated mixture weights
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

Ns = size(wc,2);
[N,L]=size(Y);
parfor n=1:N % Update pixelwise
    ii=ind0{n}; % detection of empty histograms
    if isempty(ii) % If no received photon : max the posterior relies to max the prior
        wc(n,:)=(alpha(1:Ns)-1)/(sum(alpha)-length(alpha)); % mode of the Dirichlet distribution
    else % else : max the posterior relies to max the prior
        
        
    ind=find(P(n,:)>0.01*max(P(n,:))); % Select useful value of the posterior - speed-up the estimation
    y=Ypos{n};
    p=P(n,ind)';
    wc(n,:) = Newton_Raphson(y',wc(n,:),L,p,alpha,Fc(ii,ind,:));
    end
end


