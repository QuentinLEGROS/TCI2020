function [alpha]=update_Dir_parameters(w,alpha)


% Update the parameters of the mixture weight C-Dirichlet prior
% 
% INPUT:
% w        : Current mixture weights
% alpha    : C-Dirichlet prior hyperparameter
% 
% OUTPUT:
% alpha     : Updated C-gamma prior hyperparameter for current patch 
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

% Gamma parameters
k = 1;
theta = 4;
slog=sum(log(w));
sminlog=sum(log(1-sum(w,2)));

N=size(w,1);
for rep = 1:30
    for t=1:size(w,2)
        err=10;
        alpha0=alpha(t);
        while err>1e-4
            g=slog(t) + N*(psi(sum(alpha))-psi(alpha(t))) + (k-1)/alpha(t) - (1/theta);
            gp= N*(psi(1,sum(alpha))-psi(1,alpha(t))) - (k-1)/(alpha(t)^2);
            alpha(t)=min(max(alpha(t)-g/gp,1.01),100);
            err=abs(alpha(t)-alpha0);
            alpha0=alpha(t);
        end
    end
    alpha0=alpha(end);
    err=10;
    while err>1e-4
        g=sminlog +N*(psi(sum(alpha))-psi(alpha(end))) + (k-1)/alpha(end) - (1/theta);
        gp= N*(psi(1,sum(alpha))-psi(1,alpha(end))) - (k-1)/(alpha(end)^2);
        alpha(end)=min(max(alpha(end)-g/gp,1.01),100);
        err=abs(alpha(end)-alpha0);
        alpha0=alpha(end);
    end
end




