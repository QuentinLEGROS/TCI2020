function C=compute_plik(Y,Fc,wc,ind0)


% Compute the log likelihood 
% 
% INPUT:
% Y        : Histograms of photon count
% ind0     : set of non-empty time bins
% Fct      : Impulse response functions
% wc       : current mixture weights

% OUTPUT:
% C        : Log likelihood
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.


[N,L]=size(Y);
[~,Nx,Ns]=size(Fc);
C=zeros(N,Nx);
for n=1:N
ii=ind0{n};
    if ~isempty(ii)     
        A =(1-sum(wc(n,:)))/L;
        for ain = 1:Ns
            A = A+(wc(n,ain)*Fc(ii,:,ain));
        end
        C(n,:)=sum((Y(n,ii)'*ones(1,Nx)).*log(A),1);
    end
end

