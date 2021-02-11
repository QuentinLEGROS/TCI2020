function wc = Newton_Raphson(y,wc,L,p,alpha,Ft)


% Perform second order gradient descent
% 
% INPUT:
% y           : photon counts in non-empty time bins at current pixel n
% wc          : current mixture weight at pixel n
% L           : histogram length
% p           : posterior distribution at pixel n estimated in bin where
%               photons are detected
% alpha       : Dirichlet prior parameters
% Ft          : Impulse response function estimated in the current
%               estimated depth
% 
% OUTPUT:
% wc         : estimated mixtre weight
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

%% parameters and matrix pre computation
err=10;
it=0;
Ns = size(Ft,3);
wc = transpose(wc);
H = zeros(length(wc),length(wc));
% scalar 'y' for sum over p (num of photon in each bin) and 'p' being A_ijk
scal = repmat(full(y),[1 size(Ft,2)]).*repmat(p',[size(Ft,1) 1]);
%% Main loop
while (err>0.0001 && it<5 )
    w10=wc;
    A =(1-sum(wc))/L;
    for ain = 1:Ns
        A = A+(wc(ain)*Ft(:,:,ain));
    end
    
    A = (wc(1)*Ft(:,:,1) + wc(2)*Ft(:,:,2) + wc(3)*Ft(:,:,3)) + (1-sum(wc))/L;
    P = (Ft-(1/L)) ./ repmat(A,[1,1,Ns]);
    %% Gradient
    graf = squeeze(sum(sum(scal.*P,1),2));

    %% Hessian
    H = MyHess_multi(scal,P,H);
    
    %% W-Dirichlet prior
    graf = graf + (alpha(1:Ns)'-1)./wc -repmat((alpha(end)-1)/(1-sum(wc)),[Ns,1]);
    H = H -((alpha(end)-1)/(1-sum(wc))^2);
    H = H - diag((alpha(1:Ns)'-1)./(wc.^2));

    %% NR step
    G = pinv(H)*graf;

    %% projection onto [0,1]^3;
    % Weight ensuring the update to stay in the compact [0,1]^3;
    if (sum(G)<0 && numel(find(G>0))==0)
        epsi = 0.99*min( (sum(wc)-1)/(sum(G)),1); % only ensure sum-to-one constraint of the updates
    elseif (sum(G)<0 && numel(find(G>0))>0)
        wtem = wc./G; % ensure sum-to-one constraint and positivity of the updates
        epsi = 0.99*min(min(min(wtem(G>0)), (sum(wc)-1)/(sum(G))),1);
    elseif sum(G)>0
        wtem = wc./G;
        epsi = 0.99*min(min(wtem(G>0)),1); % only ensure positivity of the updates
    else
        epsi = 1e-3;
    end
      
    %% NR update
    wc = wc - epsi*G;

    %% Compute error
    err = sum(abs(wc-w10)./min(wc,w10));
    it = it + 1;
    
end

