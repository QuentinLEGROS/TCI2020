function [W_out,P,ind0]=Estim_W(Y,Fc,Nrow,Ncol,c,reg_w,si,step_Nx,nb_label,w0)
    

% Main algorithm: estimate the mixture weights
% 
% INPUT:
% Y         : Histograms of photon count
% Fc        : Impulse response functions
% Ns        : Number of spectral component +1 for the background
% Nrow      : Number of rows
% Ncol      : Number of columns
% c         : Depth TV prior weight
% si        : Neighborhood for patch k-means
% step_Nx   : Likelihood subsampling factor
% nb_label  : Number of cluster for patch k-means   
% reg_w     : Prior model for W
%
% OUTPUT:
% W_out    : Mixture weight estimates
% P        : Depthp Posterior distribution
% ind0     : Non-empty time bins
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.


%% Initialization
[N,T]=size(Y);% number of pixels * histogram length
Nx=size(Fc,2); % number of admissible depth positions
Ns = size(Fc,3); % number of wavelength
Nz = Nx/step_Nx;
ind0=cell(N,1); % cell for non-empty time bins
Ypos=cell(N,1); % cell for photon count in the non-empty time bins
np=zeros(N,1); % cell with number of non-empty bins in each histogram

for n=1:N
    ind0{n}=find(Y(n,:)>0); % find non-empty time bins
    Ypos{n}=Y(n,ind0{n}); % find photon counts in these bins
    np(n)=numel(ind0{n}); % find number of non-empty bins
end

wc=0.01*ones(N,Ns); % Initialization current W
wold = wc; % Store previous iteration of W
iteEM = 0; % Number of EM iteration used to average W
pT=ones(1,Nz)./Nx; % Initialization uniform prior for T (first EM iterations)
alpha=1.01*ones(1,Ns+1); % Initialization W priors fixed : alpha,beta (first EM iterations)
if strcmp(reg_w,'C-Dirichlet')
    % Initialization W priors : alpha,beta (further updated)
    alphaloc = 1.01*ones(nb_label,length(alpha));
    cond=0;
end

Fct = zeros(T,Nx/step_Nx,Ns);
for i = 1:Nx/step_Nx
    Fct(:,i,:) = sum(Fc(:,(i-1)*step_Nx+1:i*step_Nx,:),2);
end

C=compute_plik(Y,Fct,wc,ind0); % log-Likelihood
C=C-max(C,[],2)*ones(1,Nz); % PAS SUR
P=exp(C);
P=P./(sum(P,2)*ones(1,Nz)); % Normalization
T0=denSampling(1:Nz,P); % Gibbs sampling
T0 = T0.*step_Nx;
    
   
W_out=zeros(N,Ns); % Initialization output 
m_compt=1; % Initialization of EM iteration Count
Stopping_EM = 0; % EM stopping criterion
err0 = 10e9;
derr = 10e9;
errold = 10e9;
%% Main iterative process
% while ((m_compt<=11 && derr>1e-3 && err0>1e-10) || m_compt<=10)
while (iteEM<5)
    disp(['EM iteration : ',num2str(m_compt)])
    %% Expectation / Stochastic Expectation step
     %TV-depth regularization
    if m_compt<3 % first iterations without TV
       C=compute_plik(Y,Fct,wc,ind0); 
       C=C+ones(N,1)*log(pT);
       C=C-max(C,[],2)*ones(1,Nz);
       P=exp(C);
       P=P./(sum(P,2)*ones(1,Nz));
    else % with TV prior on depth
        C=compute_plik(Y,Fct,wc,ind0); % Compute log likelihood P(S|T,W)
        [P,T0]=Prior_TV_depth(C,reshape(T0./step_Nx,Nrow,Ncol),c); % Generating depth samples from posterior P T
        T0 = T0.*step_Nx;
    end

    %% Maximization step
    if m_compt<3 % first iterations without prior on T and W-Dirichlet on W
        wc=update_Mstep(Y,Ypos,P,Fct,wc,alpha,ind0);
    else % With prior on T
        if  strcmp(reg_w,'C-Dirichlet') % if C-Dirichlet prior on W
            % Initialization of the clusters 
            if ((cond==0) && (m_compt==5))
                 WLab=Patch_cube_multi(reshape(wc,Nrow,Ncol,size(wc,2)),si,nb_label,2); 
                 cond = 1;
                 [wc,alphaloc]=Mstep_CDir_cube(Y,Ypos,P,Fct,wc,alphaloc,ind0,WLab,nb_label);
            elseif ((cond==0) && (m_compt<5))
                wc=update_Mstep(Y,Ypos,P,Fct,wc,alpha,ind0);
            elseif cond == 1 % If cluster already defined
                [wc,alphaloc]=Mstep_CDir_cube(Y,Ypos,P,Fct,wc,alphaloc,ind0,WLab,nb_label);
            else 
                wc=update_Mstep(Y,Ypos,P,Fct,wc,alpha,ind0);
            end
        else % W-Dirichlet prior
            wc=update_Mstep(Y,Ypos,P,Fct,wc,alpha,ind0);
        end
    end
    
%     
%     figure
%     subplot(3,2,1)
%     imagesc(reshape(w0(:,1),Nrow,Ncol))
%     title('Current r1')
%     subplot(3,2,3)
%     imagesc(reshape(w0(:,2),Nrow,Ncol))
%     title('Current r2')
%     subplot(3,2,5)
%     imagesc(reshape(w0(:,3),Nrow,Ncol))
%     title('Current r3')
%     subplot(3,2,2)
%     imagesc(reshape(wc(:,1),Nrow,Ncol))
%     title('Estim w1')
%     subplot(3,2,4)
%     imagesc(reshape(wc(:,2),Nrow,Ncol))
%     title('Estim w2')
%     subplot(3,2,6)
%     imagesc(reshape(wc(:,3),Nrow,Ncol))
%     title('Estim w3')
%     pause(0.01)
    
    err0=sum(sum(abs(wc-wold)))./min(sum(sum(abs(wc))),sum(sum(abs(wold))));
%     err0=sum(norm(wc-wold)./min(norm(wold),norm(wc))); % relative l2 norm between consecutive W estimates
    derr = norm(err0-errold);
    wold=wc; % store previous W estimate
    errold = err0;
    
    m_compt=m_compt+1;

    if (((m_compt>5) && (derr<1e-2)) || m_compt==20)  % Begin the last 5 iterations
        Stopping_EM = 1;
    end
    if Stopping_EM == 1
        W_out=W_out+wc;
        iteEM = iteEM + 1;
    end
end
W_out = W_out ./ iteEM; % mean of the (N_iter-N_bi) iterations

