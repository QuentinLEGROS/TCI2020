function WLab=Patch_cube_multi(w,si,nb_label,met)


% Compute the patches to initialize the C-Dirichlet prior model
% 
% INPUT:
% w            : current mixture weights
% si           : Neighborhood for patch k-means
% nb_label     : Number of cluster for patch k-means
% met          : for image border extension: 1=symmetric | 2=periodic
% 
% OUTPUT:
% WLab         : Clustered image indicating which pixel are in each patch
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen% , 2020. 
%        Expectation-Maximization based approach to 3D reconstruction from single-
%        waveform multispectral Lidar data. IEEE Transactions on Computational 
%        Imaging, 6, pp.1033-1043.

% Init
[Nrow,Ncol,Ns] = size(w);
N = Nrow*Ncol;
Pat = zeros(Nrow,Ncol,Ns*(2*si+1)^2);

if met == 1 % symetric
    Pat_temp = [w(:,si:-1:1,:) w w(:,end-1:-1:end-si,:)];
    Pat_temp = [Pat_temp(si:-1:1,:,:); Pat_temp; Pat_temp(end-1:-1:end-si,:,:)];
else % periodic
    Pat_temp = [w(:,end-si:end,:) w w(:,1:si,:)];
    Pat_temp = [Pat_temp(end-si:end,:,:); Pat_temp; Pat_temp(end-si:end,:,:)];
end
% Init patches

for i = si+1:Nrow+si
    for j = si+1:Ncol+si
        tem = Pat_temp(i-si:i+si,j-si:j+si,:);
        Pat(i-si,j-si,:) = tem(:);    
    end
end
Pat = reshape(Pat,N,size(Pat,3));

% Comparaison with K-mean
% Pa = sort(Pa,2);
id = kmeans(Pat,nb_label);
Idx1 = id(:);

% Output treatment
indk = cell(1,nb_label);
WLab = zeros(Nrow*Ncol,1);
for i = 1:nb_label
    [indk{i}, ~] = find(Idx1==i);
    WLab(indk{i}) = i;
end

% reshape
WLab = reshape(WLab,sqrt(length(WLab)),sqrt(length(WLab)));

