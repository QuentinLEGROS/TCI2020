clear all
close all
clc

% Generate data for the online version of the TIP and EUSIPO papers



A = zeros(150,150,3);
C = [75 75 30
    50 30 20];
C2 = [150 50 10];


temp = insertShape(A(:,:,1),'FilledCircle',C);
A(:,:,1) = A(:,:,1)+((20*temp(:,:,1)));
temp = insertShape(A(:,:,1),'FilledCircle',C2);
A(:,:,1) = A(:,:,1)+((40*temp(:,:,1)));
A = squeeze(A(:,:,1));
A = A+27;

A(A>500) = 135;
A(100:140,10:50)=95;
A = min(A,135);
imagesc(squeeze(A))



B = zeros(size(A,1)*size(A,2),3);
ind1 = find(A == 27); % back
ind2 = find(A == 135); % big circles
ind3 = find(A == 51); % small circle
ind4 = find(A == 95); % rectangle
ind5 = 1:numel(A);
ind5([ind1;ind2;ind3;ind4]) = [];
A(ind5)=27;
ind1 = sort([ind1;ind5']);

% back             big circles        small circle       rectangle
B(ind1,1) = 0.2;   B(ind2,1) = 0.4;   B(ind3,1) = 0.5;   B(ind4,1) = 0.5;   
B(ind1,2) = 0.4;   B(ind2,2) = 0.7;   B(ind3,2) = 0.5;   B(ind4,2) = 0.2;
B(ind1,3) = 0.2;   B(ind2,3) = 0.4;   B(ind3,3) = 0.3;   B(ind4,3) = 0.5;
% B(ind1,4) = 0.2;   B(ind2,4) = 0.5;   B(ind3,4) = 0.6;   B(ind4,4) = 0.6;

% figure(1)
% imagesc(A)
% 
% figure(2)
% subplot(3,1,1)
% imagesc(reshape(B(:,1),150,150))
% subplot(3,1,2)
% imagesc(reshape(B(:,2),150,150))
% subplot(3,1,3)
% imagesc(reshape(B(:,3),150,150))



T0 = A;
A0 = reshape(B,size(A,1),size(A,2),3);

save test_data.mat T0 A0

