function aug = eigDA(eigvec, org, n_aug)
%--------------------------------------------------------------------------
% This function generates augemented data by eigenfunction Data
% Augmentation. The mathematical detail is given in [1]. 
%
% eigvec   :  eigenfunctions
% org      :  real data
% n_aug    :  number of augmented data
% aug      :  augmented data
%
%
% Reference:
% [1] Huang, S.-G., Chung, M.K., Qiu, A.: Fast Mesh Data Augmentation via 
% Chebyshev Polynomial of Spectral filtering. arXiv:2010.02811, 2020.
%
%
% (C) 2020  Shih-Gu Huang    shihgu@gmail.com
%           Anqi Qiu         bieqa@nus.edu.sg
%           National University of Singapore
%
% Update history:
%     Oct 6, 2020 created by Huang
%--------------------------------------------------------------------------

n_org=size(org,2);       % number of real data
n_eig=size(eigvec,2);    % number of eigenfunctions


coeff=pinv(eigvec)*org;  % eigenfunction coefficients
res=org-eigvec*coeff;    % residual between real data and the apprximation by a subset of eigenfunctions

p=permutation(n_org, n_aug, n_eig);      % permutations

coeffaug=zeros(n_eig, n_aug);            
for j=1:n_eig
    coeffaug(j, :)=coeff(j, p(:,j));     % permutations of coefficients
end
aug=eigvec*coeffaug + res(:, p(:,end));  % augmented data (adding residual back to avoid data smoothing)

