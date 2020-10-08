function p=permutation(n_org, n_aug, n_perm)
%--------------------------------------------------------------------------
% Generate random permutations for data augmentation
%
% n_org  :  number of real data
% n_aug  :  number of augmented data
% n_perm :  number of permutations (= number of eigenfunctions or filters)
% p      :  random permutations of integers in [1, n_org] (size = n_aug x n_perm)
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

if n_aug<=n_org
    p=[];
    for i=1:n_perm
        p=[p randperm(n_org, n_aug).'];	% n_aug unique integers selected randomly from [1, n_org]
    end                                     % repeated for each eigenfunction/filter
    
else
    run=floor(n_aug/n_org)                 % each eigenfunction/filter needs mutiple permutations if n_aug>n_org
    p=[];
    for k=1:run
        tmp=[];
        for i=1:n_perm
            tmp=[tmp randperm(n_org, n_org).'];
        end
        p=[p;tmp];
    end
    
    tmp=[];
    for i=1:n_perm
        tmp=[tmp randperm(n_org, n_aug-n_org*run).'];
    end
    p=[p;tmp];
end