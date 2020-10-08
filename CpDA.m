function aug = CpDA(Delta, eigmax, coeff, org, n_aug)
%--------------------------------------------------------------------------
% This function generates augmented data using Chebyshev polynomial Data
% Augmentation (CpDA). The mathematical detail is given in [1]. 
% The recurrence relation part is modified from the Spectral Graph Wavelet 
% Transform (SGWT) [2] (https://wiki.epfl.ch/sgwt).
%
% Delta	 :	LB-operator or graph Laplacian
% eigmax :  maximun eigenvalue of LB-operator
% coeff	 :	Chebyshev expansion coefficients of spectral filters
% org    :	real data
% n_aug  :	number of augmented data
% aug    :	augmented data
%
%
% Reference:
% [1] Huang, S.-G., Chung, M.K., Qiu, A.: Fast Mesh Data Augmentation via 
% Chebyshev Polynomial of Spectral filtering. arXiv:2010.02811, 2020.
%
% [2] Hammond, D.K., Vandergheynst, P., Gribonval, R.: Wavelets on graphs via 
% spectral graph theory. Applied and Computational Harmonic Analysis 30, 129-150, 2011
%
%
% (C) 2020  Shih-Gu Huang    shihgu@gmail.com
%           Anqi Qiu         bieqa@nus.edu.sg
%           National University of Singapore
%
% Update history:
%     Oct 6, 2020 created by Huang
%--------------------------------------------------------------------------

n_org=size(org, 2);         % number of real data
n_filt=size(coeff, 2);      % number of filters
K=size(coeff, 1);           % order of Chebyshev polynomials


p=permutation(n_org, n_aug, n_filt+1);      % permutations


fmean=mean(org,1);  % signal mean over the surface
f=org-fmean;

aug=fmean(p(:,1));  % reample the signal means


% initial conditions of recurrence relation
Tf_old =0;   % T_{-1}(Delta)f = 0
Tf =f;       % T_0(Delta)f  = f


% update the augmented data by recurrence and resampling (in matrix computation)
coeffmtx=zeros(n_org,n_aug);	
for l=1:n_filt
    idx=sub2ind([n_org n_aug], p(:,l+1), [1:n_aug].');
    coeffmtx(idx)=coeffmtx(idx)+coeff(1,l);  % build coefficeint matrix for matrix computation
end
aug =aug+Tf*coeffmtx;


% recurrence relation of T_k(Delta)f for k=0
Tf_new =Delta*Tf*2/eigmax -Tf -Tf_old;
Tf_old =Tf;
Tf =Tf_new;

for k=1:K-1
    % update the augmented data by recurrence and resampling (in matrix computation)
    coeffmtx=zeros(n_org,n_aug);    
    for l=1:n_filt
        idx=sub2ind([n_org n_aug], p(:,l+1), [1:n_aug].');
        coeffmtx(idx)=coeffmtx(idx)+coeff(k+1,l); % build coefficeint matrix for matrix computation
    end
    aug =aug+Tf*coeffmtx;    
    
    % recurrence relation of T_k(Delta)f for k>0
    Tf_new =Delta*Tf*4/eigmax -Tf*2 -Tf_old;
    Tf_old =Tf;
    Tf =Tf_new;
    
end
