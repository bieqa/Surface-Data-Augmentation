function BPcoeff = BP_chebyshev(cutoff, K)
%--------------------------------------------------------------------------
% Compute Chebyshev expansion coefficients of all bandpass filters
%
% cutoff   :	cutoff frequencies in [-1,1] of all bandpass filters
% K        :	order of Chebyshev polynomials
% BPcoeff  :    Chebyshev expansion coefficients, size = K x (length(cutoff)-1)
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

cutoff=cutoff(:);


%% Chebyshev polynomials of the second kind, U_k, usded for the closed-form solution

% initial conditions of recurrence relation of U_k
Uold =zeros(size(cutoff));    % U_{-1} = 0
U =ones(size(cutoff));        % U_0  = 1
Uall=U;

% recurrence relation of U_k
for n=0:K-2
    Unew = 2 * cutoff .* U  - Uold;
    Uold =U;
    U =Unew;
    Uall=[Uall U];
end


%% Chebyshev expansion coefficients in closed form
BPcoeff=[];
for l=1:length(cutoff)-1
    a=cutoff(l);                    % lower cutoff frequency of l-th filter
    b=cutoff(l+1);                  % upper cutoff frequency of l-th filter      
    
    coeff0=(acos(a)-acos(b))/pi;    % expansion coefficient associated with T_0
%     Kall=[1:K];
    coeff1=sqrt(1-a^2).*Uall(l,1:end-1)-sqrt(1-b^2).*Uall(l+1,1:end-1);
    coeff1=coeff1./[1:K-1]*2/pi;           % expansion coefficients associated with T_1, T_2, ..., T_{K-1}
    
    BPcoeff=[BPcoeff [coeff0; coeff1.']];
end
