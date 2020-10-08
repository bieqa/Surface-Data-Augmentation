function [eigvec, eigval] = eigen(Delta, K)
%--------------------------------------------------------------------------
% Compute all or the first K eigenfunctions of LB-operator or graph Laplacian
%
% Delta  :  LB-operator or graph Laplacian
% K      :	number of eigfunctions
% eigvec :	eigfunctions
% eigval :	eigvalues 
%
%
% (C) 2020  Shih-Gu Huang    shihgu@gmail.com
%           Anqi Qiu         bieqa@nus.edu.sg
%           National University of Singapore
%
% Update history:
%     Oct 6, 2020 created by Huang
%--------------------------------------------------------------------------

nvertex=size(Delta, 1); % number of vertices 

if K==nvertex
    [eigvec, D] = eig(full(Delta));         % compute all eigenfunctions
else
    [eigvec, D] = eigs(Delta, K, 'sr');     % compute the first K eigenfunctions
end
eigval=diag(D);                         % eigenvalues

% sort eigenvalues in an ascend order, and sort eigenfunctions by the order
[~,ind]=sort(abs(eigval), 'ascend');
eigval=eigval(ind);
eigvec=eigvec(:,ind);
