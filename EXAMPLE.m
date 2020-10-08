% This example performs Laplace-Beltrami eigenfunction Data Augmentation
% (LB-eigDA) and Chebyshev polynomial Data Augmentation (C- pDA) to
% generate augmented data from simulated data on a hippocampus surface mesh
% [1] and generate Figure 2 in [1]. 
%
%
% References:
% [1] Huang, S.-G., Chung, M.K., Qiu, A.: Fast Mesh Data Augmentation via 
% Chebyshev Polynomial of Spectral filtering. arXiv:2010.02811, 2020.
%
% [2] Tan, M., Qiu, A.: Spectral Laplace-Beltrami wavelets with applications
% in medical images. IEEE Transactions on Medical Imaging 34, 1005-1017, 2015
%
%
% (C) 2020  Shih-Gu Huang    shihgu@gmail.com
%           Anqi Qiu         bieqa@nus.edu.sg
%           National University of Singapore
%
% Update history:
%     Oct 6, 2020 created by Huang
%--------------------------------------------------------------------------

%% Load hippocampus surface mesh.
load('hippocampus_l.mat')           % left hippocampus surface
nvertex=size(surf.vertices, 1)       % number of vertices



%% Simulated data on the hippocampus
n_sim0=500  % number of simulated data in Group 0
n_sim1=500  % number of simulated data in Group 1
n_aug0=500  % number of augmented data in Group 0
n_aug1=500  % number of augmented data in Group 1

sigma=0.6

% Simulated data in Group 0
sim0=normrnd(0, sigma, nvertex, n_sim0);                % normal distribution N(0,sigma^2)


% Simulated data in Group 1
sim1=normrnd(0, sigma, nvertex, n_sim1);                % normal distribution N(0,sigma^2)
idx=find(vecnorm((surf.vertices-[101 90 123]).')<3.2);  % a small patch on the hippocampus
signal=zeros(nvertex, 1);
signal(idx)=1;
sim1=sim1+signal;                                       % Add signal 1 in the small patch



%% Discretization of Laplace-Beltrami (LB) operator
addpath('./LB/')
% folder containing some functions modified from Spectral Laplace-Beltrami Wavelets [2] for computing LB-operator

Delta=LB_operator(surf);
eigmax=eigs(Delta, 1, 'lm');  % maximun eigenvalue of LB-operator



%% C-pDA method
K = 5000                            % Use K Chebyshev polynomials for filter approximation
cutoff = linspace(0,eigmax,110);    % Each bandpass filter has bandwidth 0.1 in this example
cutoff = 2*cutoff/eigmax-1;         % Shift and scale from [0, eigmax] to [-1, 1]

BPcoeff = BP_chebyshev(cutoff, K);  % Chebyshev expansion coefficients of all bandpass filters

Cp_aug0 = CpDA(Delta, eigmax, BPcoeff, sim0, n_aug0);   % augmentation for Group 0
Cp_aug1 = CpDA(Delta, eigmax, BPcoeff, sim1, n_aug1);   % augmentation for Group 1



%%  LB-eigDA method
n_eig=nvertex       % number of eigenfunctions (used all 1184 eigenfunctions in this example)

[eigvec, eigval] = eigen(Delta, n_eig);         % Compute eigenfunctions of the LB-operator

eig_aug0 = eigDA(eigvec, sim0, n_aug0);     % augmentation for Group 0
eig_aug1 = eigDA(eigvec, sim1, n_aug1);     % augmentation for Group 1



%% Show the simulated and augmented data in Group 1
surf_plot=surf;
surf_plot.vertices=surf_plot.vertices(:, [2 1 3]);       % used for visualization

figure
subplot(2,6,1); figure_trimesh(surf_plot, mean(sim1, 2), 'rywb');
set(gca,'Zdir','reverse'); view([75 25]); caxis([-1 1]); colorbar off; camlight;
title('(a)')

for i=2:6
subplot(2,6,i); figure_trimesh(surf_plot, eig_aug1(:,i), 'rywb');
set(gca,'Zdir','reverse'); view([75 25]); caxis([-1 1]); colorbar off; camlight;
if i==2, title('(b) LB-eigDA'); end
end

for i=2:6
subplot(2,6,i+6); figure_trimesh(surf_plot, Cp_aug1(:,i), 'rywb');
set(gca,'Zdir','reverse'); view([75 25]); caxis([-1 1]); colorbar off; camlight;
if i==2, title('(c) C-pDA'); end
end
colorbar




