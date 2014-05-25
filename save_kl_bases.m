%%
close all; clear all;

% Get the PDE geometry, mesh, and boundary data
pde_data = get_pde_data();

% Get KL bases
corr_length = 1; % correlation length for PDE random coefficients
filename='pde/long_corr.mat';
m = 100; % number of parameter in high-d space
[U,~] = get_kl_bases(corr_length,m,pde_data,filename);

% Get KL bases
corr_length = 0.01; % correlation length for PDE random coefficients
filename='pde/short_corr.mat';
m = 100; % number of parameter in high-d space
[U,~] = get_kl_bases(corr_length,m,pde_data,filename);
