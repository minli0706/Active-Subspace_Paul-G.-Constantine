%%
% Write and store a set of initial points for testing.
clear all; close all;
M=300; m=100;
X = randn(3*M+25,m);
save('gp/full_training.mat','X');
