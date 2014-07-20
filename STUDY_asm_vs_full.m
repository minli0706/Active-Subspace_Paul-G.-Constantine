%%
% Compare testing error between GP built on a the full 100-d space to a GP
% built on the active variables. Attempted to make the comparison fair in
% terms of work / number of function evaluations.

% This script produces Figure 5.5.

clear all; close all

% Get the PDE geometry, mesh, and boundary data
pde_data = get_pde_data();

% Get KL bases
corr_length = 1; % correlation length for PDE random coefficients

% Load the initial random study
T=load('gp/testing0.mat'); X0 = T.X; clear T;
% M is number of samples
% m is dimension of space
[M,m] = size(X0);

% Get the PDE solutions and gradients
if corr_length == 1
    filename='long_corr.mat';
    testfilename='test_long_corr.mat';
    trainfilename='train_long_corr.mat';
elseif corr_length == 0.01
    filename='short_corr.mat';
    testfilename='test_short_corr.mat';
    trainfilename='train_short_corr.mat';
else
    filename=sprintf('%0.10d.mat',randi(1e9));
end
[U,~] = get_kl_bases(corr_length,m,pde_data,filename);
[f,G] = get_pde_solutions(X0,U,pde_data,filename);

%%
T=load('gp/test_full.mat'); Xtest = T.X; clear T;
ftest = get_pde_solutions(Xtest,U,pde_data,testfilename);
T=load('gp/train_full.mat'); Xtrain = T.X; clear T;
ftrain = get_pde_solutions(Xtrain,U,pde_data,trainfilename);
[gp_mean_full,gp_var_full,parms] = gpml_regression(Xtrain,ftrain,Xtest);
errs_full = abs(gp_mean_full-ftest)./abs(ftest);
fprintf('Error: %6.4e\n',mean(errs_full));

%% Dimension of subspace
% Choose n to be 1 or 2
n = 1;
if n~=1 && n~=2, error('Error: n must be 1 or 2'); end

% Domain and design on the subspace
N=5; ubnd=3; lbnd=-3;

% ACTIVE SUBSPACES
% Getting the active subspace
[~,Sig,W] = svd(G,'econ');
lambda = (1/M)*diag(Sig).^2;
W1 = -W(:,1:n);

% Get the reduced domain samples
if n==1
    Ygrid_asm=linspace(lbnd,ubnd,N)';
    Xnew_asm = Ygrid_asm*W1';
else
    [Y1,Y2]=meshgrid(...
        linspace(lbnd,ubnd,N),...
        linspace(lbnd,ubnd,N));
    Ygrid_asm = [Y1(:) Y2(:)];
    Xnew_asm = Ygrid_asm*W1';
end

fnew_asm = get_pde_solutions(Xnew_asm,U,pde_data,[]);

% Compute errors at test sites
gamma_bnds = [var(f)/sum(lambda) 6*sqrt(m)/pi];
Ystar_asm = Xtest*W1;
[gp_mean_asm,gp_var_asm,parms] = ...
    gpml_regression(Ygrid_asm,fnew_asm,Ystar_asm,gamma_bnds,lambda);
errs_asm = abs(gp_mean_asm-ftest)./abs(ftest);
fprintf('Error: %6.4e\n',mean(errs_asm));

%%
% Plot the error histograms
close all;
figure(1)
[barz,cz] = hist([log10(errs_asm) log10(errs_full)],25);
bar(cz,barz/300,1.5);
set(gca,'FontSize',14);
xlabel('$\log_{10}\;\left[(|f-\tilde{f}|)/|f|\right]$','Interpreter','latex','FontSize',18);
xlim([cz(1)-0.1 cz(end)+0.1]);
legend('ASM','Full','Location','NorthWest');
if corr_length==0.01
    if n==1
        print(sprintf('figs/full_err_vsmall_corr_1d'),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_vsmall_corr_2d'),'-depsc2','-r300');
    end
elseif corr_length==1
    if n==1
        print(sprintf('figs/full_err_large_corr_1d'),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_large_corr_2d'),'-depsc2','-r300');
    end
else 
    if n==1
        print(sprintf('figs/full_err_%0.4d_corr_1d',randi(1000,1)),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_%0.4d_corr_2d',randi(1000,1)),'-depsc2','-r300');
    end
end

%%
% Plot comparison of errors plots.
figure(2)
plot(ftest,gp_mean_full,'bx',...
    ftest,gp_mean_asm,'ro',...
    'MarkerSize',12,'LineWidth',2);
axis square; grid on;
set(gca,'FontSize',14);
xlabel('True');
ylabel('Kriging');
xlim([min(ftest) max(ftest)]);
ylim([min(ftest) max(ftest)]);
legend('Full','ASM','Location','NorthWest');
if corr_length==0.01
    if n==1
        print(sprintf('figs/full_err_vsmall_corr_1d_comp'),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_vsmall_corr_2d_comp'),'-depsc2','-r300');
    end
elseif corr_length==1
    if n==1
        print(sprintf('figs/full_err_large_corr_1d_comp'),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_large_corr_2d_comp'),'-depsc2','-r300');
    end
else 
    if n==1
        print(sprintf('figs/full_err_%0.4d_corr_1d_comp',randi(1000,1)),'-depsc2','-r300');
    else
        print(sprintf('figs/full_err_%0.4d_corr_2d_comp',randi(1000,1)),'-depsc2','-r300');
    end
end



