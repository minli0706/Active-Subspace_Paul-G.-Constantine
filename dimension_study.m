%%
% PLOTS
% Four 1d plots: 
% (1) corr 0.1, ASM, projected 300 samples, 5 design points, mean and 2sigma
% (2) corr 0.1, KL, projected 300 samples, 5 design points, mean and 2sigma
% (3) corr 1, ASM, projected 300 samples, 5 design points, mean and 2sigma
% (4) corr 1, KL, projected 300 samples, 5 design points, mean and 2sigma
% Singular value plots
% KL singular vals, ASM with corr 0.1, ASM with corr 1
% Four 2d plots:
% (1) corr 0.1, ASM, projected 300 samples, 25 design points, mean and 2sigma
% (2) corr 0.1, KL, projected 300 samples, 25 design points, mean and 2sigma
% (3) corr 1, ASM, projected 300 samples, 25 design points, mean and 2sigma
% (4) corr 1, KL, projected 300 samples, 25 design points, mean and 2sigma


clear all; close all

% Get the PDE geometry, mesh, and boundary data
pde_data = get_pde_data();

% Get KL bases
corr_length = 1; % correlation length for PDE random coefficients

% Load the initial random study
T=load('testing0.mat'); X = T.X; clear T;
% M is number of samples
% m is dimension of space
[M,m] = size(X);

% Get the PDE solutions and gradients
if corr_length == 1
    filename='experiment1.mat';
elseif corr_length == 0.01
    filename='experiment2.mat';
else
    filename=sprintf('%0.10d.mat',randi(1e9));
end
[U,~] = get_kl_bases(corr_length,m,pde_data,filename);
[f,G] = get_pde_solutions(X,U,pde_data,filename);


%% ACTIVE SUBSPACES
% Domain and design on the subspace
N=5; ubnd=3; lbnd=-3;

% Getting the active subspace
for n=1:5

    [~,Sig,W] = svd(G,'econ');
    lambda = (1/M)*diag(Sig).^2;
    W1 = -W(:,1:n);

    % Get the reduced domain samples
    yy = linspace(lbnd,ubnd,N)';
    Ygrid = 1;
    for i=1:n
        Ygrid = [kron(Ygrid,ones(N,1)) kron(ones(size(Ygrid,1),1),yy)];
    end
    Ygrid = Ygrid(:,2:end);
    Xnew_asm = Ygrid*W1';
    
    if corr_length == 1
        fnew_asm = get_pde_solutions(Xnew_asm,U,pde_data,sprintf('train_n%d_c1.mat',n));
    elseif corr_length == 0.01
        fnew_asm = get_pde_solutions(Xnew_asm,U,pde_data,sprintf('train_n%d_c2.mat',n));
    end
    

    % Compute errors at test sites
    gamma_bnds = [var(f)/sum(lambda) 6*sqrt(m)/pi];
    Ystar_asm = X*W1;
    [gp_mean_asm,gp_var_asm,parms] = ...
        gpml_regression(Ygrid,fnew_asm,Ystar_asm,gamma_bnds,lambda);
    errs_asm(:,n) = abs(gp_mean_asm-f)./abs(f);
    fprintf('Error: %6.4e\n',mean(errs_asm));

end
%%
% Plot the error histograms
close all;
figure(1)
[barz,cz] = hist([log10(errs_asm) log10(errs_sens)],25);
bar(cz,barz/300,1.5);
set(gca,'FontSize',14);
xlabel('$\log_{10}\;\left[(|f-\tilde{f}|)/|f|\right]$','Interpreter','latex','FontSize',18);
xlim([cz(1)-0.1 cz(end)+0.1]);
legend('ASM','Sens','Location','NorthWest');
if corr_length==0.01
    if n==1
        print(sprintf('figs/err_vsmall_corr_1d'),'-depsc2','-r300');
    else
        print(sprintf('figs/err_vsmall_corr_2d'),'-depsc2','-r300');
    end
elseif corr_length==1
    if n==1
        print(sprintf('figs/err_large_corr_1d'),'-depsc2','-r300');
    else
        print(sprintf('figs/err_large_corr_2d'),'-depsc2','-r300');
    end
else 
    if n==1
        print(sprintf('figs/err_%0.4d_corr_1d',randi(1000,1)),'-depsc2','-r300');
    else
        print(sprintf('figs/err_%0.4d_corr_2d',randi(1000,1)),'-depsc2','-r300');
    end
end




