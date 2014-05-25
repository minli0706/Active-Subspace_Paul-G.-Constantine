%%
% Compare GPs built on active subspace variables to GPs built on most
% important coordinates as determined by a local sensitivity analysis at
% the origin.

% This script produces Figures 5.2, 5.3, and 5.4.

clear all; close all

% Get the PDE geometry, mesh, and boundary data
pde_data = get_pde_data();

% Get KL bases
corr_length = 1; % correlation length for PDE random coefficients

% Load the initial random study
T=load('gp/testing0.mat'); X = T.X; clear T;
% M is number of samples
% m is dimension of space
[M,m] = size(X);

% Get the PDE solutions and gradients
if corr_length == 1
    filename='long_corr.mat';
elseif corr_length == 0.01
    filename='short_corr.mat';
else
    filename=sprintf('%0.10d.mat',randi(1e9));
end
[U,~] = get_kl_bases(corr_length,m,pde_data,filename);
[f,G] = get_pde_solutions(X,U,pde_data,filename);

%% Dimension of subspace
% Choose n to be 1 or 2
n = 2;
if n~=1 && n~=2, error('Error: n must be 1 or 2'); end

% Domain and design on the subspace
N=5; ubnd=3; lbnd=-3;


%% ACTIVE SUBSPACES
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
Ystar_asm = X*W1;
[gp_mean_asm,gp_var_asm,parms] = ...
    gpml_regression(Ygrid_asm,fnew_asm,Ystar_asm,gamma_bnds,lambda);
errs_asm = abs(gp_mean_asm-f)./abs(f);
fprintf('Error: %6.4e\n',mean(errs_asm));

% Get points for plots
if n==1
    Yte = linspace(lbnd,ubnd,501)';
    [gp_mean_asm_yte,gp_var_asm_yte,parms] = ...
        gpml_regression(Ygrid_asm,fnew_asm,Yte,gamma_bnds,lambda);
else
    [Y1,Y2] = meshgrid(linspace(lbnd,ubnd,101),linspace(lbnd,ubnd,101));
    Yte = [Y1(:) Y2(:)];
    [gp_mean_asm_yte,gp_var_asm_yte,parms] = ...
        gpml_regression(Ygrid_asm,fnew_asm,Yte,gamma_bnds,lambda);
end

%% Plots
close all;
if n==1
    figure(1); hold on;
    plot(Ystar_asm,f,'b.',...
        Yte,gp_mean_asm_yte,'k-',...
        Ygrid_asm,fnew_asm,'ko',...
        'LineWidth',1,'MarkerFace','k','MarkerSize',11);
    set(gca,'FontSize',16);
    axis square
    xlim([lbnd ubnd]);
    legend('test','GP mean','train','Location','NorthWest');
    xlabel('y','FontSize',18);
    ylabel('f','FontSize',18);
    fill([Yte;flipud(Yte)],...
        [gp_mean_asm_yte+2*sqrt(gp_var_asm_yte); flipud(gp_mean_asm_yte-2*sqrt(gp_var_asm_yte))],...
        'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
    hold off;
    if corr_length==0.01
        print(sprintf('figs/gp_asm_vsmall_corr_1d'),'-depsc','-r300');
    elseif corr_length==1
        print(sprintf('figs/gp_asm_large_corr_1d'),'-depsc','-r300');
    else
        print(sprintf('figs/gp_asm_%0.4d_corr_1d',randi(1000,1)),'-depsc','-r300');
    end
elseif n==2
    figure(1); hold on
    surf(Y1,Y2,reshape(gp_mean_asm_yte,101,101),'EdgeColor','None','FaceAlpha',0.5);
    colormap('Hot'); colorbar;
    plot3(Ystar_asm(:,1),Ystar_asm(:,2),f,'b.','MarkerSize',20);
    plot3(Ygrid_asm(:,1),Ygrid_asm(:,2),fnew_asm,'ko','MarkerFace','k','MarkerSize',11);
    set(gca,'FontSize',14);
    axis square; grid on; view(-37.5-120,10);
    xlim([lbnd ubnd]); ylim([lbnd ubnd]);
    xlabel('y_1','Interpreter','tex','FontSize',18);
    ylabel('y_2','Interpreter','tex','FontSize',18);
    zlabel('f','Interpreter','tex','FontSize',18);
    %surf(Y1,Y2,reshape(gp_mean_asm_yte+2*sqrt(gp_var_asm_yte),101,101),...
    %    'EdgeColor','None','FaceAlpha',0.1);
    %surf(Y1,Y2,reshape(gp_mean_asm_yte-2*sqrt(gp_var_asm_yte),101,101),...
    %    'EdgeColor','None','FaceAlpha',0.1);
    hold off;
    if corr_length==0.01
        print(sprintf('figs/gp_asm_vsmall_corr_2d'),'-depsc2','-r300');
    elseif corr_length==1
        print(sprintf('figs/gp_asm_large_corr_2d'),'-depsc2','-r300');
    else
        print(sprintf('figs/gp_asm_%0.4d_corr_2d',randi(1000,1)),'-depsc2','-r300');
    end
end

%% 
% Use the n (= 1 or 2) most important coordinates as determined by a local
% sensitivity analysis at the origin.

ind = get_local_sensitivity_indices(X,U,pde_data);
fprintf('Sensitivity indices: %d, %d\n',ind(1), ind(2));

if n==1
    Ygrid_sens = linspace(lbnd,ubnd,N)';
    Xnew_sens = zeros(N,m);
    Xnew_sens(:,ind(1)) = linspace(lbnd,ubnd,N)';
elseif n==2
    [Y1,Y2]=meshgrid(...
        linspace(lbnd,ubnd,N),...
        linspace(lbnd,ubnd,N));
    Ygrid_sens = [Y1(:) Y2(:)];
    Xnew_sens = zeros(N^2,m);
    Xnew_sens(:,ind(1)) = Y1(:);
    Xnew_sens(:,ind(2)) = Y2(:);
end

fnew_sens = get_pde_solutions(Xnew_sens,U,pde_data,[]);

% Train the kriging surface
if n==1
    Ystar_sens = X(:,ind(1));
    Yte = linspace(lbnd,ubnd,501)';

else
    Ystar_sens = [X(:,ind(1)) X(:,ind(2))];
    [Y1,Y2] = meshgrid(...
        linspace(lbnd,ubnd,101),...
        linspace(lbnd,ubnd,101));
    Yte = [Y1(:) Y2(:)];
end
[gp_mean_sens,gp_var_sens,parms] = gpml_regression(Ygrid_sens,fnew_sens, ...
    [Ystar_sens; Yte]);
errs_sens = abs(gp_mean_sens(1:M)-f)./abs(f);
gp_mean_sens_yte = gp_mean_sens(M+1:end);
gp_var_sens_yte = gp_var_sens(M+1:end);

fprintf('Error: %6.4e\n',mean(errs_sens));


%% Plots
% Plot response surface on low-d subspace
close all;
if n==1
    figure(1); hold on;
    plot(Ystar_sens,f,'b.',...
        Yte,gp_mean_sens_yte,'k-',...
        Ygrid_sens,fnew_sens,'ko',...
        'LineWidth',1,'MarkerFace','k','MarkerSize',11);
    set(gca,'FontSize',14);
    axis square
    xlim([lbnd ubnd]);
    legend('test','GP mean','train','Location','NorthWest');
    xlabel(sprintf('x_%d',ind(1)),'Interpreter','tex','FontSize',18);
    ylabel('f','Interpreter','tex','FontSize',18);
    fill([Yte;flipud(Yte)],...
        [gp_mean_sens_yte+2*sqrt(gp_var_sens_yte); flipud(gp_mean_sens_yte-2*sqrt(gp_var_sens_yte))],...
        'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
    hold off;
    if corr_length==0.01
        print(sprintf('figs/gp_sens_vsmall_corr_1d'),'-depsc2','-r300');
    elseif corr_length==1
        print(sprintf('figs/gp_sens_large_corr_1d'),'-depsc2','-r300');
    else
        print(sprintf('figs/gp_sens_%0.4d_corr_1d',randi(1000,1)),'-depsc2','-r300');
    end
elseif n==2
    figure(1); hold on
    surf(Y1,Y2,reshape(gp_mean_sens_yte,101,101),'EdgeColor','None','FaceAlpha',0.5);
    colormap('Hot'); colorbar;
    plot3(Ystar_sens(:,1),Ystar_sens(:,2),f,'b.','MarkerSize',20);
    plot3(Ygrid_sens(:,1),Ygrid_sens(:,2),fnew_sens,'ko','MarkerFace','k','MarkerSize',11);
    set(gca,'FontSize',14);
    axis square; grid on; view(-37.5-120,10);
    xlim([lbnd ubnd]); ylim([lbnd ubnd]);
    xlabel(sprintf('x_%d',ind(1)),'Interpreter','tex','FontSize',18);
    ylabel(sprintf('x_%d',ind(2)),'Interpreter','tex','FontSize',18);
    zlabel('f','Interpreter','tex','FontSize',18);
    %surf(Y1,Y2,reshape(gp_mean_kl_yte+2*sqrt(gp_var_kl_yte),101,101),...
    %    'EdgeColor','None','FaceAlpha',0.1);
    %surf(Y1,Y2,reshape(gp_mean_kl_yte-2*sqrt(gp_var_kl_yte),101,101),...
    %    'EdgeColor','None','FaceAlpha',0.1);
    hold off;
    if corr_length==0.01
        print(sprintf('figs/gp_sens_vsmall_corr_2d'),'-depsc2','-r300');
    elseif corr_length==1
        print(sprintf('figs/gp_sens_large_corr_2d'),'-depsc2','-r300');
    else
        print(sprintf('figs/gp_sens_%0.4d_corr_2d',randi(1000,1)),'-depsc2','-r300');
    end
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




