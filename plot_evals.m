%%
% Plot the singular values
close all; clear all;
pde_data = get_pde_data();

T=load('randos.mat'); X = T.X; clear T;
[M,m] = size(X);

filename='experiment1.mat';
[U1,kl_sv1] = get_kl_bases(1,m,pde_data,filename);
[f1,G1] = get_pde_solutions(X,U1,pde_data,filename);
[~,Sig1,W1] = svd(G1,'econ');
lambda1 = (1/M)*diag(Sig1).^2;
W1 = -W1(:,1:2);

[kl_sv1(1:5).^2/kl_sv1(1)^2 lambda1(1:5)/lambda1(1)]

filename='experiment2.mat';
[U0,kl_sv0] = get_kl_bases(0.01,m,pde_data,filename);
[f0,G0] = get_pde_solutions(X,U0,pde_data,filename);
[~,Sig0,W0] = svd(G0,'econ');
lambda0 = (1/M)*diag(Sig0).^2;
W0 = -W0(:,1:2);

[kl_sv0(1:5).^2/kl_sv0(1)^2 lambda0(1:5)/lambda0(1)]


%%
close all;
figure(1);
plot(1:m,W0(:,1),'bx',...
    1:m,W1(:,1),'ro',...
    'MarkerSize',12,'LineWidth',2);
set(gca,'FontSize',14);
grid on; axis square; xlim([0 m+1]);
xlabel('Index');
legend('\beta=0.01','\beta=1');
print(sprintf('figs/asm_eig1'),'-depsc2','-r300');

figure(2);
plot(1:m,W0(:,2),'bx',...
    1:m,W1(:,2),'ro',...
    'MarkerSize',12,'LineWidth',2);
set(gca,'FontSize',14);
grid on; axis square; xlim([0 m+1]);
xlabel('Index');
legend('\beta=0.01','\beta=1');
print(sprintf('figs/asm_eig2'),'-depsc2','-r300');
