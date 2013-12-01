%%
close all; clear all;

[X,Y] = meshgrid(linspace(-1,1,101),linspace(-1,1,101));
Z = exp(0.7*X+0.3*Y);
surf(X,Y,Z,'EdgeColor','none'); view(2);
set(gca,'FontSize',14);
xlabel('x_1'); ylabel('x_2');
axis square; xlim([-1,1]); ylim([-1,1]);
set(gca,'XTick',[-1 0 1]); set(gca,'YTick',[-1 0 1]);
colormap('Hot'); colorbar;
set(gca,'FontSize',14);
print('figs/expfun','-depsc2','-r300');