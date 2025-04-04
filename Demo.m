%% 
clearvars
close all
clc

%% Dimension and ensemble size
n = 100;
ne = 20;

%% Make a covariance matrix
C     = getCov(n,[1 8], "Satellite");
[u,l] = getSVD(C);
X     = getSamples(ne,u,l);
Cml   = cov(X');

%% Use NICE
[Cov_NICE,Corr_NICE] = NICE(X,X,1);

figure
subplot(131)
imagesc(C)
colormap(brewermap(100,'RdBu'))
clim([-1 1])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'FontSize',20)
title('True covariance')

subplot(132)
imagesc(Cml)
colormap(brewermap(100,'RdBu'))
clim([-1 1])
set(gca,'FontSize',20)
set(gca,'XTick',[])
set(gca,'YTick',[])
title('Sample covariance')

subplot(133)
imagesc(Cov_NICE)
colormap(brewermap(100,'RdBu'))
clim([-1 1])
set(gca,'FontSize',20)
set(gca,'XTick',[])
set(gca,'YTick',[])
title('NICE')

set(gcf,'Color','w')
f = gcf;
f.Position = [680 605 1159 273];
