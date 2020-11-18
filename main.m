
close all;clear all;clc;
%% generating ris profile
config = GetConfig();
config.risPhaseMethod = "random"; % in the set {"random","directional"}
rng(100)
crb = CRBcal(config);
crb = crb.PEBcalc; % calculating PEB

%%   plots

% --------------- random, PEB
figure 
[X,Y] = meshgrid(crb.config.xvec,crb.config.yvec);
crbR = reshape(crb.PEB,size(X));
contourf(X,Y,crbR,10.^[-4:.5:0],'--')
set(gca,'ColorScale','log')
xlabel('x (m)');ylabel('y(m)');
title('PEB (m)');
cl = colorbar;


%---------------------------- random SE---------
crb.config.risPhaseMethod = "random";
Rate = crb.getRate;
figure 
[X,Y] = meshgrid(crb.config.xvec,crb.config.yvec);
crbR = reshape(Rate.rate,size(X));
contourf(X,Y,crbR,10.^[-3,-2,-1.5,-1,-.5,-.2,0],'--')
set(gca,'ColorScale','log')
xlabel('x (m)');ylabel('y(m)');
title('SE (bits/sec/Hz)');
cl = colorbar;


%------------------------------------------directional SE--------
crb.config.risPhaseMethod = "directional";
Rate = crb.getRate;
figure 
[X,Y] = meshgrid(crb.config.xvec,crb.config.yvec);
crbR = reshape(min(Rate.rate,15),size(X));
contourf(X,Y,crbR,[0:.2:1,1:10],'--')
hold on
xlabel('x (m)');ylabel('y(m)');
title('Rate (bits/sec/Hz)');
cl = colorbar;



