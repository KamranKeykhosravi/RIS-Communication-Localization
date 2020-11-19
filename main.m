
close all;clear all;clc;
%% generating ris profile
config = GetConfig();
rng(100)
crb = CRBcal(config);
crb = crb.PEBcalc; % calculating PEB
Rate = crb.getRate;

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
figure 
crbR = reshape(Rate.rate_rand,size(X));
contourf(X,Y,crbR,[.03,.1,.2,.5,1],'--')
xlabel('x (m)');ylabel('y(m)');
title('SE (bits/sec/Hz)');
cl = colorbar;


%------------------------------------------directional SE--------
figure 
crbR = reshape(min(Rate.rate_dir,15),size(X));
contourf(X,Y,crbR,[0:.2:1,1:10],'--')
hold on
xlabel('x (m)');ylabel('y(m)');
title('Rate (bits/sec/Hz)');
cl = colorbar;



