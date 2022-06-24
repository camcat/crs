function [] = plot_outputs(outname);
%function [] = plot_outputs(outname);

%This file contains the background rate:
l=load([outname '_LogLikelihood.dat']);
r0=l(3);

%Forecast in space and time:
fx=load([outname '_foremap.dat']);
ft=load([outname '_forecast.dat']);

%Background rate (r0) and tot number of events (n0) per cell:
r0=r0/length(fx);
T=range(ft(:,1))*(length(ft)+1)/length(ft);
n0=r0*T;

%The 6th column contains forecasted no. of events. Here I normalize it by background:
n_norm=fx(:,6)/n0;

% Plot spatial forecast.
% plot log10(x), with x=normalized number of events.
ops=struct('cloud',1, 'field',n_norm, 'fun', @(x) log10(x),'clim',[-3 3]);
m=mapme(fx,ops);
hcb=colorbar;
ylabel(hcb,'log_{10}(Normalized no. of events)');

% Plot temporal forecast: cumulative no. of events.
figure
plot(ft(:,1), cumsum(ft(:,2)),'LineWidth',2)
xlabel('Days since IssueTime')
ylabel('Cumulative no. of earthquakes')

