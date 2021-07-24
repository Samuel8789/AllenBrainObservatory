%first addpath
clear;
close all;
%%

addpath(genpath('deconvolution_dep'));
DDFF=double(dff_traces);
Selectcell=DDFF(1,:);
figure; subplot(3,3,1);
cnt=1;
for i=[1 20 30 45 60 150 168 89 200]
    subplot(3,3,cnt);
    cnt=cnt+1;
    plot(TimeStmps,DDFF(i,:));
    drawnow()
end

%%
tau_rise = 2;
tau_decay = 10;

[g2,h2] = tau_c2d(tau_rise,tau_decay,1);   % generate discrete time constants

[ca_foopsi,cb,c1,~,~,spikes_foopsi] = constrained_foopsi(Selectcell,[],[],[]); %g2

%%


params.p = 2;
%params.g = g2;
params.sp = spikes_foopsi;   % pass results of foopsi for initialization (if not, they are computed)
params.c = ca_foopsi;
params.b = cb;
params.c1 = c1;
%params.sn = sg;
params.marg = 0;

SAMPLES = cont_ca_sampler(Selectcell,params);    %% MCMC   
plot_continuous_samples(SAMPLES,Selectcell(:));
M = plot_marginals(SAMPLES.ss,T,spikeRaster);