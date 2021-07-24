%% Data

raster=Events>1;

%% Parameters

samples_per_second=30; % In Hz
smoothing_frames=3; % in frame number
fixed_width_frames=20; % in frame nnumber
division_window=20; % in ms

samples_per_minute = 60*samples_per_second;
period=round(1000*1/samples_per_second);  % in ms
smooth_filter=smoothing_frames*period; % In integer multiples of period in ms
fixed_width=fixed_width_frames*period; % In integer multiples of period in ms

%% Firing Rates 
bin_window=200 % in ms
binning_step_frames=1
bin_step=binning_step_frames*period % in ms, it has to be at least teh same as 
bin=bin_window*samples_per_second/1000;
step=round(bin_step*samples_per_second/1000);

[c,f]=size(raster);
norm_frequencies=zeros([c f]);
frequencies=zeros([c f]);
for i=1:c
    n=f;
    n_sum=floor((n-bin)/step+1);
    spike_rate=zeros(1,n);
    for j=1:n_sum
        % Get spike rate
        ini=(j-1)*step+1;
        fin=(j-1)*step+bin;
        spike_rate(i,ini:(ini+step-1))=sum(raster(i,ini:fin));
    end

    frequencies(i,:)=spike_rate(i,:);
    SubAver_frequencies(i,:)=spike_rate(i,:)-mean(spike_rate(i,:));
    if(std(spike_rate(i,:)))
       norm_frequencies(i,:)=SubAver_frequencies(i,:)/std(spike_rate(i,:));
    end
end
 
%% Coactivity calculation
% Coactivity measures the sum of spikes for each frame. The objective is to
% determine which frames have highr ocactivity than expected by chance.

coactivity=sum(raster)';
smooth_coactivity=smooth(coactivity,smoothing_frames);
smooth_coactivity=smooth(smooth_coactivity,smoothing_frames);

% Remove Oscillations
percentage=0.01;
removed=smooth(coactivity,percentage,'loess'); % quadratic fit
smoothed=coactivity-removed;

% Normalize Coactivity
zscore_window=300;   % In Frames
F = length(coactivity);
    %With Window
    n_final = round(F/zscore_window);
    for i = 1:n_final
        inicio = zscore_window*(i-1)+1;
        fin = inicio+zscore_window-1;
        if fin>F
            fin = F;
        end
        SubMeanCoactivityBinned(inicio:fin) = coactivity(inicio:fin)-mean(coactivity(inicio:fin));
        zscoreCoactivityBinned(inicio:fin) = (coactivity(inicio:fin)-mean(coactivity(inicio:fin)))/std(coactivity(inicio:fin));
       
    end
    
    %Full Raster
    SubMeanCoactivityAll = coactivity-mean(coactivity);
    zscoreCoactivityAll = (coactivity-mean(coactivity))/std(coactivity);

%% Test the coactivity threshold

threshold=1.5 % In which units



%% Vector calculation



    
    
    
    
    
    
    
    
    
%% Plotting

Plot_Raster(raster)
Plot_Coactivity(zscoreCoactivityBinned,'Coactivity',threshold,samples_per_second);

