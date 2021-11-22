%% data preprocessing
fs=1000;
PEAK_FREQUENCY_INTERVAL = [8 14];
HILBERTWIN=128;
A=readmatrix('eyeclose.csv');
A=A(:,3).';
T = table(A,'RowNames', {'subj1'},'VariableNames',{'data'});
row_index='subj1';
clear A;
% filter design method for phastimate (order and peak frequency is variable)
design_phastimate_filter = @(ord, freq, fs) designfilt('bandpassfir', ...
    'FilterOrder', ord, 'CutoffFrequency1', freq-1, 'CutoffFrequency2',...
    freq+1, 'SampleRate', fs, 'DesignMethod', 'window');
[epochs,time]= create_epochs_overlapping(T{'subj1','data'},fs);
[peak_frequency, peak_SNR] = estimate_SNR(epochs, fs, PEAK_FREQUENCY_INTERVAL);
%% a family of equivalelnt  filters to find the true phase
filter_objects = {};

for ord = [2 3 4 5] % FIR - windowed sinc
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
end
for ord = [3 4 5] % FIR - least squares (equiripple is similar)
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-4, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
end
for ord = [4 8 12] % IIR - butterworth
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, 'HalfPowerFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'butter')};
end
for ord = [4 6 8] % IIR - chebychev I
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, 'PassbandFrequency2', peak_frequency+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
end
for attenuation = [10 20] % IIR - elliptic
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+2, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
end    

[truephase_mean, truephase_variance, trueamp_mean, trueamp_variance] = truephase(epochs, filter_objects);


T=addvars(T,nan(height(T),size(epochs,1),size(epochs,2)),nan(height(T),size(epochs,1),size(epochs,2)),...
    nan(height(T),size(epochs,1),size(epochs,2)),nan(height(T),size(epochs,1),size(epochs,2)),...
    'NewVariableNames',{'epochs_truephase_mean','epochs_truephase_ang_var','epochs_trueamp_mean','epochs_trueamp_var'});

T(row_index,:).epochs_truephase_mean = reshape(truephase_mean,[1,size(truephase_mean)]);
T(row_index,:).epochs_truephase_ang_var = reshape(truephase_variance,[1,size(truephase_variance)]);

T(row_index,:).epochs_trueamp_mean = reshape(trueamp_mean,[1,size(trueamp_mean)]);
T(row_index,:).epochs_trueamp_var = reshape(trueamp_variance,[1,size(trueamp_variance)]);
clear truephase_mean truephase_variance trueamp_mean trueamp_variance;
%% genetic algorithm
filter_order_range = 100:250;
filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
for ord = filter_order_range
    filter_objects_by_order{ord} = design_phastimate_filter(ord, peak_frequency, fs);
end

bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
bounds_window = [400 750];
bounds_edge = [30 120];
bounds_ar_order = [5 60];

% the includemask allows optimizing for a subset of epochs
% it makes sense to exclude epochs that would also be excluded by the
% real-time system, e.g. if artifacts are detected so as to not optimize
% for noisy epochs that wouldn't result in a stimulus anyway

% subselect according to truephase variance
%includemask = T(row_index,:).epochs_truephase_angdev <= quantile(T(row_index,:).epochs_truephase_angdev, 0.5);

% subselect according to true amplitude
[epochs_midphase_mean, ~, epochs_midamp_mean, ~] = phastimate_truephase(epochs, filter_objects);
includemask = epochs_midamp_mean >= quantile(epochs_midamp_mean, 0.5);

[optimal_parameters, ga_output] = phastimate_optimize(epochs(1:ceil(end/2),includemask), ...
    epochs_midphase_mean(includemask), filter_objects_by_order, bounds_filter_order,...
    bounds_window, bounds_edge, bounds_ar_order, HILBERTWIN);
save('optimal_parameters.mat','optimal_parameters')
clear filter_objects_by_order filter_objects;
%% phase error vs time
ang_diff = @(x, y) angle(exp(1i*x)./exp(1i*y));
D = design_phastimate_filter(optimal_parameters.filter_order, peak_frequency, fs);
stop_time=1*fs;
plt_X=(0:stop_time)/fs;
plt_Y=zeros(size(plt_X));
for t=0:stop_time
    [estphase, estamp] = predict_phase(epochs(((-optimal_parameters.window_length+1):0)+ceil(end/2),:), ...
        D, optimal_parameters.edge, optimal_parameters.ar_order, [HILBERTWIN,stop_time],t);
    truephase=T{row_index,'epochs_truephase_mean'}(1,ceil(end/2)+t,:);
    estphase=reshape(estphase,1,[]);
    truephase=reshape(truephase,1,[]);
    phase_error=rad2deg(ang_diff(estphase,truephase));
    
    phase_error=mean(abs(phase_error),'all');
    plt_Y(1,t+1)=phase_error;
    if t<10
        disp(size(estphase));
        disp(size(truephase));
        fprintf('t=%4d,\nphase error=',t);
        disp(phase_error);
    end
    
end
figure;
plot(plt_X,plt_Y);