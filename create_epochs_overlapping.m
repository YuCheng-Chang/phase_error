function [epochs, time] = create_epochs_overlapping(signal, fs)
% disp('create_epochs_overlapping');
% disp('signal');
% disp(size(signal));
nr_seconds = 2;
nr_samples = nr_seconds * fs;
epoch_overlap = 0.75*nr_samples; %2048 * 3.5;
[epochs, ~] = buffer(signal, nr_samples, epoch_overlap);
% y=buffer(x,n); The number of columns in y is ceil(L/(n-p)).L is the
% length of x

% remove some starting samples
epochs = epochs(:, ceil(nr_samples/(nr_samples-epoch_overlap)):end);

% create timeaxis in milliseconds
diff_t = 1000/fs; 
time = 0:diff_t:(nr_samples-1)*diff_t;
time = time-time(floor(numel(time)/2));

% linear detrending
epochs = detrend(epochs, 'constant');
disp('epochs');
disp(size(epochs));

    
end