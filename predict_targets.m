function [prediction,target_intervals] = predict_targets(epochs,target_val,predict_len,ar_ord)
%FIND_PEAK_INTERVAL Summary of this function goes here
%   Detailed explanation goes here
%   target_val can be 1, 2, 3, 4, indicating 0, 90, 180, 270 deg
%   respectively
%   epochs: time x channel
targets=nan(size(epochs));
for i=1:size(targets,2)
    targets(:,i)=find_target_phase(epochs(:,i).');
end

targets=targets==target_val;
%chnnel number
M=size(epochs,2);
prediction=false(predict_len,M);
coefficients=nan(ar_ord,M);
target_intervals=cell(1,M);
n_intervals=nan(1,M);
for i=1:M
    target_arr=targets(:,i);
    target_interval=diff(find(target_arr));
    n_intervals(1,i)=numel(target_interval);
    target_intervals{1,i}=target_interval;
%     fprintf('channel %d\n',i);
%     disp(target_interval);
    a = aryule(target_interval, ar_ord);
    coefficients(:,i) = -1 * flip(a(:, 2:end));
end

min_n_intervals=min(n_intervals);
target_intervals=cellfun(@(x)x(end-min_n_intervals+1:end,1),target_intervals,...
    'UniformOutput',false);
target_intervals=cell2mat(target_intervals);%min_n_interval x M
% disp(size(coefficients));
%   predicted location of next target


next_loc=zeros(1,M);
while median(next_loc)<predict_len%stop when half of the channels exceed the prediction length
    next_interval=sum(coefficients .*target_intervals(end-ar_ord+1:end,:),1);
%     print_arr(target_intervals(end-ar_ord+1:end,:));
    mask=next_interval<0;
    next_interval(mask)=mean(target_intervals(end-ar_ord+1:end,mask),1);
    target_intervals(end+1,:)=next_interval;
    
    next_loc=next_loc+next_interval;
    
        
%     fprintf('next interval=')
%     print_arr(next_interval);
%     fprintf('next location=');
%     print_arr(next_loc);
    for i=1:M
        if round(next_loc(i))<=predict_len
            try
                prediction(round(next_loc(i)),i)=true;
            catch ME
                fprintf('index=%d\n',round(next_loc(i)));
                rethrow(ME);
            end
        end
    end
end 
end

