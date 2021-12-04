function [targets] = find_target_phase(EEG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% EEG: 1 x time points
% target has 5 type of elememnt 1=>0 degree 2=.90 degree 3=>180 degree
% 4=>270 degree 0=> otherwise
targets=zeros(size(EEG));
[~,loc2]=findpeaks(EEG);
% disp('loc2');
% disp(loc2);
[~,loc4]=findpeaks(-EEG);
% disp('loc4');
% disp(loc4);
peak_loc=sort(horzcat(loc2,loc4));
grad=gradient(EEG);
[~,loc1]=findpeaks(grad);
[~,loc3]=findpeaks(-grad);
for i=loc1
    targets(i)=1;
end
for i=loc2
    targets(i)=2;
end
for i=loc3
    targets(i)=3;
end
for i=loc4
    targets(i)=4;
end
targets=remove_incorrect_phase(targets,peak_loc);
end
function targets=remove_incorrect_phase(targets,peak_loc)              
    for i=1:numel(peak_loc)-1
        peak_type=targets(peak_loc(i));
        for j=peak_loc(i)+1:peak_loc(i+1)-1
            %3 must follow 2, 1 must follow 4. Otherwise, remove that
            %target.
            if targets(j)~=mod(peak_type+1,4)
                targets(j)=0;
            end
        end
    end
    %remove the incorrect phases after the last peak
    peak_type=targets(peak_loc(end));
    for i=peak_loc(end)+1:numel(targets)
        if targets(i)~=mod(peak_type+1,4)
            targets(i)=0;
        end
    end
    %remove the incorrect phases before the first peak
    peak_type=targets(peak_loc(1));
    for i=peak_loc(1)-1:-1:1
        %1 must go before 2; 3 must go before 4.  Otherwise, remove that
        %target.
        if targets(i)~=mod(peak_type-1,4)
            targets(i)=0;
        end
    end
        
end

