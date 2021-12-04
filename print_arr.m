function print_arr(matrix)
%PRINT_ARR Summary of this function goes here
%   Detailed explanation goes here
% matrix = magic(4) % example matrix
[mrows, ncols] = size(matrix);
outputstr = ['%' num2str(mrows) 'i '] ;% template for the string, you put your datatype here
outputstr = repmat(outputstr, 1, ncols); % replicate it to match the number of columns
% disp(outputstr);
outputstr = [outputstr '\n']; % add a new line if you want
fprintf(outputstr, matrix.'); % write it
end

