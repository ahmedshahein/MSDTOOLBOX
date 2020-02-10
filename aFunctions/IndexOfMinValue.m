function index = IndexOfMinValue(x)

%
% index = IndexOfMinValue(x)
%
% This function export the index of the minimum value in a vector
% Example:
%   x = [5 8 3 4 9 13 4];
%   index = index_of_min_value(x);
%   index = 3

for i = 1 : length(x),
    if min(x) == x(i),
        index = i;
    end
end