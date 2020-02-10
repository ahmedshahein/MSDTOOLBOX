function y = SPT(x)

%
% y = SPT(x)
%
% This function used to convert any decimal number to its crossponding
% Signed Power of Two 'SPT'. We can also use Matlab commands, such as; s =
% HIST(x, max(x)), or TABULATE(x). But both of them for integer numbers
% only. It works for vectors and matrices.
%   Example:
%   y = SPT(53);
%   y = 64
%

[rows columns] = size(x);

if rows == 1,   
    for i = 1 : length(x),
        if x(i) > 0,
            y(i) = 2^round(log2(x(i)));
        elseif x(i) < 0,
            tmp = 2^round(log2(abs(x(i))));
            y(i) = -tmp;
        else
            y(i) = 0;
        end
    end
else
    for j = 1 : rows,
        x_row = x(j, :);
        for k = 1 : length(x_row),
            if x_row(k) > 0,
                y(j,k) = 2^round(log2(x_row(k)));
            elseif x_row(k) < 0,
                tmp = 2^round(log2(abs(x_row(k))));
                y(j,k) = -tmp;
            else
                y(j,k) = 0;
            end 
        end
    end
end