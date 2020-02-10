function matrix = PPDMatrix(coeff, M)

%
% matrix = PPDMatrix(coeff, M)
%
% This function used to form the polyphase matrix of coefficients according
% to the filter coefficients and the decimation factor for the same stage.
%
%   coeff:  Filter coefficients
%   M:      Decimationfactor for the same stage
%
% Example:
% coeff  = [1 3 5 7 9 11 13 15 31 31 15 13 11 9 7 5 3 1];
% M      = 3;
% matrix = PPDMatrix(coeff, M);
%        = [ [1     7    13    31    11     5];
%            [3     9    15    15     9     3];
%            [5    11    31    13     7     1] ]


rows = M;
columns = ceil(length(coeff)/rows);

index = rows * columns;

coeffVector = [zeros(1,(index-length(coeff))) coeff];

matrix = zeros(rows, columns);

for i = columns : -1: 1,
    for j = rows : -1: 1,
        matrix(j, i) = coeffVector(index);
        index = index - 1;
    end
end