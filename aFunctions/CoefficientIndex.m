function indecies = CoefficientIndex(OriginalCoefficients, MixedCoefficients)

%
% indecies = CoefficientIndex(OriginalCoefficients, MixedCoefficients)
%
% This fucntion used to export the indecies of the filter that has been
% changed during the rounding process to be able to estimate the saving in
% adders afterwards.
% Example:
% OriginalCoefficients  = [3 4 7 8 9 8 7 4 3]
% MixedCoefficients     = [2 4 6 8 9 8 6 4 2]
% indecies  = CoefficientIndex(OriginalCoefficients, MixedCoefficients)
%           = [1 3 7 9]

indecies = [];
j = 1;

for i = 1 : length(OriginalCoefficients),
    if OriginalCoefficients(i) ~= MixedCoefficients(i),
        indecies(j) = i;
        j = j + 1;
    else
        i = i + 1;
    end
end

    