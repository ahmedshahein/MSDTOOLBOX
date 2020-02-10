function y = RoundFraction(x)

%
% y = RoundFraction(x)
%
% This function round the fractions to nearest integer fraction.
%   Example:
%   y = round_fraction(0.156);
%   Y = 0.2

if x >=1,
    y = round(x);
else
    for i = 1 : 10,
        tmp = x * 10^i;
        if tmp < 1,
            i = i + 1;
        else
            y = round(tmp)/(10^i);
            return
        end
    end
end