function y = LOG2_Width(x)

%
% y = LOG2_Width(x)
%
% This function is used to estimate the required bit-width for multipliers
% or adders.
% Example:
% x = [12589 6598 45 32]; 
% y = LOG2_Width(x);
% y = [14 13 6 5]
% where x is filter coefficients and y is the constant multiplier
% bit-width.

for i = 1 : length(x),
    y(i) = LOG2(x(i));
end

%%%%%%%%%%%%%%%%
% Sub-function %
%%%%%%%%%%%%%%%%
function y = LOG2(x)
if x ==0,
    y = 1;
else
     v=x;
     y=0;
     while v>1,
         y=y+1;
         v=v/2;
     end
end