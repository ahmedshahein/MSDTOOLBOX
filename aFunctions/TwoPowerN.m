function y = TwoPowerN(number)

% y = TwoPowerN(number)
% The name comes from 2^n. The functions output is the maximal power of 2 which
% fits in number
%
% This function was developed by Niklas Lotze.
% lotze@imtek.uni-freiburg.de
%

y = power(2,floor(log2(number)));