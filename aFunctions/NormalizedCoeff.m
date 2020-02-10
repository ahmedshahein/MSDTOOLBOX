function ncoeff = NormalizedCoeff(coeff, n, c)

% ncoeff = NormalizedCoeff(coeff, n, c)
% This function used to convert floating point format coefficients to
% signed integers in 2's complement format, for hardware implementation.

%       coeff   : floating point coeff., which is the output from remez commands
%       n       : number of bits used to represent the data
%       c       : single coeff.
%
%       ncoeff  : scaled coefficients

if nargin == 2,
    max_value = 2^(n-1)-1;
    ncoeff = round(coeff.*max_value);
else
    max_value = 2^(n-1)-1;
    ncoeff = round(c*max_value);
end