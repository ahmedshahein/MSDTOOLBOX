function [zerocoeff] = ZeroOutCoeff(coeff)

%
% [zerocoeff] = ZeroOutCoeff(coeff)
%
% This function to set every second coefficient to 0 except 
% the middle one useful for the design of halfband filters
% length of coeff. has to be uneven
%
% This function was developed by Niklas Lotze.
% lotze@imtek.uni-freiburg.de
%

middlecoeff = (length(coeff)+1)/2;
% check if middle coefficient has even index
if floor(middlecoeff/2) == middlecoeff/2
    coeff(2:2:middlecoeff-1) = 0;
    coeff(middlecoeff+2:2:length(coeff)) = 0;
else
    coeff(1:2:middlecoeff-1) = 0;
    coeff(middlecoeff+2:2:length(coeff)) = 0;
end
zerocoeff=coeff;