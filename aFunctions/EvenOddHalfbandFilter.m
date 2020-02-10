function hb = EvenOddHalfbandFilter(coeff, type)

%
% hb = EvenOddHalfbandFilter(coeff, type)
%
% This function assert 0 every second coefficient according to filter
% type, i.e. even or odd.
%   coeff:  filter coefficient
%   type:   filter type, either even or odd
%
% Example:
% coeff = [1 2e-12 3 7e-16 9 7e-16 3 2e-12 1];
% type  = 'odd'
% hb = EvenOddHalfbandFilter(coeff, type)
%    = [1 0 3 0 9 0 3 0 1]

for i = 1 : length(coeff),
    if type == 'odd',
        if mod(i,2) == 1,
            hb(i) = 0;
        else
            hb(i) = coeff(i);
        end
    else
        if mod(i,2) == 0,
            hb(i) = 0;
        else
            hb(i) = coeff(i);
        end      
    end
end