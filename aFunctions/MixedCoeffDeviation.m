function [mixed_coeff_matrix count]= MixedCoeffDeviation(coeff, Fs, Q, N_vector, deviation_values)

% [mixed_coeff_matrix count]= MixedCoeffDeviation(coeff, Fs, Q, N_vector,
% deviation_values)
% This function used to export a matrix of normalized filter coefficients.
% Those coefficients are the mixed coefficients, mixed between normal
% mormalized coefficients and rounded coefficients. Rounding strategy is
% upon to deviation in coefficients.
%
%   coeff:              Filter coeffieicnts
%   Q:                  Quantization factor for this filter stage
%   N_vector:           Vector holding the number of bits (N) for Signed Power of Two
%                       values (SPT/NPT)
%   deviation_values:   Vector holds a different values for acceptable range
%                       of deviation in coefficients. Deviation = Normalized(coeff.) - NPT(coeff.)
%
%   mixed_coeff_matrix: This is a table holding:
%                       Column(1)   Column(2)           Column(3)       Column(4,:)
%                       N           Deviation value     count           Mixed coeff.

mixed_coeff_matrix = [];
k = 1;

for i = 1 : length(N_vector),
    N = N_vector(i);
    for j = 1 : length(deviation_values),
        deviation_value = deviation_values(j);
        [deviation_in_coeff norm_coeff qcoeff mixed_coeff(k,:) count(k)] = FilterCoeffDeviation(coeff, Q, N, deviation_value, Fs);        
        
        mixed_coeff_matrix(k,:) = [N deviation_value count(k) mixed_coeff(k,:)];
       
        k = k + 1;
    end
end

% End