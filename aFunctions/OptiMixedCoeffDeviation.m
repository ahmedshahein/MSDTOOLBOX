function [mixed_coeff_matrix count] = OptiMixedCoeffDeviation(coeff, Q, Nb, deviation_values)

%
% [mixed_coeff_matrix count] = OptiMixedCoeffDeviation(coeff, Q, deviation_values)
%
% The function is just a loop iteration for the function titled
% "OptiMixedCoeff".

mixed_coeff_matrix = [];

    for j = 1 : length(deviation_values),
        deviation_value = deviation_values(j);
        [deviation_in_coeff norm_coeff mixed_coeff(j,:) count(j)] = OptiMixedCoeff(coeff, Q, Nb, deviation_value);        
        
        mixed_coeff_matrix(j,:) = [deviation_value count(j) mixed_coeff(j,:)];
    end
