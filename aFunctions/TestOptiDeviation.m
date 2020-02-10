function coeff_matrix = TestOptiDeviation(filter_coefficients, filter_lengths, Fs, OSR, Fsignal, K, M, Q, Nb, sdm_data, deviation_values)

%
% coeff_matrix = TestOptiDeviation(filter_coefficients, filter_lengths, Fs, OSR, Fsignal, K, M, Q, sdm_data, deviation_values)
%
% This function used to test the effect of mixed coefficients due to
% acceptable deviation value in the coefficients. The test procedure as
% follow:
%   1. Calculate the IBN for this filter stage without any rounding in the
%      coefficients.
%   2. Create a matrix of new filter coefficietns holding a mixed coefficients
%      due to the acceptable deviation range and the number of bits in SPT
%      numbers.
%   3. Test the effect of each new filter coeffiecients 'mixed one' in the IBN
%      after this stage.
%   4. Export a table contains:
%       Column(1)   Column(2)           Column(3)                   Column(4)
%       N           Deviation Value     Count of rounded coeff.     IBN
%
%   filter_coefficients:    Matrix of filter coefficients exported from 'decimation_filters' function, 
%                           which represents the coefficients at each stage
%   filter_lengths:         Vector of filter lengths exported from 'decimation_filters' function, which 
%                           represents the filter length at each stage 
%   sdm_data:               Array holds the Delta-Sigma Modulator data
%   Fs:                     Sampling frequency
%   OSR:                    Over Sampling Ratio
%   Fsignal:                Signal frequency
%   K:                      Number of decimation stages
%   M:                      Vector holds the decimation factor at each stage
%

depth = length(deviation_values)+1;

filter_downsample = struct('filter_data', [], 'downsample_data', []);

for i = 1 : K,
    if i == 1,
        filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], sdm_data);
        filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
    else
        filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], filter_downsample(i-1).downsample_data);
        filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
    end
end

for j = 1 : K,
    if j == 1,
        coeff_matrix(1+(j-1)*depth:j*depth,:) = test_opti_deviation_single_stage(filter_coefficients(j,:), Fs, OSR, Fsignal, M(j), Q(j), Nb, sdm_data, deviation_values(j,:), prod(M(1:j)));
    else
        coeff_matrix(1+(j-1)*depth:j*depth,:) = test_opti_deviation_single_stage(filter_coefficients(j,:), Fs, OSR, Fsignal, M(j), Q(j), Nb, filter_downsample(j-1).downsample_data, deviation_values(j,:), (prod(M(1:j))));
    end
end


%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function coeff_matrix= test_opti_deviation_single_stage(coeff, Fs, OSR, Fsignal, M, Q, Nb, data, deviation_values, Fs_stage)

% This function used to test the effect of mixed coefficients due to
% acceptable deviation value in the coefficients. The test procedure as
% follow:
%   1. Calculate the IBN for this filter stage without any rounding in the
%      coefficients.
%   2. Create a matrix of new filter coefficietns holding a mixed coefficients
%      due to the acceptable deviation range and the number of bits in SPT
%      numbers.
%   3. Test the effect of each new filter coeffiecients 'mixed one' in the IBN
%      after this stage.
%   4. Export a table contains:
%       Column(1)   Column(2)           Column(3)                   Column(4)
%       N           Deviation Value     Count of rounded coeff.     IBN

[IBN_ideal IBN_ideal_norm] = TestIBNforNormCoeffMatrix(Fs, OSR, Fsignal, M, coeff, Q, data, Fs_stage);

[mixed_coeff_matrix count] = OptiMixedCoeffDeviation(coeff, Q, Nb, deviation_values);

[IBN_without_norm IBN_norm] = TestIBNforNormCoeffMatrix(Fs, OSR, Fsignal, M, mixed_coeff_matrix(:,3:length(mixed_coeff_matrix)), Q, data, Fs_stage);

table = [mixed_coeff_matrix(:,1) reshape(count,length(count),1) reshape(IBN_norm,length(IBN_norm),1) mixed_coeff_matrix(:,3:length(mixed_coeff_matrix))];

ideal_filter = [nan, nan, IBN_ideal_norm, zeros(1,length(coeff))];

coeff_matrix = [ideal_filter; table];

fprintf('\n');
fprintf('Deviation -- Count -- IBN\n');
for i = 1 : length(deviation_values)+1,    
    fprintf('%2.1f         %3.0f       %3.3f\n', coeff_matrix(i,1:3));
end
