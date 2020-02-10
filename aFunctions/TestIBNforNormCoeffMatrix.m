function [IBN_without_norm IBN_norm] = TestIBNforNormCoeffMatrix(Fs, OSR, Fsignal, M, coeff_matrix, Q, data, Fs_stage)

%
% [IBN_without_norm IBN_norm] = TestIBNforNormCoeffMatrix(Fs, OSR, Fsignal, M, coeff_matrix, Q, data, Fs_stage)
%
% This functin used to estimate the IBN after single decimation stage and
% export the IBN for the same filter with and without normalized
% coefficients. Normalized coefficients means that the coefficients are
% represented in fixed point format using Q-bits.
%
%   Fs:                     Sampling frequency 
%   OSR:                    Over Sampling Rationes
%   M:                      Decimation factor at this stage
%   coeff:                  Filter coefficients 
%   Q:                      Quantization factor
%   data:                   Inout data for this stage
%

[rows columns] = size(coeff_matrix);

for i = 1 : rows,
    coeff_tmp = coeff_matrix(i, :);
    [IBN_without_norm(i) IBN_norm(i)] = test_IBN_for_norm_coeff(Fs, OSR, Fsignal, M, coeff_tmp, Q, data, Fs_stage);
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function [IBN_without_norm IBN_norm] = test_IBN_for_norm_coeff(Fs, OSR, Fsignal, M, coeff, Q, data, Fs_stage)

% Calculating the baseband frequency
Fb = Fs/(2*OSR);

if sum(coeff) < 2,
    % Filtering and downsampling using coefficients in floating point format
    filtered_data       = filter(coeff, [1], data);
    downsampled_data    = downsample(filtered_data, M);

    % Calculating IBN for coefficients in floating point format
    % IBN_without_norm = LPIBN(downsampled_data, Fsignal, Fs_stage/M, Fb, 0);
    tmpdata1                         = plotFunction(downsampled_data,'OSR',OSR/Fs_stage,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/Fs_stage);
    IBN_without_norm = tmpdata1.P_IBN_dB;

    % Representing coefficients in fixed point format
    coeff_norm = NormalizedCoeff(coeff, Q);

    % Filtering and downsampling using coefficients in fixedg point format
    filtered_data_norm      = filter(coeff_norm, [1], data);
    downsampled_data_norm   = downsample(filtered_data_norm, M);

    % Calculating IBN for coefficients in floating point format
    % IBN_norm = LPIBN(downsampled_data_norm./(sum(coeff_norm)), Fsignal, Fs_stage/M, Fb, 0);
    tmpdata2                         = plotFunction(downsampled_data_norm./(sum((coeff_norm))),'OSR',OSR/Fs_stage,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/Fs_stage);
    IBN_norm = tmpdata2.P_IBN_dB;
else
    IBN_without_norm = nan;
    % Filtering and downsampling using coefficients in fixedg point format
    filtered_data_norm      = filter(coeff, [1], data);
    downsampled_data_norm   = downsample(filtered_data_norm, M);

    % Calculating IBN for coefficients in floating point format
    % IBN_norm = LPIBN(downsampled_data_norm./(sum(coeff)), Fsignal, Fs_stage/M, Fb, 0);
    tmpdata3                         = plotFunction(downsampled_data_norm./(sum((coeff))),'OSR',OSR/Fs_stage,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/Fs_stage);
    IBN_norm = tmpdata3.P_IBN_dB;    
end

% End