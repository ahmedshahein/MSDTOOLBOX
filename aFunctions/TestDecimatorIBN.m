function IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, filter_coefficients, filter_lengths, sdm_data)

%
% IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, filter_coefficients, filter_lengths, sdm_data)
%
% This function tests each decimation stage after quantization and
% optimizing the filter coefficients. It exports the IBN after each
% deciamtion stage.
%
%   Fs:                     Sampling frequency
%   OSR:                    Over Sampling Ratio
%   Fsignal:                Signal frequency
%   K:                      Number of decimation stages
%   M:                      Vector holds the decimation factor at each stage
%   filter_coefficients:    Matrix of filter coefficients exported from 'decimation_filters' function, 
%                           which represents the coefficients at each stage
%   filter_lengths:         Vector of filter lengths exported from 'decimation_filters' function, which 
%                           represents the filter length at each stage 
%   sdm_data:               Array holds the Delta-Sigma Modulator data
%

% Calculating the baseband frequency
% Fb = Fs/(2*OSR);

filter_downsample = struct('filter_data', [], 'downsample_data', []);

for i = 1 : K,
    if i == 1,
        filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], sdm_data);
        filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
        % IBN_norm(i) = LPIBN((filter_downsample(i).downsample_data)./(sum(filter_coefficients(i,:))), Fsignal, Fs/M(i), Fb, 0);        
        tmpdata1(i)                         = plotFunction((filter_downsample(i).downsample_data)./(sum(filter_coefficients(i,:))),'OSR',OSR/prod(M(1:i)),'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/(M(i)));
        IBN_norm(i) = tmpdata1(i).P_IBN_dB;
    else
        filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], (filter_downsample(i-1).downsample_data)./sum(filter_coefficients(i-1,:)));
        filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
        % IBN_norm(i) = LPIBN((filter_downsample(i).downsample_data)./(sum(filter_coefficients(i,:))), Fsignal, Fs/prod(M(1:i)), Fb, 0);  
        tmpdata2(i)                         = plotFunction((filter_downsample(i).downsample_data)./(sum(filter_coefficients(i,:))),'OSR',OSR/prod(M(1:i)),'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/(prod(M(1:i))));
        IBN_norm(i) = tmpdata2(i).P_IBN_dB;
    end
end