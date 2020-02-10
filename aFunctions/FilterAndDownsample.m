function [deci_data IBN Sig] = FilterAndDownsample(filter_coefficients, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, plot_psd, export_IBN, print_IBN, print_Sig)

%
% deci_data IBN Sig] = FilterAndDownsample(filter_coefficients, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, plot_psd, export_IBN, print_IBN, print_Sig)
%
% This function filter and downsdample the output patern from Delta-Sigma Modulator due to the number 
% of decimation stages and its coefficients. Plots the PSD before and after
% the decimation process. Export the estimated IBN before and after
% deciamtion process and the Signal amplitude after decimation.
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
%   plot_psd:               Optional flag for plotting the PSD before and after 
%                           decimation, 1 -> Plot or 0 -> don't Plot
%   export_IBN:             Optinal flag to export vector holding the IBN
%                           before and after decimation, 1 -> export, 0 ->
%                           don't export, [IBN_before_decimation  IBN_after_decimation]
%   print_IBN:              Optinal flag to print the IBn before and after
%                           decimation in the command window, where
%                           1 -> print or 0 -> don't print
%   print_Sig:              Optinal flag to print the Signal peak value
%                           before and after decimation
%                           1 -> print or 0 -> don't print
%
%   deci_data:              Filtered output data
%   IBN:                    IBN before and after decimation respectively
%   Sig:                    Signal power before and after decimation respectively


filter_downsample = struct('filter_data', [], 'downsample_data', []);


for i = 1 : K,
%     if max(abs(filter_coefficients(i,:))) < 1,
        if i == 1,
            filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], sdm_data);
            filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
        else
            filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], filter_downsample(i-1).downsample_data);
            filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
        end
%     else
%         if i == 1,
%             filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], sdm_data./sum(filter_coefficients(i,1:filter_lengths(i))));
%             filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
%         else
%             filter_downsample(i).filter_data = filter(filter_coefficients(i,1:filter_lengths(i)), [1], filter_downsample(i-1).downsample_data./sum(filter_coefficients(i,1:filter_lengths(i))));
%             filter_downsample(i).downsample_data = downsample(filter_downsample(i).filter_data, M(i));
%         end  
%     end
end

deci_data = filter_downsample(K).downsample_data;

% Estimate the IBN before and after decimation
if (nargin == 10 | nargin == 11 | nargin == 12) & export_IBN == 1,
    Fb=Fs/(2*OSR);         
    sdmdata                         = plotFunction(sdm_data,'OSR',OSR,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs);
    IBN_DSM = sdmdata.P_IBN_dB;
    Sig_DSM = sdmdata.P_signal_dB;
    decidata                        = plotFunction(deci_data,'OSR',1,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'plot_fft',false,'fs',Fs/(prod(M)));
    IBN_Decimator = decidata.P_IBN_dB;
    Sig_Decimator = decidata.P_signal_dB;
    IBN = [IBN_DSM IBN_Decimator];
    Sig = [Sig_DSM Sig_Decimator];
else
    IBN = [nan nan];
end

if print_IBN == 1,
        fprintf('\n');
        fprintf('The IBN before decimation = %f\n The IBN after decimation = %f\n', IBN);
        if nargin == 12 & print_Sig == 1,
            fprintf('The Signal Peak before decimation = %f\n The Signal Peak after decimation = %f\n', Sig);
            fprintf('\n');
        end
        fprintf('\n');
end

% Plot PSD before and after decimation
% if (nargin == 9 | nargin == 10 | nargin == 11 | nargin == 12) & plot_psd == 1,
%     FIG = figure('Name', 'Power Spectral Density', 'NumberTitle' , 'off');
%     Np1 = TwoPowerN(length(sdm_data));
%     [spectr_psd1, freq_psd1] = GetSpectrum(sdm_data, Np1, Fs);
%     if min(spectr_psd1) < 1e-20,
%         semilogx(freq_psd1, 10*log10(spectr_psd1*Fs/Np1));
%     else
%         semilogx(freq_psd1, 10*log10(spectr_psd1));
%     end
%     hold on
%     Np = TwoPowerN(length(filter_downsample(K).downsample_data));
%     [spectr_psd, freq_psd] = GetSpectrum(filter_downsample(K).downsample_data, Np, Fs/prod(M));
%     if min(spectr_psd) < 1e-20,
%         semilogx(freq_psd, 10*log10(spectr_psd*(Fs/prod(M))/Np), 'r');
%     else
%         semilogx(freq_psd, 10*log10(spectr_psd), 'r');
%     end
%     axis auto
%     grid on
%     xlabel('Frequency - Hz');
%     ylabel('PSD - dB');
%     title(['IBN = ', num2str(IBN), ', Sig = ', num2str(Sig)]);
% end


% Plot PSD before and after decimation
if (nargin == 9 | nargin == 10 | nargin == 11 | nargin == 12) & plot_psd == 1,
    FIG = figure('Name', 'Power Spectral Density', 'NumberTitle' , 'off');
    plotFunction(sdm_data,'OSR',OSR,'fsig',Fsignal,'stats',true,'no_nz_bins',20,'fs',Fs);
    hold on
    %Np = TwoPowerN(length(filter_downsample(K).downsample_data));
    %[spectr_psd, freq_psd] = GetSpectrum(filter_downsample(K).downsample_data, Np, Fs/prod(M));
    plotFunction(filter_downsample(K).downsample_data,'OSR',OSR/prod(M),'fsig',Fsignal,'stats',true,'no_nz_bins',20,'fs',Fs/prod(M));
    axis auto
    grid on
    xlabel('Frequency - Hz');
    ylabel('PSD - dB');
    title(['IBN = ', num2str(IBN), ', Sig = ', num2str(Sig)]);
end
