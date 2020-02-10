function TestQuantizedFilters(quantized_filter_coefficients, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, q, plot_psd, export_IBN, print_IBN, print_Sig)

%
% TestQuantizedFilters(quantized_filter_coefficients, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, q, plot_psd, export_IBN, print_IBN, print_Sig)
%
% This function test the quantization effect of the filter coefficients, by
% examining the its effect on IBn and Signal peak.
%
%   quantized_filter_coefficients:      Matrix of filter coefficients exported from 'decimator_coefficient_sensitivity' function, 
%                                       which represents the coefficients at each stage
%   filter_lengths:                     Vector of filter lengths exported from 'decimation_filters' function, which 
%                                       represents the filter length at each stage 
%   sdm_data:                           Array holds the Delta-Sigma Modulator data
%   Fs:                                 Sampling frequency
%   OSR:                                Over Sampling Ratio
%   Fsignal:                            Signal frequency
%   K:                                  Number of decimation stages
%   M:                                  Vector holds the decimation factor at each stage
%   plot_psd:                           Optional flag for plotting the PSD before and after 
%                                       decimation, 1 -> Plot or 0 -> don't Plot
%   export_IBN:                         Optinal flag to export vector holding the IBN
%                                       before and after decimation, 1 -> export, 0 ->
%                                       don't export, [IBN_before_decimation  IBN_after_decimation]
%   print_IBN:                          Optinal flag to print the IBn before and after
%                                       decimation in the command window, where
%                                       1 -> print or 0 -> don't print
%   print_Sig:                          Optinal flag to print the Signal peak value
%                                       before and after decimation
%                                       1 -> print or 0 -> don't print
%

% Filter and downsample the DSm output using the quantized filters in the
% previous step and pltting the PSD for each trial
for i = 1 : length(q),
    filter_coefficients_temporary = quantized_filter_coefficients((i-1)*K+1:i*K,:);
    [deci_data(:,i) IBN(i,:) Sig(i,:)] = FilterAndDownsample(filter_coefficients_temporary, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, plot_psd, export_IBN, print_IBN, print_Sig);
end

% Calculate the difference in IBN and Signal peak before and after
% quntization
for j = 1 : length(IBN),
    diff_IBN(j) = IBN(j,1) - IBN(j,2);
    diff_Sig(j) = Sig(j,1) - Sig(j,2);
end

[rows columns] = size(q);

if rows == 1,
    
    % Form a matrix holding the 
    % Quantization bit width                                    -- 1st column
    % IBN before quantization                                   -- 2nd column
    % IBN after quantization                                    -- 3rd column
    % Difference in IBN before and after quantization           -- 4th column
    % Signal peak before quantization                           -- 5th column
    % Signal peak after quantization                            -- 6th column
    % Difference in Signal peak before and after quantization   -- 7th column

    quatization_effect_on_IBN_Sig           = zeros(length(IBN), 7);
    quatization_effect_on_IBN_Sig(:,1)      = q;
    quatization_effect_on_IBN_Sig(:,2:3)    = IBN;
    quatization_effect_on_IBN_Sig(:,4)      = diff_IBN;
    quatization_effect_on_IBN_Sig(:,5:6)    = Sig;
    quatization_effect_on_IBN_Sig(:,7)      = diff_Sig;

    % Print the table in the command window
    fprintf('Q -- IBN before decimation -- IBN after decimation -- Difference in IBN -- Sig after decimation -- Sig before decimation -- Difference in Sig \n');
    for k = 1 :length(IBN),
        fprintf('%2.0f   %3.4f                 %3.4f                %3.4f              %3.4f                 %3.4f                  %2.4f \n', quatization_effect_on_IBN_Sig(k,:));
    end

else

    q_tmp = flipud(rot90(q));
    quatization_effect_on_IBN_Sig               = zeros(length(IBN), K+6);
    quatization_effect_on_IBN_Sig(:,1:K)        = q_tmp;
    quatization_effect_on_IBN_Sig(:,K+1:K+2)    = IBN;
    quatization_effect_on_IBN_Sig(:,K+3)        = diff_IBN;
    quatization_effect_on_IBN_Sig(:,K+4:K+5)    = Sig;
    quatization_effect_on_IBN_Sig(:,K+6)        = diff_Sig;
    if K == 2,
        % Print the table in the command window
        fprintf('Q 1st Stage -- Q 2nd Stage -- IBN before decimation -- IBN after decimation -- Difference in IBN -- Sig after decimation -- Sig before decimation -- Difference in Sig \n');
        for k = 1 :length(IBN),
            fprintf('%2.0f             %2.0f             %3.4f                 %3.4f                %3.4f              %3.4f                 %3.4f                  %2.4f \n', quatization_effect_on_IBN_Sig(k,:));
        end
    elseif K == 3,
        % Print the table in the command window
        fprintf('Q 1st Stage -- Q 2nd Stage -- Q 3nd Stage -- IBN before decimation -- IBN after decimation -- Difference in IBN -- Sig after decimation -- Sig before decimation -- Difference in Sig \n');
        for k = 1 :length(IBN),
            fprintf('%2.0f             %2.0f             %2.0f             %3.4f                 %3.4f                %3.4f              %3.4f                 %3.4f                  %2.4f \n', quatization_effect_on_IBN_Sig(k,:));
        end 
    elseif K == 4,
        % Print the table in the command window
        fprintf('Q 1st Stage -- Q 2nd Stage -- Q 2nd Stage -- Q 2nd Stage -- IBN before decimation -- IBN after decimation -- Difference in IBN -- Sig after decimation -- Sig before decimation -- Difference in Sig \n');
        for k = 1 :length(IBN),
            fprintf('%2.0f             %2.0f             %2.0f             %2.0f             %3.4f                 %3.4f                %3.4f              %3.4f                 %3.4f                  %2.4f \n', quatization_effect_on_IBN_Sig(k,:));
        end 
    end
end