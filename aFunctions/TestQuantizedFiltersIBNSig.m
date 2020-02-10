function  quatization_effect_on_IBN_Sig = TestQuantizedFiltersIBNSig(quantized_filter_coefficients, filter_lengths, quantization_coefficients, sdm_data, Fs, OSR, Fsignal, K, M, q, plot_psd, export_IBN, print_IBN, print_Sig, IBN_penalty, Sig_penalty) %, auto_manual)

%
% quatization_effect_on_IBN_Sig = TestQuantizedFiltersIBNSig(quantized_filter_coefficients, filter_lengths, quantization_coefficients, sdm_data, Fs, OSR, Fsignal, K, M, q, plot_psd, export_IBN, print_IBN, print_Sig, IBN_penalty, Sig_penalty)
%
% This function test the quantization effect of the filter coefficients, by
% examining the its effect on IBN and Signal peak. This function might
% exhaust some time during calculations since the number of iterations
% depend on the length of the q matrix, e.g. 
% K = 3;
% q = [[15 14 13 12 11 10];
%      [15 14 13 12 11 10];
%      [15 14 13 12 11 10]];
% The number of iterations here will be 6^3 = length(q)^K.
% The testing procudure will be as follow;
% [15 15 15] 15-bit for the first filter, 15-bit for the second filter, 15-bit for the third filter 
% [15 15 14] 15-bit for the first filter, 15-bit for the second filter, 14-bit for the third filter 
% [15 15 13] 15-bit for the first filter, 15-bit for the second filter, 13-bit for the third filter 
%     .
%     .
%     .
% [15 14 15] 15-bit for the first filter, 14-bit for the second filter, 15-bit for the third filter 
% [15 14 14] 15-bit for the first filter, 14-bit for the second filter, 14-bit for the third filter 
% [15 14 13]
%     .
%     .
%     .
% and so on.
%
% The function will print on the command window a table holding the
% quantization bit width for each stage and the penalty in IBN and Signal
% peak, if the designer constrain the acceptable penalty in IBN and Signal,
% the table will hold only the values satisfy those two constarins.
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
%   IBN_penalty:                        Acceptable penalty in IBN after deciamtion, 
%                                       IBN_penalty = IBN'before decimation' - IBN'after deciamtion'
%                                       e.g.
%                                       -70 - (-66) = -4dB
%   Sig_penalty                         Acceptable penalty in the signal peak,
%                                       Sig_penalty = Signal peak'before deciamtion'-Signal peak'after deciamtion'
%                                       e.g.
%                                       -3dB - (-3.2) = 0.2dB
%   auto_manual                         Flag to calculate the optimized quantization factors, wether automatically or manually.
%                                       1 -> Auto, 0 -> Manual
%                                       Auto means that the function will export the
%                                       quantization factor correspondent to the
%                                       minimum IBN.
%                                       Manual means that the function will export a
%                                       table of the quantization factors that
%                                       satisfy the acceptable penalty range, and the
%                                       designer has to choose the q factor manually.
%

% Filter and downsample the DSm output using the quantized filters in the
% previous step and pltting the PSD for each trial
[rows columns] = size(quantized_filter_coefficients);
for i = 1 : rows/K,
    filter_coefficients_temporary = quantized_filter_coefficients((i-1)*K+1:i*K,:);
    [deci_data(:,i) IBN(i,:) Sig(i,:)] = FilterAndDownsample(filter_coefficients_temporary, filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, plot_psd, export_IBN, print_IBN, print_Sig);
    q_coeff_K_stages(i,1:K) = quantization_coefficients((i-1)*K+1:i*K);
end

% Calculate the difference in IBN and Signal peak before and after
% quntization
for j = 1 : size(IBN, 1), %length(IBN),
    diff_IBN(j) = IBN(j,1) - IBN(j,2);
    diff_Sig(j) = Sig(j,1) - Sig(j,2);
end

if nargin > 14,
    for l = 1 : size(IBN, 1), %length(IBN),
        if diff_IBN(l) >= IBN_penalty & diff_IBN(l) <= 0,
            if diff_Sig(l) >= -Sig_penalty & diff_Sig(l) <= Sig_penalty, % -Sig_penalty%%%%
            %if diff_Sig(l) <= Sig_penalty,
                accepted_diff_IBN(l)    = diff_IBN(l);
                accepted_diff_Sig(l)    = diff_Sig(l);
                accepted_index(l)       = l;
            else
                accepted_diff_IBN(l)    = 0;
                accepted_diff_Sig(l)    = 0;
                accepted_index(l)       = 0;
            end
        else
            accepted_diff_IBN(l)    = 0;
            accepted_diff_Sig(l)    = 0;
            accepted_index(l)       = 0;
        end
    end
end

[qrows qcolumns] = size(q);

if qrows == 1,
    
    % Form a matrix holding the 
    % Quantization bit width                                    -- 1st column
    % Difference in IBN before and after quantization           -- 2nd column
    % Difference in Signal peak before and after quantization   -- 3rd column

    quatization_effect_on_IBN_Sig           = zeros(length(IBN), 3);
    quatization_effect_on_IBN_Sig(:,1)      = q;
    quatization_effect_on_IBN_Sig(:,2)      = diff_IBN;
    quatization_effect_on_IBN_Sig(:,3)      = diff_Sig;

    % Print the table in the command window
    fprintf('Q -- IBN before decimation -- IBN after decimation -- Difference in IBN -- Sig after decimation -- Sig before decimation -- Difference in Sig \n');
    for k = 1 :length(IBN),
        fprintf('%2.0f   %3.4f                 %3.4f                %3.4f              %3.4f                 %3.4f                  %2.4f \n', quatization_effect_on_IBN_Sig(k,:));
    end

else
    if nargin > 14, % the penalty number are given as constrint
        q_tmp = q_coeff_K_stages;
        quatization_effect_on_IBN_Sig               = zeros(length(IBN), K+2);
        quatization_effect_on_IBN_Sig(:,1:K)        = q_tmp;
        quatization_effect_on_IBN_Sig(:,K+1)        = accepted_diff_IBN;
        quatization_effect_on_IBN_Sig(:,K+2)        = accepted_diff_Sig;
    else
        q_tmp = q_coeff_K_stages;
        quatization_effect_on_IBN_Sig               = zeros(length(IBN), K+2);
        quatization_effect_on_IBN_Sig(:,1:K)        = q_tmp;
        quatization_effect_on_IBN_Sig(:,K+1)        = diff_IBN;
        quatization_effect_on_IBN_Sig(:,K+2)        = diff_Sig;
    end
    
    if K == 2,
        fprintf('Q 1st Stage -- Q 2nd Stage -- Difference in IBN -- Difference in Sig \n');
    elseif K == 3, %& quatization_effect_on_IBN_Sig(k,K+1) ~= 0,
        fprintf('Q 1st Stage -- Q 2nd Stage -- Q 3rd Stage -- Difference in IBN -- Difference in Sig \n');
    elseif K == 4,% & quatization_effect_on_IBN_Sig(k,K+1) ~= 0,
        fprintf('Q 1st Stage -- Q 2nd Stage -- Q 3rd Stage -- Q 4th Stage -- Difference in IBN -- Difference in Sig \n');
    end
    
    for k = 1 : length(IBN),
        if K == 2 & quatization_effect_on_IBN_Sig(k,K+1) ~= 0,
            fprintf('%2.0f             %2.0f             %3.4f              %3.4f \n', quatization_effect_on_IBN_Sig(k,:));
        elseif K == 3 & quatization_effect_on_IBN_Sig(k,K+1) ~= 0,
            fprintf('%2.0f             %2.0f             %2.0f             %3.4f              %3.4f \n', quatization_effect_on_IBN_Sig(k,:));
        elseif K == 4 & quatization_effect_on_IBN_Sig(k,K+1) ~= 0,
            fprintf('%2.0f             %2.0f             %2.0f             %2.0f             %3.4f               %3.4f \n', quatization_effect_on_IBN_Sig(k,:));
        end
    end

end

% End