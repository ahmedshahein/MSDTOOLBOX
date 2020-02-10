function EstimateOptiDeviation(filter_coefficients, filter_lengths, Fs, K, Q, Nb, deviation_value)

%
% EstimateOptiDeviation(filter_coefficients, filter_lengths, Fs, K, Q, deviation_value)
%
% This function estimates the effect of filter coefficent optimzation at
% each decimation stage according to the deviation criteria. It exports a
% plot, in which each row of subplots represent a single deciamtion stage
% while the total number of subplot rows represent the total number of
% decimation stages. In each row there are two subplots, the first subplot
% represents the deviation in filter coefficeints while the second one
% represents the filter response before and after optimization.
%
%   filter_coefficients:    Matrix of filter coefficients exported from 'decimation_filters' function, 
%                           which represents the coefficients at each stage
%   filter_lengths:         Vector of filter lengths exported from 'decimation_filters' function, which 
%                           represents the filter length at each stage 
%   Fs:                     Sampling frequency
%   K:                      Number of decimation stages
%   Q:                      Quantization bit width
%   deviation_value:        Any arbitrary value for the accepted value of
%                           deviation factor between coefficient and its rounded value, i.e.
%                           coeff           = [53 67 98 28]
%                           rounded_coeff   = [48 64 96 32]
%                           deviation       = [5 3 2 4]
%                           if the deviation_value = 3, then only two factor will be rounded,
%                           coeff(2) and coeff(3)

color = ['b', 'm', 'r', 'g', 'k', 'c', 'y', 'd'];

if K==2,
    p = [[1 2];[3 4]];
elseif K == 3,
    p = [[1 2];[3 4];[5 6]];
elseif K == 4,
    p = [[1 2];[3 4];[5 6];[7 8]];
else
    fprintf('This is not a suitable number of decimation stages');
end    

FIG = figure('Name', 'Test Deviation Effect', 'NumberTitle' , 'off');
for i = 1 : K,
    [deviation_in_coeff(i,:) norm_coeff(i,:) mixed_coeff(i,:) count(i,:)] = OptiMixedCoeff(filter_coefficients(i,:), Q(i), Nb, deviation_value);

    subplot(K,2,p(i));
    plot(deviation_in_coeff(i,1:filter_lengths(i)));
    grid on
    xlabel('Filter Coefficients');
    title(['Deviation in Coefficients for Stage - ' num2str(i)]);

    subplot(K,2,p(i+K));
    filters = [DenormalizeCoeff(norm_coeff(i,1:filter_lengths(i))); DenormalizeCoeff(mixed_coeff(i,1:filter_lengths(i)))];
    for j = 1:2,
        [H F]=freqz(filters(j,:), [1], 2^10, Fs);
        Mag = 20*log10(abs(H));
        plot(F, Mag, color(j));
        hold on
        grid on
        xlabel('Frequency');
        ylabel('Magnitude');
        title(['Filter Response for Stage - ', num2str(i)]);
        legend('Normalized Coeff.', 'Mixed Coeff.');
    end    
end
