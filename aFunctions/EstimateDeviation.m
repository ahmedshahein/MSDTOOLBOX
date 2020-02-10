function EstimateDeviation(filter_coefficients, filter_lengths, Fs, K, Q, N, deviation_value)

%
% EstimateDeviation(filter_coefficients, filter_lengths, Fs, K, Q, N, deviation_value)
%
% This function export plots of:
%   Coefficients Deviation v.s. Coefficients
%   &
%   Filter frequency resoponse for normal filter and mixed coeffieicients
%   filter
%
%   filter_coeffieicients:  Filter coefficients for each deciamtion stage                           
%   filter_lengths:         Filter length for each deciamtion stage                           
%   Fs:                     Sampling frequency
%   Q:                      Quantization factor for each deciamtion stage
%   N:                      Number of bits for SPT rounding, N = [1 2]
%   deviation_value:        Random initial value for acceptable deviation
%                           value

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
    [deviation_in_coeff(i,:) norm_coeff(i,:) qcoeff(i,:) mixed_coeff(i,:) count(i,:)] = FilterCoeffDeviation(filter_coefficients(i,:), Q(i), N, deviation_value, Fs);   

    
    subplot(K,2,p(i));
    plot(deviation_in_coeff(i,1:filter_lengths(i)));
    grid on
    xlabel('Filter Coefficeints');
    ylabel('Deviation in Coefficients');
    title('Deviation in Filter Coefficients');

    subplot(K,2,p(i+K));
    filters = [norm_coeff(i,1:filter_lengths(i)); qcoeff(i,1:filter_lengths(i)); mixed_coeff(i,1:filter_lengths(i))];
    for j = 1:3,
        [H F]=freqz(filters(j,:), [1], 2^10, Fs);
        Mag = 20*log10(abs(H));
        plot(F, Mag, color(j));
        hold on
        grid on
        xlabel('Frequency');
        ylabel('Magnitude');
        legend('Normalized Coeff.','NPT Coeff.','Mixed Coeff.');
    end
    
end

% End