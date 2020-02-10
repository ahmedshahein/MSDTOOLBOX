function [qcoeff Sn] = FilterCoefficientSensitivity(coeff, q, Fs, axis_ranges, plot_freq_response, plot_Sn)

%
% [qcoeff Sn] = FilterCoefficientSensitivity(coeff, q, Fs, axis_ranges, plot_freq_response, plot_Sn)
%
% This function calculate the quantized filter coefficients due to the 
% different bit widths specificed in the vector 'q'. It exports matrix
% holding the quantized filter coefficients and matrix holds the senstivity
% of the filter coefficients. It is for single filter only.
%
%   coeff:              Ideal filter coefficients
%   q:                  vector holds different bit widths for quantizing the filter
%                       coefficients, it can holds from 1 till 6 values
%   Fs:                 Sampling frequency
%   axis_ranges:        Set a specific axes range
%   plot_freq_response: Flag to enable (1) or disable (0) ploting the filter
%                       frequency response
%   plot_Sn:            Flag to enable (1) or disable (0) ploting the
%                       coefficient sensitivity

% Initialization for internal variables
qcoeff=[];
Hq=[];
Fq=[];
Magq=[];
Sn=[];
sort_Sn=[];

% Calculate the quantized coefficients due to the input quantized bit
% widths in vector n
for i = 1 : length(q),
    qcoeff(i, :) = QuantizedCoeff(q(i), 'fixed', coeff);
end

% The frequency response of the ideal input filter
[H F]=freqz(coeff, [1], 2^10, Fs);

% The frequency response of the quantized filters
for j = 1 : length(q),
    [Hq(j, :) Fq(j, :)]=freqz(qcoeff(j, :), [1], 2^10, Fs);
end

% The magnitude of the ideal input filter
Mag = 20*log10(abs(H));

% The magnitude of the quantized filters
for k = 1 : length(q),
    Magq(k, :) = 20*log10(abs(Hq(k, :)));
end

% Color vector to be used in loop plotting
color = ['r', 'k', 'm', 'g', 'c', 'y', 'd', '*', '+'];

% Plot the frequency response of the ideal and quantized filters
if nargin == 5 | nargin == 6 & plot_freq_response == 1,
    FIG = figure('Name', 'Filter Frequency Response', 'NumberTitle' , 'off');
    plot(F, Mag, 'b');
    hold on
    for l = 1 : length(q),
        plot(F, Magq(l, :), color(l));
    end
    grid on
    xlabel('Frequency-Hz');
    ylabel('Amplitude-dB');
    if axis_ranges ~= 0,
        axis(axis_ranges);
    else
        axis auto;
    end

    legend('Ideal');
    if length(q) == 1,
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit']);
    elseif length(q) == 2, 
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit']);
    elseif length(q) == 3, 
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'], ['Quantized ' num2str(q(3)) '-bit']);
    elseif length(q) == 4, 
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'], ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit']);
    elseif length(q) == 5, 
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'],  ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit'], ['Quantized ' num2str(q(5)) '-bit']);
    elseif length(q) == 6, 
        legend('Ideal', ['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'],  ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit'], ['Quantized ' num2str(q(5)) '-bit'], ['Quantized ' num2str(q(6)) '-bit']);
    end
end

% Estimate the sensitvity of each quantized filter compared to the ideal
% filter
for m = 1 : length(q),
    qcoeff_tmp = qcoeff(m,:);
    for o = 1 : length(coeff),
        Sn(m, o) =sqrt((qcoeff_tmp(o) - coeff(o))^2);
    end
end

% Sort the sensitvity in ascending form
for n = 1 : length(q),
    sort_Sn(n, :) = sort(Sn(n, :));
end

% Plot the sensitvity for the quantized coefficeints
if plot_Sn == 1,
    FIG = figure('Name', 'Coefficients Sensitvity', 'NumberTitle' , 'off');
    for p = 1: length(q),
        plot(sort_Sn(p, :), (color(p))); 
        hold on
    end
    grid on
    xlabel('Filter Length');
    ylabel('Filter Coefficient Sensitvity')

    if length(q) == 1,
        legend(['Quantized ' num2str(q(1)) '-bit']);
    elseif length(q) == 2, 
        legend(['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit']);
    elseif length(q) == 3, 
        legend(['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'], ['Quantized ' num2str(q(3)) '-bit']);
    elseif length(q) == 4, 
        legend(['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'], ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit']);
    elseif length(q) == 5, 
        legend(['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'],  ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit'], ['Quantized ' num2str(q(5)) '-bit']);
    elseif length(q) == 6, 
        legend(['Quantized ' num2str(q(1)) '-bit'], ['Quantized ' num2str(q(2)) '-bit'],  ['Quantized ' num2str(q(3)) '-bit'], ['Quantized ' num2str(q(4)) '-bit'], ['Quantized ' num2str(q(5)) '-bit'], ['Quantized ' num2str(q(6)) '-bit']);
    end
end