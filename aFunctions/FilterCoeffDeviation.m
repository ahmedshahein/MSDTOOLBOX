function [deviation_in_coeff norm_coeff nptcoeff mixed_coeff count] = FilterCoeffDeviation(coeff, Q, N, deviation_value, Fs, plot_deviation, plot_freq)

%
% [deviation_in_coeff norm_coeff nptcoeff mixed_coeff count] = FilterCoeffDeviation(coeff, Q, N, deviation_value, Fs, plot_deviation, plot_freq)
%
% This function used to calculate the deviation in filter coefficients
% between the normalized value and the rounded value wether it is SPT or
% NPT. It exports the deviation values betweent the normalized and rounded
% coefficients, the normalized coeff., the rounded coefficients, and the
% mixed coeff. which is a mix between the normal and rounded coeff. due to
% the acceptable deviation value. A plot of the frequency response is
% plotted the 3 filter responses.
%
%   coeff:              Filter coefficients
%   Q:                  Quantization factor
%   N:                  Roundeing factor
%                           N = 1 -> SPT  
%                           N = 2 -> NPT for N = 2
%                           N = 3 -> NPT for N = 3
%   deviation_value:    Accepted deviation value 
%   Fs:                 Sampling frequency
%   plot_deviation:     Flag to plot the deviation, 1 -> plot
%   plot_frequency:     Flag to plot the frequency response for each
%                       new rounded filter
 
count = 0;

norm_coeff = NormalizedCoeff(coeff,Q); 

if N == 1,
    SPT_filter = SPT(NormalizedCoeff(coeff, Q)); 
    for j = 1 : length(coeff),
        deviation_in_coeff(j) = norm_coeff(j) - SPT_filter(j);
        if abs(deviation_in_coeff(j)) <= deviation_value & coeff(j) ~= 0,
            mixed_coeff(j) = SPT(NormalizedCoeff(coeff, Q, coeff(j))); 
            count = count + 1;
        else
            mixed_coeff(j) = NormalizedCoeff(coeff, Q, coeff(j));
        end
    end
    nptcoeff = SPT_filter;
elseif N == 2,
    N2PT_filter = NPT(NormalizedCoeff(coeff, Q),N,Q); 
    for k = 1 : length(coeff),
        deviation_in_coeff(k) = norm_coeff(k) - N2PT_filter(k);
        if abs(deviation_in_coeff(k)) <= deviation_value & coeff(k) ~= 0,
            mixed_coeff(k) = N2PT_filter(k); 
            count = count + 1;
        else
            mixed_coeff(k) = NormalizedCoeff(coeff, Q, coeff(k));
        end
    end
    nptcoeff = N2PT_filter;
elseif N == 3,
    N3PT_filter = NPT(NormalizedCoeff(coeff, Q),N,Q); 
    for l = 1 : length(coeff),
        deviation_in_coeff(l) = norm_coeff(l) - N3PT_filter(l);
        if abs(deviation_in_coeff(l)) <= deviation_value & coeff(l) ~= 0,
            mixed_coeff(l) = N3PT_filter(l); 
            count = count + 1;
        else
            mixed_coeff(l) = NormalizedCoeff(coeff, Q, coeff(l));
        end
    end
    nptcoeff = N3PT_filter;
end 

deviation_in_coeff = abs(deviation_in_coeff);

if  nargin > 5 & plot_deviation == 1, 
    FIG = figure('Name', 'Deviation in Filter Coefficients', 'NumberTitle' , 'off');
    plot(abs(deviation_in_coeff));
    grid on
    xlabel('Filter Coefficeints');
    ylabel('Deviation in Coefficients');
    title('Deviation in Filter Coefficients');
end

if nargin > 6 & plot_freq == 1,
    FIG = figure('Name', 'Filters Frequency Responses', 'NumberTitle' , 'off');
    filters = [norm_coeff; mixed_coeff];
    PlotFreqResponse(filters, Fs)
    title('Filter Response');
    legend('Normalized Coeff.','NPT Coeff.','Mixed Coeff.');
    hold on
end

% End