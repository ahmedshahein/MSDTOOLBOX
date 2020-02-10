function PlotPSD(data, Fs, color, Title)

%
% PlotPSD(data, Fs, color, Title)
%
% This function used to plot the PSD of lowpass sigma delta modulator 
% patterns.
%
%   data:   Pattern to be ploted
%   Fs:     Sampling frequency
%   color:  Color of the graph 'optional'
%   Title:  Title of the graph 'optional'

if nargin == 2,
    Np = TwoPowerN(length(data));
    [spectr_psd, freq_psd] = GetSpectrum(data, Np, Fs);
    semilogx(freq_psd, 10*log10(spectr_psd*Fs/Np));
else
    Np = TwoPowerN(length(data));
    [spectr_psd, freq_psd] = GetSpectrum(data, Np, Fs);
    semilogx(freq_psd, 10*log10(spectr_psd*Fs/Np), color);
end    

if nargin == 4 & Title ~= 0,
    title([Title]);
end
    
grid on
xlabel('Frequency - Hz');
ylabel('PSD - dB');