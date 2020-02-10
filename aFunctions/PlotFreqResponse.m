function PlotFreqResponse(coeff, Fs, colour)

%
% PlotFreqResponse(coeff, Fs, color)
%
% This function used to plot filter frequency response, it is useful for looping 
%
% P.S. Plotted responses are denormalized, i.e. if scaled coefficeints are
% used the plotted response will have maximum gain of 0.

[rows columns] = size(coeff);
color = ['b', 'm', 'r', 'g', 'k', 'c', 'y', 'd'];

if rows == 1,
    if nargin == 2,
        [H F]=freqz(coeff./sum(coeff), [1], 2^10, Fs);
        Mag = 20*log10(abs(H));
        plot(F, Mag);
        grid on
        xlabel('Frequency');
        ylabel('Magnitude');
    else
        [H F]=freqz(coeff./sum(coeff), [1], 2^10, Fs);
        Mag = 20*log10(abs(H));
        plot(F, Mag, colour);
        grid on
        xlabel('Frequency');
        ylabel('Magnitude');
    end
else
    if nargin == 2,
        for i = 1 : rows,
            [H F]=freqz(coeff(i,:)./sum(coeff(i,:)), [1], 2^10, Fs);
            Mag = 20*log10(abs(H));
            plot(F, Mag, color(i));
            hold on
            grid on
            xlabel('Frequency');
            ylabel('Magnitude');
        end
    else
        for i = 1 : rows,
            [H F]=freqz(coeff(i,:)./sum(coeff(i,:)), [1], 2^10, Fs);
            Mag = 20*log10(abs(H));
            plot(F, Mag, color(i), colour);
            hold on
            grid on
            xlabel('Frequency');
            ylabel('Magnitude');
        end
    end
end