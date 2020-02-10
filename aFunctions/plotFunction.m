function data = plotFunction(varargin)

%data = plotFunction(yout,parameter)
%parameter.name = value
%
%data = plotFunction(yout,...,'name',value,...)
%
%Optionen f�r die Simulationsparameter
%Np  -------------------------------------- Die gesamte Punkte des PSDs
%segment_length  -------------------------- segment_length = Np
%fs  -------------------------------------- Samplefrequenz
%OSR  ------------------------------------- Oversampling Rate
%fsig  ------------------------------------ Signal Frequenzen [1 Hz 5 Hz ...]
%FS  -------------------------------------- Full Scale des Quantisierers [2]
%
%Optionen f�r das Windowing der FFT
%overlap_no  ------------------------------ Anzahl der �berlappungen
%segment_no  ------------------------------ Anzahl der Segmente [1]
%overlap_string  -------------------------- Prozentualer Anteil der �berlappung [50]
%scale  ----------------------------------- Sin Scalin 'sin' oder 'NBW' Scaling [sin,NBW]
%
%Optionen f�r die Berechnung
%no_nz_bins  ------------------------------ Anzahl der Bins zur Berechnung des Signales [5]
%f0 --------------------------------------- Mittenfrequenz [0 Hz]
%
%Optionen f�r die Darstellung
%plot_fft --------------------------------- Ausgabe der FFT [true|flase]
%new_fig ---------------------------------- Erstellen eines neuen Fensters [true|false]
%label  ----------------------------------- Erstellen der Labels [true|false]
%title ------------------------------------ Bild�berschrift [true|false]
%normfs ----------------------------------- Frequenz auf fs normieren [true|false]
%stats ------------------------------------ Erstellt ein Textfeld mit IBN [false|true]
%
%Falls eine Option nicht definiert wird, werden die Standardeinstellungen verwendet.
%
% Fr�here Optionen (nicht mehr unterst�tzt)
%f_start  --------------------------------- Startfrequenz f�r das IBN [1]
%
% This function is developed by:
%   Copyright 2007, 2008 Alexander Buhmann, Michael Maurer, Matthias Keller
%
% For any questions or permission for reuse contact:
%   michael.maurer@imtek.de
%
%   This file is part of DISCO Toolbox.

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Number of output arguments
nout    = nargout;
error(nargoutchk(0,5,nout))

% Checking arguments
data    = ChkArgs(varargin,nin);

% Calculate spectrum
data    = PowerSpec(data);

% Calculate signal power
data    = SignalPower(data);

% Calculate IBN
data    = EstimatingIBN(data);

% Plot function
if data.plot_fft
    plotPSD(data)
end

if data.info.verbose
    fprintf(data.info.txt)
end

function data = PowerSpec(data)
% Extract parameters
yout            = data.yout;
segment_length  = data.segment_length;
overlap_no      = data.overlap_no;
fs              = data.fs;

%W=blackmanharris(segment_length);
W=blackman(segment_length,'periodic');
%W=blackman(segment_length);

% Anzahl der Bins neben dem Hauptbin (einseitig), auf die Signalenergie auf Grund des verwendeten Fensters verteilt wird
% berechne Einseitiges Leistungsdichtespektrum in dBFS (normiert auf
% Leistung eines Sinussignals mit halber FS Amplitude multipliziert mit W(O) (1-Norm, FFT von verwendetem Window f��r f=0)) oder in dB (normiert
% auf Leistung des verwendeten Windows)
% ACHTUNG: Leistung des Window ~ Np (Spektrum um Np nach unten geshiftet, weil FFT^2 ~ Np^2 --> Multiplikation mit 1/Np notwendig)
%
% Window Skalierung (Matlab Default) --> Leistung von weissem Rauschen unabh?ngig von
% verwendetem Fenster sowie Np direkt im Spektrum ablesbar, sofern Rauschleistung konstant
% (NICHT der Fall f��r Band-Limited-White-Noise (Rauschleitung pro Bin = Noise_Power_Value*fs --> gr?sseres fs liefert mehr Rauschen)
% Um gesamte Rauschleitung konstant zu halten: Noise_Power_Value*1/fs
[ys,f]          = pwelch(yout,W,overlap_no,segment_length,fs);
data.sp_psd     = ys;
data.Np         = 2*length(ys); %Two sided spectrum
data.f          = f;
data.W          = W;

function data = SignalPower(data)
% Extract parameters
sp_psd          = data.sp_psd;
bins_signal     = data.bins_signal;
fs              = data.fs;
Np              = data.Np;
FS              = data.FS;

% Leistung im Hauptbin
P_signal        = sum(sp_psd(bins_signal))*fs/Np;
data.P_signal   = P_signal;

%Berechne Gesamtleistung in dB
data.P_signal_dBFS = 10*log10(P_signal/((FS/2)^2/2));
data.P_signal_dB = 10*log10(P_signal);


function data = EstimatingIBN(data)
% Extract parameters
sp_psd_noise    = data.sp_psd;
bins_signal     = data.bins_signal;
bins_fband      = data.bins_fband;
fs              = data.fs;
Np              = data.Np;
FS              = data.FS;

% Berechne aufsummiertes Rauschen f�r jeden Bin
try
    IBN_psd         = sp_psd_noise(bins_fband);
    IBN             = cumsum(IBN_psd)*fs/(Np-length(bins_signal));
    data.IBN        = IBN;
    data.IBN_psd    = IBN_psd;
    data.P_IBN_dBFS = 10*log10(IBN(end)/((FS/2)^2/2));
    data.P_IBN_dB   = 10*log10(IBN(end));
    data.SNDR_dB    = data.P_signal_dBFS-data.P_IBN_dBFS;
catch
    data.IBN        = [];
    data.IBN_psd    = [];
    data.P_IBN_dBFS = -Inf;
    data.P_IBN_dB   = -Inf;
    data.SNDR_dB    = +Inf;
end



function plotPSD(data)
sp_psd          = data.sp_psd;
IBN             = data.IBN;
new_fig         = data.new_fig;
normfs          = data.normfs;
scale           = data.scale;
bins_fband      = data.bins_fband;
bins_signal     = data.bins_signal;
str_title       = data.str_title;
flabel          = data.flabel;
f               = data.f;
fs              = data.fs;
Np              = data.Np;
FS              = data.FS;
stats           = data.stats;
W               = data.W;

% New figure
if new_fig
    fprintf('Creating new figure\n')
    figure;
    clf;
end

% normalize fs
if normfs
    f = f / fs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSD Plot mit dBFS/bin %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIER FEHLT NOCH NORMIERUNG AUF FS BEI JEDEM SEMILOGX mit sp_psd--> ...*1/(FS/2)^2/2
if strcmpi(scale,'sin')
    % Spektrum in dBFS/bin OHNE Signal
    semilogx(f,10*log10(sp_psd*fs/Np/((FS/2)^2/2)),'r');
    hold on;
    % Zeige aufsummiertes Rauschen
    semilogx(f(bins_fband),10*log10(IBN/((FS/2)^2/2))); % --> Normierung hier ebenfalls n�tig ??
    % Signalpeak hervorheben
    semilogx(f(bins_signal),10*log10(sp_psd(bins_signal)*fs/Np/((FS/2)^2/2)),'g*-');
    grid on;
    ylabel('PSD [dBFS/bin]');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PSD Plot mit dBFS/NBW %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % --> s. Schreier, liefert PSD, bei dem Signalleistung in dBFS direkt im Hauptbin ablesbar ist, Nebenbins werden NICHT ben�tigt
% % Grund: Skalierung der FFT durch verwendetes Fenster und somit der Signalamplitude bei Signalfrequenz wird aufgehoben
if strcmpi(scale,'nbw')
    % Spektrum in dBFS/NBW OHNE Signal
    semilogx(f,10*log10(sp_psd*fs*(W'*W)/(norm(W,1)^2*(FS/2)^2/2)),'r');
    hold on;
    % Zeige aufsummiertes Rauschen
    semilogx(f(bins_fband),10*log10(IBN/((FS/2)^2/2))); % --> Normierung hier ebenfalls n�tig ??
    % Signalpeak hervorheben
    semilogx(f(bins_signal),10*log10(sp_psd(bins_signal)*fs*(W'*W)/(norm(W,1)^2*(FS/2)^2/2)),'g*-');
    grid on;
    ylabel('PSD [dBFS/NBW]');
end

% Plot labels & text
if normfs
    xlabel('f [f_s]');
    xlim([1/Np 1/2]);
    fs = 1;
else
    xlabel('f [Hz]');
    xlim([fs/Np fs/2]);
end
title(str_title)

if stats
    text(1.5*fs/Np, -100, {['IBN = ',num2str(data.P_IBN_dBFS,'%6.2f'),' [dBFS];'],['P_{sig} = ',num2str(data.P_signal_dBFS,'%6.2f'),' [dBFS]']}, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
end

if flabel
    if (segment_no == 1)
        text(1.5*fs/Np, -10, {'Blackman-Harris window', ['2\^' num2str(floor(log2(Np))) ' pt FFT'], [num2str(segment_no) ' segment, no overlap' ]}, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
    else
        text(1.5*fs/Np, -10, {'Blackman-Harris window', ['2\^' num2str(floor(log2(Np))) ' pt FFT'], [num2str(segment_no) ' segments, ' num2str(overlap_string) '% overlap']}, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
    end
end

function data = ChkArgs(args,ni)
if ni == 0
    fprintf('\ndata = plotFunction(yout,parameter)\nparameter.name = value\n\ndata = plotFunction(yout,...,name,value,...)\n\nOptionen f�r die Simulationsparameter\nNp  -------------------------------------- Die gesamte Punkte des PSDs\nsegment_length  -------------------------- segment_length = Np\nfs  -------------------------------------- Samplefrequenz\nOSR  ------------------------------------- Oversampling Rate \nfsig  ------------------------------------ Signal Frequenz\nFS  -------------------------------------- Full Scale des Quantisierers\n\nOptionen f�r das Windowing der FFT\noverlap_no  ------------------------------ Anzahl der �berlappungen\nsegment_no  ------------------------------ Anzahl der Segmente\noverlap_string  -------------------------- Prozentualer Anteil der �berlappung\nscale  ----------------------------------- Sin Scalin sin oder NBW Scaling\n\nOptionen f�r die Berechnung\nf_start  --------------------------------- Startfrequenz f�r das IBN\nno_bins  --------------------------------- Anzahl der Bins zur Berechnung des Signales\n\nOptionen f�r die Darstellung\nplot_fft --------------------------------- Ausgabe der FFT [True Flase]\nnew_fig ---------------------------------- Erstellen eines neuen Fensters [True False]\nlabel  ----------------------------------- Erstellen der Labels [True False]\ntitle ------------------------------------ Bild�berschrift\nnormfs ----------------------------------- Frequenz auf fs normieren [True False]\n\nFalls eine Option nicht definiert wird, werden die Standardeinstellungen verwendet.\n')
    error('Syntax error::at least one argument required');
elseif ~isa(args{1},'numeric')
    error('Syntax error:First argument is not numeric')
elseif ~isvector(args{1})
    error('Syntax error::First argument is not an array')
else
    data.yout = args{1};
end

if ni==2
    if isa(args{2},'struct')
        data        = args{2};
        data.yout   = args{1};
    else
        error('Syntax error:command called with two arguments, where non of them is a structure')
    end
else
    % Check syntax {yout,'str_1',var_1,'str_2',var_2,...}
    for i = 2:2:ni-1
        if ~isa(args{i},'char')
            error('Syntax error::string assumed')
        end
    end
    data = ChkString(data,args,'fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'normfs','logical','norm fs is not logical');
    data = ChkString(data,args,'scale','char','PSD scaling is not a string');
    data = ChkString(data,args,'Np','numeric','number of points is not an integer');
    data = ChkString(data,args,'segment_no','numeric','segment number is not an integer');
    data = ChkString(data,args,'OSR','numeric','OSR is not numeric');
    data = ChkString(data,args,'fsig','numeric','signal frequency is not a number');
    data = ChkString(data,args,'FS','numeric','full scale is not numeric');
    data = ChkString(data,args,'new_fig','logical','boolean expression assumed');
    data = ChkString(data,args,'plot_fft','logical','boolean expression assumed');
    data = ChkString(data,args,'f0','numeric','center frequency is not numeric');
    data = ChkString(data,args,'stats','logical','stats is not logical');
    data = ChkString(data,args,'flabel','logical','flabel is not logical');
    data = ChkString(data,args,'no_nz_bins','numeric','bins is not numeric');
    data = ChkString(data,args,'verbose','logical','verbose is not true or false');
    data = ChkString(data,args,'f_start','numeric','sample frequency is not numeric');
    
end

data = ChkStruct(data);


function data = ChkString(data,args,string,type,err_msg)
pos = find(strcmp(args,string),1);
if isempty(pos)
elseif ~isa(args{pos+1},type) || (strcmp('numeric',type) && (any(isnan(args{pos+1})) || any(isinf(args{pos+1}))))
    error(['Syntax error::',err_msg])
else
    %        data=setfield(data,string,args{pos+1});
    data.(string)=args{pos+1};
end


function data = ChkStruct(data)
data.info.txt        = '\n';
if ~isfield(data,'normfs')
    data.normfs = 0;
end

if ~isfield(data,'scale')
    data.scale = 'sin';
end

if ~isfield(data,'Np')
    data.Np = length(data.yout);
end

if ~isfield(data,'segment_no')
    data.segment_no = 1;
elseif mod(data.segment_no,2)
    data.info.txt=[data.info.txt,'Warning:signal is not located on a bin\n'];
end

if ~isfield(data,'fs')
    data.fs = 1;
end

if ~isfield(data,'OSR')
    data.OSR = 50;
end

if ~isfield(data,'fsig')
    data.fsig = data.fs/(8*2*data.OSR);
end

if ~isfield(data,'FS')
    data.FS = 2;
end

if isfield(data,'f_start') && isfield(data,'f0') && data.f0~=0
    error('Start frequency is not supported for bandpass-systems');
end

% if isfield(data,'f_start') == 1
%     data.bin_start   = ceil(data.f_start*data.Np/(data.fs*data.segment_no));
% else
%     data.bin_start = 1;
% end

if ~isfield(data,'overlapRatio')
    data.overlapRatio = 0.5;
end

if ~isfield(data,'plot_fft')
    data.plot_fft     = true;
end

if ~isfield(data,'no_nz_bins')
    data.no_nz_bins = 4;
end

if ~isfield(data,'new_fig')
    data.new_fig     = false;
end

if ~isfield(data,'flabel')
    data.flabel      = false;
end

if ~isfield(data,'str_title')
    data.str_title   = '';
end

if ~isfield(data,'f0')
    data.f0   = 0;
end

if ~isfield(data,'stats')
    data.stats      = false;
end

if ~isfield(data,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.verbose;
    data = rmfield(data,'verbose');
end


% Bestimmung der Bins
fsigs                = abs(data.fsig-data.fs*round(data.fsig/data.fs));
bin_signal           = ceil(fsigs*data.Np./(data.fs*data.segment_no))+1;
bin_f0               = ceil(data.f0*data.Np/(data.fs*data.segment_no))+1;
if bin_f0 == 1
    bin_fband        = ceil(data.Np/(data.segment_no*2*data.OSR))+1;
else
    bin_fband        = ceil(data.Np/(data.segment_no*4*data.OSR))+1;
end

if isfield(data,'f_start')
    bin_fstart      = ceil(data.f_start*data.Np/(data.fs*data.segment_no));
    bin_f0          = ceil((bin_fband+bin_fstart)/2);
    bin_fband       = bin_f0-bin_fstart;
end
 

Nsig                 = length(bin_signal);
bin_signal_min       = max((bin_signal-data.no_nz_bins),1);
bin_signal_max       = min((bin_signal+data.no_nz_bins),floor(data.Np/2));
bins_signal          = cell(Nsig,1);
for i=1:Nsig
    bins_signal{i} = bin_signal_min(i):bin_signal_max(i);
end
data.bins_signal     = unique(cell2mat(bins_signal));
data.bins_fband      = max((bin_f0-bin_fband),1):min((bin_f0+bin_fband),floor(data.Np/2));
data.bins_fband      = setdiff(data.bins_fband,data.bins_signal);
data.segment_length  = floor(data.Np/data.segment_no);
data.overlap_no      = floor(data.Np/data.segment_no*data.overlapRatio);
data.overlap_string  = floor(data.overlapRatio * 100);
