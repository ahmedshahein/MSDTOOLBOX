% This is a complete design example, it is divided to a set of cells, each
% cell demonstrate a certain design step. It is recommended to run each
% cell indvidually and analyze the output on the command window and then
% import any required data from it and then proceed to the following cell.
% In this example there the optimization criteria using coefficient deviation 
% has been introduced for the coefficient.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-1: Load the input data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Workspace and command window intialization
close all;
clear;
clc;

% Load a 3nd Order LPDSM 3-bit pattern for testing purpose
sdm_data    = load('iPatterns/OSR48_9600kHz_3rdOrder_3bit.txt');

% Design input parameters
Fs          = 9.6e6;        % Sampling frequency
OSR         = 48;           % Over Sampling Frequency
Fsignal     = 20e3;          % Input signal frequency

delta_F     = 0.15;         % Transition bandwidth = (cutoff band-pass band)/cutoff band
rp_tune     = [0.00001 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1];        % Passband ripples
rc_tune     = [0.00001 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1];        % Cutoff suppression

Stages      = [2 3 4];      % Number of decimation stages to be tunned, instead of entering K and M
                            % Here it will test 2, 3 and 4 decimation
                            % stages

IBN_penalty     = -4;       % Acceptable penalty in IBN after deciamtion, 
                            % penalty = IBN'before decimation' - IBN'after deciamtion'
                            % e.g.
                            % -70 - (-66) = -4dB
Sig_penalty     = 0.5;      % Acceptable penalty in the signal peak,
                            % Sig_penalty = Signal peak'before deciamtion'-Signal peak'after deciamtion'
                            % e.g.
                            % -3dB - (-3.2) = 0.2dB
                            
Filter_Type = 'mb';         % Flag to constrain the filter structure, wether to be a 
                            % reqular 'rg' FIR filter structure or multi-band 'mb' 
                            % FIR filter structure
mb_type     = 'WB';         % Type of multiband filter whether narrowband 'NB' or wideband 'WB'

Pass_Stop   = false;        % Flag to determine that base-band frequency is constrinted as
                            % passband or stopband frequency corner.
                            % false -> Pass
                            % true  -> Stop
                            
q           = [[20 18 16 1 12 10 8]; [20 18 16 1 12 10 8]; [20 18 16 1 12 10 8]];
                            % This matrix holds the quantization bit widths
                            % for each decimation stages, since , in this
                            % example, the first row for the first stage
                            % and the second row for the second decimation
                            % stage, the number of columns suppose to be
                            % equivelent to the number of decimation stages

POT                  = false;   % Flag used for the generation of the decimation factor, to constrain it to a 
                            % Power-Of-Two 'POT' factors or not                            
plot_filter_response = true;   % Flag to plot the frequency response of each decimation stage
                            % 1 -> Plot, 0 -> don't plot
plot_psd             = true;   % Flag to plot the PSD before and after deciamtion, where
                            % 1 -> Plot, 0 -> don't plot
export_IBN           = true;   % Flag to export the IBN before and after decimation, where
                            % 1 -> Export, 0 -< don't export
                            % It is exported in a vector format IBN = [b a]
                            % b element is the IBN before decimation
                            % while a after decimation
print_IBN            = true;   % Flag to print the IBN before and after decimation on the command window, where
                            % 1 -> print, 0 -> don't print1
plot_RT              = true;   % Flag to plot the computation effort bar diagram, where
                            % 1 -< Plot, 0 -> don't plot
plot_freq_response   = false;   % Flag to plot the filter frequency response with quantization effect, where
                            % 1 -> Plot, 0 -> don't plot
plot_Sn              = false;   % Flag to plot the filter coefficient sensitvity, where
                            % 0 -> don't plot, 1 -> plot 
print_Sig            = 0;   % Flag to print the signal peak before and after decimation proces, where
                            % 1 -> Print, 0 -> don't print
PlotEachStage        = 1;   % 1 -> Bar plot for power consumption in each decimation stage
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-2: Tuning K and M for minimum computational effort RT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rp = 0.01;
rc = 0.001;

% Tunning the number of decimation stages 'K' and the deciamtion factor at
% each stage 'M', to satsfy minimum computational effort 'RT' after we had
% tune thr ripples factor.
dataTuneK = TuneK('Fs', Fs, 'OSR', OSR, 'delta_F', delta_F, 'rp', rp, 'rc', rc, 'Stages', Stages, 'POT', POT, 'plot_RT', plot_RT);

% Pass K and M to the design functions
K   = length(dataTuneK.min_stage);
M   = dataTuneK.min_stage;

Print_KM(K, M);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-3: Tune the passband ripples and stopband attenuation factors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Tune ripples 
% In this step the objective is to find the optimal passband ripples and
% cutoff suppresion values to satisfy the accepted penalty value by the
% designer. Here you need to setup an intial value for the number of
% decimation stages 'K' and the decimation value for each stage 'M', the
% easiest values is to set K = 2 and M = [OSR/2 2].
dataTuneR = TuneR('sdm_data', sdm_data, 'rpb_tune', rp_tune, 'rsb_tune',rc_tune, 'Fs',  Fs, 'OSR', OSR, 'delta_F', delta_F, 'Fsignal', Fsignal, 'K', K, 'M', M, 'Filter_Type', Filter_Type, 'mb_type', mb_type, 'Pass_Stop', Pass_Stop, 'IBN_penalty', IBN_penalty);

% 
fprintf('rp  --  rc  --   differencr in IBN\n');
for i = 1 : size(dataTuneR.accepted_ripples, 1),    
    fprintf('%2.4f  %3.4f   %3.4f      \n', dataTuneR.accepted_ripples(i,:));
end
% Analyzing the output from the previous step, we get the following table;
% The 3rd column represent the difference between IBN before and after
% decimation and which in the range of the accepted penalty value
% constrained by the designer. I will choose;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-4: Design the decimator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rp = 0.00001; 
rc = 0.000001; 

% Design the decimation stages due to the input parameters and constrains
dataDecimationFilters = DecimationFilters('Fs', Fs, 'OSR', OSR, 'K', K, 'M', M, 'delta_F', delta_F, 'rp', rp, 'rc', rc, 'Filter_Type', Filter_Type, 'mb_type', mb_type, 'Pass_Stop', Pass_Stop, 'plot_filter_response', plot_filter_response);

% Filter and downsample the LPDSM data using the designed decimation stages
% in the previous step
[deci_data IBN] = FilterAndDownsample(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, sdm_data, Fs, OSR, Fsignal, K, M, plot_psd, export_IBN, print_IBN);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-5: Tune coefficeint quantization factor % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Quantized decimation filters, the output from this step will be the
% quantized filter coefficients for each decimation stage and the
% coefficient sensitvity
dataDecimatorQuantizationCoefficientSensitivity = DecimatorQuantizationCoefficientSensitivity('filter_coefficients', dataDecimationFilters.filter_coefficients, 'filter_lengths', dataDecimationFilters.filter_lengths, 'Q', q, 'Fs', Fs, 'plot_freq_response', plot_freq_response, 'plot_Sn', plot_Sn);


% This function test the quantization effect of the filter coefficients, by
% examining the its effect on IBn and Signal peak. It exports the
% quantization bitwidth for each stage and its effect on the IBN and Signal
% peak
TestQuantizedFiltersIBNSig(dataDecimatorQuantizationCoefficientSensitivity.quantized_filter_coefficients, dataDecimationFilters.filter_lengths, dataDecimatorQuantizationCoefficientSensitivity.quantization_coefficients, sdm_data, Fs, OSR, Fsignal, K, M, q, 0, export_IBN, 0, print_Sig, IBN_penalty, Sig_penalty);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-6: SPT/NPT rounding mechanism %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Analyzing the output table on the command window, we can obserev is those
% 2 values satisfy the penalties range and support optimal bit width for
% coefficient quantization.

Q = [16 14 12]; 
Nb = 2;
dataNormalizedCoefficients = NormalizedCoefficients('filter_coefficients', dataDecimationFilters.filter_coefficients, 'filter_lengths', dataDecimationFilters.filter_lengths, 'Q', Q, 'K', K);
EstimateOptiDeviation(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, Fs, K, Q, Nb, 20)
N_vector = [1 2];
deviation_values = [[0 2 3 6 8 10 11 12 13 14]; [4 10 11 12 13 14 16 18 20 28]; [4 10 11 12 13 14 16 18 20 28]];
coeff_matrix = TestOptiDeviation(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, Fs, OSR, Fsignal, K, M, Q, Nb, sdm_data, deviation_values);

mixed_filter_coefficients = [coeff_matrix(2,4:size(coeff_matrix,2)); coeff_matrix(13, 4:size(coeff_matrix,2)); coeff_matrix(24, 4:size(coeff_matrix,2))];
IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, mixed_filter_coefficients, dataDecimationFilters.filter_lengths, sdm_data)
%%