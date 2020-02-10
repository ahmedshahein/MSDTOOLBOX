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

% Load a 2nd Order LPDSM 3-bit pattern for testing purpose
sdm_data    = load('iPatterns/OSR24_3bit_960kHz.txt');

% Design input parameters
Fs          = 960e3;        % Sampling frequency
OSR         = 24;           % Over Sampling Frequency
Fsignal     = 5e3;          % Input signal frequency

delta_F     = 0.15;         % Transition bandwidth = (cutoff band-pass band)/cutoff band
rpb_tune    = [0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01];        % Passband ripples
rsb_tune    = [0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01];        % Cutoff suppression

Stages      = [2 3 4];      % Number of decimation stages to be tunned, instead of entering K and M
                            % Here it will test 2, 3 and 4 decimation
                            % stages

IBN_penalty     = -4;       % Acceptable penalty in IBN after deciamtion, 
                            % penalty = IBN'before decimation' - IBN'after deciamtion'
                            % e.g.
                            % -70 - (-66) = -4dB
Sig_penalty     = 0.6;      % Acceptable penalty in the signal peak,
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
                            
q           = [[20 18 16 14 12 10 8]; [20 18 16 14 12 10 8]; [20 18 16 14 12 10 8]];
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
plot_freq_response   = false;   % Flag to plot the filter frequency response with quantization effect, where
                            % 1 -> Plot, 0 -> don't plot
plot_RT              = true;   % Flag to plot the computation effort bar diagram, where
                            % 1 -< Plot, 0 -> don't plot
plot_Sn              = false;   % Flag to plot the filter coefficient sensitvity, where
                            % 0 -> don't plot, 1 -> plot 
print_Sig            = true;   % Flag to print the signal peak before and after decimation proces, where
                            % 1 -> Print, 0 -> don't print
PlotEachStage        = 1;   % 1 -> Bar plot for power consumption in each decimation stage
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-2: Tuning K and M for minimum computational effort RT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rpb = 0.01;
rsb = 0.001;

% Tunning the number of decimation stages 'K' and the deciamtion factor at
% each stage 'MK', to satsfy minimum computational effort 'RT' after we had
% tune thr ripples factor.
dataTuneK = TuneK('Fs', Fs, 'OSR', OSR, 'delta_F', delta_F, 'rp', rpb, 'rc', rsb, 'Stages', Stages, 'POT', POT, 'plot_RT', plot_RT);

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
% decimation stages 'K' and the deciamtion value for each stage 'M', the
% easiest values is to set K = 2 and M = [OSR/2 2].
dataTuneR = TuneR('sdm_data', sdm_data, 'rpb_tune', rpb_tune, 'rsb_tune',rsb_tune, 'Fs',  Fs, 'OSR', OSR, 'delta_F', delta_F, 'Fsignal', Fsignal, 'K', K, 'M', M, 'Filter_Type', Filter_Type, 'mb_type', mb_type, 'Pass_Stop', Pass_Stop, 'IBN_penalty', IBN_penalty);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-4: Design the decimator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rpb = 0.0001;
rsb = 0.00001;

% Design the decimation stages due to the input parameters and constrains [filter_coefficients, filter_lengths]
dataDecimationFilters = DecimationFilters('Fs', Fs, 'OSR', OSR, 'K', K, 'M', M, 'delta_F', delta_F, 'rp', rpb, 'rc', rsb, 'Filter_Type', Filter_Type, 'mb_type', mb_type, 'Pass_Stop', Pass_Stop, 'plot_filter_response', plot_filter_response);

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

Q = [18 14 8];
Nb=2;
dataNormalizedCoefficients = NormalizedCoefficients('filter_coefficients', dataDecimationFilters.filter_coefficients, 'filter_lengths', dataDecimationFilters.filter_lengths, 'Q', Q, 'K', K);
EstimateOptiDeviation(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, Fs, K, Q, Nb, 4)
N_vector = [1 2];
deviation_values = [[6 7 8 9 10 12 13 14 15]; [6 7 8 9 10 11 12 13 20]; [0.7 0.8 1 1.5 2 4 5 10 20]];
coeff_matrix = TestOptiDeviation(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, Fs, OSR, Fsignal, K, M, Q, Nb, sdm_data, deviation_values);

mixed_filter_coefficients = [coeff_matrix(6,4:length(coeff_matrix)); coeff_matrix(12, 4:length(coeff_matrix)); coeff_matrix(28, 4:length(coeff_matrix))];
IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, mixed_filter_coefficients, dataDecimationFilters.filter_lengths, sdm_data)

% Print out the total number of approximated/rounded coefficient for each
% stage
fprintf('Total number of approximated coefficients in the first stage = %d\n', coeff_matrix(6,2));
fprintf('Total number of approximated coefficients in the second stage = %d\n', coeff_matrix(12,2));
fprintf('Total number of approximated coefficients in the third stage = %d\n', coeff_matrix(28,2));

% Export the fixed point filter coefficients into text files
% dlmwrite('./Exported_Coefficients/First_Stage_Coeffs_960kHz.txt', mixed_filter_coefficients(1,:), ';');
% dlmwrite('./Exported_Coefficients/Second_Stage_Coeffs_960kHz.txt', mixed_filter_coefficients(2,:), ';');
% dlmwrite('./Exported_Coefficients/Third_Stage_Coeffs_960kHz.txt', mixed_filter_coefficients(3,:), ';');

%%