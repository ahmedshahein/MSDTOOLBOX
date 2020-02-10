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
gyro_data    = load('iPatterns/I_Gyroscope_103720_2_20');

% Forming complex multiplier for downmixing BPSDM output signal to baseband
% n=[0:length(gyro_data)-1];
% multi = exp(i*2*pi*(0.125)*n); %exp(pi*n/4 + 2*pi*0.00454);% cos(2*pi*0.45)*cos(pi*n/4)-sin(2*pi*0.45)*sin(pi*n/4); % %
% multi=reshape(multi,length(multi),1);

% Mixing the BPSDM data to baseband, the output of this stage is complex
% signal (I/Q signals)
%sdm_data=real(gyro_data.*multi);%real(gyro_data.*multi);
%sdm_data    = I_bitstream(:,2);
sdm_data    = gyro_data.I;
% Design input parameters
Fs          = 103720;        % Sampling frequency
OSR         = 1000;           % Over Sampling Frequency
Fsignal     = 0; %12954;          % Input signal frequency

delta_F     = 0.072;         % Transition bandwidth = (cutoff band-pass band)/cutoff band
rp_tune     = [0.0001 0.0005 0.001 0.005 0.01 0.1];        % Passband ripples
rc_tune     = [0.0001 0.0005 0.001 0.005 0.01 0.1];        % Cutoff suppression

Stages      = [2 3 4];      % Number of decimation stages to be tunned, instead of entering K and M
                            % Here it will test 2, 3 and 4 decimation
                            % stages

IBN_penalty     = -15;       % Acceptable penalty in IBN after deciamtion, 
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
mb_type     = 'NB';         % Type of multiband filter whether narrowband 'NB' or wideband 'WB'

Pass_Stop   = false;        % Flag to determine that base-band frequency is constrinted as
                            % passband or stopband frequency corner.
                            % false -> Pass
                            % true  -> Stop
                            
q           = [[18 16 15 14 13 12 10]; [18 16 14 12 11 10 9]; [18 16 14 12 11 10 9]];
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
print_Sig            = 1;   % Flag to print the signal peak before and after decimation proces, where
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
% each stage 'MK', to satsfy minimum computational effort 'RT' after we had
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
% decimation stages 'K' and the deciamtion value for each stage 'M', the
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell-4: Design the decimator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rpb = 0.1;  %0.1; for BW 100Hz
rsb = 0.01; %0.0001;

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

Q = [12 8 6];
Nb=2;
dataNormalizedCoefficients = NormalizedCoefficients('filter_coefficients', dataDecimationFilters.filter_coefficients, 'filter_lengths', dataDecimationFilters.filter_lengths, 'Q', Q, 'K', K);
IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, dataNormalizedCoefficients.normalized_filter_coeff, dataDecimationFilters.filter_lengths, sdm_data)
% IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, normalized_filter_coeff, filter_lengths, sdm_data)
% EstimateOptiDeviation(filter_coefficients, filter_lengths, Fs, K, Q, Nb, 4)
% N_vector = [1 2];
% deviation_values = [[6 7 8 9 10 12 13 14 15]; [1 2 3 4 5 6 7 8 9]; [ 1 1.5 2 4 5 10 15 20 25]];
% coeff_matrix = TestOptiDeviation(filter_coefficients, filter_lengths, Fs, OSR, Fsignal, K, M, Q, Nb, sdm_data, deviation_values);
% 
% mixed_filter_coefficients = [coeff_matrix(3,4:length(coeff_matrix)); coeff_matrix(13, 4:length(coeff_matrix)); coeff_matrix(28, 4:length(coeff_matrix))];
% IBN_norm = TestDecimatorIBN(Fs, OSR, Fsignal, K, M, mixed_filter_coefficients, filter_lengths, sdm_data)

%%