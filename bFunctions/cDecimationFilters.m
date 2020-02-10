function data = cDecimationFilters(varargin)

%
% data = cDecimationFilters(varargin)
% This is a wrapper function for 'DecimationFilters' +
% 'CascadedCICCompensationFilter' function. The main purpose is to choose
% between a normal cascaded decimation filter (CIC = false) and cascaded
% decimation filter with the first stage is a CIC filter (CIC = true).
% The function inputs are given as:
%
%   Fs:                     Sampling frequency 
%   OSR:                    Over Sampling Ration
%   K:                      Number of decimation stages
%   M:                      Decimation vector for decimation factor at each stage
%   delta_F:                Decimator transtion bandwidth
%   rp:                     Pass-band ripples
%   rc:                     Cutoff-band ripples
%   Filter_Type:            Regular FIR or Multi-band FIR -> 'rg' or 'mb'
%   plot_filter_response:   Flag to plot filter frequency response or not -> 1 or 0
%   mb_type:                Define the type of the multi-band filter,
%                           whether narrowband 'NB' or wideband 'WB'
%   mb_taps_lengths:        Vector holds the filter lengths for multi-band filters 
%   CIC:                    Flag to set the first decimation stage to CIC
%                           filter (CIC = true) or not (CIC = false)
%   oSDM:                   Sigma-Delta modulator order, it is important
%                           for CIC filter order calculation
%   sdm_data:               Input stimuli data from prior SDM
%
%   filter_coefficients:    Decimator filter coefficients for eac decimation stage, it is 
%                           in a form of matrix, each row represents a decimation stages
%   filter_lengths:         Vector hold the length of each decimation filter stage
%   DataOut:                Output filtered data from the CIC filter
%   Hcascaded:              CIC+Compensation transfer function

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

if data.CIC == false,
    data = DecimationFilters('Fs', data.Fs, 'OSR', data.OSR, 'K', data.K, 'M', data.M, 'delta_F', data.delta_F, 'rp', data.rp, 'rc', data.rc, 'Filter_Type', data.Filter_Type, 'mb_type', data.mb_type, 'Pass_Stop', data.Pass_Stop, 'plot_filter_response', data.plot_filter_response);
else
    [data.DataOut data.Hcascaded] = CascadedCICCompensationFilter(1, data.M(1), data.oSDM+1, data.Fs, data.sdm_data);
    K = data.K-1;
    M = data.M(2:length(data.M));
    Fs = data.Fs/data.M(1);
    OSR = data.OSR/data.M(1);
    dataDecimationFilters = DecimationFilters('Fs', Fs, 'OSR', OSR, 'K', K, 'M', M, 'delta_F', data.delta_F, 'rp', data.rp, 'rc', data.rc, 'Filter_Type', data.Filter_Type, 'mb_type', data.mb_type, 'Pass_Stop', data.Pass_Stop, 'plot_filter_response', data.plot_filter_response);
    data.filter_coefficients    = dataDecimationFilters.filter_coefficients;
    data.filter_lengths         = dataDecimationFilters.filter_lengths;
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%  
% This function is implemented by: 
%   Micheal Maurer 
% and modified by: 
%   Ahmed Shahein

function data = ChkArgs(args,ni)

if ni == 0
    error('Syntax error::at least one argument required');
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
    % Check syntax {'str_1',var_1,'str_2',var_2,...}
    for i = 1:2:ni-1
        if ~isa(args{i},'char')
            error('Syntax error::string assumed')
        end
    end
       
    data = ChkString(data,args,'Fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'OSR','numeric','OSR is not numeric');
    data = ChkString(data,args,'K','numeric','number of decimation stages is not a number');
    data = ChkString(data,args,'M','numeric','decimation factor is not a number');  
    data = ChkString(data,args,'delta_F','numeric','transition band is not numeric');
    data = ChkString(data,args,'rp','numeric','passband ripples is not numeric');
    data = ChkString(data,args,'rc','numeric','stopband ripples is not numeric');   
    data = ChkString(data,args,'mb_taps_lengths','numeric','');
    
    data = ChkString(data,args,'sdm_data','numeric','sigma delta modulator data is not numeric');
    data = ChkString(data,args,'oSDM','numeric','sigma delta modulator order is not numeric');
    
    data = ChkString(data,args,'Filter_Type','char','filter architecture type is not a string');
    data = ChkString(data,args,'mb_type','char','multi-band filter type is not a string');
    
    data = ChkString(data,args,'Pass_Stop','logical','boolean expression assumed');
    data = ChkString(data,args,'plot_filter_response','logical','boolean expression assumed');
    
    data = ChkString(data,args,'CIC','logical','boolean expression assumed');
    
end

data = ChkStruct(data);

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
% This function is implemented by: 
%   Micheal Maurer 
% and modified by: 
%   Ahmed Shahein

function data = ChkString(data,args,string,type,err_msg)

pos = find(strcmp(args,string),1);
if isempty(pos)
elseif ~isa(args{pos+1},type) || (strcmp('numeric',type) && (any(isnan(args{pos+1})) || any(isinf(args{pos+1}))))
    error(['Syntax error::',err_msg])
else
    data.(string)=args{pos+1};
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
% This function is implemented by: 
%   Micheal Maurer 
% and modified by: 
%   Ahmed Shahein

function data = ChkStruct(data)

data.info.txt        = '\n';
if ~isfield(data,'Fs')
    data.Fs = 1;
end

if ~isfield(data,'OSR')
    data.OSR = 1;
end

if ~isfield(data,'delta_F')
    data.delta_F = 0.15;
end

if ~isfield(data,'K')
    error('Decimation stages is not define');
end

if ~isfield(data,'M')
     error('Decimation factor is not define');
end

if ~isfield(data,'oSDM')
     error('Sigma-Delta Modulator order is not define');
end

if ~isfield(data,'sdm_data')
     error('SDM bit-stream is not define');
end

if ~isfield(data,'rp')
    data.rp = 0.01;
end

if ~isfield(data,'rc')
    data.rc = 0.001;
end

if ~isfield(data,'Filter_Type')
    data.Filter_Type ='rg';
end

if ~isfield(data,'mb_type')
    data.mb_type ='WB';
end

if ~isfield(data,'Pass_Stop')
    data.Pass_Stop = true;
end

if ~isfield(data,'CIC')
    data.CIC = false;
end

if ~isfield(data,'plot_filter_response')
    data.print_IBN = false;
end

if ~isfield(data,'print_Sig')
    data.print_Sig = false;
end

if ~isfield(data,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.verbose;
    data = rmfield(data,'verbose');
end    

