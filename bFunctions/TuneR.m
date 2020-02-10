function data = TuneR(varargin)

%
% accepted_ripples = TuneR(rp_tune, rc_tune, Fs, OSR, delta_F, Fsignal, K, M, Filter_Type, sdm_data, IBN_penalty, print_IBN, print_Sig) 
%
% This fnction used for the ripples tunning purpose.
%
%   rp_tune                 Passband ripples factor vector
%   rc_tune                 Cutoff band ripples factor vector
%   Fs:                     Sampling frequency
%   OSR:                    Over Sampling Ratio
%   delta_F:                Transition bandwidth = (fc-fp)/fc
%   K:                     Number of decimation stages
%   M:                     Decimation vector for decimation factor at each stage
%   Filter_Type:            Regular FIR or Multi-band FIR -> 'rg' or 'mb'
%   sdm_data:               Array holds the Delta-Sigma Modulator data
%
%   diff:                   The variations in IBN before and after
%                           decimation using the filters with the selected
%                           ripples factor in the 3rd columns, while the
%                           first column holds the passband ripples value
%                           and the second column hold the cutoff band
%                           ripples value
%
% P.S. This is wrapper function, i.e. to mask the inputs without order

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

%filter_accepted_ripples=[];
k=1;

if nargin > 12,
    printIBN = data.print_IBN;
    printSig = data.print_Sig;
else
    printIBN = 0;
    printSig = 0;
end

for i = 1 : length(data.rpb_tune),
    rp = data.rpb_tune(i);
    for j = 1 : length(data.rsb_tune),
        rc = data.rsb_tune(j);
        dataDecimationFilters = DecimationFilters('Fs', data.Fs, 'OSR', data.OSR, 'K', data.K, 'M', data.M, 'delta_F', data.delta_F, 'rp', rp, 'rc', rc, 'Filter_Type', data.Filter_Type, 'mb_type', data.mb_type, 'Pass_Stop', data.Pass_Stop, 'plot_filter_response', false);
        [deci_data IBN(k,:) Sig(k,:)] = FilterAndDownsample(dataDecimationFilters.filter_coefficients, dataDecimationFilters.filter_lengths, data.sdm_data, data.Fs, data.OSR, data.Fsignal, data.K, data.M, 0, 1, printIBN, printSig);
        ripples(k,:)=[rp rc];
        k=k+1;
    end
end

for l = 1 : size(IBN,1),
    diff(l,3) = IBN(l, 1) - IBN(l,2); % before - after;
    diff(l,1:2) = ripples(l,:);
end 

n = 1;
data.accepted_ripples = [];
for m = 1 : size(IBN,1), 
    if diff(m,3) >= data.IBN_penalty & diff(m,3) < 0,
        data.accepted_ripples(n,:) = diff(m,:);
        n = n + 1;
    end
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
    % Check syntax {sdm_data, 'str_1',var_1,'str_2',var_2,...}
    for i = 1:2:ni-1
        if ~isa(args{i},'char')
            error('Syntax error::string assumed')
        end
    end
    
    data = ChkString(data,args,'sdm_data','numeric','sdm bit-stream is not numeric');
    
    data = ChkString(data,args,'rpb_tune','numeric','passband ripples is not numeric');
    data = ChkString(data,args,'rsb_tune','numeric','stopband ripples is not numeric');
    
    data = ChkString(data,args,'Fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'OSR','numeric','OSR is not numeric');
    data = ChkString(data,args,'delta_F','numeric','transition band is not numeric');
    data = ChkString(data,args,'Fsignal','numeric','signal frequency is not a number');
    data = ChkString(data,args,'K','numeric','number of decimation stages is not a number');
    data = ChkString(data,args,'M','numeric','decimation factor is not a number');
    
    data = ChkString(data,args,'Filter_Type','char','filter architecture type is not a string');
    data = ChkString(data,args,'mb_type','char','multi-band filter type is not a string');
    
    data = ChkString(data,args,'IBN_penalty','numeric','penalty in IBN is not a number');
    
    data = ChkString(data,args,'Pass_Stop','logical','boolean expression assumed');
    data = ChkString(data,args,'print_IBN','logical','boolean expression assumed');
    data = ChkString(data,args,'print_Sig','logical','boolean expression assumed');
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

if ~isfield(data,'rpb_tune')
    data.rpb_tune = [0.1 0.01 0.001 0.0001];
end

if ~isfield(data,'rsb_tune')
    data.rsb_tune = [0.1 0.01 0.001 0.0001];
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

if ~isfield(data,'print_IBN')
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
% EOF