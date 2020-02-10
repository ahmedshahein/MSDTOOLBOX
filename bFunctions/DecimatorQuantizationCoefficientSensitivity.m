function data = DecimatorQuantizationCoefficientSensitivity(varargin)

%
% [quantized_filter_coefficients Sn quantization_coefficients] = DecimatorQuantizationCoefficientSensitivity(filter_coefficients, filter_lengths, q, Fs, plot_freq_response, plot_Sn)
%
% This function calculate the quantized filter coefficients due to the 
% different bit widths specificed in the vector 'q'. It exports matrix
% holding the quantized filter coefficients and matrix holds the senstivity
% of the filter coefficients. It is for multistage decimator filter.
%
%   filter_coefficients:   Decimator filter coefficients for eac decimation stage, it is 
%                          in a form of matrix, each row represents a decimation stages
%   filter_lengths:        Vector hold the length of each decimation filter stagents
%   q:                     vector holds different bit widths for quantizing the filter
%                          coefficients, it can holds from 1 till 6 values
%   Fs:                    Sampling frequency
%
% P.S. This is wrapper function, i.e. to mask the inputs without order

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

% Intialization for some variablles in the function
quantized_coeff_temporary=[];
data.quantized_filter_coefficients=[];
data.quantization_coefficients=[];
[rows q_length] = size(data.Q);
K = length(data.filter_lengths);
m=1;

for a = 1 : K,
   filter_coeff_tmp = data.filter_coefficients(a, 1:data.filter_lengths(a));
   [quantized_coeff_temporary((a-1)*q_length+1:a*q_length,1:data.filter_lengths(a)) Sn_temporary((a-1)*q_length+1:a*q_length,1:data.filter_lengths(a))] = FilterCoefficientSensitivity(filter_coeff_tmp, data.Q(a,:), data.Fs, 0, data.plot_freq_response, data.plot_Sn);
end

q_vector=zeros(1,q_length*K);

for b = 1 : K,
    q_vector(q_length*(b-1)+1 : b*q_length) = data.Q(b,1:q_length);
end

if K == 2,
    for i = 1 : q_length,
        for j = 1 : q_length,
                index((m-1)*K+1:m*K) = [i j+q_length];
                m=m+1;
        end
    end
elseif K == 3,
    for i = 1 : q_length,
        for j = 1 : q_length,
            for k = 1 : q_length,
                index((m-1)*K+1:m*K) = [i j+q_length k+(2*q_length)];
                m=m+1;
            end
        end
    end
elseif K == 4,
    for i = 1 : q_length,
        for j = 1 : q_length,
            for k = 1 : q_length,
                for l = 1 : q_length,
                    index((m-1)*K+1:m*K) = [i j+q_length k+(2*q_length) l+(3*q_length)];
                    m=m+1;
                end
            end
        end
    end
end

for w = 1 : length(index),
    data.quantized_filter_coefficients(w,:)=quantized_coeff_temporary(index(w),:);
    data.Sn(w,:) = Sn_temporary(index(w),:); 
    data.quantization_coefficients(w) = q_vector(index(w));
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
       

    data = ChkString(data,args,'filter_coefficients','numeric','filter coefficeints is not numeric');
    data = ChkString(data,args,'filter_lengths','numeric','filter lengths is not numeric');
    data = ChkString(data,args,'Fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'Q','numeric','quantization bit-width is not a number');
    
    data = ChkString(data,args,'plot_freq_response','logical','boolean expression assumed');
    data = ChkString(data,args,'plot_Sn','logical','boolean expression assumed');
    
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
elseif ~isa(args{pos+1},type),
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
if ~isfield(data,'filter_coefficients')
    error('Filter coefficients is not define');
end

if ~isfield(data,'filter_lengths')
    error('Filter lengths is not define');
end

if ~isfield(data,'Fs')
    data.Fs = 1;
end

if ~isfield(data,'Q')
    error('Quantization bit-width is not define');
end

if ~isfield(data,'plot_freq_response')
    data.plot_freq_response = false;
end

if ~isfield(data,'plot_Sn')
    data.plot_Sn = false;
end

if ~isfield(data,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.verbose;
    data = rmfield(data,'verbose');
end    