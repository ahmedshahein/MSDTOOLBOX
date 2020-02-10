function  data = NormalizedCoefficients(varargin)

%
% normalized_filter_coeff = NormalizedCoefficients(filter_coefficients, filter_lengths, n, K)
%
% This function is just for looping the 'normalized_coeff' function for
% each deciamtion stage from the multistage deciamtion filters, i.e. K
% times.
%
% P.S. This is wrapper function, i.e. to mask the inputs without order

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

for i = 1 : data.K,
    data.normalized_filter_coeff(i,1:data.filter_lengths(i)) = NormalizedCoeff(data.filter_coefficients(i,1:data.filter_lengths(i)), data.Q(i));
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
       
    data = ChkString(data,args,'filter_coefficients','numeric','filter coefficeints is not numeric');
    data = ChkString(data,args,'filter_lengths','numeric','filter lengths is not numeric');
    data = ChkString(data,args,'K','numeric','number of decimation stages is not a number');
    data = ChkString(data,args,'Q','numeric','quantization bit-width is not a number');  
    
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
elseif ~isa(args{pos+1},type) 
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

if ~isfield(data,'K')
    error('Number of decimation stages is not define');
end

if ~isfield(data,'Q')
    error('Quantization bit-width is not define');
end

if ~isfield(data,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.verbose;
    data = rmfield(data,'verbose');
end