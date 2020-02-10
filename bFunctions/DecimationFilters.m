function data = DecimationFilters(varargin)

%
% [filter_coefficients, filter_lengths] = DecimationFilters(Fs, OSR, K, M, delta_F, rp, rc, Filter_Type, plot_filter_response, mb_taps_lengths)
%
% This function design and export LOWPASS decimator filters coefficients, for both
% regular and multi-band structures.
%
%   Fs:                    Sampling frequency 
%   OSR:                   Over Sampling Ration
%   K:                     Number of decimation stages
%   M:                     Decimation vector for decimation factor at each stage
%   delta_F:               Decimator transtion bandwidth
%   rp:                    Pass-band ripples
%   rc:                    Cutoff-band ripples
%   Filter_Type:           Regular FIR or Multi-band FIR -> 'rg' or 'mb'
%   plot_filter_response:  Flag to plot filter frequency response or not -> 1 or 0
%   mb_type:               Define the type of the multi-band filter,
%                          whether narrowband 'NB' or wideband 'WB'
%   mb_taps_lengths:       Vector holds the filter lengths for multi-band filters 
%
%   filter_coefficients:   Decimator filter coefficients for eac decimation stage, it is 
%                          in a form of matrix, each row represents a decimation stages
%   filter_lengths:        Vector hold the length of each decimation filter stage
%
% P.S. This is wrapper function, i.e. to mask the inputs without order

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

% Intialization for some variablles in the function
    wi = [];
    filter_lengths_rg = [];
    filter_coefficients_rg = [];
    
% Turn off the warning of Divide by zero

%     1/0;
%     [msg, msgid] = lastwarn;
%     s = warning('off', msgid);

% Some default values in the case of the absence of some parameters or for
% intial testing

    if nargin == 4,
        data.delta_F         = 0.1;
        data.rp              = 0.001;
        data.rc              = 0.00001;
        data.Filter_Type     = 'rg';
    end

% Design decimator filters in a regular FIR structure
    
    Fb = data.Fs/(2*data.OSR);
    
    if data.Pass_Stop == false,
        % If this flag is set to zero/false
        % The passband frequency corner is set to 
        % baseband frequency, i.e.
        % Fp = Fb = Fs/(2OSR)
        for i = 1 : data.K,
            if i < data.K,
                Fp(i) = Fb; 
                Fc(i) = data.Fs/prod(data.M(1:i)) - Fb/(1-data.delta_F); 
                Rp(i) = data.rp/data.K;
                Rc(i) = data.rc;
            else
                Fp(i) = Fb; 
                Fc(i) = Fb/(1-data.delta_F);
                Rp(i) = data.rp/data.K;
                Rc(i) = data.rc;  
            end
        end
        
    else
        % If this flag is set to zero/false
        % The stopband frequency corner is set to 
        % baseband frequency, i.e.
        % Fs = Fb = Fs/(2OSR)        
        for i = 1 : data.K,
            Fp(i) = Fb*(1-data.delta_F); 
            Fc(i) = data.Fs/prod(data.M(1:i)) - Fb; 
            Rp(i) = data.rp/data.K;
            Rc(i) = data.rc;
        end    
        
    end

    filter = struct('n', [], 'fo', [], 'ao', [], 'w', []);

    filter_stages_rg=struct('b', []);

    for j = 1 : data.K,
        if j == 1,
            [filter(j).n, filter(j).fo, filter(j).ao, filter(j).w] = remezord([Fp(j) Fc(j)], [1 0], [Rp(j) Rc(j)], data.Fs);
            if nargin == 11 & mb_taps_lengths ~= 0,            
                filter_stages_rg(j).b  = remez(mb_taps_lengths(j), filter(j).fo, filter(j).ao, filter(j).w);
            else
                filter_stages_rg(j).b  = remez(filter(j).n, filter(j).fo, filter(j).ao, filter(j).w);
            end
        else
            [filter(j).n, filter(j).fo, filter(j).ao, filter(j).w] = remezord([Fp(j) Fc(j)], [1 0], [Rp(j) Rc(j)], data.Fs/prod(data.M(1:j-1)));
            if nargin == 11 & mb_taps_lengths ~= 0,            
                filter_stages_rg(j).b  = remez(mb_taps_lengths(j), filter(j).fo, filter(j).ao, filter(j).w);
            else
                filter_stages_rg(j).b  = remez(filter(j).n, filter(j).fo, filter(j).ao, filter(j).w);
            end
        end
            filter_lengths_rg(j) = length(filter_stages_rg(j).b);
            filter_coefficients_rg(j, 1:filter_lengths_rg(j))=filter_stages_rg(j).b;
    end

% Design decimator filters in a multi-band FIR structure
if nargin >= 10 & data.mb_type ~= 0,
    [f, a, w, lengths] = MultiBandFilters(data.K, data.M, data.OSR, data.mb_type);
else
    [f, a, w, lengths] = MultiBandFilters(data.K, data.M, data.OSR, 'WB');
end

% Designing Half-band stage
    N = 25;    % Filter Length
    P = (N+1)/2; % Filter Peak
    Filter_Coeff_Length = [1:N];
    Filter_Coeff_Length = Filter_Coeff_Length - P;
    HB_Sinc_Filter = sin(pi*Filter_Coeff_Length/2)./(pi*Filter_Coeff_Length);
    HB_Sinc_Filter(P)=.5;
    thewin=window(@hamming, N);
    HB_Window = HB_Sinc_Filter .* thewin';
    HB_Filter = HB_Window/sum(HB_Window);
    HB_stage = ZeroOutCoeff(HB_Filter);


    if nargin == 11 & mb_taps_lengths ~= 0,
        n = mb_taps_lengths;
    else
        n = filter_lengths_rg;
    end
    
    filter_stages_mb=struct('b', []);

    for o = 1 : data.K,
       if o == data.K,
            filter_lengths_mb(o) = length(HB_stage);
            filter_coefficients_mb(o, 1:filter_lengths_mb(o))=HB_stage;
       else
           if mod(lengths(o), 2) == 0,
                filter_stages_mb(o).b = remez(n(o), f(o,1:lengths(o)), a(o,1:lengths(o)), w(o,1:lengths(o)/2));
                filter_lengths_mb(o) = length(filter_stages_mb(o).b);
                filter_coefficients_mb(o, 1:filter_lengths_mb(o))=filter_stages_mb(o).b;
           else
                filter_stages_mb(o).b = remez(n(o), f(o,1:lengths(o)-1), a(o,1:lengths(o)-1), w(o,1:(lengths(o)-1)/2));
                filter_lengths_mb(o) = length(filter_stages_mb(o).b);
                filter_coefficients_mb(o, 1:filter_lengths_mb(o))=filter_stages_mb(o).b;
           end
       end
    end
    
    if data.K == 1 & data.Filter_Type == 'mb',
        fprintf('WARNING: You cannot use Multi-band structure for single decimation stage \n');
    end

% Output filter length and coefficients due to the the specified structure

     if data.Filter_Type == 'rg',
         data.filter_coefficients    = filter_coefficients_rg;
         data.filter_lengths         = filter_lengths_rg;     
     elseif data.Filter_Type == 'mb' & data.K ~= 1,
         data.filter_coefficients   = filter_coefficients_mb;
         data.filter_lengths         = filter_lengths_mb; 
     else    
        data.filter_lengths = nan;
        data.filter_coefficients = nan;
        fprintf('Not a Recognized Filter Type \n');   
     end

% Plotting decimator filters frequency response option

    if nargin ~= 4 & data.Filter_Type == 'mb' & data.K == 1,
        return
    else
        if data.plot_filter_response == 1,
            filter_response=struct('h', [], 'f', []);
            FIG = figure('Name', 'Filter Frequency Response', 'NumberTitle' , 'off');
            for p = 1 : data.K,
                if p == 1,
                    [filter_response(p).h, filter_response(p).f] = freqz(data.filter_coefficients(p,1:data.filter_lengths(p)),[1],2^14,data.Fs);
                    subplot(data.K,1,p)
                    plot(filter_response(p).f, 20*log10(abs(filter_response(p).h)));
                    xlabel('Frequency - Hz');
                    ylabel('Amplitude - dB');
                    grid on;
                    title(['Filter Response for Stage - ', num2str(p)]);
                else
                    [filter_response(p).h, filter_response(p).f] = freqz(data.filter_coefficients(p,1:data.filter_lengths(p)),[1],2^14,data.Fs/(prod(data.M(1:p-1))));
                    subplot(data.K,1,p)
                    plot(filter_response(p).f, 20*log10(abs(filter_response(p).h)));
                    xlabel('Frequency - Hz');
                    ylabel('Amplitude - dB');
                    grid on;
                    title(['Filter Response for Stage - ', num2str(p)]);
                end
            end
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
       
    data = ChkString(data,args,'Fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'OSR','numeric','OSR is not numeric');
    data = ChkString(data,args,'K','numeric','number of decimation stages is not a number');
    data = ChkString(data,args,'M','numeric','decimation factor is not a number');  
    data = ChkString(data,args,'delta_F','numeric','transition band is not numeric');
    data = ChkString(data,args,'rp','numeric','passband ripples is not numeric');
    data = ChkString(data,args,'rc','numeric','stopband ripples is not numeric');   
    data = ChkString(data,args,'mb_taps_lengths','numeric','');
    
    data = ChkString(data,args,'Filter_Type','char','filter architecture type is not a string');
    data = ChkString(data,args,'mb_type','char','multi-band filter type is not a string');
    
    data = ChkString(data,args,'Pass_Stop','logical','boolean expression assumed');
    data = ChkString(data,args,'plot_filter_response','logical','boolean expression assumed');
    
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