function data = TuneK(varargin)

%
% [min_RT min_stage] = TuneK(Fs, OSR, delta_F, rp, rc, Stages, plot_RT, POT) 
%
% Thsi function estimate the computational effort due to number of
% decimation stages. It exports the min computational effort, the
% crossponding decimation value for each stage as a vector, and a matrix
% holding the computational effort and the crossponding decimation factor
% at each stage for all number of stage, i.e. the first column will hold
% the computational effort value while the the rest will holds the
% decimation value at each stage, since second column holds the value of
% the first decimation stage, the third column holds the value of the
% second decimation stage and so on.
% 
%   Fs:                     Sampling frequency
%   OSR:                    Over Sampling Ratio
%   delta_F:                Transition bandwidth = (fc-fp)/fc
%   K:                      Number of decimation stages
%   M:                      Vector holds the decimation factor at each stage
%   rp:                     Pass-band ripples
%   rc:                     Cutoff-band ripples
%   Stages:                 Vector holds the number of stages to be tested
%   RT_Stages:              Matrix holds the computation effort and its
%                           crossponding decimation vector
%
% P.S. This is wrapper function, i.e. to mask the inputs without order

% Check & assigning arguments
% Number of input arguments
nin     = nargin;
% Checking arguments
data    = ChkArgs(varargin,nin);

% Internal variables intialization
l=0;
RT=[];
depth=[];
width=[];

v = Factorize(data.OSR);

if length(v) <= 2,
    fprintf('This decimation factor cannot hold more than 2 stages\n');
    fprintf('Please change Stages to 2 only\n');
    data.min_RT = nan;
    data.min_stage = nan;
end    
    decimation_matrices = struct('matrix', [], 'sizes', []);

    for i = 1 : length(data.Stages),
        if data.POT == 1,
            decimation_matrices(i).matrix = DecimationMatrix(data.OSR,data.Stages(i));
        else
            decimation_matrices(i).matrix = DecimationMatrixNonPOT(data.OSR,data.Stages(i));
        end
        [depth(i) width(i)] = size(decimation_matrices(i).matrix);
        decimation_matrices(i).sizes = depth(i);
    end

    computationaleffort = struct('RT', [], 'Stages', []);
    l = 0;
    for j = 1 : length(data.Stages),
        matrix = decimation_matrices(j).matrix;
         for k = 1 : depth(j),
             l=l+1;
             M = matrix(k,:);
             RT(l) = ComputationalEffort(data.Fs, data.OSR, data.delta_F, data.Stages(j), M, data.rp, data.rc);
             computationaleffort(l).RT = RT(l);
             computationaleffort(l).Stages = M;
         end
    end

    data.min_RT = min(RT);
    for o = 1 : length(RT),
        if RT(o) == data.min_RT,
            data.min_stage = computationaleffort(o).Stages;
        end
    end

    x=[1:sum(depth)];

    color = ['b', 'r', 'k', 'g', 'c', 'y', 'm'];
    StagesVector = [];
    
    if data.plot_RT == 1,        
        if sum(x) > 1,
          FIG = figure('Name', 'Computational Effort', 'NumberTitle' , 'off');
            for i = 1 : length(computationaleffort),    
                if length(computationaleffort(i).Stages) == 2,
                    clr = color(1);
                    StagesVector(i) = 2;
                elseif length(computationaleffort(i).Stages) == 3,
                    clr = color(2);
                    StagesVector(i) = 3;
                elseif length(computationaleffort(i).Stages) == 4,
                    clr = color(3);
                    StagesVector(i) = 4;
                end
                    plot2stages = plot(0,0,color(1));
                    plot3stages = plot(0,0,color(2));
                    plot4stages = plot(0,0,color(3));
                    bar(x(i), computationaleffort(i).RT,0.3,clr);
                    hold on
                    set(gca, 'XTickLabel', 'Stages')
            end
                ylabel('RT - MADS');
                set(gca, 'XTick', [1:sum(x)], 'XTicklabel', []);
                legend([plot2stages,plot3stages,plot4stages], '2 Stages', '3 Stages', '4 Stages',2)
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
    % Check syntax {'str_1',var_1,'str_2',var_2,...}
    for i = 1:2:ni-1
        if ~isa(args{i},'char')
            error('Syntax error::string assumed')
        end
    end
    data = ChkString(data,args,'Fs','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'OSR','numeric','OSR is not numeric');
    data = ChkString(data,args,'delta_F','numeric','transition band is not numeric');
    data = ChkString(data,args,'rp','numeric','passband ripples is not numeric');
    data = ChkString(data,args,'rc','numeric','stopband ripples is not numeric');
    data = ChkString(data,args,'Stages','numeric','norm fs is not vector of numeric');
    data = ChkString(data,args,'plot_RT','logical','boolean expression assumed');
    data = ChkString(data,args,'POT','logical','boolean expression assumed');
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

if ~isfield(data,'rp')
    data.rp = 0.01;
end

if ~isfield(data,'rc')
    data.rc = 0.001;
end

if ~isfield(data,'Stages'),
    error('Decimation stages is not define');
end

if ~isfield(data,'plot_RT')
    data.plot_RT = true;
end

if ~isfield(data,'POT')
    data.POT = true;
end

if ~isfield(data,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.verbose;
    data = rmfield(data,'verbose');
end