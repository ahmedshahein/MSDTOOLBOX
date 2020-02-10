function [MultiAdderData Power] = EstimateMultiAddersPowerConsumption(filter_coefficients,filter_lengths,K,M,Wi,MultType,PlotEachStage)

%
% [MultiAdderData Power] = EstimateMultiAddersPowerConsumption(filter_coefficients,filter_lengths,K,M,Wi,MultType)
%
% This function estimates the number of fulladders for the parallel
% multi-adder tree for the transposed FIR polyphase decimator structures.
%
%   filter_coefficients:    Filter coefficients in each decimation stage
%   filter_lengths:         Lengths of decimation filters for each stage
%   K:                      Number of decimation stages
%   M:                      Deciamtion factor at each decimation stage
%   Wi:                     Input bit width at each decimation stage
%   MultType:               Multipliers types in each decimation stage,
%                           whether 'USigned' or 'Signed '
%   PlotEachStage:          Flag to plot the power conumption for each
%                           indivdual stage
%
%   MultiAdderData:         Internal distribution of adder trees and busses
%   Power:                  Estimated power, which also reflects the number
%                           of gates used by the MOA

MatrixSize = [];
for i = 1 :K,
    MatrixSize(i) = round(filter_lengths(i)/M(i));
end

MultiAdderInfos = struct('Power',[],'Adders',[],'Buses',[],'AddersList',[],'PolyPhaseMatrix',[]);

for i = 1 : K,
    MultiAdderData(i) = MultiAdderInfos;
end
    
for i = 1 : K,
    [MultiAdderData(i).Power MultiAdderData(i).Adders MultiAdderData(i).Buses MultiAdderData(i).AddersList MultiAdderData(i).PolyPhaseMatrix] = EstimateMultiAdders(filter_coefficients(i,1:filter_lengths(i)), M(i), Wi(i), MultType(i)); 
end

for i = 1 : K,
    PowerMatrixSize(i,:) = size(MultiAdderData(1,i).Power);
end

MaxSize = max(max(PowerMatrixSize));

Power = zeros(K,MaxSize);

for i = 1 : K,
    Power(i,1:PowerMatrixSize(i,2)) = MultiAdderData(1,i).Power;
    if PlotEachStage == 1,
        FIG = figure('Name', 'Power Consumption in Adders', 'NumberTitle' , 'off');
        bar(Power(i,1:PowerMatrixSize(i,2)))
    end
end

FIG = figure('Name', 'Power Consumption in Stages', 'NumberTitle' , 'off');
for i = 1 : K,
    bar(i, sum(Power(i,:)),0.6);
    hold on

end
ylabel('Power Consumption - Metric')
xlabel('Stages')
set(gca, 'XTick', [1:K]);


%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function [Power Adders Buses AddersList PolyPhaseMatrix] = EstimateMultiAdders(Coefficients, M, Wi, MultType)

% [Power Adders Buses AddersList] = EstimateMultiAdders(Coefficients, M, Wi, MultType)
% This function to estimate the total number of adders in multiadder trees
% for a polyphase FIR filter.
%   Coefficients:   Filter coefficients of the stage
%   M:              Decimation factor of the stage
%   Wi:             Input bit width of the stage
%   MultType:       The input type of the stage, wehther 'Signed ' or 'USigned'
%

% Multipliers output upper limits
MultUpperLimits = MultiplierUpperLimits(Wi, Coefficients, MultType);
% Multipliers output lower limits
MultLowerLimits = MultiplierLowerLimits(Wi, Coefficients, MultType);

% Multipliers\Coefficients in a polyphase matirx format
PolyPhaseMatrix = PPDMatrix(Coefficients, M);
% Multipliers\Coefficients upper limits in a polyphase matirx format
PolyPhaseUpperMatrix = PPDMatrix(MultUpperLimits,M);
% Multipliers\Coefficients lower limits in a polyphase matirx format
PolyPhaseLowerMatrix = PPDMatrix(MultLowerLimits,M);

% Multiadders input upper limits
MAUpperLimits = MultiAdderLimits(MultUpperLimits, M);
% Multiadders input lower limits
MALowerLimits = MultiAdderLimits(MultLowerLimits, M);

% The size of the polyphase matrix;
% rows      = M
% columns   = number of multiadders
[rows columns] = size(PolyPhaseMatrix);
NumberOfMultiAdders = columns;
NumberOfInputsPerMultiAdder = ExcludeZeroBuses(PolyPhaseMatrix); %[rows (rows+1)*ones(1,columns)];

% Format the multiplier outputs in an appropriate format for multiadder
% inputs
MAT = struct('MAlim', []);
for i = 1 : columns,
    if i == 1,
        for j = 1 : M,
            if PolyPhaseMatrix(j,i)~=0,
                MAT(j,i).MAlim = [PolyPhaseLowerMatrix(j,i) PolyPhaseUpperMatrix(j,i)];
                BI(j) = j; 
                FirstMABusesIndex = VectorWithoutZeros(BI);                 
            end
        end
    else
        for j = 1 : M+1,
            if j == M+1,
                MAT(j,i).MAlim = [MALowerLimits(i-1) MAUpperLimits(i-1)];
            else
                MAT(j,i).MAlim = [PolyPhaseLowerMatrix(j,i) PolyPhaseUpperMatrix(j,i)];
            end
        end
    end
end

%HB
% FlagHB = CheckForHBFilters(Coefficients);
% if M==2 && FlagHB==1,
%     MAT = [MAT(1,:); flipud(MAT(2:3,:))];
%     %AdderElement    = struct('aDepth', cell(28,2), 'bDepth', [], 'oDepth', [], 'aWidth', [], 'bWidth', [], 'oWidth', [], 'Number', []);
% end
AdderElement    = struct('aDepth', [], 'bDepth', [], 'oDepth', [], 'aWidth', [], 'bWidth', [], 'oWidth', [], 'Number', []);

% Get the multiadder info's, such as;
% Adders depth, organization in each depth, recoreds
Adders = zeros(columns,LOG2(rows+1));
Buses  = zeros(columns,LOG2(rows+1));
for i = NumberOfMultiAdders : -1 : 1,
    lim = [];
    if i == 1,
       if NumberOfInputsPerMultiAdder(i)~=0,
            for j = 1 : NumberOfInputsPerMultiAdder(i),
                lim(j,:)=MAT(FirstMABusesIndex(j),i).MAlim;
            end
       else
            lim(j,:) = [0 0];
       end
        [AdderElement1(i,:) a b] = ParallelMultiAdderTree(lim);
        Adders(i,1:length(a))=a;
        Buses(i,1:length(b))=b;
    else
        for j = 1 : NumberOfInputsPerMultiAdder(i),
            lim(j,:)=MAT(j,i).MAlim;
        end
        [a b c] = ParallelMultiAdderTree(lim);
        AdderElement(i-1,1:length(a))=a;
        Adders(i,1:length(b))=b;
        Buses(i,1:length(c))=c;
        %[AdderElement(i-1,:) Adders(i,:) Buses(i,:)] = ParallelMultiAdderTree(lim);
    end
end

% Preparing the adders out widths for power estimation
for i = 1 : NumberOfMultiAdders-1,
    for j = 1 : NumberOfInputsPerMultiAdder(i+1)-1,
        AddersOutputWidth(i,j) = AdderElement(i,j).oWidth;
    end
end

for i = 1 : length(AdderElement1),
    FirstAdderOutputWidth(i) = AdderElement1(i).oWidth;
end

for i = 1 : NumberOfMultiAdders,
    if i == 1,
        Power(i) = sum(FirstAdderOutputWidth);
    else
        Power(i) = sum(AddersOutputWidth(i-1,:));
    end
end

% Just to export the multiadder tree in a readable format
EmptyAdderElement = struct('aDepth', [], 'bDepth', [], 'oDepth', [], 'aWidth', [], 'bWidth', [], 'oWidth', [], 'Number', []);
for i = 1 : size(AdderElement,2)-size(AdderElement1,2),
    Filler(i) = EmptyAdderElement;
end

if FirstAdderOutputWidth == 0,
    for i = 1 : NumberOfMultiAdders-1,
        AddersList(i,:) = AdderElement(i,:);
    end
else
    for i = 1 : NumberOfMultiAdders,
        if i == 1,
            AddersList(i,:) = [AdderElement1 Filler];
        else
            AddersList(i,:) = AdderElement(i-1,:);
        end
    end
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function MultUpperLimits = MultiplierUpperLimits(Wi, Coefficients, MultType)

% MultUpperLimits = MultiplierUpperLimits(Wi, Coefficients, MultType)
% This function estimate the multiplier upper limit by simply multiplying
% the constant by the maximum value for the input bitwidth.
%
%   Wi:             Input bitwidth
%   Coefficient:    Constant or single filter coefficient
%   MultType:       Input type, wether
%                   Signed      -> 'Signed '
%                   Unsigned    -> 'USigned'
%

for i = 1 : length(Coefficients),
    MultUpperLimits(i) = MultiplierUpperLimit(Wi,Coefficients(i), MultType);
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function MultUpperLimit = MultiplierUpperLimit(Wi,Coefficient, MultType)

if strcmp(MultType, 'Signed '), %MultType == 'Signed ',
    if Coefficient > 0,
        MultUpperLimit = (2^(Wi-1)-1)*Coefficient;
    else
        MultUpperLimit = (-2^(Wi-1))*Coefficient;
    end
else
    if Coefficient > 0,
        MultUpperLimit = (2^(Wi)-1)*Coefficient;
    else
        MultUpperLimit = 0;
    end
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function MultLowerLimits = MultiplierLowerLimits(Wi, Coefficients, MultType)

% MultLowerLimits = MultiplierLowerLimits(Wi, Coefficients, MultType)
% This function estimate the multiplier lower limit by simply multiplying
% the constant by the maximum value for the input bitwidth.
%
%   Wi:             Input bitwidth
%   Coefficient:    Constant or single filter coefficient
%   MultType:       Input type, wether
%                   Signed      -> 'Signed '
%                   Unsigned    -> 'USigned'
%

for i = 1 : length(Coefficients),
    MultLowerLimits(i) = MultiplierLowerLimit(Wi,Coefficients(i), MultType);
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function MultLowerLimit = MultiplierLowerLimit(Wi,Coefficient, MultType)

if strcmp(MultType, 'Signed '), %MultType == 'Signed ',
    if Coefficient > 0,
        MultLowerLimit = (-2^(Wi-1))*Coefficient;
    else
        MultLowerLimit = (2^(Wi-1)-1)*Coefficient;
    end
else
    if Coefficient > 0,
        MultLowerLimit = 0;
    else
        MultLowerLimit = (2^(Wi)-1)*Coefficient;
    end
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function MALimits = MultiAdderLimits(MultLimits, M)

% MALimits = MultiAdderLimits(MultLimits, M)
% This function estimate the input limits for parallel multi-adder tree.
%
%   MultLimits: Multiplier output limits, upper or lower
%   M:          Decimation factor for the stage
%
Matrix = PPDMatrix(MultLimits, M);

for i = 1 : size(Matrix, 2),
    if i == 1,
        MALimits(i) = sum(Matrix(:,i));
    else
        MALimits(i) = sum(Matrix(:,i))+MALimits(i-1);
    end
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function Bus = ExcludeZeroBuses(PolyPhaseMatrix)

[rows columns] = size(PolyPhaseMatrix);

for i = 1 : columns,
    z = 1;
    for j = 1 : rows,
        if PolyPhaseMatrix(j,i)==0,
            Nothing = 1;
        else
            B(z,i) = 1;
            z = z + 1;
        end
    end
end

for i = 1 : columns,
    if i == 1,
        Bus(i) = sum(B(:,i));
    else
        Bus(i) = sum(B(:,i))+1;
    end
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function [AdderElement Adders Buses] = ParallelMultiAdderTree(MA_InputLimits)

% Defining structured variable for bus and adder elements
BusElement      = struct('Depth', [], 'UpperLimit', [], 'LowerLimit', [], 'Width', [], 'Number', []);
AdderElement    = struct('aDepth', [], 'bDepth', [], 'oDepth', [], 'aWidth', [], 'bWidth', [], 'oWidth', [], 'Number', []);

    [r c] = size(MA_InputLimits);
    
    % Total number of adders in the multi-adder tree
    NumberOfAdders = r-1;

    % Maximum depth of the multi-adder tree
    AddersDepth = LOG2(r);

    % Number of adders in each single depth
    Wi = r;
    if Wi == 1,
        Adders = 0;
    else
        for i = 1 : AddersDepth,
            if i == AddersDepth,
                Adders(i) = 1;
            else
                if mod(Wi,2) == 0,
                    Wi = Wi/2;
                    Adders(i) = Wi;
                else
                    ExtraAdder = rem(Wi,2);
                    Wi = (Wi-1)/2;
                    Adders(i) = Wi;
                    if ExtraAdder == 1,
                        Wi = Wi + 1;
                    end
                end
            end
        end
    end

    % Number of buses in each single depth
    if Adders == 0,
        Buses = 1;
    else
        for i = 1 : length(Adders),
            if i == 1,
                Buses(i) = length(MA_InputLimits);
            else
                Buses(i) = Adders(i-1);
            end
        end
    end

    % Bus elements assignments
    BusIncrement = 1;
    for i = 1 : length(Buses), 
        BusIncrement = 1;
        for j = 1 : Buses(i),
            if i == 1,
                BusElement(i,j).Depth           = 1;
                BusElement(i,j).UpperLimit      = MA_InputLimits(j,2);
                BusElement(i,j).LowerLimit      = MA_InputLimits(j,1);
                BusElement(i,j).Width           = MaxWidth(MA_InputLimits(j,2), (-1*MA_InputLimits(j,1)));
                BusElement(i,j).Number          = i+j;
            elseif BusElement(i-1,BusIncrement).Depth == BusElement(i-1,BusIncrement+1).Depth,
                BusElement(i,j).Depth           = i;
                BusElement(i,j).UpperLimit      = BusElement(i-1,BusIncrement).UpperLimit + BusElement(i-1,BusIncrement+1).UpperLimit;
                BusElement(i,j).LowerLimit      = BusElement(i-1,BusIncrement).LowerLimit + BusElement(i-1,BusIncrement+1).LowerLimit;
                BusElement(i,j).Width           = MaxWidth((BusElement(i-1,BusIncrement).UpperLimit + BusElement(i-1,BusIncrement+1).UpperLimit), -1*(BusElement(i-1,BusIncrement).LowerLimit + BusElement(i-1,BusIncrement+1).LowerLimit));
                BusElement(i,j).Number          = i+j;
                BusIncrement                    = BusIncrement + 2;
            else
                BusElement(i,j).Depth         = i;
                BusElement(i,j).UpperLimit    = BusElement(i-1,BusIncrement).UpperLimit + BusElement(1,length(MA_InputLimits)).UpperLimit;
                BusElement(i,j).LowerLimit    = BusElement(i-1,BusIncrement).LowerLimit + BusElement(1,length(MA_InputLimits)).LowerLimit;
                BusElement(i,j).Width         = MaxWidth((BusElement(i-1,BusIncrement).UpperLimit + BusElement(1,length(MA_InputLimits)).UpperLimit), -1*(BusElement(i-1,BusIncrement).LowerLimit + BusElement(1,length(MA_InputLimits)).LowerLimit));
                BusElement(i,j).Number        = i+j;
                BusIncrement = BusIncrement + 2; 
            end
        end
    end

    % Transforming the bus elements from matrix format to vector format for
    % easier manuplation
    k = 1;
    for i = 1 : length(Buses),
        for j = 1 : Buses(i),
            BusesVector(k) = BusElement(i,j);
            k = k + 1;
        end    
    end

    % Creating the multi-adder tree from the buses elements
    if Adders == 0,
        AdderElement(i).aDepth  = 0;
        AdderElement(i).bDepth  = 0;
        AdderElement(i).oDepth  = 0;
        AdderElement(i).aWidth  = 0;
        AdderElement(i).bWidth  = 0;
        AdderElement(i).oWidth  = 0;
        AdderElement(i).Number  = 0;        

    else
        AdderIndex = 1;
        for i = 1 : sum(Adders),
            AdderElement(i).aDepth  = BusesVector(AdderIndex).Depth;
            AdderElement(i).bDepth  = BusesVector(AdderIndex+1).Depth;
            AdderElement(i).oDepth  = max(BusesVector(AdderIndex+1).Depth) + 1;
            AdderElement(i).aWidth  = BusesVector(AdderIndex).Width;
            AdderElement(i).bWidth  = BusesVector(AdderIndex+1).Width;
            AdderElement(i).oWidth  = MaxBusWidth(BusesVector(AdderIndex).Width, BusesVector(AdderIndex+1).Width);
            AdderElement(i).Number  = i;
            AdderIndex              = AdderIndex + 2;
        end
    end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%    
function FlagHB = CheckForHBFilters(x)

if mod(length(x),2)==0
    N = length(x)/2;
else
    N = length(x+1)/2;
end

count = 0;
for i = 1 : N,
    if x(i) == 0 && x(i+2) == 0,
        count = count + 1;
    end
end

if count > 2,
    FlagHB = 1;
else
    FlagHB = 0;
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function w = Width(Coefficient)

count = 0;
if Coefficient > 0,
    Value = Coefficient;
else
    Value = -1*Coefficient;
end

while Value > 1,
    count = count + 1;
    Value = Value/2;
end

w = count + 1;

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function W = MaxWidth(A, B)

if A == 0 && B == 0,
    W = 0;
elseif Width(A) == Width(B),
    W = Width(A) + 1;
else
    W = max(Width(A), Width(B));
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function W = MaxBusWidth(A, B)

if A == 0 && B == 0,
    W = 0;
elseif abs(A-B) >= 1,
    W = max(A,B);
else
    W = max(A,B) + 1;
end

%%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
function y = LOG2(x)

 v=x;
 y=0;
 while v>1,
     y=y+1;
     v=v/2;
 end
 
 %%%%%%%%%%%%%%%
% Subfunction %
%%%%%%%%%%%%%%%
 function V = VectorWithoutZeros(Vector)

% V = VectorWithoutZeros(Vector)
% This function used to remove and zeros from the input vector to have
% vector with only intgers larger than zeros, which can be used as index as
% an example.
%
%   Vector: Input vector holding zeros and intgers.
%
%   Example:
%   Vector  = [0 0 3 0 5]
%   V       = [3 5]
%

j = 1;

for i = 1 : length(Vector),
    if Vector(i) ~= 0,
        V(j) = Vector(i);
        j=j+1;
    end
end