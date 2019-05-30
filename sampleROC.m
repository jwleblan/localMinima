% sampleROC - Computes the sample ROC curve given data and a detector
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   [PD,PFA,S,testGroups] = sampleROC(data,labels,trainFun,KTo1)
%   [PD,PFA,S] = sampleROC(data,labels,detFun)
%   [PD,PFA] = sampleROC(S,labels)
%
% PURPOSE:
%   This function trains a detector, and uses it to compute a sample ROC curve
%   based on data and its associated labels.  A round-robin technique is used to
%   ensure that all of the data is used for both training and testing.  If a
%   detector has already been trained, this function can be used to produce an
%   ROC curve assuming all the data is for training.  If the detector has
%   already been run, this function can be used to compute PD and PFA.
%  
% INPUT:
%   data        - [P N] Data matrix to train and test over.  This may also be
%                 any more general object that can be linearly indexed as
%                 follows someSubset = data(inds) where inds is an Nx1 logical
%                 array and N = length(data)
%   labels      - [1 N] logical array that is true under H1
%   detFun      - A function that implements a detector with the following
%                 prototype:
%                   [phi,S] = detFun(dataSubset)
%                 where phi is an Lx1 logical array that is true under H1, and
%                 S is an Lx1 floating-point array of test statistics;
%   trainFun    - A function with the following prototype:
%                   detFun = trainFun(dataSubset,labelSubset)
%                 where dataSubset is used for training, and detFun is defined
%                 above.
%   Kto1        - Positive integer that indicates the proportion of training
%                 data for each part of testing data.  For example, the standard
%                 80/20 training-testing rule is represented by Kto1 = 4.
%                   Default: 4
%   S           - [1 N] Test statistic array associated with data
% 
% OUTPUT:
%   PD          - [1 N] array of sample probabilities of detection associated
%                 with PFA
%   PFA         - [1 N] array of sample probabilities of false alarm associated
%                 with PD
%   S           - [1 N] Test statistic array associated with data
%   testGroups  - [1 N] Array if indicies between 1 and Kto1+1 that gives
%                 the testing group associated with data and S
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function [PD,PFA,S,testGroups] = sampleROC(data,labels,varargin)

% Default values
N = numel(labels);
testGroups = [];

% Parse Inputs
[data,labels,detFun,trainFun,KTo1,maxDataVals] = parseInputs(data,labels,...
    varargin{:});
trainFunProvided = ~isempty(trainFun);
detFunProvided = ~isempty(detFun);

% Get the test groups
if trainFunProvided
    testGroups = kFoldPartition(labels,KTo1+1);
    [groupEnds,groupInds] = sort(testGroups);
    [~,groupEnds] = unique(groupEnds,'last');
    groupEnds = [0 ; groupEnds];
    
    % Test against each group
    parfor k = 1:(KTo1+1)
        testMask = groupInds(groupEnds(k)+1:groupEnds(k+1));
        trainMask = groupInds([1:groupEnds(k) , groupEnds(k+1)+1:N]);
        
        % Train the detector
        detFun = trainFun(data(trainMask),labels(trainMask));
        
        % Compute the test statistic
        [~,SCell{k}] = detFun(data(testMask));
    end
    
    % Convert S back into a correctly sorted array
    S = zeros(1,N);
    if any(cellfun(@(x)size(x,1),SCell) > 1)    % We have vectors
        S(groupInds) = cat(1,SCell{:});
    else                                        % We have arrays
        S(groupInds) = cat(2,SCell{:});
    end
    clear SCell
    
else
    if detFunProvided       % We got detFun instead
        if isempty(maxDataVals)         % We don't need to worry about memory
            [~,S] = detFun(data(1:N));
            if size(S,1) > 1
                S = S.';
            end
            
        else                            % Perform blocking
            % Determine the groups
            groupEnds = 0:maxDataVals:N;
            if groupEnds(end)~=N
                groupEnds(end) = N;
            end
            
            % Generate the test statistics
            parfor k = 1:numel(groupEnds)-1
                [~,SCell{k}] = detFun(data(groupEnds(k)+1:groupEnds(k+1)));
            end
            
            % Convert S back into a correctly sorted array
            if any(cellfun(@(x)size(x,1),SCell) > 1)    % We have vectors
                S = cat(1,SCell{:}).';
            else                                        % We have arrays
                S = cat(2,SCell{:});
            end
            clear SCell
            
        end
        
    else                    % We got S and labels
        S = data;
    end
end


% Compute PD and PFA
[~,ind] = sort(S,'descend');
numH1 = sum(labels);
numH0 = N-numH1;
PD = cumsum(labels(ind))/numH1;
PFA = cumsum(~labels(ind))/numH0;

% The test function needn't be monotonically increasing so test for ROC curves
% below the chance line and invert the rule if necessary
[PFAUnique,ind] = unique(PFA,'last');   % Be conservative in flipping
if trapz(PFAUnique,PD(ind)) < .5
    temp = PD;
    PD = PFA;
    PFA = temp;
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   [data,labels,detFun,trainFun,KTo1,maxDataVals] = parseInputs(data,labels,...
%       varargin)
%
% PURPOSE:
%   This function parses the inputs and provides a consistent interface.
%  
% INPUT:
%   See sampleROC.m
% 
% OUTPUT:
%   data        - Linearly indexable object
%   labels      - [1 N] logical array of labels
%   detFun      - Detector function or [] if training is requried
%   trainFun    - Training function or [] if a detector is provided
%   KTo1        - Integer giving proportion of training to testing data
%   maxDataVals - Maximum number of data values that can be simultaneously
%                 requested within the memory limits
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function [data,labels,detFun,trainFun,KTo1,maxDataVals] = parseInputs(data,labels,...
    varargin)

% Default Values
memLimit = 1E9;     % In bytes
N = numel(labels);
detFun = [];
trainFun = [];
KTo1 = [];
maxDataVals = [];

% Parse out the two calling syntaxes
switch nargin
    case 2
        % S and labels provided
        
    case 3
        detFun = varargin{1};
        
        % Check for test statistics being provided
        numDetOutputs = nargout(detFun);
        if any(numDetOutputs==[0 1])
            error([mfilename,':invalidDetector'],...
                'Detector function must provide a test statistic');
        end
    case 4
        [trainFun,KTo1] = varargin{1:2};
        if isempty(KTo1)
            KTo1 = 4;   % Default to 80/20
        end
        
    otherwise
        error([mfilename,':invalidInput'],'Unexpected number of inputs');
end

% Cast matrix data to a linearly indexable object
if ismatrix(data) && size(data,1) > 1
    data = dataMatrix2LinearIndexObj(data);
else
    % It is already a linear indexible object, determine maxDataVals
    dataElement = data(1);
    dataElementInfo = whos('dataElement');
    maxDataVals = floor(memLimit/dataElementInfo.bytes);
    
    if maxDataVals >= N
        % We don't need to worry about this
        maxDataVals = [];
    end
end

% Be nice and help people with arrays vs vectors
labels = logical(labels(:)).';

% If KTo1 is too large, issue a warning and reduce it
if ~isempty(KTo1) && KTo1>=N
    warning([mfilename,':KTo1TooLarge'],...
        'Data does not support KTo1.  It is being automatically reduced.');
    KTo1 = N-1;
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   dataObj = dataMatrix2LinearIndexObj(data)
%
% PURPOSE:
%   Helper function to cast data matricies to linearly indexable objects.
%  
% INPUT:
%   data        - [P N] data matrix
% 
% OUTPUT:
%   dataObj     - Object that accepts a [1 K] array of linear indices and
%                 returns a [P K] matrix of data values.
%-------------------------------------------------------------------------------
%}
function dataObj = dataMatrix2LinearIndexObj(data)
dataObj = @(ind)data(:,ind);
clear('data');              % Don't allow the anonymous function to bind this
end
