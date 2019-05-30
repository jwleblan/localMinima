% kFoldPartition - Partition data for k to 1 cross validation
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   inds = kFoldPartition(labels,K)
%
% PURPOSE:
%   This function computes the test groups for k to 1 cross validation based on
%   stratefied data.  The data stratification is described by labels, a
%   partition of the total number of samples into a finite number of (G) groups.
%   This function provides randomized and stratified test-group assignment by
%   returning a vector of indicies the size of labels that assigns each member
%   of the sample population to one of the K test groups.
%   Specifically, this function ensures
%       1) All testing partitions maintain the stratification proportions.  The
%          sum of the absolute differences between the number of elements from
%          each stratification (sub-group described by labels) is minimized in
%          the groupings provided by inds.
%       2) The difference between the sizes of the final test groups is
%          minimized
%       3) Within each stratification, allocation to the testing partitions is
%          randomized.
%  
% INPUT:
%   labels      - Array of logicals (samples stratified into two groups) or
%                 integers (G groups) that gives the population stratification.
%   K           - Positive integer that indicates the number of test groups. 
%                 For example, the standard 80/20 training/testing rule is 
%                 represented by K = 5 (4 to 1 = 5 total groups).  This is
%                 sometimes called a 5-fold partition.
%                   Default: 5
% 
% OUTPUT:
%   inds        - Random partition of indicies associated with the labels.  Each
%                 member of of inds is an natural number in [1 G].
%
% NOTES:
%   This function alters the state of the random number generator.  Fixing the
%   seed (see rng.m) will made the output deterministic.
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function inds = kFoldPartition(labels,K)

% Parse Inputs
if nargin<2
    K = 5;
end

% Default Values
strataMissingWarning = false;
N = numel(labels);
inds = zeros(size(labels));

% Check for an inadequate number of labels
if K > N
    error([mfilename,':insufficientData'],...
        'Insufficient data for %g-fold partition',K);
end

% Determine the stratification
if islogical(labels)
    G = 2;
    gVals = [false true];
else
    gVals = unique(labels);
    G = numel(gVals);
end

% Iterate allocating the population to test groups
startInd = 0;       % Index to start test-group assignment in [0 G-1]
for g = 1:G
    % Get group inidicies
    gInds = labels==gVals(g);
    gNum = sum(gInds(:));
    
    % Warn if each test-group will not contain at least one member from each
    % population strata
    if ~strataMissingWarning && gNum<K
        warning([mfilename,':insufficientData'],...
            ['Insufficient data to represent each population strata in all ',...
            'test groups']);
    end
    
    % Randomize the indices and allocate them to bins 1:K in sequence.  This is
    % achieved by allocating them to bins in sequence and then randomizing the
    % how the bins are mapped onto the population strata
    bins = mod(startInd:startInd+gNum-1,K)+1;
    inds(gInds) = bins(randperm(gNum));
    
    % Update startInd
    startInd = mod(startInd+gNum,K);
end

end