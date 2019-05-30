% LBFGSOptions - Default options for LBFGS.m
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   result = LBFGSOptions;
%   LBFGSOptions(opt)
%
% PURPOSE:
%   This function returns a default options structure.  If an options structure
%   is passed in, it is tested for validity.
%  
% INPUT:
%   opt         - An existing options structure
% 
% OUTPUT:
%   result      - Default options structure.  If an option struct is passed in
%                 for testing, and it can be fixed, then result will be the
%                 corrected version.
% TODO:
%   Improve options verification
%
% Copyright (C) 2013-2019 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function result = LBFGSOptions(opt)

result = struct(...
        'maxIter',100,...
        'maxEval',250,...
        'maxTime',inf,...
        'relGradTol',1e-8,...
        'gradTol',0,...
        'termFun',[],...
        'initStep',[],...
        'precond',[],...
        'precondInv',[],...
        'blocks',{{}},...
        'precondX0',false,...
        'histLength',17,...
        'adaptiveH0',true,...
        'method','',...
        'lineSearchOpt',lineSearchMTOptions(),...
        'verbose',false,...
        'debug',false,...
        'interruptible',false);
result.lineSearchOpt.maxEval = 7;
result.lineSearchOpt.termCondition = 3;
result.lineSearchOpt.fTol = 1E-4;
result.lineSearchOpt.gTol = .9;
result.lineSearchOpt.verbose = false;
    
if(nargin)
    result = verifyOptions(result,opt);
end
end

function opt = verifyOptions(default,opt)

% Ensure that all the field names are the same
defaultFields = fieldnames(default);
optFields = fieldnames(opt);
msg = 'Unexpected Fields Found';
try
    assert(all(strcmp(sort(defaultFields),sort(optFields))), msg);
catch
    error(msg);
end

% Ensure that the fields are of the correct type
for i = 1:numel(optFields)
    fName = optFields{i};
    fClass = class(default.(fName));
    if( ~strcmp(fClass , class(opt.(fName))) )
        % Attempt to make the cast
        try
            opt.(fName) = cast(opt.(fName) , fClass);
            warning('Automatic option casting occured');
        catch
            error('Unexpected option class');
        end
    end
end

% Ensure histLength is integer and NOT inf (MATLAB bug in mod)

end