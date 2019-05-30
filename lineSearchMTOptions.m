% lineSearchMTOptions - Default options structure lineSearchMT.m
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   result = lineSearchMTOptions
%   lineSearchMTOptions(opt)
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
function result = lineSearchMTOptions(opt)

result = struct(...
        'maxEval',100,...
        'maxTime',inf,...
        'termCondition',1,...
        'minStep',0,...
        'maxStep',inf,...
        'fTol',1E-3,...
        'gTol',.1,...
        'xTol',0,...
        'xTolRel',0,...
        'xTrapL',1.1,...
        'xTrapU',4,...
        'delta',.75,...
        'minIntervalShrink',.8,...
        'minCubicShrink',1E-4,...
        'verbose',true,...
        'interruptible',false);
    
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

end