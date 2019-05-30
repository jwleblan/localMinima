% lineSearchMT - More Thuente based line search
%{
%-------------------------------------------------------------------------------
% SYNTAX:
%   [alpha, x, f, g, G, exitFlag, numEval, output] = lineSearchMT(...
%       fun, x0, d0, <flags>);
%   [alpha, x, f, g, G, exitFlag, numEval, output] = lineSearchMT(...
%       fun, x0, d0, f0, g0, <flags>);
%
% DESCRIPTION: 
%   This function performs a 1-dimensional search (line-search) using the
%   technique first described by Moré and Thuente.  Where possible, additional
%   improvments from the literature have been made, and the code is structured
%   to support any line-search of the same basic structure.  By default the
%   strong Wolfe termination conditions are used.  In the event that the
%   convergence conditions cannot be met, this function guarantees monotonic 
%   progress.
%   
% INPUT: 
%   fun     - Function handle to the objective function.  This function must
%             have the following prototype:
%               [f,g] = fun(x)
%             where f is a scalar and g is a vector (Nx1).
%
%             Sometimes this can seem inconvenient when the objective functions
%             implementation has additional arguments, or x is not the 1st
%             argument.  This can easily be remedied as follows
%               fun = @(x) myOtherFun(var1,x,var2,var3)
%             where here we have defined fun as a function pointer to
%             myOtherFun.  var1, var2, and var3 will have their values captured
%             when fun is created and these values won't change during fun's
%             existance.
%             
%   x0      - Starting point for the search (Nx1)
%   d0      - Initial search step (Nx1)
%             (NOT necessarily a unit vector)
%   f0      - Objective function value at x0
%   g0      - Gradient of the objective function at x0 (Nx1)
%   flags   - Any of the following flags may be passed in any order.  A
%             structure of the form returned by lineSearchOptions.m may also be
%             used.  Options structures are more efficient, but they should
%             first be verified as follows
%               >> opt = lineSearchMTOptions;        % Get defaults
%                   ... Set your options here ...
%               >> lineSearchMTOptions(opt);         % Verify the options
%
%   --- Termination ---
%       - maxEval           - Maximum number of function evals.  If this is set
%                             to 1 and f0 and g0 are not provided, then 1 eval 
%                             "Away from x0" will be made.
%                             Default: 100
%       - maxTime           - Maximum amount of optimization (CPU) time in
%                             seconds
%                             Default: inf
%       - termCondition     - Termination condition
%                               1 - Strong Wolfe conditions
%                               2 - Armijio and Goldstein lenient Wolfe
%                               3 - Xie and Schlick conditions
%                             Default: 1
%
%                             The strong Wolfe conditions are the most commonly
%                             used condition, and are also the most strict.  The
%                             degree to which they ensure an approximate
%                             minimizer is driven by gTol and fTol.  This
%                             condition ensures that the stopping point lies in
%                             in a convex region of the function.
%                             1)  f(x) <= f0 + fTol*f'(x0)*alpha
%                             2)  abs(g(x)) <= gTol*abs(f'(x0))
%
%                             The Armijio and Goldstein conditions are a common
%                             relaxation of the Wolfe conditions which accept a
%                             broader range of gradients.  This condition also
%                             ensures the stopping point lies in a convex region.
%                             1)  f(x) <= f0 + fTol*f'(x0)*alpha
%                             2)  g(x) >= -gTol*abs(f'(x0))
%                             
%                             The Xie and Schlick conditions allow stopping
%                             points in non-convex regions.  This is the most
%                             lenient set of conditions, and are generally more
%                             efficient when used in conjunction with
%                             large-scale problems that are expensive to
%                             evaluate.
%                             1)  f(x) <= f0 + fTol*f'(x0)*alpha
%                             2)  g(x) >= -gTol*abs(f'(x0))
%                                   or
%                                 g(x) <= -(2-gTol)*abs(f'(x0))
%       - minStep           - Minimum allowed step length
%                             Default: 0
%       - maxStep           - Maximum allowed step length
%                           - Default: inf
%       - fTol              - Function decrease tolerance. 0 < fTol < gTol
%                             In practice fTol is typically less than .5 because
%                             when the objective function is quadratic (the goal
%                             of most preconditioners) then we have
%                             f(x0+alpha*d0) = f(x0) + .5*alpha*f'(x0)
%                             where x0+alpha*d0 is the minimizer.  As fTol gets
%                             smaller the measure of the set of acceptable step
%                             lengths increases.
%                             Default: 1E-3
%       - gTol              - Derivative tolerance.  0 < fTol < gTol < 1.  In
%                             practice gTol is tyically .1 when linsearch is
%                             being used alone, and .9 when used in conjunction
%                             with an encompasing optimization scheme.  As gTol
%                             increases, so to does the measure of acceptable
%                             step lengths.
%                             Default: .1
%       - xTol              - Interval length tolerance.  Termination will
%                             occur if the local minimizer is determined up to
%                             this amount.  This is not a relative quantity as
%                             was originally proposed.
%                             Default: 0
%       - xTolRel           - Relative interval tolerance.  This is what was
%                             originally designated xTol.  Termination will
%                             occur if the current search interval is reduced by
%                             a fraction that is smaller than this value.  In
%                             practice 0 < xTol < .5
%                             Default: 0
%       - interruptible     - Use the "interruptible" utility to allow the user
%                             to terminate the search
%                             Default: false
%
%   --- Search ---
%       - xTrapL            - Lower extrapolation limit.  When safeguarded steps
%                             are needed, this is used to determine the minimium
%                             extrapolation distance. 1 < xTrapL < xTrapU
%                             Default: 1.1
%       - xTrapU            - Upper extrapolation limit.  When safeguarded steps
%                             are needed, this is used to determine the maximum
%                             extrapolation distance. xTrapU > xTrapL
%                             Default: 4
%       -delta              - This constant limits extrapolation steps when the
%                             minimizer has been bracketed, but the remaining
%                             search interval contains a non-convex segment.
%                             When this occurs the trial step length becomes
%                             delta times the remaining search interval or the
%                             extrapolation depending on which is more
%                             conservative.  0 < delta < 1.  In practice
%                             .5 < delta < 1
%                             Default: .75
%       - minIntervalShrink - When the minimimum has been bracketed, we require
%                             that each step reduce the bracket width by this
%                             fraction.  As this value increases, we become
%                             increasingly trusting of the local function
%                             approximations.  As this value decreases, we tend
%                             toward bisection steps. .5 < minIntervalShrink < 1
%                             Default: .8
%       - minCubicShrink    - Limits how quickly a cubic interpolant may shrink
%                             the search interval.  A cubic can become a poor
%                             approximation to the actual underlying function
%                             causing very small steps to be taken.  This is the
%                             minimum fraction of the existing interval length
%                             that a cubic interpolant can shrink to.
%                             (Addresses the Xie and Schlick concern)
%                             Default: 1E-4
%
%   --- Output ---
%       -verbose            - If passed, verbose output is provided.
%                             Default: false
%       
% RETURN VALUES: 
%   alpha       - Step length that satifies the termination conditions
%   x           - x0 + alpha * d0
%   f           - Function value at x
%   g           - Gradient at x
%   G           - Directional derivative at x
%   exitFlag    - Flag giving the exiting conditons
%                   0 - Termination condition was met
%                   1 - maxEval was reached
%                   2 - maxTime was reached
%                   3 - xTol was satisfied
%                   4 - xTolRel was satisfied
%                   5 - minStep was reached
%                   6 - maxStep was reached
%                   7 - round-off error prevented progress
%                   8 - User interrupt
%   numEval     - Number of function evaluations made
%   output      - A detailed structure describing what happened.  This is useful
%                 for debugging. 
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
% TODO:
%   Address the "negAlpha" question for dealing with local minima
%
% Copyright (C) 2013-2019 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function [alpha,xk,fk,gk,Gk,exitFlag,numEval,output] = lineSearchMT(fun,x0,d,varargin)

% Parse Inputs
[d,config,opt,output] = parseInputArguments(fun,x0,d,nargout,varargin{:});
alpha = 1;
updateMethod = 6;
methodArray ={'bisect','interp','interp','extrap','extrap','bisect','init'};

while(1)
    %%%%%%%%%%%%%%%%
    %%% Evaluate %%%
    %%%%%%%%%%%%%%%%
    xk = x0 + alpha*d;
    [fk,gk] = fun(xk);
    config.numEval = config.numEval + 1;
    Gk = gk.'*d;
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Update Display %%%
    %%%%%%%%%%%%%%%%%%%%%%
    if opt.verbose
        fprintf(' Eval: %-3g |  step = %-8.2g |  f = %-11.5g |  g = %-11.5g  | %s\n',...
            config.numEval,alpha,fk,Gk,methodArray{updateMethod+1});
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Check For Termination %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [termFlag,exitFlag] = checkTermination(alpha,fk,Gk,config,opt);
    if termFlag
        % Ensure monotonically non-increasing
        if exitFlag > 0 && ( fk>config.fL || ~isfinite(fk) )
            alpha = config.alphaL;
            xk = x0 + alpha*d;
            fk = config.fL;
            gk = config.gL;
        end
        break
    end
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get Safeguarded Step %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [alphaNew,config,updateMethod] = getSafeguardedAlpha(alpha,fk,gk,Gk,...
        config,opt);
    
    
    %%%%%%%%%%%%%%
    %%% Update %%%
    %%%%%%%%%%%%%%    
    if config.buildOutput
        output.iter(end+1) = struct(...
            'alpha',alpha,...
            'fk',fk,...
            'Gk',Gk,...
            'alphaL',config.alphaL,...
            'fL',config.fL,...
            'GL',config.gL,...
            'alphaU',config.alphaU,...
            'fU',config.fU,...
            'GU',config.GU,...
            'bracketed',config.bracketed,...
            'stage',config.stage,...
            'updateMethod',updateMethod);
    end
    alpha = alphaNew;
    
end

%%% Cleanup %%%
numEval = config.numEval;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUB-FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{   
% SYNTAX:
%   [d,config,opt,output] = parseInputArguments(fun,x0,d,numArgOut,...
%       varargin{:})
%
% DESCRIPTION: 
%   This function parses the input arguments.
%
% INPUTS:
%   See lineSearchMT.m
%
% RETURN VALUES:
%   d       - Unit search direction
%   config  - Configuration structure
%       .T0             - Time we started
%       .numEval        - Number of function evaluations used
%       .f0             - Initial objective function value
%       .G0             - Initial directional derivative
%       .alphaL         - "Lower" alpha.  Best alpha in the bracket so far.
%       .fL             - Function value at alphaL (smaller than fU by definition)
%       .gL             - Gradient at alphaL
%       .GL             - Directional derivative at alphaL
%       .alphaU         - Upper end of the bracketing window
%       .fU             - function value at alphaU
%       .GU             - Directional derivative at alphaU
%       .stage          - In stage 1 we use a surrogage objective function, and
%                         in stage 2 we use the actual function
%       .bracketed      - True if a solution is bracketed between alphaL and
%                         alphaU
%       .dUnit          - d/norm(d)
%       .fTolG0         - fTol*G0
%       .gTolG0         - gTol*G0
%       .minSubStep     - Smallest step for extrapolation.  When bracketed, this
%                         becomes the larger of alphaL and alphaU
%       .maxSubStep     - Largest step for extrapolation.  When bracketed, this
%                         becomes the larger of alphaL and alphaU
%       .builtOutput    - True if output should be populated
%       .negAlpha       - True if G0 was positive
%       
%       
%   opt     - Options structure
%   output  - Output structure for 1st iterate or empty
%   
%}
function [d,config,opt,output] = parseInputArguments(fun,x0,d,...
    numArgOut,varargin)

%%%%%%%%%%%%%%%%%%%%%%
%%% Default Values %%%
%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;               % Do time sensative stuff first
opt = lineSearchMTOptions;
config = struct(...
    'T0',T0,...
    'numEval',0,...
    'f0',[],...
    'G0',[],...
    'alphaL',0,...
    'fL',[],...
    'gL',[],...
    'GL',[],...
    'alphaU',[],...
    'fU',[],...
    'GU',[],...
    'stage',1,...
    'bracketed',false,...
    'fTolG0',[],...
    'gTolG0',[],...
    'minSubStep',[],...
    'maxSubStep',[],...
    'buildOutput',numArgOut>6,...
    'negAlpha',false);  % True if the user started us pointing the wrong way
numArgs = numel(varargin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse the variable arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We expect the numeric types to be f0 then g0, followed by strings (flags),
% or an options structure
numericCount = 0;
i = 1;
while i <= numArgs
    arg = varargin{i};
    
    if isnumeric(arg)                   % We found dim or axis
        i = i+1;
        switch numericCount
            case 0
                config.f0 = arg;
            case 1
                config.gL = arg;        % This is holding g0
            otherwise
                error('Unexpected numeric argument');
        end
        numericCount = numericCount + 1;
    
    else                                % We're done with numerics onto flags
        break;
    end
end

while i <= numArgs
    arg = varargin{i};
    i = i+1;
    
    switch class(arg)
        case 'char'
            if arg(1)=='-'
                compArg = lower(arg(2:end));
            else
                compArg = lower(arg);
            end
            switch compArg
                case 'maxeval'
                    opt.maxEval = varargin{i};
                    i = i+1;
                case 'maxTime'
                    opt.maxTime = varargin{i};
                    i = i+1;                    
                case 'termCondition'                    
                    opt.termCondition = varargin{i};
                    i = i+1;
                case 'minStep'                    
                    opt.minStep = varargin{i};
                    i = i+1;
                case 'maxStep'                    
                    opt.maxStep = varargin{i};
                    i = i+1;
                case 'ftol'
                    opt.fTol = varargin{i};
                    i = i+1;
                case 'gtol'
                    opt.gTol = varargin{i};
                    i = i+1;                    
                case 'xtol'
                    opt.xTol = varargin{i};
                    i = i+1;                    
                case 'xtolrel'
                    opt.xTolRel = varargin{i};
                    i = i+1;
                case 'xtrapl'
                    opt.xTrapL = varargin{i};
                    i = i+1;                    
                case 'xtrapu'
                    opt.xTrapU = varargin{i};
                    i = i+1;                    
                case 'delta'
                    opt.delta = varargin{i};
                    i = i+1;
                case 'minintervalshrink'
                    opt.minIntervalShrink = varargin{i};
                    i = i+1;                    
                case 'mincubicshrink'
                    opt.minCubicShrink = varargin{i};
                    i = i+1;                                     
                case 'verbose'
                    opt.verbose = true;
                otherwise
                    error('Unexpected flag');
            end  
            
        case 'struct'
            % Apply the options without checking them for correctness
            %opt = lineSearchMTOptions(arg);         % <-- This would check them
            opt = arg;
            
        otherwise
            error('Unexpected argument class');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get Initial Values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure we have f and g for x0
if isempty(config.f0) || isempty(config.gL)
    [config.f0,config.gL] = fun(x0);
    config.numEval = config.numEval + 1;
end
config.G0 = config.gL.'*d;
if config.G0 > 0
    config.G0 = -config.G0;
    d = -d;
    config.negAlpha = true;
end


%%%%%%%%%%%%%%%%%%%%%
%%% Update Config %%%
%%%%%%%%%%%%%%%%%%%%%
config.fL = config.f0;
config.GL = config.G0;
config.alphaU = opt.maxStep;
config.fTolG0 = opt.fTol*config.G0;
config.gTolG0 = opt.gTol*config.G0;
config.minSubStep = max([opt.xTrapL+1,opt.minStep]);
config.maxSubStep = min([opt.xTrapU+1,opt.maxStep]);


%%%%%%%%%%%%%%%%%%%%
%%% Setup Output %%%
%%%%%%%%%%%%%%%%%%%%
if config.buildOutput
    output = struct(...
        'opt',opt,...
        'numEval',[],...
        'cpuTime',[],...
        'x0',x0,...
        'd',d,...
        'f0',config.f0,...
        'G0',config.G0,...
        'iter',struct(...
            'alpha',{},...
            'fk',{},...
            'Gk',{},...
            'alphaL',{},...
            'fL',{},...
            'GL',{},...
            'alphaU',{},...
            'fU',{},...
            'GU',{},...
            'bracketed',{},...
            'stage',{},...
            'updateMethod',{}) );
else
    output = struct([]);
end


%%%%%%%%%%%%%%%%%%%%%%
%%% Update Display %%%
%%%%%%%%%%%%%%%%%%%%%%
if opt.verbose
    fprintf(' Eval: %-3g |  step = %-8.2g |  f = %-11.5g |  g = %-11.5g  | %s\n',...
            config.numEval,0,config.f0,config.G0,'init');
end

end


%{   
% SYNTAX:
%   [result,exitFlag] = checkTermination(alpha,fk,Gk,config,opt)
%
% DESCRIPTION: 
%   This function checks for termination
%
% RETURN VALUES:
%   result      - True if we should terminiate and false otherwise
%   exitFlag    - Exit value (see "exitFlag" in lineSearchMT)
%   
%}
function [result,exitFlag] = checkTermination(alpha,fk,Gk,config,opt)

%%%%%%%%%%%%%%%%%%%%%%
%%% Default Values %%%
%%%%%%%%%%%%%%%%%%%%%%
result = true;          % Be optimistic
exitFlag = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for non-finite return values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isinf(fk) && fk<0
    % We are done by any reasonable measure
    exitFlag = 0;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Primary Objective %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fThresh = config.f0 + config.fTolG0*alpha;
switch opt.termCondition
    case 1  % Strong Wolfe
        if (fk <= fThresh) && (abs(Gk) <= abs(config.gTolG0))
            exitFlag = 0;
            return
        end
        
    case 2  % Armijio and Goldstein
        if (fk <= fThresh) && (Gk >= -abs(config.gTolG0))
            exitFlag = 0;
            return
        end
        
    case 3  % Xie and Schlick
        if (fk <= fThresh) && ...
                (Gk >= -abs(config.gTolG0) || Gk <= (opt.gTol-2)*abs(config.G0) )
            exitFlag = 0;
            return
        end
        
    otherwise
        error('Unhandled termination condition');
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for timeout %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if config.numEval >= opt.maxEval
    exitFlag = 1;
    return
end

if toc(config.T0) >= opt.maxTime
    exitFlag = 2;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check interval conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if config.bracketed            % TODO: Test this part
    width = config.maxSubStep-config.minSubStep;
    % xTol
    if width < opt.xTol
       exitFlag = 3;
       return
    end
    
    % xTolRel
    if width < opt.xTolRel*config.maxSubStep
        exitFlag = 4;
        return
    end
    
    % Stagnated
    if alpha<=config.minSubStep || alpha>=config.maxSubStep
        exitFlag = 7;
        return
    end
    
else
    if alpha<=opt.minStep
        exitFlag = 5;
        return
    elseif alpha>=opt.maxStep
        exitFlag = 6;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%%% User Interrupt %%%
%%%%%%%%%%%%%%%%%%%%%%
if opt.interruptible && interruptible()
    exitFlag = 8;
    return
end

result = false;
end


%{   
% SYNTAX:
%   [result,config,updateMethod] = getSafeguardedAlpha(alpha,fk,gk,Gk,config,opt)
%
% DESCRIPTION: 
%   Determines a safeguarded alpha and updates the search interval according to
%   the new information.
%
% RETURN VALUES:
%   result          - A safeguarded step about alpha
%   config          - Updated config structure
%   updateMethod    - How we made the update (1-4)
%   
%}
function [result,config,updateMethod] = getSafeguardedAlpha(alpha,fk,gk,Gk,...
    config,opt)

%%%%%%%%%%%%%%%%%%%%%%
%%% Default Values %%%
%%%%%%%%%%%%%%%%%%%%%%
alphaL = config.alphaL;
alphaU = config.alphaU;
bracketed = config.bracketed;
fThresh = config.f0 + config.fTolG0*alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update algorithm stage %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fk <= fThresh) && (Gk > 0)
    config.stage = 2;
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup Surrogate %%%
%%%%%%%%%%%%%%%%%%%%%%%
useSurrogate = config.stage==1 && fk<=config.fL && fk>fThresh;
if useSurrogate
    fL = config.fL - fThresh;
    GL = config.GL - config.fTolG0;
    
    fU = config.fU - fThresh;
    GU = config.GU - config.fTolG0;
    
    fk = fk - fThresh;
    Gk = Gk - config.fTolG0;
    
else
    fL = config.fL;
    GL = config.GL;
    
    fU = config.fU;
    GU = config.GU;
    
end

if ~isfinite(fk) || ~isfinite(Gk)
    %%%%%%%%%%%%%%
    %%% Case 0 %%%
    %%%%%%%%%%%%%%
    % "Bad" return values are denying us information
    % Bisect toward the best known point
    result = (alphaL + alpha)/2;
    updateMethod = 0;
    bracketed = true;
    
elseif fk > fL          % TODO:  These need to be looked at with "negAlpha"
    %%%%%%%%%%%%%%
    %%% Case 1 %%%
    %%%%%%%%%%%%%%
    % We have a function increase so a local minimizer has been "bracketed"
    
    % Get a quadratic interpolate
    alphaQ = quadraticMinimizer(alphaL,fL,GL,alpha,fk,NaN);
    
    % Get a cubic interpolate
    alphaC = cubicMinimizer(alphaL,fL,GL,alpha,fk,Gk);
    
    if isnan(alphaC)
        result = alphaQ;
        
    else
        % Prevent overshrinking
        if (alphaC-alphaL)/(alpha-alphaL) < opt.minCubicShrink
            alphaC = alphaL + opt.minCubicShrink*(alpha-alphaL);
        end
    
        % Choose between the steps
        if abs(alphaC-alphaL) < abs(alphaQ-alphaL)
            result = alphaC;
        else
            result = (alphaC+alphaQ)/2;
        end
    end
        
    updateMethod = 1;
    bracketed = true;
    
elseif GL*Gk <= 0
    %%%%%%%%%%%%%%%
    %%% Case II %%%
    %%%%%%%%%%%%%%%
    % We have a derivative sign change so the minimizer has been "bracketed"
    
    % Get a cubic interpolate
    alphaC = cubicMinimizer(alphaL,fL,GL,alpha,fk,Gk);
    
    % Get a quadratic interpolate (secant)
    alphaS = quadraticMinimizer(alpha,fk,Gk,alphaL,fL,NaN);
    
    %Choose the best step
    if abs(alphaC-alpha) >= abs(alphaS-alpha)
        result = alphaC;
    else
        result = alphaS;
    end
    
    updateMethod = 2;
    bracketed = true;
    
elseif abs(Gk) <= abs(GL)
    %%%%%%%%%%%%%%%%
    %%% Case III %%%
    %%%%%%%%%%%%%%%%
    % A cubic minimizer may not exist or may be in the wrong direction.  We will
    % only use a cubic if it's minima is beyond alpha, or tends to inifinity in
    % the direction of the step
    
    % Get a quadratic interpolate (secant)
    alphaS = quadraticMinimizer(alphaL,fL,GL,alpha,NaN,Gk);
    
    [alphaC,cubicCoeff] = cubicMinimizer(alphaL,fL,GL,alpha,fk,Gk);
    
    temp = alphaC-alphaL;
    
    % If alphaC is in the right direction and beyond alpha
    if temp/(alpha-alphaL)>1
        
        % If the cubic tends to inifinity in the cubic steps direction
        if temp*cubicCoeff > 0
            % Keep the cubic step
        else
            alphaC = config.maxSubStep;
        end
    else
        alphaC = alphaS;
    end
    
    % Choose the best step
    if config.bracketed
        if abs(alphaC-alpha) < abs(alphaS-alpha)
            result = alphaC;
        else
            result = alphaS;
        end
        
        if alpha > alphaL
            result = min([result,alpha + opt.delta*(alphaU-alpha)]);
        else
            result = max([result,alpha + opt.delta*(alphaU-alpha)]);
        end
        
    else
        if abs(alphaC - alpha) < abs(alphaS - alpha)
            result = alphaC;
        else
            result = alphaS;
        end
        result = min([result,config.maxSubStep]);
        result = max([result,config.minSubStep]);
        
    end
    updateMethod = 3;
    
else
    %%%%%%%%%%%%%%%
    %%% Case IV %%%
    %%%%%%%%%%%%%%%
    % Objective function decreased with derivatives of the same sign, and the
    % magnitude of the derivative did not decrease.  We need to extrapolate.
    if config.bracketed
        result = cubicMinimizer(alpha,fk,Gk,alphaU,fU,GU);
        
    elseif alpha>alphaL
        result = config.maxSubStep;
        
    else
        result = config.minSubStep;
        
    end
    updateMethod = 4;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Keep the step valid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = max([result, opt.minStep]);
result = min([result, opt.maxStep]);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update Endpoints %%%
%%%%%%%%%%%%%%%%%%%%%%%%
if updateMethod == 0
    % Do not update into a "bad" point
    
elseif updateMethod == 1
    config.alphaU = alpha;
    config.fU = fk + useSurrogate*fThresh;
    config.GU = Gk + useSurrogate*config.fTolG0;

else
    if updateMethod == 2
        config.alphaU = alphaL;
        config.fU = config.fL;
        config.GU = config.GL;
        
    end
    
    config.alphaL = alpha;
    config.fL = fk + useSurrogate*fThresh;
    config.gL = gk;
    config.GL = Gk + useSurrogate*config.fTolG0;
    
end

% We sort these now to simiplify the code that follows
alphaL = min([config.alphaL config.alphaU]);
alphaU = max([config.alphaL config.alphaU]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bisect if necessary %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if config.bracketed         % If we we're bracketed before (and we are now)
    oldWidth = config.maxSubStep-config.minSubStep;
    maxAllowedWidth = opt.minIntervalShrink*oldWidth;
    [widthSort,widthInd] = sort([result-alphaL alphaU-result]);
    
    if widthSort(2)>maxAllowedWidth && widthSort(1)<maxAllowedWidth;
        updateMethod = 5;
        % config.maxSubStep-maxAllowedWidth to config.minSubStep+maxAllowedWidth
        currentWidth = alphaU-alphaL;
        maxAllowedWidth = opt.minIntervalShrink*currentWidth;
        if widthInd(2)==1
            result = max([alphaL+maxAllowedWidth alphaU-maxAllowedWidth]);
        else
            result = min([alphaL+maxAllowedWidth alphaU-maxAllowedWidth]);
        end
    end
end
config.bracketed = bracketed;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Substep Limits %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if bracketed
    config.minSubStep = alphaL;
    config.maxSubStep = alphaU;
else
    temp = result - config.alphaL;
    config.minSubStep = result + opt.xTrapL*temp;
    config.maxSubStep = result + opt.xTrapU*temp;
end

end
