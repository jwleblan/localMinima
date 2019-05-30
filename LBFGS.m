% LBFGS - Limited memory quasi-Newton search
%{
%-------------------------------------------------------------------------------
% SYNTAX:
%   [x,f,g,exitFlag,output] = LBFGS(fun,x0,<flags>)
%
% DESCRIPTION: 
%   This function performs a limited-memory quasi-newton search.  The function
%   guarantees a monotonically non-increasing search, and tries to identify a
%   local minima.
%
% INPUT: 
%   fun     - Function handle to the objective function.  This function must
%             have the following prototype:
%               [f,g] = fun(x)
%             where f is a scalar and g is a vector (Nx1).  A function pointer
%             can be used if additional input arguments are needed.
%   x0      - Starting point for the search (Nx1)  This must be a double or a
%             single.
%   flags   - Any of the following flags may be passed in any order.  A
%             structure of the form returned by lineSearchOptions.m may also be
%             used.  Options structures are more efficient, but they should
%             first be verified as follows
%               >> opt = LBFGSOptions;               % Get defaults
%                   ... Set your options here ...
%               >> LBFGSOptions(opt);                % Verify the options
%
%   --- Termination ---
%       - maxIter           - Maximum number of search iterations.  An iteration
%                             is defined as a direction update followed by a
%                             line search.  A termination check is made, and if
%                             convergence has not occured, a quasi-newton update
%                             is made followed by another iteration.
%                             Default: 100
%       - maxEval           - Maximum number of function evals.  If this is set
%                             to 1, then 1 eval "Away from x0" will be made.
%                             Default: 250
%       - maxTime           - Maximum amount of optimization (CPU) time in
%                             seconds
%                             Default: inf
%       - relGradTol        - Relative gradient tolerance until convergence
%                             Default: 1e-8
%       - gradTol           - Absolute gradient tolerance until convergence
%                             Default: 0
%       - interruptible     - Use the "interruptible" utility to allow the user
%                             to terminate the search
%                             Default: false
%       - termFun           - Custom defined function pointer that can be used
%                             to terminate the search.  The prototype must be
%                               exitFlag = termFun(x,f,g,xm1,fm1,gm1,...
%                                   lineExitFlag,restartDirection,config,opt,...
%                                   removePrecondFun)
%                             where removePrecondFun is a function that removes
%                             any applied preconditioner and has the prototype
%                                   [x,g] = removePrecond(x,g,config)
%                             Default: []
%
%   --- Search ---
%       - initStep          - Initial step.  This may be specified a few ways:
%                             (1) If left empty, the initial step will be half
%                                 the initial gradient.  This is optimal when
%                                 the objective function is a quadratic and a 
%                                 good choice when it is well conditioned
%                             (2) If a scalar, the initial step will be the
%                                 length of the step taken along the steepest 
%                                 descent direction.
%                             (3) Vector; the expected deviation of each of the
%                                 parameters.
%                             Default: []
%       - histLength        - Length of the step history used to estimate the
%                             action of the inverse hessian
%                             Default: 17
%       - method            - How the search directions will be established
%                              'lbfgs'  - Limited memory BFGS
%                              'bfgs'   - Full bfgs
%                              ''       - Automatic
%                             Default: ''
%       - adaptiveH0        - When lbfgs is being used, this causes an adaptive
%                             estimate of the initial hessian to be used rather
%                             than the identity.  In practice this tends to
%                             improve performance on poorly conditioned problems
%                             but sometimes causes a degredation of performance
%                             when the conditioning is rather good.  Thus we
%                             default this to on with the suggestion that
%                             advanced users look to disable it when good
%                             preconditioners are known.
%                             Default: true
%       - precond           - Preconditioner.  If a vector is passed it is
%                             interpreted as a diagonal (but evaluated more
%                             efficiently).  A rectangular matrix may also be
%                             passed.  Typically a preconditioner seeks to
%                             define a surrogate objective function whos hessian
%                             is approximately the identity matrix.  The
%                             preconditioner is defined through y = A*x,
%                             where x is in the unpreconditioned space and y is
%                             in the preconditioned space.
%
%                             An offset in the unconditioned space may also be
%                             used by passing precond as a cell array {A,b}.
%                             This causes the preconditioner to be defined as
%                             y = A(x-b)
%
%                             The transformation from the unpreconditoned
%                             gradient to the preconditioned one involves the
%                             pseudo-inverse of A.  This is hard to compute from
%                             A because it is a O(n^3) calculation and so
%                             providing it is generally quite helpful.  See
%                             processPrecond.m for the code we use to make this
%                             transformation.
%                             Default: []
%       - precondInv        - Moore-Penrose pseudo-inverse of precond.  Like
%                             precond, this may be a cell array to provide b.
%                             Default: []
%       - blocks            - [k numBlocks] Variable block structure.  This
%                             cell array of indicies into x0 specifies the
%                             block problems block structure.  If a block
%                             structure is specified, LBFGS will assume the
%                             hessian of the objective function is zero
%                             everywhere except within the blocks.  If precondX0
%                             is true, then the blocks are being specified in
%                             the preconditioned space.
%
%                             If k is 2, then the blocks are allowed to overlap
%                             and the 2nd row provides indicies into the first
%                             that points to which indicies should also be used
%                             for direction updates.  As such, if the 2nd row
%                             is omitted it is interpreted as "1:numInd"
%                             Default: {}
%       - precondX0         - Flag indicating whether x0 is being provided in
%                             the preconditioned space.  By default it is
%                             assumed to be in the original space, and all
%                             results are returned in the original space.  If
%                             true, {x,f,g} will be left in the preconditioned
%                             space.
%                             Default: false
%       - lineSearchOpt     - Options for lineSearchMT.  These are defaulted to
%                             reasonable parameters for LBFGS.
%
%   --- Output ---
%       - verbose           - If passed, verbose output is provided.
%                             Default: false
%       - debug             - If true, output will contain per-iteration info
%                             Default: false
%       
% OUTPUT: 
%   x           - Approximate local minimizer of fun.
%   f           - Function value at x
%   g           - Gradient at x
%   exitFlag    - Flag giving the exiting conditons
%                   0 - Convergence occured
%                   1 - maxEval was reached
%                   2 - maxTime was reached
%                   3 - maxIter was reached
%                   4 - Round-off error prevents progress
%                   5 - User interrupt
%                   6 - Non-finite return values
%                   7 - Custom termination function
%   output      - A detailed structure describing what happened.  This is useful
%                 for debugging. 
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
% TODO:
%   Complete test coverage
%
% Copyright (C) 2013-2019 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function [xk,fk,gk,exitFlag,output] = LBFGS(fun,xk,varargin)

% Parse Inputs
[xk,fk,gk,config,opt,fun,output] = parseInputArguments(fun,xk,...
    nargout,varargin{:});
xkm1 = [];
fkm1 = [];
gkm1 = [];
gkNorm = norm(gk);
alpha = 1;
restartDirection = false;

while(1)
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Update Display %%%
    %%%%%%%%%%%%%%%%%%%%%%
    updateDisplay(fk,gkNorm,alpha,restartDirection,config,opt);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get search direction %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [d,restartDirection] = getSearchDirection(xk,fk,gk,xkm1,fkm1,gkm1,config,opt);
    if restartDirection
        config.searchHistory = resetSearchHistory(config.updateMethod);
        config.restartFlag = true;
    end
    
    
    %%%%%%%%%%%%%%%%%%%
    %%% Line Search %%%
    %%%%%%%%%%%%%%%%%%%
    %config.searchHistory.length
    lineOpt = opt.lineSearchOpt;
    if restartDirection
        if isinf(opt.maxEval)
            lineOpt.maxEval = 50;           % Try VERY hard
        else
            lineOpt.maxEval = opt.maxEval-config.numEval;
        end
    end
    lineOpt.maxTime = opt.maxTime - toc(config.T0);
    
    % Make sure explicit line search settings are respected
    lineOpt.maxEval = min(lineOpt.maxEval, opt.lineSearchOpt.maxEval);
    lineOpt.maxTime = min(lineOpt.maxTime, opt.lineSearchOpt.maxTime);
    
    if config.buildOutput && opt.debug
        [alpha,xkp1,fkp1,gkp1,Gkp1,exitFlag,numEval,lineOutput] = lineSearchMT(...
            fun,xk,d,fk,gk,lineOpt);
        
        output.iter(end+1) = struct(...
            'x',xkp1,...
            'f',fkp1,...
            'g',gkp1,...
            'd',d,...
            'alpha',alpha,...
            'lineIter',lineOutput.iter,...
            'restartDirection',restartDirection);
    else
        [alpha,xkp1,fkp1,gkp1,Gkp1,exitFlag,numEval] = lineSearchMT(...
            fun,xk,d,fk,gk,lineOpt);
    end
    if fkp1>=fk || exitFlag
        % We don't want to perform restarts unless absolutely necessary
        restartNeeded = true;
    else
        restartNeeded = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Update Iteration Info %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    config.numIter = config.numIter+1;
    config.numEval = config.numEval+numEval;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Check For Termination %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gkp1Norm = norm(gkp1);
    [termFlag,exitFlag] = checkTermination(xkp1,fkp1,gkp1,xk,fk,gk,gkp1Norm,...
        exitFlag,restartDirection,config,opt);
    if termFlag
        % Ensure monotonically non-increasing
        %warning('Need code');
        xk = xkp1;
        fk = fkp1;
        gk = gkp1;
        gkNorm = gkp1Norm;
        break
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Update the search history %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if restartNeeded
        config.searchHistory = resetSearchHistory(config.updateMethod);
        config.restartFlag = true;
    else
        config.searchHistory = updateSearchHistory(xkp1-xk,gkp1-gk,config,opt);
        config.restartFlag = false;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Accept the new step %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    xkm1 = xk;
    fkm1 = fk;
    gkm1 = gk;
    xk = xkp1;
    fk = fkp1;
    gk = gkp1;
    gkNorm = gkp1Norm;
    
end


%%%%%%%%%%%%%%%%%%%%%%
%%% Update Display %%%
%%%%%%%%%%%%%%%%%%%%%%
updateDisplay(fk,gkNorm,alpha,restartDirection,config,opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove Preconditioner %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xk,gk] = removePrecond(xk,gk,config);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Final Output Update %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if config.buildOutput
    output.numEval = config.numEval;
    output.numIter = config.numIter;
    output.cpuTime = toc(config.T0);
    output.config = config;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUB-FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{   
% SYNTAX:
%   [x0,f,g,config,opt,fun,output] = parseInputArguments(fun,x0,...
%       numArgOut,varargin)
%
% DESCRIPTION: 
%   This function parses the input arguments.
%
% INPUTS:
%   See LBFGS.m
%
% OUTPUT:
%   x0      - x0 in the preconditioned space
%   f       - Function value at x0
%   g       - Gradient at x0
%   config  - Configuration structure
%       .T0             - Time we started
%       .numIter        - Number of LBFGS iterations used
%       .numEval        - Number of function evaluations used
%       .f0             - Initial objective function value
%       .g0             - Initial gradient
%       .PA             - Preconditioner "A" matrix
%       .PAInv          - Preconditioner inverse "A" matrix
%       .Pb             - Preconditioner offset
%       .blocks         - Blocking structure
%       .initStep       - Initial step information
%       .updateMethod   - How directions are selected
%                           1 - lbfgs
%                           2 - bfgs
%       .restartFlag    - True if a restart is needed
%       .searchHistory  - Search history
%       .gradTol        - Gradient tolerance needed for convergence
%       .builtOutput    - True if output should be populated     
%       
%   opt     - Options structure
%   fun     - Objective function possibly preconditioned
%   output  - Output structure for 1st iterate or empty
%   
%}
function [x0,f,g,config,opt,fun,output] = parseInputArguments(fun,x0,...
    numArgOut,varargin)

%%%%%%%%%%%%%%%%%%%%%%
%%% Default Values %%%
%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;               % Do time sensative stuff first
opt = LBFGSOptions;
numArgs = numel(varargin);
maxSVDSize = 3000;      % Biggest size we'll try to SVD
maxBFGSDim = 2000;      % LBFGS is automatically selected for larger problems 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse the variable arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We expect funArgs as a cell array followed by flags or an options struct
i = 1;
while( i <= numArgs )
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
                case 'maxiter'
                    opt.maxIter = varargin{i};
                    i = i+1;
                case 'maxeval'
                    opt.maxEval = varargin{i};
                    i = i+1;
                case 'maxtime'
                    opt.maxTime = varargin{i};
                    i = i+1;                    
                case 'relgradtol'                    
                    opt.relGradTol = varargin{i};
                    i = i+1;
                case 'gradtol'                    
                    opt.gradTol = varargin{i};
                    i = i+1;
                case 'initstep'                    
                    opt.initStep = varargin{i};
                    i = i+1;
                case 'precond'
                    opt.precond = varargin{i};
                    i = i+1;
                case 'precondinv'
                    opt.precondInv = varargin{i};
                    i = i+1;
                case 'precondx0'
                    opt.precondX0 = true;
                case 'blocks'
                    opt.blocks = varargin{i};
                    i = i+1;
                case 'histlength'
                    opt.histLength = varargin{i};
                    i = i+1;
                case 'inithessinvfun'
                    opt.initHessInvFun = varargin{i};
                    i = i+1;
                case 'linesearchopt'
                    opt.lineSearchOpt = varargin{i};
                    i = i+1;
                case 'verbose'
                    opt.verbose = true;
                otherwise
                    error('Unexpected flag: %s',arg);
            end
            
        case 'struct'
            % Apply the options without checking them for correctness
            %opt = LBFGSOptions(arg);         % <-- This would check them
            opt = arg;
            
        otherwise
            error('Unexpected argument class');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Propagate interuptible %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.interruptible
    opt.lineSearchOpt.interruptible = true;
    
    % Purge any latent calls
    interruptible();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply Preconditioning %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fun,x0,initStep,PA,PAInv,Pb,blocks] = applyPrecond(fun,x0,opt,maxSVDSize);
numBlocks = size(blocks,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Function Eval %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f,g] = fun(x0);
%config.numEval = config.numEval + 1;           %Already did this


%%%%%%%%%%%%%%%%%%%%
%%% Setup Config %%%
%%%%%%%%%%%%%%%%%%%%
% Determine update method
switch lower(opt.method)
    case 'lbfgs'
        updateMethod(1:max(numBlocks,1)) = 1;
        
    case 'bfgs'
        updateMethod(1:max(numBlocks,1)) = 2;
        
    case ''
        if isempty(blocks)              % Global update setting
            if numel(x0)<maxBFGSDim
                updateMethod = 2;
            else
                updateMethod = 1;
            end
        else                            % Blockwise update setting
            updateMethod = cellfun(@numel,blocks(1,:));
            temp = updateMethod<maxBFGSDim;
            updateMethod(temp) = 2;
            updateMethod(~temp) = 1;
        end
        
    otherwise
        error('Unexpected update method');
end
config = struct(...     % Blocks will force the use of {}
    'T0',{T0},...
    'numIter',{0},...
    'numEval',{1},...
    'f0',{f},...
    'g0',{g},...
    'PA',{PA},...
    'PAInv',{PAInv},...
    'Pb',{Pb},...
    'blocks',{blocks},...
    'initStep',{initStep},...
    'updateMethod',{updateMethod},...
    'restartFlag',{false},...
    'searchHistory',{resetSearchHistory(updateMethod)},...
    'gradTol',{max([opt.gradTol, opt.relGradTol*norm(g)])},...
    'buildOutput',{numArgOut>4});


%%%%%%%%%%%%%%%%%%%%
%%% Setup Output %%%
%%%%%%%%%%%%%%%%%%%%
if config.buildOutput
    output = struct(...
        'opt',opt,...
        'numIter',[],...
        'numEval',[],...
        'cpuTime',[],...
        'x0',x0,...
        'f0',config.f0,...
        'g0',config.g0,...
        'iter',struct(...
            'x',{},...
            'f',{},...
            'g',{},...
            'd',{},...
            'alpha',{},...
            'lineIter',{},...
            'restartDirection',{}) );
else
    output = struct([]);
end

end

%{   
% SYNTAX:
%   [d, restartDirection] = getSearchDirection(x,f,g,xm1,fm1,gm1,config,opt)
%
% DESCRIPTION: 
%   This function computes a limited-memory search direction based on the
%   current point as well as the stored search history.
%
% INPUTS:
%   See LBFGS.m
%
% OUTPUT:
%   d           - Search direction.  This is an approximation to the negative of
%                 the inverse hessian times the gradient (a Newton step)
%   restartDirection - True if a restart direction is used
%   
%}
function [d,restartDirection] = getSearchDirection(x,f,g,xm1,fm1,gm1,config,opt)

% Default Values
restartDirection = false;

if config.numIter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Apply Hessian Updates %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config.restartFlag
        % Need a restart "initStep"
        d = restartInitStep(x,f,g,xm1,fm1,gm1,config);
        restartDirection = true;
        
    else
        d = zeros(numel(x),1);
        numBlocks = numel(config.searchHistory);
        
        for block = 1:numBlocks
            % Pointer to aid readability
            SH = config.searchHistory{block};
            
            switch config.updateMethod(block)
                case 1  % lbfgs
                    % Get the block components
                    if isempty(config.blocks)
                        gb = g;
                        gm1b = gm1;
                        dxb = x-xm1;
                    else
                        blockInds = config.blocks{1,block};
                        gb = g(blockInds);
                        gm1b = gm1(blockInds);
                        dxb = x(blockInds)-xm1(blockInds);
                    end
                    
                    % Indexing needed to step from oldest history to newest
                    historyIndex = [opt.histLength - (SH.length - SH.index)+1:...
                        opt.histLength, 1:SH.index];
                    alpha = zeros(SH.length,1);
                    
                    % Project the direction vector backwards
                    db = gb;
                    for i = historyIndex(end:-1:1)
                        alpha(i) = SH.rho(i)*(SH.s(:,i)'*db);
                        db = db - alpha(i)*SH.y(:,i);
                    end
                    
                    % Apply the initial inverse hessian approximation
                    if opt.adaptiveH0
                        dgb = gb-gm1b;
                        temp = norm(dgb)^2;
                        if temp>eps
                            % 'abs' because of possible neg. curv.
                            gamma = abs(dgb'*dxb)/temp;
                            db = gamma*db;
                        end
                    end
                    
                    % Project the direction vector forwards
                    for i = historyIndex
                        beta = SH.rho(i)*(SH.y(:,i)'*db);
                        db = db + SH.s(:,i)*(alpha(i)-beta);
                    end
                    
                    if isempty(config.blocks)
                        d = -db;
                    else
                        temp = config.blocks{2,block};
                        if isempty(temp)
                            d(blockInds) = -db;
                        else
                            d(blockInds(temp)) = -db(temp);
                        end
                    end
                    
                case 2 % bfgs
                    if isempty(config.blocks)
                        d = -(SH.H*g);
                    else
                        blockInds = config.blocks{1,block};
                        temp = config.blocks{2,block};
                        if isempty(temp)
                            d(blockInds) = -(SH.H*g(blockInds));
                        else
                            db = SH.H*g(blockInds);
                            d(blockInds(temp)) = -db(temp);
                        end
                    end
                    
                    
                otherwise
                    error('Unexpected updated method');
            end
            
        end
        
        % Ensure everything is well defined
        if ~all(isfinite(d))
            restartDirection = true;
            d = restartInitStep(x,f,g,xm1,fm1,gm1,config);
        end
        
    end
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get 1st search direction %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numElements = numel(config.initStep);
    switch numElements
        case 0
            d = g/(-2);    % Best step for isotropic quadratic objective function
        case 1
            gNorm = norm(g);
            if gNorm
                d = g*(-config.initStep/gNorm);
            else
                d = g;
            end
        case numel(x)
%             % Projection of expected deviation onto steepest descent step
%             % (trust the direction)
%             gHat = g/norm(g);
%             d = (-(abs(gHat)'*config.initStep))*gHat;

            % Minimal expected risk step (trust only the sign)
            d = -sign(g).*config.initStep;
            
        otherwise
            error('Not yet implemented');
    end
end

end


%{   
% SYNTAX:
%   d = restartInitStep(x,f,g,xm1,fm1,gm1,config)
%
% DESCRIPTION: 
%   This function computes a initial step after a restart
%
% INPUTS:
%   See LBFGS.m
%
% OUTPUT:
%   d           - Search direction
%   
%}
function d = restartInitStep(x,f,g,xm1,fm1,gm1,config)
initStep = config.initStep;
numElements = numel(initStep);
switch numElements
    case 0
        d = -g/2;    % Best step for isotropic quadratic objective function
    case 1
        gNorm = norm(g);
        if gNorm/2 > initStep
            d = g*(-initStep/norm(g));
        else
            d = g/(-2);
        end
    case numel(x)
        % Define some things for readability
        gHat = g/norm(g);
%         n = numel(x);
%         fFull = f(ones(n,1));
%         fm1Full = fm1(ones(n,1));
%         
%         xq = quadraticMinimizer(x,fFull,g,xm1,fm1Full,NaN(n,1));
%         [xc a] = cubicMinimizer(x,fFull,g,xm1,fm1Full,gm1);
%         xc(g.*a >= 0) = NaN;
%         
%         initStep = min([abs([xq-x xc-x]),initStep],[],2);
%         d = (-(abs(gHat)'*initStep))*gHat;
        
        %d = (-(abs(gHat)'*config.initStep))*gHat;
        
        % Minimal expected risk step (trust only the sign)
        d = -sign(g).*config.initStep;
end

end


%{   
% SYNTAX:
%   [result,exitFlag] = checkTermination(fk,fkp1,Gkp1,lineExitFlag,...
%                                        restartDirection,config,opt)
%
% DESCRIPTION: 
%   This function checks for termination
%
% OUTPUT:
%   result      - True if we should terminiate and false otherwise
%   exitFlag    - Exit value (see "exitFlag" in LBFGS.m)
%   
%} 
function [result,exitFlag] = checkTermination(x,f,g,xm1,fm1,gm1,G,lineExitFlag,...
    restartDirection,config,opt)


%%%%%%%%%%%%%%%%%%%%%%
%%% Default Values %%%
%%%%%%%%%%%%%%%%%%%%%%
result = true;          % Be optimistic
exitFlag = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for non-finit return values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfinite(f) || ~isfinite(G)
    if f < 0
        exitFlag = 0;
    else
        exitFlag = 6;
    end
    return
end


%%%%%%%%%%%%%%%%
%%% Gradient %%%
%%%%%%%%%%%%%%%%
if abs(G) < config.gradTol
    exitFlag = 0;
    return
end


%%%%%%%%%%%%%%%
%%% Timeout %%%
%%%%%%%%%%%%%%%
if config.numEval >= opt.maxEval
    exitFlag = 1;
    return
end

if config.numIter >= opt.maxIter
    exitFlag = 3;
    return
end

if toc(config.T0) >= opt.maxTime
    exitFlag = 2;
    return
end


%%%%%%%%%%%%%%%%%%%%%%
%%% User Interrupt %%%
%%%%%%%%%%%%%%%%%%%%%%
% This needs to happen before the stagnation check
if lineExitFlag==8 || (opt.interruptible && interruptible)
    exitFlag = 5;
    return
end


%%%%%%%%%%%%%%%%%%
%%% Stagnation %%%
%%%%%%%%%%%%%%%%%%
if restartDirection && f>=fm1 
    exitFlag = 4;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Custom Termination Fun. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(opt.termFun,'function_handle') && ...
    opt.termFun(x,f,g,xm1,fm1,gm1,lineExitFlag,restartDirection,config,opt,...
        @removePrecond)
    exitFlag = 7;
    return
end

result = false;

end


%{   
% SYNTAX:
%   searchHistory = updateSearchHistory(dx,dg,config,opt)
%
% DESCRIPTION: 
%   This function checks for termination
%
% INPUTS:
%   dx              - x_{k+1} - x_k
%   dg              - g_{k+1} - g_k
%   config          - Config struct
%   opt             - Options struct
%
% OUTPUT:
%   searchHistory   - Updated search history
%   
%} 
function SH = updateSearchHistory(dx,dg,config,opt)
% For readability
SH = config.searchHistory;
numBlocks = numel(SH);

for block = 1:numBlocks
    % Get the block indicies
    if isempty(config.blocks)
        dxb = dx;
        dgb = dg;
    else
        blockInds = config.blocks{1,block};
        dxb = dx(blockInds);
        dgb = dg(blockInds);
    end
    
    %Check for convexity (We assume continuous differentiability)
    rho = 1/(dgb'*dxb);
    
    if rho > 0
        switch config.updateMethod(block)
            case 1  % lbfgs
                %Get the new index
                index = mod(SH{block}.index,opt.histLength)+1;
                
                %Add the new information
                SH{block}.index = index;
                SH{block}.length = min(SH{block}.length+1,opt.histLength);
                SH{block}.s(:,index) = dxb;
                SH{block}.y(:,index) = dgb;
                SH{block}.rho(index) = rho;
                
            case 2 % bfgs
                temp = (rho*dxb)';
                V = (-dgb)*temp;
                V(1:numel(dxb)+1:end) = diag(V)+1;   %diag(V)+=1
                if SH{block}.length==0 && opt.adaptiveH0
                    SH{block}.H = SH{block}.H/(norm(dgb)^2*rho);
                end
                SH{block}.H = V'*SH{block}.H*V + dxb*temp;
                
                SH{block}.length = SH{block}.length + 1;
                
            otherwise
                error('Unexpected update method');
        end
    else
        % Assume the previous curvature (we're making progress)
    end
    
end



end


%{   
% SYNTAX:
%   result = resetSearchHistory(method)
%
% DESCRIPTION: 
%   Resets the search history to its initial (empty) state
%
% INPUT:
%   method      - Scalar describing the history update method or an array of
%                 methods corresponding to the block structure.
%                 See config.updateMethod
%
% OUTPUT:
%   result      - Cell aray of empty search history structures
%      LBFGS structure
%       .length         - Current search history length
%       .index          - Index to the latest history information.  All buffers
%                         are modulo (FILO) buffers of vertical vectors
%       .s              - x_k - x_{k-1} (modulo buffer of vectors)
%       .y              - g_k - g_{k-1} (modulo buffer of vectors)
%       .rho            - 1/s_k'*y_k  (modulo buffer of scalars)
%
%      BFGS structure
%       .length         - Current search history length
%       .H              - Inverse Hessian approximation
%} 
function result = resetSearchHistory(method)

for i = 1:numel(method)
    switch method(i)
        case 1
            result{i} = struct(...
                'length',0,...
                'index',0,...
                's',[],...
                'y',[],...
                'rho',[]);
        case 2
            result{i} = struct(...
                'length',0,...
                'H',1);
        otherwise
            error('Unexpected search update method');
    end
end

end


%{   
% SYNTAX:
%   updateDisplay(fk,gkNorm,alpha,restartDirection,config,opt)
%
% DESCRIPTION: 
%   Updates the display with the status of the solver
%
% OUTPUT:
%   May print output to the screen
%}
function updateDisplay(fk,gkNorm,alpha,restartDirection,config,opt)

if ~restartDirection || config.numIter==0
    msg = '';
elseif config.restartFlag
    msg = 'restart';
else
    msg = 'neg curv';
end
if opt.verbose
    fprintf(' Iter: %-3g |  f = %-11.5g |  g = %-11.5g | step = %-8.1g %s\n',...
        config.numIter,fk,gkNorm,alpha,msg);
end
end



%{   
% SYNTAX:
%   [f,g] = precondFun(x,fun,PInv,b)
%
% DESCRIPTION: 
%   Evaluates fun at x after first unpreconditioning it (moving it out of the
%   preconditioned space.  The function values and gradients are then
%   transformed into the preconditioned space.
%}
function [f,g] = precondFun(x,fun,PInv,b)
    if isvector(PInv) && ~isscalar(x)
        if isempty(b)
            [f,g] = fun(PInv.*x);
        else
            [f,g] = fun(PInv.*x+b);
        end
        g = PInv.*g;
    else
        if isempty(b)
            [f,g] = fun(PInv*x);
        else
            [f,g] = fun(PInv*x+b);
        end
        g = (g.'*PInv).';
    end
end


%{   
% SYNTAX:
%   [f,g] = precondFun2(x,fun,P,b)
%
% DESCRIPTION: 
%   Evaluates fun at x after first preconditioning it (moving it into the
%   unpreconditioned space.  The function values and gradients are then
%   transformed into the preconditioned space.  This function is like precondFun
%   except the linear system involving P is solved each time.
%}
function [f,g] = precondFun2(x,fun,P,b)
if isempty(b)
    [f,g] = fun(P\x);
else
    [f,g] = fun(P\x+b);
end
g = (g.'/P).';
end


%{   
% SYNTAX:
%   [fun,x0,initStep,PA,PAInv,Pb,blocks] = applyPrecond(fun,x0,opt,maxSVDSize)
%
% DESCRIPTION: 
%   This function uses the options structure to apply any preconditioners
%}
function [fun,x0,initStep,PA,PAInv,Pb,blocks] = applyPrecond(fun,x0,opt,maxSVDSize)

% Default Values
initStep = opt.initStep;
precondX0 = opt.precondX0;
blocks = opt.blocks;

% Decompose the preconditioner and its inverse
if iscell(opt.precond)
    PA = opt.precond{1};
    Pb = opt.precond{2};
else
    PA = opt.precond;
    Pb = [];
end
if iscell(opt.precondInv)
    PAInv = opt.precondInv{1};
    if isempty(Pb)
        Pb = opt.precondInv{2};
    elseif Pb ~= opt.precondInv{2} % This should get checked in LBFGSOptions.m
        error('The preconditioner and its inverse are inconsistent');
    end
else
    PAInv = opt.precondInv;
end

% If an explicit preconditioner was provided
if ~isempty(PA) || ~isempty(PAInv)
    
    % Process the preconditioner
    [PA,PAInv,Pb,x0] = processPrecond(PA,PAInv,Pb,maxSVDSize,x0,precondX0);
    if isempty(PAInv)
        fun = @(x)precondFun2(x,fun,PA,Pb);     % hard way
    else
        fun = @(x)precondFun(x,fun,PAInv,Pb);   % easy way
    end
    
    % Check for preconditioner compatibility with the blocking
    blocks = checkPrecondBlockingCompatibility(PA,PAInv,opt.blocks,precondX0);
    
    % Adjust the initStep for the preconditioner
    if ~precondX0
        if numel(initStep)>1
            % We need to move the initstep vector into the preconditioned space
            if size(PA,2)==0
                initStep = PAInv\initStep;
            elseif size(PA,2)==1
                initStep = PA.*initStep;
            else
                initStep = PA*initStep;
            end
        elseif numel(initStep)==1
            % ???
        end
    end

elseif ~isempty(blocks)
    
    % Gracefully fix blocks that have been transposed
    % (why is this not in options?)
    [k,numBlocks] = size(blocks);
    if numBlocks==1 && k>2
        blocks = blocks.';
        temp = numBlocks;
        numBlocks = k;
        k = temp;
    elseif k==1
        blocks{2,end} = [];
    end
    
    % Look for block collisions or gaps
    x = zeros(numel(x0),1);
    for blockInd = 1:numBlocks
        if isempty(blocks{2,blockInd})
            x(blocks{1,blockInd}) = x(blocks{1,blockInd})+1;
        else
            temp = blocks{1,blockInd}(blocks{2,blockInd});
            x(temp) = x(temp)+1;
        end
    end
    
    if any(x>1)
        error('Block collisions.  Blocks may not overlap.');
    end
    
    x = find(x==0); % Indentify unallocated dimensions
    if ~isempty(x)
        warning('Automatically completing block coverage');
        blocks{1,end+1} = x;
    end
    
end

end


%{   
% SYNTAX:
%   [x,g] = removePrecond(x,g,config)
%
% DESCRIPTION: 
%   This function removes the preconditioner if it has been applied.
%}
function [x,g] = removePrecond(x,g,config)

% For readability
PInv = config.PAInv;
P = config.PA;
b = config.Pb;
PNotEmpty = ~isempty(P);

% Short circuit on empty (user preconditioned out everything)
if isempty(x)
    x = b;
    g = zeros(size(b));
    return
end


if isempty(PInv)
    if PNotEmpty
        % P will always be a full matrix (otherwise PInv would have been formed)
        x = P\x;
        g = (g.'*P).';
    end
else
    PInvVector = isvector(PInv) && ~isscalar(x);
    if PInvVector
        x = PInv.*x;
    else
        x = PInv*x;
    end
    
    if PNotEmpty
        if size(P,2)==1
            g = P.*g;
        else
            g = (g.'*P).';
        end
    else
        % PInv will always be a full matrix (otherwise P would have been formed)
        g = (g'/PInv).';
    end
end
if ~isempty(b)
    x = x+b;
end

end


%{   
% SYNTAX:
%   checkPrecondBlockingCompatibility(PA,PAInv,blocks)
%
% DESCRIPTION: 
%   This function checks the compatibility of the preconditioner with the
%   blocking, and if possible, updates the blocking structure to be compatible
%   with the preconditioner.
%
%   The code is principally interested in evaluating the user-provided objective
%   function in the unpreconditioned space.  Thus, we really want to have PInv
%   so we can use some form of "fun(PInv.*x)".  If PInv is not provided, and
%   could not be efficiently found, we'll do things the hard way "fun(P\x)"
%}
function blocks = checkPrecondBlockingCompatibility(PA,PAInv,blocks,precondX0)

% If there are no blocks or their already in the right space we're done
if isempty(blocks) || precondX0
    return
end

% Default Values
[k,numBlocks] = size(blocks);
if numBlocks==1 && k>2
    blocks = blocks.';
    temp = numBlocks;
    numBlocks = k;
    k = temp;
elseif k==1
    blocks{2,end} = [];
end
PAEmpty = isempty(PA);
PAInvEmpty = isempty(PAInv);

% Determine the dimensions of the problem and try to exit early
if PAInvEmpty           % We'll be doing things the hard way w/ PA
    if PAEmpty || isvector(PA)
        % Blocks all stay the same (only diagonal scaling is occuring)
        return
    end
    
    [M,N] = size(PA);       %  Number of dimensions in the two spaces (M<=N)
    
else                    % Using PInv
    if isvector(PAInv)
        % Blocking not effected
        return
    end
    
    [N,M] = size(PAInv);    %  Number of dimensions in the two spaces (M<=N)
end


% Project energy back through the blocks to see where it lands
if PAEmpty
    projectFun = @(x)find(PAInv\x~=0);
else
    projectFun = @(x)find(PA*x~=0);
end
   
for blockInd = 1:numBlocks
    % Map the direction indicies to their direct preconditioned indicies
    if ~isempty(blocks{2,blockInd})
        temp = blocks{2,blockInd};
        blockSize = numel(temp);
        x = zeros(N,1);
        x(blocks{1,blockInd}(temp)) = randn(blockSize,1);
        blocks{2,blockInd} = projectFun(x);
    end
    
    % Map the hessian indicies to their preconditioned values
    blockSize = numel(blocks{1,blockInd});
    x = zeros(N,1);
    x(blocks{1,blockInd}) = randn(blockSize,1);
    blocks{1,blockInd} = projectFun(x);
    
    % Update the subindexing
    if ~isempty(blocks{2,blockInd})
        [~,blocks{2,blockInd}] = intersect(...
            blocks{1,blockInd},...
            blocks{2,blockInd});
    end
    
end


% Look for block collisions or gaps
x = zeros(M,1);
for blockInd = 1:numBlocks
    if isempty(blocks{2,blockInd})
        x(blocks{1,blockInd}) = x(blocks{1,blockInd})+1;
    else
        temp = blocks{1,blockInd}(blocks{2,blockInd});
        x(temp) = x(temp)+1;
    end
end

if any(x>1)
    error('Preconditioner is inconsistent with blocking strategy');
end

x = find(x==0); % Indentify unallocated dimensions
if ~isempty(x)
    warning('Assigning block for preconditioner near-null-space');
    blocks{1,end+1} = x;
end


end
