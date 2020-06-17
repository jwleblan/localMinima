% localMinimaDetectionPoster - Even simpler example of local minima detection
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   localMinimaDetectionPoster()
%
% PURPOSE:
%   This function demonstrates the key ideas behind the improved test for
%   convergence to a local (rather than global) optima on a simple example.
%
%   Specifically, this code is intended to demonstrate the key ideas in a manner
%   that is very easy to follow.  The example based on detecting the frequency
%   of a sinusoid in noise.  The key features of this example are 
%   1) The problem is very easy to understand
%   2) An ML estimator is used, and the search problem is chosen to ensure
%      multiple local minima exist
%   3) The problem is low-dimensional, simplifiying visualization
%  
% INPUT:
%   NONE
% 
% OUTPUT:
%   Plots are to illustrate some key results, but the code is meant to be read
%   and understood.
%
% NOTES:
%   This code makes use of lambda functions (functions that return functions),
%   which leads me to perform some scope encapsulation.  If you don't know why
%   older versions of Matlab benefit from scope encapsulation just ingore it.
%   Hopefully the use of lambda functions will make it clearer where
%   problem-dependent things are happening.
%
%   "numRelaxationDims" determines what test is compared with Biernacki's work.
%   When 0, a one-sided Biernacki test is used.  For location families this
%   adjustment leads to a dominating test.  For high integers, this is the
%   number of non-parametrically determined relaxation dims to be used in
%   conjunction with a one-sided Biernacki test.
%
%   Whether one searches the entire relaxed space, or only the relaxation
%   dimensions appended to the canonical repsresentaion, is controlled by the
%   variable 'searchThetaPrimeOnly' in the subfunction leblancTest.
%-------------------------------------------------------------------------------
%}
function localMinimaDetectionPoster()

% Default Values
N = 100;
thetaTrue = 3*pi;
sigmaN = 1;
numTrials = 5000;			% 5000 was used for publication (~2hrs)
numRelaxationDims = 1;      % Must be 1 to get the 3D plot (2 and 3 are better)
rng(0);                     % For total reproducability

% Setup the forward model
xx = linspace(0,1,N).';
dFun = getForwardModel(xx,sigmaN);

% Identify the local minima for this problem
dFunNF = @(theta)dFun(theta,false);     % Noise-free data function
d = dFunNF(thetaTrue);                  % Get noise free data
Lambda = getLambda(dFunNF, d, sigmaN);  % Get the negative log-likelihood
theta0 = pi/10;                         % Initial estimate
LBFGSOpts = LBFGSOptions();
LBFGSOpts.initStep = .1;
[thetaLocalMin,f,g,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);

% Plot the noise-free data and local minima
hFig = figure(1);
plot(...
    xx,d,'-',...
    xx,dFun(thetaLocalMin,false),'--');
xlabel('x');
ylabel('E[d]');
hLegend = legend({...
    '$\mu: \theta = \theta_0$',...
    '$\mu: \theta = \hat\theta$' });
set(hLegend,'Interpreter','latex');
title('Expected signal at two stationary points');
prepareFigure();

% Plot signal realizations at two stationary points
% This plot will show why the problem is hard
hFig = figure(2);
plot(...
    xx,dFun(thetaTrue),'o',...
    xx,dFun(thetaLocalMin),'x');
ylim([-4 4]);
xlabel('x');
ylabel('d');
hLegend = legend({...
    '$d: \theta = \theta_0$',...
    '$d: \theta = \hat\theta$' });
set(hLegend,'Interpreter','latex');
title('Signal realizations at two stationary points');
prepareFigure();


% Plot Figures 1 and 2 on the same plot
figure(3);
plot(...
    xx,dFun(thetaTrue),'bo',...
    xx,dFun(thetaLocalMin),'rx',...
    xx,d,'b-',...
    xx,dFun(thetaLocalMin,false),'r--');
ylim([-4 4]);
xlabel('x');
hLegend = legend({...
    '$d: \theta = \theta_0$',...
    '$d: \theta = \hat\theta$',...
    '$\mu: \theta = \theta_0$',...
    '$\mu: \theta = \hat\theta$'});
set(hLegend,'Interpreter','latex');
title('Signal at two stationary points');
prepareFigure();


% Plot the objective function
figure(4);
tt = linspace(0,4*pi,500);
ll = arrayfun(Lambda,tt);
plot(tt,ll);
xlim([0 4*pi]);
temp = ylim;
hold on
h = plot(...
    [thetaLocalMin thetaLocalMin],temp,'--r',...
    [thetaTrue thetaTrue],temp,'--r');
hold off;
axis tight
title('Negative Log-Likelihood');
prepareFigure();


% Find a helpful nonparametric embedding
nominalParameters = linspace(0,4*pi,50);
startingPoints = linspace(0,4*pi,10);
[B,s] = findRelaxationBasis(@(theta)dFun(theta,false),...% Need noise-free model
    nominalParameters,...                               % Nominal true solutions
    startingPoints,...                                  % Nominal starting points
    numRelaxationDims,sigmaN);                                 

% Plot the relaxation basis
if numRelaxationDims
    figure(3);
    plot(xx,B);
    xlabel('x');
    ylabel('r');
    title('Embedding basis determined by Algorithm 1');
    grid on
    prepareFigure();
end

% Plot the low-dimensional relaxation objective space
if numRelaxationDims==1
    % Construct an objective function with the non-parametric embedding
    dFunB = getForwardModel(xx,sigmaN,B(:,1));
    
    % Make some plots that help clarify how the test will work
    plotRelaxationObjectiveSpace(xx,dFunB,thetaLocalMin,thetaTrue,sigmaN);
end

% Note: Right here is where you can change the model from naive to the
% algorithmically selected relaxation.  This code was used for lots of plots and
% diagnostics so it has a lot of switches.  By changing the forward model from
% dFun to dFunB you are using the "basis" version.  You can use any basis you
% want, but in the paper I show one way to find "good" ones.  
%
% Thanks for checking out the code, and if you have ideas on how to either get
% better relaxation bases to prove the ones I have are are in some sense the 
% best, let me know.  I've started a proof along these lines, but chose to
% abandon it in favor of a simpler paper because it was clear from our reviewers
% that we were struggling to communicate even the basic ideas -- writing is hard.  

% dFun = dFunB;


% Setup the tests.  The argument to leblancTest is the number
testArray = {@biernackiTest, leblancTest(numRelaxationDims)};
if numRelaxationDims>0
    testNameArray = {'Biernacki',sprintf('leblanc(%g)',numRelaxationDims)};
else
    testNameArray = {'Biernacki [4]','Proposed one-sided'};
end


% Note: This is how I configured the script to merge what was once 3 different
% plots.  Originally, the paper threaded the simple example through all of the
% theoretical discussion.  Ultimately, we decided to condense all of the example
% discussion into a single section/plot.  Below is how I configured the script
% to do this.  I ran the final curve separately and added it 
% (see above w/ dFun = dFunB)

if false
    testArray = {@biernackiTest, leblancTest(0), leblancTest(1),...
        leblancTest(3)};
    testNameArray = {'Biernacki [4]','Proposed one-sided','Proposed k=1',...
        'Proposed k=3'};
end


disp('Running Monte Simulations...');
for testInd = 1:numel(testArray)
    currentTest = testArray{testInd};
    testName = testNameArray{testInd};
    
    % Generate data under H0 (Validate the null distributin)
    theta0 = thetaTrue;
    thetaHatArray = zeros(1,numTrials);
    for i = numTrials:-1:1
        % Generate some real data and find a the global minima
        d = dFun(thetaTrue);
        Lambda = getLambda(@(x)dFun(x,false),d,sigmaN);
        [thetaHat,f,g,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);
        
        % Test for convergence to a local minima
        thetaHatArray(i) = thetaHat;
        [phi0(i),s0(i)] = currentTest(d,thetaHat,dFun,@getLambda);
    end
    
    % This next line just gives more detailed diagnostics
    %showMonteCarloPlots(thetaHatArray,s0,thetaTrue,thetaLocalMin,sprintf('H0 (%s)',testName));
    
    
    % Generate data under H1
    theta0 = thetaLocalMin;
    for i = numTrials:-1:1
        % Generate some real data and find a the local minima
        d = dFun(thetaTrue);
        Lambda = getLambda(dFun,d,sigmaN);
        [thetaHat,f,g,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);
        
        % Test for convergence to a local minima
        thetaHatArray(i) = thetaHat;
        [phi1(i),s1(i)] = currentTest(d,thetaHat,dFun,@getLambda);
    end
    
    % This next line just gives more detailed diagnostics
    %showMonteCarloPlots(thetaHatArray,s1,thetaTrue,thetaLocalMin,sprintf('H1 (%s)',testName));
    
    % Generate a ROC curve
    [PD(:,testInd),PFA(:,testInd)] = sampleROC(...
        [s0,s1],...
        [false(1,numTrials),true(1,numTrials)]);
    figure();
    plot(PFA(:,testInd),PD(:,testInd));
    xlabel('PFA');
    ylabel('PD');
    title(sprintf('Global Optimum Detection Performance (%s)',...
        testName));
    prepareFigure();
    
end

% Show both ROC Curves
figure();
plot(PFA,PD);
xlabel('PFA');
ylabel('PD');
title('Global Maximum Detection Performance');
grid on
legend(testNameArray);
prepareFigure;

figure();
semilogx(PFA,PD);
xlim([10/numTrials 1]);
xlabel('log_{10} PFA');
ylabel('PD');
title('Global Maximum Detection Performance');
grid on
legend(testNameArray);
prepareFigure;

end



%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   dFun = getForwardModel(xx,sigmaN)
%   dFun = getForwardModel(xx,sigmaN,B)
%   [d,dd] = dFun(theta,<addNoise>)
%
% PURPOSE:
%   This function returns a function that implements the forward model.  The
%   resulting function easily be evaluated with and without noise, and
%   optionally evaluates the gradient w.r.t. the parameters.
%
% INPUT:
%   xx          - [N 1] Ordinates where the data is measured
%   sigmaN      - Standard deviation of the normally distributed read noise
%   B           - [N P-1] Basis function for the relaxation
%   theta       - [P 1] Model parameterization
%   addNoise    - True if noise is added and false otherwise.  Each time noise
%                 is requested the state of the random number generator is
%                 altered.  See rng.m for more information.
%                   Default: true
% 
% OUTPUT:
%   dFun        - Forward model function with the prototype
%                   [d,dd] = dHat(theta,<addNoise>)
%                 where d is the data and dd its gradient w.r.t. theta
%   d           - [N 1] Data vector
%   dd          - [N P] Gradient of d w.r.t. theta
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function dFun = getForwardModel(xx,sigmaN,B)

% This function basically provides scope encapsulation
switch nargin
    case 2
        dFun = @(theta,varargin) evaluateForwardModel(xx,sigmaN,theta,...
            varargin{:});
        
    case 3
        dFun = @(theta,varargin) evaluateForwardModelBasis(xx,sigmaN,B,theta,...
            varargin{:});
        
    otherwise
        error('Unexpected number of arguments');
end

end


% This is an implementation of the forward model
% It uses a variable length parameterization of theta of the form
%   sin(theta(1)*x + theta(2)*x^2 + ...) 
function [d,dd] = evaluateForwardModel(xx,sigmaN,theta,addNoise)

% Parse Inputs
if nargin<4
    addNoise = true;
end

% Get the argument to sine
thetaLen = numel(theta);
dArg = [xx,bsxfun(@power,xx,2:thetaLen)];
arg = dArg*theta;

% Generage d
d = sin(arg);
if addNoise
    d = d + sigmaN*randn(size(xx));
end

% Generate dd
if nargout>1
    dd = bsxfun(@times,dArg,cos(arg));
end

end


% This is the implementation of the forward model with relaxation
% It uses a variable length parameterization of theta of the form
%   sin(theta(1)*x) + B*theta(2:end)
function [d,dd] = evaluateForwardModelBasis(xx,sigmaN,B,theta,addNoise)

% Parse Inputs
if nargin<5
    addNoise = true;
end

% Generage d
thetaLen = numel(theta);
d = sin(theta(1)*xx);
if thetaLen > 1
    d = d + B*theta(2:thetaLen);
end

if addNoise
    d = d + sigmaN*randn(size(xx));    
end

% Generate dd
if nargout>1
    dd = xx.*cos(theta(1)*xx);
    
    if thetaLen > 1
        dd = [dd B];
    end
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   Lambda = getLambda(dHatFun, d, sigmaN)
%   [f,g] = Lambda(theta)
%
% PURPOSE:
%   This function provides Lambda; negative log-liklihood function in terms of
%   the parameters theta.  The problem parameterization is abstracted away, and
%   instead a general forward model is used.  The correct way to think of
%   Lambda is as follows
%       [f,g] = lambda(dHat(theta) ; d,sigmaN)      % A function of theta
%
%   Most tests for convergence to a local minima exploit assymptotic normality,
%   so the affine dependency of Lambda on sigmaN will not matter.  I'm including
%   it here to not confuse those less familiar with how these tests play out.
%
% INPUT:
%   dFun        - A function that provides a noise-free data estimate with the 
%                 prototype
%                   [dHat,ddHat] = dHatFun(theta)
%                 where dHat is the predicted data and ddHat is its gradient
%                 w.r.t. theta.
%   d           - [N 1] Measured data vector
%   theta       - [P 1] The data parameterization
% 
% OUTPUT:
%   Lambda      - A function with the prototype
%                   [f,g] = Lambda(theta)
%                 with f the function value and g its gradient w.r.t. theta
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%-------------------------------------------------------------------------------
%}
function Lambda = getLambda(dHatFun,d,sigmaN)

% This function is providing scope encapsulation
switch nargin(dHatFun)
    case 1
        Lambda = @(theta)evaluateLambda(dHatFun,d,sigmaN,theta);
    
    case {-2,2}
        % We got a non-conforming function that probably has <addNoise>)
        Lambda = @(theta)evaluateLambda(@(x)dHatFun(x,false),d,sigmaN,theta);
        
    otherwise
        error('Non-conforming forward model: dHatFun');
end

end

% Implementation of Lambda
function [f,g] = evaluateLambda(dFun,d,sigmaN,theta)

% Default Values
[dHat,ddHat] = dFun(theta);
n = numel(d);

% Compute f
res = dHat-d;
f = norm(res).^2;       % norm(res).^2 would be sufficient for most problems
f = f/(2*sigmaN^2) + n*log(sigmaN) + n/2*log(2*pi);

g = ddHat.'*(2*res);    % gradient when we ignore the affine pieces
g = g/(2*sigmaN^2);

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   B = findRelaxationBasis(dFun,nominalParameters,startingPoints,...
%                           numRelaxationDims,sigmaN)
%
% PURPOSE:
%   This function identifies a non-parametric basis which maximizes the power of
%   of a problem-specific Rao test designed to identify convergence to a local
%   minima.
%
% INPUT:
%   dFun                - A function that provides a noise-free data estimate
%                         with the prototype
%                           [dHat,ddHat] = dHatFun(theta)
%                         where dHat is the predicted data and ddHat is its
%                         gradient w.r.t. theta.
%   nominalParameters   - [P K] Array of nominal parameters
%   startingPoints      - [P L] Array of starting points
%   numRelaxationDims   - Number of relaxation dimensions for the embedding
%   sigmaN              - Standard deviation of the noise
% 
% OUTPUT:
%   B       - [N numRelaxationDims] Relaxation basis
%   s       - [numRelaxationDims 1] Singular values associated with B
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%-------------------------------------------------------------------------------
%}
function [B,s] = findRelaxationBasis(dFun,nominalParameters,startingPoints,...
    numRelaxationDims,sigmaN)

% Default Values
B = [];
[P,numNominalParams] = size(nominalParameters);
N = numel(dFun(startingPoints(:,1)));
numStartingPoints = size(startingPoints,2);

LBFGSOpts = LBFGSOptions();
LBFGSOpts.initStep = .1;

% Threshold used to determine convergence to the correct solution
thresh = 1E-3;          % This is problem dependent, but usually easy to get

% Record the lagrange multipliers for local minima
lagrange = zeros(N,numNominalParams*numStartingPoints);
J = 0;
for i = 1:(numNominalParams*numStartingPoints)
    
    % Get the parmater and starting point
    [iParam,iStart] = ind2subVec([numNominalParams numStartingPoints], i);
    thetaTrue = nominalParameters(:,iParam);
    theta0 = startingPoints(:,iStart);
    
    % Solve the noise-free problem
    d = dFun(thetaTrue);                % Get noise-free data
    Lambda = getLambda(dFun,d,sigmaN);  % Get objective function
    [thetaHat,f,g,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);
    
    % Record the score
    if norm(thetaHat-thetaTrue)>thresh
        J = J+1;
        lagrange(:,J) = (dFun(thetaHat)-d)/sigmaN;  
    end
end

% Free the excess memory
lagrange(:,J+1:(numNominalParams*numStartingPoints)) = [];

% Identify the basis
[B,s] = svds(lagrange,numRelaxationDims,'L');

% Provide the singular values if requested
if nargout>1
    s = diag(s);
end

end

%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   plotRelaxationObjectiveSpace(xx,dFun,thetaLocal,thetaTrue,sigmaN)
%
% PURPOSE:
%   This function plots the expected log-likelihood over the relaxed space.
%  
% INPUT:
%   xx          - Ordinates of the data
%   dFun        - Data generation function for the embedding with the prototype
%                   [dHat,ddHat] = dHatFun(theta,<addNoise>)
%   thetaLocal  - Local optima in the canonical space
%   thetaTrue   - Global optima in the canonical space
%   sigmaN      - Noise std. dev.
% 
% OUTPUT:
%   3D plot of the local minima structure under relaxation
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%-------------------------------------------------------------------------------
%}
function plotRelaxationObjectiveSpace(xx,dFun,thetaLocal,thetaTrue,sigmaN)

% Default values
p = numel(xx);                      % We are in the embedded space
muTrue = dFun(thetaTrue,false);     % True expectation

% Non-centrality parameter of the distribution of the log-likelihood
lambda = @(theta)sum( (muTrue-dFun(theta,false)).^2/sigmaN.^2 );

% Negative expectation of the log-likelihood
EmEll = @(theta) (p+lambda(theta))/2 + p*log(sigmaN) + p*log(2*pi)/2;

% Negative expectation of the log-likelihood under thetaHat
% This is a costant for location-family members and is put here to remind the
% reader of this fact.
EmEllHat = p/2 + p*log(sigmaN) + p*log(2*pi)/2;

% Determine how much of the relaxed dimension we need to insert to reach a minima
x0 = 1;
[xMin,fVal] = fminsearch(@(x)EmEll([thetaLocal,x]),x0);
relaxLim = [0 2*abs(xMin)];       % Should start at 0 and go positive for plotting
relaxSign = sign(xMin);

% Determine limits for the natural space
delta = abs(thetaTrue - thetaLocal);
if thetaLocal<thetaTrue
    canonLim = [thetaLocal-.03*delta thetaTrue+.2*delta];
else
    canonLim = [thetaTrue-.03*delta thetaLocal+.2*delta];
end

% Evaluate the surface
numGridPoints = 100;
theta1Array = linspace(canonLim(1),canonLim(2),numGridPoints);
theta2Array = relaxSign * linspace(relaxLim(1),relaxLim(2),numGridPoints);
EmEllArray = zeros(numGridPoints);
for i = 1:numGridPoints
    for j = 1:numGridPoints
        EmEllArray(i,j) = EmEll([theta1Array(i),theta2Array(j)]);
    end
end

% Here we determine the expectation of the min as opposed to the min of the
% expection as given below.  For this problem the two are essentially the same,
% so I will use the faster of the two calculations...  Feel free to run this and
% check for yourself. ~Joel
% Determine the expected minima over the relaxed dimension.
if false
    x0 = 0;                                 % Start with no relaxation
    numTrials = 1000;
    for i = 1:numGridPoints
        theta1 = theta1Array(i);
        for trialInd = numTrials:-1:1
            d = dFun(thetaTrue,true);                       % Noisy realization
            Lambda = getLambda(@(x)dFun(x,false),d,sigmaN); % Log likelihood
            [xMin,fVal(trialInd)] = fminsearch(@(x)Lambda([theta1,x]),x0);
        end
        EMinHat(i) = mean(fVal);
    end
end

% Here I am measuring the the very small offset lost due to fitting the noise in
% the neighborhood of the true solution
if false
    x0 = [thetaTrue;0];
    opts = optimset('fminsearch');
    opts.MaxFunEvals = 2000;
    opts.MaxIter = 1000;
    numTrials = 100;
    for trialInd = numTrials:-1:1
        d = dFun(thetaTrue,true);
        Lambda = getLambda(@(x)dFun(x,false),d,sigmaN);
        [xMin,fMin] = fminsearch(@(x)Lambda(x),x0,opts);
        fVal(trialInd) = EmEll([xMin(1),0])-fMin;
    end
    delta = mean(fVal);
else
    delta = 3.01;   % <-- Value from code in "if false" statement
end

% Show the plot
hFig = figure();
h = surfc(relaxSign*theta2Array,theta1Array,EmEllArray,...
    'LineStyle','none');
hAxes = gca;
grid off
set(hAxes,'XDir','reverse');

levelList = get(h(2),'LevelList');
set(h(2),'LevelList',linspace(levelList(1),levelList(end),10));
set(h(2),'LineWidth',2);

view(82,40);
caxis([130 210]);   % Colors
set(hAxes,'Color','none');
prepareFigure();

[minCurve,minCurveInd] = min(EmEllArray,[],2);
minTheta2 = relaxSign*theta2Array(minCurveInd);

hold on
h2 = plot3(zeros(1,numGridPoints),theta1Array,EmEllArray(:,1),'-b');
set(h2,'LineWidth',3);
h3 = plot3(zeros(1,numGridPoints),theta1Array,minCurve,'-g');
set(h3,'LineWidth',3);
h4 = plot3(minTheta2,theta1Array,minCurve,'--g');
set(h4,'LineWidth',3);
hold off

% I used "parts" of this plot to place a transparent plane in front for improved
% visibility
xLimits = xlim();
yLimits = ylim();
zLimits = zlim();


end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   [phi,s] = biernackiTest(d,thetaHat,dFun,Lambda,<alpha>)
%
% PURPOSE:
%   Biernacki test for convergence to a local minima.
%  
% INPUT:
%   d           - Data associated with thetaHat
%   thetaHat    - Parameter estimate associated with a local minima
%   dFun        - A data generation function with the prototype
%                   d2 = dFun(theta,<addNoise>);
%   getLambda   - See the getLambda sub-function
%   alpha       - Type 1 error rate of the detector
%                   Default: 0.01
% 
% OUTPUT:
%   phi         - Detector output.  True if a local minima is detected and false
%                 otherwise 
%   s           - Test statistic
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function [phi,s] = biernackiTest(d,thetaHat,dFun,getLambda,alpha)

% Default values
if nargin<5
    alpha = .01;
end
numTrials = 50;
sigmaN = 1;     % <- Function is invariant to this

% Get the log-likelihood value at the estimate
Lambda = getLambda(dFun,d,sigmaN);
f = Lambda(thetaHat);

% Run the Monte-Carlo trials to estimate the mean and variance
for i = numTrials:-1:1
    % Generate data assuming the thetaHat is the true parameterization
    dTest = dFun(thetaHat);
    
    % Compute the log-likelihood assuming dTest
    Lambda = getLambda(dFun,dTest,sigmaN);
    fTest(i) = Lambda(thetaHat);
end

% Compute the test statistic
s = (f-mean(fTest))^2/var(fTest);

% Compute the decision statistic
s = gammainc(s/2,.5,'upper');
phi = s<alpha;

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   testFun = leblancTest(numDims);
%   [phi,s] = testFun(d,thetaHat,dFun,getLambda,numDims,<alpha>)
%
% PURPOSE:
%   LeBlanc test for convergence to a local minima.  Here we assume a location
%   family and incorporate the one-sided Biernacki test (numDims == 0).
%  
% INPUT:
%   
%   d           - Data associated with thetaHat
%   thetaHat    - Parameter estimate associated with a local minima
%   dFun        - A data generation function with the prototype
%                   d2 = dFun(theta,<addNoise>);
%   getLambda   - See the getLambda sub-function
%   alpha       - Type 1 error rate of the detector
%                   Default: 0.01
% 
% OUTPUT:
%   phi         - Detector output.  True if a local minima is detected and false
%                 otherwise 
%   s           - Test statistic
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function [phi,s,s2] = leblancTest(d,thetaHat,dFun,getLambda,numDims,alpha)

% Allow the number of relaxation dimensions to be returned programatically
if nargin==1
    numDims = d;
    phi = @(d,thetaHat,dFun,getLambda,varargin)leblancTest(...
        d,thetaHat,dFun,getLambda,numDims,varargin{:});
    return
end

% Default values
if nargin<6
    alpha = .01;
end
numTrials = 50;
sigmaN = 1;     % <- Function is invariant to this

% Thrm. 2 yields a nicer theoretical result if this restriction is in place
% Default: false
searchThetaPrimeOnly = false;   

% Get the log-likelihood value at the estimate
Lambda = getLambda(dFun,d,sigmaN);
f = Lambda(thetaHat);

% Run the Monte-Carlo trials to estimate the mean and variance
rngState = rng();               % We might need to replay this data
for i = numTrials:-1:1
    % Generate data assuming the thetaHat is the true parameterization
    dTest = dFun(thetaHat);
    
    % Compute the log-likelihood assuming dTest
    Lambda = getLambda(dFun,dTest,sigmaN);
    fTest(i) = Lambda(thetaHat);
end

%%% Perform the right-hand side of the test %%%
% This is a one-sided variation of Biernacki's test
fTestMu = mean(fTest);
fTestSigma = std(fTest);
s = (f-fTestMu)/fTestSigma;
s = .5 - .5*erf(s/sqrt(2));
phi = s<alpha/2;

% Exit early we have already decided against H0 and user isn't specifically
% interested in the 2nd half of the test
if phi && nargout<3  
    return
end

% Exit early if we are just running a one-sided Biernacki test
if numDims==0
    s2 = s;
    return
end

%%% Perform the left-hand side of the test %%%
LBFGSOpts = LBFGSOptions();
LBFGSOpts.initStep = .1;
theta0 = [thetaHat ; zeros(numDims,1)];

if searchThetaPrimeOnly
    % Precondition out dimensions we don't want to search
    LBFGSOpts.precond = zeros(size(theta0));
    LBFGSOpts.precond(end-numDims+1:end) = 1;
    warning('off','processPrecond:dimRemoved');
end

% Re-run the Monte-Carlo trials to get the log-likelihood under relaxation
rng(rngState);
for i = numTrials:-1:1
    % Generate data assuming the thetaHat is the true parameterization
    dTest = dFun(thetaHat);
    
    % Compute the log-likelihood assuming dTest
    Lambda = getLambda(dFun,dTest,sigmaN);
    [thetaHat2,fTestRelax(i),g2,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);
end
gap = fTest-fTestRelax;
gapMu = mean(gap);
gapSigma = std(gap);

% Solve under he relaxation
Lambda = getLambda(dFun,d,sigmaN);
[thetaHat2,f2,g2,exitFlag] = LBFGS(Lambda,theta0,LBFGSOpts);
 
s2 = ((f2-f)+gapMu)/gapSigma;
s2 = .5 + .5*erf(s2/sqrt(2));
phi2 = s2<alpha/2;

% Report w.r.t. the min of s and s2
if s2<s
    phi = phi2;
    temp = s;
    s = s2;
    s2 = temp;
end

end

%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   showMonteCarloPlots(thetaHatArray,s,thetaTrue,thetaLocalMin,HStr)
%
% PURPOSE:
%   This function helps generate consistent plots.
%
% INPUT:
%   thetaHatArray   - Array of parameter estimates
%   s               - Array of test statistics
%   thetaTrue       - True parameter
%   thetaLocalMin   - Parameter at local min
%   HStr            - "H0" or "H1"
%
% OUTPUT:
%   result      - Subreferenced object
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function showMonteCarloPlots(thetaHatArray,s,thetaTrue,thetaLocalMin,HStr)

figure;
hist(thetaHatArray,100);
yLimits = ylim;
hold on
h = plot(thetaTrue([1 1]),yLimits,'-r');
h = plot(thetaLocalMin([1 1]),yLimits,'--r');
hold off
xlabel('Theta');
ylabel('Count');
title(sprintf('Histogram of parameter estimates under %s',HStr));
legend('Histogram','True Solution','Local Minima');
prepareFigure();

figure;
hist(s,100);
xlabel('Test Statistic');
ylabel('Count');
title(sprintf('Test statistic under %s',HStr));
prepareFigure();


end