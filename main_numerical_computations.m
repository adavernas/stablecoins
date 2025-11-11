clearvars
set(0,'DefaultFigureWindowStyle','docked')
addpath('./files/')
addpath('./code/')
dbstop if error

%% PARAMETERS
% Model parameters (names must match parfun expectations)
r   = 0.06;      % Interest rate
xmu = 0.05;      % Drift parameter
sig = 0.1;       % Volatility

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 6;
lambda = 0.1;

muk = 0.055;

delt0 = (r + sig^2/2 - exp(-1)*r)*1.15;
varphi0 = 0.99;

% Parameter validation
if r - xmu + lambda/(xi+1) < 0
    warning('ParameterCheck:NegativePhi', 'phi is negative! This may cause numerical issues.');
end

if r + xmu - sig^2 - lambda/(xi-1) < 0
    warning('ParameterCheck:NegativeKappa', 'kappa is negative! This may cause numerical issues.');
end

gmin = -1;

% Solver tolerances
odeTol = 1e-8;
tolalb = 1e-6;
tolESa = 1e-3;
tolaub = 1e-2;
display_plot = 'off';
display_iter = 'off';

% Parallel computation option
use_parallel = false;  % Set to true to enable parallel processing

%% END
writepar(mfilename)
par = parfun;
options = optimset('tolF', 1e-20, 'tolX', 1e-20);

% Setup parallel pool
ParallelUtils.setupParallelPool(use_parallel);

% Create parameter grids
[lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid] = ParameterGridUtils.createAllGrids();
numLambda = length(lambdaGrid);
numAlpha = length(alphaGrid);

% Pre-compute file existence cache
fprintf('Pre-computing file existence cache...\n');
fileExistenceCache = FileCacheUtils.createFileExistenceCache(...
    lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid);

% Main computation loop
for lambdaIdx = 1:numLambda
    currentLambda = lambdaGrid(lambdaIdx);
    
    % Solve initial system to get starting values
    par.lambda = currentLambda;
    initialGuess = [amax; amax/2; 0];
    solutionVector = fsolve(@(x) Fsyst(x, par, 'case', 'jump commit'), initialGuess, options);
    [~, initialSolution] = Fsyst(solutionVector, par, 'case', 'jump commit');
    
    if strcmp(display_iter, 'on')
        disp(sum(abs(initialSolution.F)))
    end
    
    % Set initial values for this lambda
    defaultLowerBound = 0.02;
    defaultUpperBound = initialSolution.aub;
    defaultEpSa = initialSolution.EpSa;
    
    % Select varphi grid based on lambda value
    varphiGrid = ParameterGridUtils.selectVarphiGrid(currentLambda, varphiGridLow, varphiGridHigh);
    numVarphi = length(varphiGrid);
    
    % Process all combinations for this lambda
    processLambdaCombinations(lambdaIdx, lambdaGrid, varphiGrid, alphaGrid, par, ...
        defaultLowerBound, defaultUpperBound, defaultEpSa, tolaub, ...
        fileExistenceCache, use_parallel);
end

eval('save ''./files/savefile_grid''')
fprintf('Computation complete.\n')

%% Helper Functions

function processLambdaCombinations(lambdaIdx, lambdaGrid, varphiGrid, alphaGrid, par, ...
    defaultLowerBound, defaultUpperBound, defaultEpSa, toleranceAub, fileCache, useParallel)
    % PROCESSLAMBDACOMBINATIONS Process all varphi/alpha combinations for one lambda
    %
    % Optimized to minimize file I/O and class instantiation
    % Supports parallel processing over varphi values
    
    numVarphi = length(varphiGrid);
    numAlpha = length(alphaGrid);
    currentLambda = lambdaGrid(lambdaIdx);
    
    % Determine if we should use parallel processing
    shouldUseParallel = useParallel && numVarphi > 1;
    
    if shouldUseParallel
        % Parallel loop over varphi values (each worker processes one varphi)
        parfor varphiIdx = 1:numVarphi
            % Each worker needs its own class instance
            solver = CommitCollateralClass();
            solver.alb00  = defaultLowerBound;
            solver.aub00  = defaultUpperBound;
            solver.EpSa00 = defaultEpSa;
            
            % Process all alpha values for this varphi
            for alphaIdx = 1:numAlpha
                processSingleCombination(solver, lambdaIdx, lambdaGrid, varphiGrid, varphiIdx, ...
                    alphaGrid, alphaIdx, par, currentLambda, toleranceAub, fileCache);
            end
        end
    else
        % Sequential processing - reuse single class instance
        solver = CommitCollateralClass();
        solver.alb00  = defaultLowerBound;
        solver.aub00  = defaultUpperBound;
        solver.EpSa00 = defaultEpSa;
        
        % Process all combinations
        for varphiIdx = 1:numVarphi
            for alphaIdx = 1:numAlpha
                processSingleCombination(solver, lambdaIdx, lambdaGrid, varphiGrid, varphiIdx, ...
                    alphaGrid, alphaIdx, par, currentLambda, toleranceAub, fileCache);
            end
        end
    end
end

function processSingleCombination(solver, lambdaIdx, lambdaGrid, varphiGrid, varphiIdx, ...
    alphaGrid, alphaIdx, par, currentLambda, toleranceAub, fileCache)
    % PROCESSSINGLECOMBINATION Process a single varphi/alpha combination
    %
    % Extracted to support both parallel and sequential execution
    
    filename = getSaveFileName(alphaGrid(alphaIdx), varphiGrid(varphiIdx), currentLambda);
    
    % Check cache instead of file system
    if fileCache.isKey(filename) && fileCache(filename)
        % File exists - check if valid
        try
            loadedSolution = load(filename, '-mat');
            
            % Quick validation without loading full structure
            if isfield(loadedSolution, 'par') && isfield(loadedSolution, 'fstar') && ...
                    isfield(loadedSolution, 'aub')
                % Check if recomputation needed
                needsRecompute = FileCacheUtils.checkIfRecomputeNeeded(...
                    loadedSolution, currentLambda, toleranceAub);
                
                if needsRecompute
                    solver.optInit(lambdaGrid, lambdaIdx, varphiGrid, varphiIdx, ...
                        alphaGrid, alphaIdx, fileCache);
                    solver.optFunLambda(alphaGrid(alphaIdx), par, lambdaGrid, lambdaIdx, ...
                        varphiGrid, varphiIdx);
                end
            else
                % Invalid file structure - recompute
                solver.optInit(lambdaGrid, lambdaIdx, varphiGrid, varphiIdx, ...
                    alphaGrid, alphaIdx, fileCache);
                solver.optFunLambda(alphaGrid(alphaIdx), par, lambdaGrid, lambdaIdx, ...
                    varphiGrid, varphiIdx);
            end
        catch
            % File corrupted or unreadable - recompute
            solver.optInit(lambdaGrid, lambdaIdx, varphiGrid, varphiIdx, ...
                alphaGrid, alphaIdx, fileCache);
            solver.optFunLambda(alphaGrid(alphaIdx), par, lambdaGrid, lambdaIdx, ...
                varphiGrid, varphiIdx);
        end
    else
        % File doesn't exist - compute solution
        solver.optInit(lambdaGrid, lambdaIdx, varphiGrid, varphiIdx, ...
            alphaGrid, alphaIdx, fileCache);
        solver.optFunLambda(alphaGrid(alphaIdx), par, lambdaGrid, lambdaIdx, ...
            varphiGrid, varphiIdx);
    end
end
