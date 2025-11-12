% ============================================================================
% MAIN NUMERICAL COMPUTATIONS
% ============================================================================
% This script performs parameter sweeps and computes solutions for the
% stablecoin model across a grid of lambda, varphi, and alpha values.
% ============================================================================

clearvars
set(0, 'DefaultFigureWindowStyle', 'docked')
addpath('./files/')
addpath('./code/')
dbstop if error

% ============================================================================
% PARAMETER SETUP
% ============================================================================
setup_parameters();

% Script-specific configuration
use_parallel = false;  % Enable/disable parallel processing

% ============================================================================
% INITIALIZATION
% ============================================================================
writepar(mfilename)
par = parfun;
options = optimset('tolF', 1e-20, 'tolX', 1e-20);

% ============================================================================
% SETUP
% ============================================================================
ParallelUtils.setupParallelPool(use_parallel);

% Create parameter grids
[lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid] = ...
    ParameterGridUtils.createAllGrids();
numLambda = length(lambdaGrid);
numAlpha = length(alphaGrid);

% Pre-compute file existence cache (speeds up file checking)
fprintf('Pre-computing file existence cache...\n');
fileExistenceCache = FileCacheUtils.createFileExistenceCache(...
    lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid);

% ============================================================================
% MAIN COMPUTATION LOOP
% ============================================================================
% Iterate over all lambda values and compute solutions for each
% varphi/alpha combination
% ============================================================================
for lambdaIdx = 1:numLambda
    currentLambda = lambdaGrid(lambdaIdx);
    
    % Solve initial system to get starting values for this lambda
    par.lambda = currentLambda;
    initialGuess = [amax; amax/2; 0];
    solutionVector = fsolve(@(x) Fsyst(x, par, 'case', 'jump commit'), ...
        initialGuess, options);
    [~, initialSolution] = Fsyst(solutionVector, par, 'case', 'jump commit');
    
    if strcmp(display_iter, 'on')
        fprintf('  Lambda %.3f: Residual = %.2e\n', ...
            currentLambda, sum(abs(initialSolution.F)));
    end
    
    % Set initial values for optimization
    defaultLowerBound = 0.02;
    defaultUpperBound = initialSolution.aub;
    defaultEpSa = initialSolution.EpSa;
    
    % Select appropriate varphi grid based on lambda value
    varphiGrid = ParameterGridUtils.selectVarphiGrid(...
        currentLambda, varphiGridLow, varphiGridHigh);
    
    % Process all varphi/alpha combinations for this lambda
    processLambdaCombinations(lambdaIdx, lambdaGrid, varphiGrid, alphaGrid, ...
        par, defaultLowerBound, defaultUpperBound, defaultEpSa, tolaub, ...
        fileExistenceCache, use_parallel);
end

% ============================================================================
% SAVE RESULTS
% ============================================================================
eval('save ''./files/savefile_grid''')
fprintf('Computation complete.\n')

% ============================================================================
% HELPER FUNCTIONS
% ============================================================================

function processLambdaCombinations(lambdaIdx, lambdaGrid, varphiGrid, alphaGrid, ...
    par, defaultLowerBound, defaultUpperBound, defaultEpSa, toleranceAub, ...
    fileCache, useParallel)
    % Process all varphi/alpha combinations for one lambda value
    %
    % This function handles both parallel and sequential processing.
    % In parallel mode, each worker processes one varphi value with all
    % alpha values. In sequential mode, a single solver instance is reused.
    
    numVarphi = length(varphiGrid);
    numAlpha = length(alphaGrid);
    currentLambda = lambdaGrid(lambdaIdx);
    
    % Decide whether to use parallel processing
    shouldUseParallel = useParallel && numVarphi > 1;
    
    if shouldUseParallel
        % PARALLEL MODE: Each worker processes one varphi value
        parfor varphiIdx = 1:numVarphi
            % Each worker needs its own solver instance
            solver = CommitCollateralClass();
            solver.alb00  = defaultLowerBound;
            solver.aub00  = defaultUpperBound;
            solver.EpSa00 = defaultEpSa;
            
            % Process all alpha values for this varphi
            for alphaIdx = 1:numAlpha
                processSingleCombination(solver, lambdaIdx, lambdaGrid, ...
                    varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, ...
                    currentLambda, toleranceAub, fileCache);
            end
        end
    else
        % SEQUENTIAL MODE: Reuse single solver instance
        solver = CommitCollateralClass();
        solver.alb00  = defaultLowerBound;
        solver.aub00  = defaultUpperBound;
        solver.EpSa00 = defaultEpSa;
        
        % Process all combinations sequentially
        for varphiIdx = 1:numVarphi
            for alphaIdx = 1:numAlpha
                processSingleCombination(solver, lambdaIdx, lambdaGrid, ...
                    varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, ...
                    currentLambda, toleranceAub, fileCache);
            end
        end
    end
end

function processSingleCombination(solver, lambdaIdx, lambdaGrid, varphiGrid, ...
    varphiIdx, alphaGrid, alphaIdx, par, currentLambda, toleranceAub, fileCache)
    % Process a single varphi/alpha combination
    %
    % This function checks if a solution file exists and is valid.
    % If not, it computes a new solution.
    
    % Generate filename for this parameter combination
    filename = getSaveFileName(alphaGrid(alphaIdx), varphiGrid(varphiIdx), ...
        currentLambda);
    
    % Check if file exists (using cache for speed)
    fileExists = fileCache.isKey(filename) && fileCache(filename);
    
    if fileExists
        % File exists - validate it
        try
            loadedSolution = load(filename, '-mat');
            
            % Check if solution structure is valid
            hasRequiredFields = isfield(loadedSolution, 'par') && ...
                isfield(loadedSolution, 'fstar') && ...
                isfield(loadedSolution, 'aub');
            
            if hasRequiredFields
                % Check if recomputation is needed (e.g., lambda changed)
                needsRecompute = FileCacheUtils.checkIfRecomputeNeeded(...
                    loadedSolution, currentLambda, toleranceAub);
                
                if needsRecompute
                    computeNewSolution(solver, lambdaIdx, lambdaGrid, ...
                        varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, fileCache);
                end
            else
                % Invalid file structure - recompute
                computeNewSolution(solver, lambdaIdx, lambdaGrid, ...
                    varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, fileCache);
            end
        catch
            % File corrupted or unreadable - recompute
            computeNewSolution(solver, lambdaIdx, lambdaGrid, ...
                varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, fileCache);
        end
    else
        % File doesn't exist - compute new solution
        computeNewSolution(solver, lambdaIdx, lambdaGrid, ...
            varphiGrid, varphiIdx, alphaGrid, alphaIdx, par, fileCache);
    end
end

function computeNewSolution(solver, lambdaIdx, lambdaGrid, varphiGrid, ...
    varphiIdx, alphaGrid, alphaIdx, par, fileCache)
    % Compute a new solution for the given parameter combination
    
    solver.optInit(lambdaGrid, lambdaIdx, varphiGrid, varphiIdx, ...
        alphaGrid, alphaIdx, fileCache);
    solver.optFunLambda(alphaGrid(alphaIdx), par, lambdaGrid, lambdaIdx, ...
        varphiGrid, varphiIdx);
end
