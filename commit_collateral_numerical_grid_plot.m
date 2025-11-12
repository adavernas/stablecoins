% ============================================================================
% COMMIT COLLATERAL NUMERICAL GRID PLOT
% ============================================================================
% This script loads or computes optimal solutions across a parameter grid
% and generates plots showing the relationship between lambda, varphi, and
% the optimal solutions.
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
use_parallel = true; % Enable/disable parallel processing

% Configuration options
loadsolution = 'on';  % 'on' to load saved solutions, 'off' to recompute
savegraph = 'off';   % 'on' to save graphs, 'off' to only display
display = 'off';     % 'on'/'all' to show detailed plots, 'off' to suppress

% ============================================================================
% INITIALIZATION
% ============================================================================
writepar(mfilename)
par = parfun;
options = optimset('tolF', 1e-20, 'tolX', 1e-20);

% Setup parallel pool
ParallelUtils.setupParallelPool(use_parallel);

% Initialize figures
initializeFigures(10);

if strcmp(loadsolution, 'off')
    % Compute solutions from scratch
    fprintf('Computing optimal solutions...\n');
    [lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid] = ParameterGridUtils.createAllGrids();
    
    % Pre-compute file existence cache
    fprintf('Pre-computing file existence cache...\n');
    fileCache = FileCacheUtils.createFileExistenceCache(...
        lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid);
    
    [optimalVarphiVec, optimalFstarVec, optimalAstarVec, referenceFunction, ...
        referenceFstarVec, referenceVarphiVec] = computeOptimalSolutions(...
        lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid, par, ...
        options, display, fileCache, use_parallel);
    
    save('commit_collateral_solution', 'lambdaGrid', 'optimalVarphiVec', 'optimalFstarVec', ...
        'optimalAstarVec', 'referenceFunction', 'referenceFstarVec', 'referenceVarphiVec')
    fprintf('Solutions computed and saved.\n');
else
    % Load precomputed solutions
    load('commit_collateral_solution')
    % Rename for compatibility with old variable names
    if exist('lvec', 'var')
        lambdaGrid = lvec;
    end
    if exist('vstarvec_lambda', 'var')
        optimalVarphiVec = vstarvec_lambda;
        optimalFstarVec = fstarvec_lambda;
        optimalAstarVec = astarvec_lambda;
    end
end

% Process and plot results
if sum(~isnan(optimalFstarVec)) > 4
    if ~exist('referenceFunction', 'var')
        referenceFunction = [];
        referenceFstarVec = [];
        referenceVarphiVec = [];
    end
    solution = processAndPlotResults(lambdaGrid, optimalVarphiVec, optimalFstarVec, ...
        optimalAstarVec, amax, par, options, savegraph, ...
        referenceFunction, referenceFstarVec, referenceVarphiVec);
end

%% Helper Functions

function initializeFigures(numFigures)
    % INITIALIZEFIGURES Clear and prepare figure windows
    
    for figNum = 1:numFigures
        figure(figNum); clf(figNum);
    end
end

function [optimalVarphiVec, optimalFstarVec, optimalAstarVec, referenceFunction, ...
    referenceFstarVec, referenceVarphiVec] = computeOptimalSolutions(...
    lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid, par, options, ...
    display, fileCache, use_parallel)
    % COMPUTEOPTIMALSOLUTIONS Compute optimal solutions for each lambda
    %
    % Supports parallel processing over lambda values
    
    numLambda = length(lambdaGrid);
    optimalVarphiVec = NaN(1, numLambda);
    optimalFstarVec = NaN(1, numLambda);
    optimalAstarVec = NaN(1, numLambda);
    
    referenceLambdaIdx = find(abs(lambdaGrid - 0.10) == min(abs(lambdaGrid - 0.10)));
    referenceFunction = [];
    referenceFstarVec = [];
    referenceVarphiVec = [];
    
    % Determine if we should use parallel processing
    shouldUseParallel = use_parallel && numLambda > 1;
    
    if shouldUseParallel
        % Parallel loop over lambda values
        tempOptimalVarphi = NaN(1, numLambda);
        tempOptimalFstar = NaN(1, numLambda);
        tempOptimalAstar = NaN(1, numLambda);
        
        parfor lambdaIdx = 1:numLambda
            [optimalVarphi, optimalFstar, optimalAstar] = processSingleLambda(...
                lambdaIdx, lambdaGrid, varphiGridLow, varphiGridHigh, ...
                alphaGrid, par, options, displayMode, fileCache);
            tempOptimalVarphi(lambdaIdx) = optimalVarphi;
            tempOptimalFstar(lambdaIdx) = optimalFstar;
            tempOptimalAstar(lambdaIdx) = optimalAstar;
        end
        
        % Copy results back
        optimalVarphiVec = tempOptimalVarphi;
        optimalFstarVec = tempOptimalFstar;
        optimalAstarVec = tempOptimalAstar;
        
        % Find reference function for referenceLambdaIdx (must be done sequentially)
        if ~isempty(referenceLambdaIdx)
            varphiGridRef = ParameterGridUtils.selectVarphiGrid(...
                lambdaGrid(referenceLambdaIdx), varphiGridLow, varphiGridHigh);
            [optimalAstarVecForRef, optimalFstarVecForRef, ~, ~, ~] = ...
                loadSolutionsForLambda(lambdaGrid(referenceLambdaIdx), varphiGridRef, ...
                alphaGrid, options, display, fileCache);
            [~, ~, ~, fittedFunction] = findOptimalVarphi(...
                varphiGridRef, optimalAstarVecForRef, optimalFstarVecForRef, ...
                referenceLambdaIdx, referenceLambdaIdx, lambdaGrid);
            if ~isempty(fittedFunction)
                referenceFunction = fittedFunction;
                referenceFstarVec = optimalFstarVecForRef;
                referenceVarphiVec = varphiGridRef;
            end
        end
    else
        % Sequential processing
        for lambdaIdx = 1:numLambda
            if mod(lambdaIdx, 10) == 0
                fprintf('  Processing lambda %d/%d...\n', lambdaIdx, numLambda);
            end
            
            % Select varphi grid based on lambda
            varphiGrid = ParameterGridUtils.selectVarphiGrid(...
                lambdaGrid(lambdaIdx), varphiGridLow, varphiGridHigh);
            
            % Load and process solutions for this lambda
            [optimalAstarVecForLambda, optimalFstarVecForLambda, ~, ~, ~] = ...
                loadSolutionsForLambda(lambdaGrid(lambdaIdx), varphiGrid, alphaGrid, ...
                options, display, fileCache);
            
            % Find optimal varphi and alpha for this lambda
            [optimalVarphi, optimalFstar, optimalAstar, fittedFunction] = findOptimalVarphi(...
                varphiGrid, optimalAstarVecForLambda, optimalFstarVecForLambda, ...
                lambdaIdx, referenceLambdaIdx, lambdaGrid);
            
            optimalVarphiVec(lambdaIdx) = optimalVarphi;
            optimalFstarVec(lambdaIdx) = optimalFstar;
            optimalAstarVec(lambdaIdx) = optimalAstar;
            
            % Store reference function for referenceLambdaIdx
            if lambdaIdx == referenceLambdaIdx && ~isempty(fittedFunction)
                referenceFunction = fittedFunction;
                referenceFstarVec = optimalFstarVecForLambda;
                referenceVarphiVec = varphiGrid;
            end
        end
    end
end

function [optimalVarphi, optimalFstar, optimalAstar] = processSingleLambda(...
    lambdaIdx, lambdaGrid, varphiGridLow, varphiGridHigh, alphaGrid, par, ...
    options, display, fileCache)
    % PROCESSSINGLELAMBDA Process a single lambda value (for parallel execution)
    %
    % Returns optimal values for this lambda
    
    if mod(lambdaIdx, 10) == 0
        fprintf('  Processing lambda %d/%d...\n', lambdaIdx, length(lambdaGrid));
    end
    
    % Select varphi grid based on lambda
    varphiGrid = ParameterGridUtils.selectVarphiGrid(...
        lambdaGrid(lambdaIdx), varphiGridLow, varphiGridHigh);
    
    % Load and process solutions for this lambda
    [optimalAstarVec, optimalFstarVec, ~, ~, ~] = ...
        loadSolutionsForLambda(lambdaGrid(lambdaIdx), varphiGrid, alphaGrid, ...
        options, display, fileCache);
    
    % Find optimal varphi and alpha for this lambda
    referenceLambdaIdx = find(abs(lambdaGrid - 0.10) == min(abs(lambdaGrid - 0.10)));
    [optimalVarphi, optimalFstar, optimalAstar, ~] = findOptimalVarphi(...
        varphiGrid, optimalAstarVec, optimalFstarVec, lambdaIdx, referenceLambdaIdx, lambdaGrid);
end

function [optimalAstarVec, optimalFstarVec, fstarMatrix, astarMatrix, aubMatrix] = ...
    loadSolutionsForLambda(lambdaValue, varphiGrid, alphaGrid, options, display, fileCache)
    % LOADSOLUTIONSFORLAMBDA Load and validate solutions for a given lambda
    %
    % OPTIMIZED: Uses file cache, skips expensive aub recomputation when possible
    
    numVarphi = length(varphiGrid);
    numAlpha = length(alphaGrid);
    
    optimalAstarVec = NaN(1, numVarphi);
    optimalFstarVec = NaN(1, numVarphi);
    fstarMatrix = NaN(numVarphi, numAlpha);
    astarMatrix = NaN(numVarphi, numAlpha);
    aubMatrix = NaN(numVarphi, numAlpha);
    
    for varphiIdx = 1:numVarphi
        for alphaIdx = 1:numAlpha
            filename = getSaveFileName(alphaGrid(alphaIdx), varphiGrid(varphiIdx), lambdaValue);
            
            % Check cache first (much faster than file system)
            if fileCache.isKey(filename) && fileCache(filename)
                try
                    loadedSolution = load(filename, '-mat');
                    
                    % Quick validation without expensive aub recomputation
                    if FileCacheUtils.validateSolution(loadedSolution, lambdaValue, 1e-2)
                        fstarMatrix(varphiIdx, alphaIdx) = loadedSolution.fstar;
                        astarMatrix(varphiIdx, alphaIdx) = loadedSolution.astar;
                        aubMatrix(varphiIdx, alphaIdx) = loadedSolution.aub;
                        
                        % Display if requested
                        if strcmp(display, 'on') || strcmp(display, 'all')
                            displaySolutionDetails(loadedSolution);
                        end
                    end
                catch
                    % File corrupted - skip
                end
            end
        end
        
        % Find optimal alpha for this varphi
        [optimalAstarVec(varphiIdx), optimalFstarVec(varphiIdx)] = findOptimalAlpha(...
            alphaGrid, fstarMatrix(varphiIdx, :), astarMatrix(varphiIdx, :), ...
            aubMatrix(varphiIdx, :), display);
    end
end

function displaySolutionDetails(loadedSolution)
    % DISPLAYSOLUTIONDETAILS Display detailed solution information
    
    alphaVector = linspace(0, 2, 100);
    
    % Plot e(a)
    figure(1); hold on;
    plot(alphaVector, loadedSolution.efun(alphaVector), 'b')
    plot(loadedSolution.aub, loadedSolution.efun(loadedSolution.aub), '*k')
    plot(loadedSolution.astar, loadedSolution.efun(loadedSolution.astar), '*r')
    if isfield(loadedSolution, 'astar_')
        plot(loadedSolution.astar_, loadedSolution.efun(loadedSolution.astar_), 'or')
    end
    title('$e(a)$', 'Interpreter', 'LateX')
    
    % Plot p(a)
    figure(2); hold on;
    plot(alphaVector, loadedSolution.pfun(alphaVector), 'b')
    plot(loadedSolution.aub, loadedSolution.pfun(loadedSolution.aub), '*k')
    plot(loadedSolution.astar, loadedSolution.pfun(loadedSolution.astar), '*r')
    if isfield(loadedSolution, 'astar_')
        plot(loadedSolution.astar_, loadedSolution.pfun(loadedSolution.astar_), 'or')
    end
    title('$p(a)$', 'Interpreter', 'LateX')
    
    % Plot E[p(Sa)]
    figure(3); hold on;
    plot(alphaVector, loadedSolution.EpSa(alphaVector), 'b')
    plot(loadedSolution.aub, loadedSolution.EpSa(loadedSolution.aub), '*k')
    plot(loadedSolution.astar, loadedSolution.EpSa(loadedSolution.astar), '*r')
    if isfield(loadedSolution, 'astar_')
        plot(loadedSolution.astar_, loadedSolution.EpSa(loadedSolution.astar_), 'or')
    end
    ylim([0.0 1])
    title('$E[p(Sa)]$', 'Interpreter', 'LateX')
    
    % Plot eb(a)
    figure(4); hold on;
    if isfield(loadedSolution, 'eb') && isfield(loadedSolution, 'astar')
        plot(alphaVector, abs(loadedSolution.eb(alphaVector, alphaVector, ...
            max(alphaVector, loadedSolution.astar))), 'b')
        plot(loadedSolution.aub, loadedSolution.eb(loadedSolution.aub, loadedSolution.aub, ...
            max(loadedSolution.aub, loadedSolution.astar)), '*k')
    end
    title('$eb(a)$', 'Interpreter', 'LateX')
    
    % Display numerical values
    if isfield(loadedSolution, 'estar') && isfield(loadedSolution, 'astar_') && ...
            isfield(loadedSolution, 'par')
        fstarComputed = (loadedSolution.estar(loadedSolution.aub, loadedSolution.astar_) + ...
            1 - loadedSolution.par.varphi0) / loadedSolution.astar_;
        disp(['E[p(Sastar)]: ', num2str(loadedSolution.EpSa(loadedSolution.astar_)), ...
            ' lambda: ', num2str(loadedSolution.par.lambda), ...
            ' varphi: ', num2str(loadedSolution.par.varphi0), ...
            ' fstar: ', num2str(loadedSolution.fstar), ...
            ' astar: ', num2str(loadedSolution.astar), ...
            ' aub: ', num2str(loadedSolution.aub)]);
    end
end

function [optimalAlpha, optimalFstar] = findOptimalAlpha(alphaGrid, fstarVector, ...
    astarVector, aubVector, display)
    % FINDOPTIMALALPHA Find optimal alpha for given varphi
    
    optimalAlpha = NaN;
    optimalFstar = NaN;
    
    % Find valid solutions (fstar > 0, not NaN, astar >= aub)
    validIndices = and(fstarVector > 0, and(~isnan(fstarVector), astarVector >= aubVector));
    
    if sum(validIndices) <= 4
        return;  % Need at least 4 points for polynomial fit
    end
    
    % Find first and last valid indices
    firstValidIdx = min(find(validIndices == 1, 1, 'first'), ...
        find(validIndices == 0, 1, 'last'));
    if isempty(firstValidIdx)
        firstValidIdx = 1;
    end
    
    lastValidIdx = find(validIndices == 1, 1, 'last');
    
    % Fit polynomial and find maximum
    polynomialCoeffs = polyfit(alphaGrid(validIndices), fstarVector(validIndices), 3);
    fittedFunction = @(a) polyval(polynomialCoeffs, a);
    
    optimalAlpha = fminbnd(@(a) -fittedFunction(a), 1, alphaGrid(lastValidIdx));
    optimalFstar = fittedFunction(optimalAlpha);
    
    % Display if requested
    if strcmp(display, 'on') || strcmp(display, 'all')
        alphaVectorFine = linspace(min(alphaGrid(validIndices)), max(alphaGrid(validIndices)), 200);
        figure(5); hold on;
        plot(alphaGrid, fstarVector, 'ob')
        plot(alphaGrid(validIndices), fstarVector(validIndices), 'ok')
        plot(alphaVectorFine, fittedFunction(alphaVectorFine))
        plot(optimalAlpha, optimalFstar, '*r')
        title('$f^\star(a^\star;\varphi,\lambda)$', 'Interpreter', 'LateX')
    end
end

function [optimalVarphi, optimalFstar, optimalAstar, fittedFunction] = ...
    findOptimalVarphi(varphiGrid, optimalAstarVec, optimalFstarVec, lambdaIdx, ...
    referenceLambdaIdx, lambdaGrid)
    % FINDOPTIMALVARPHI Find optimal varphi for given lambda
    
    optimalVarphi = NaN;
    optimalFstar = NaN;
    optimalAstar = NaN;
    fittedFunction = [];
    
    % Determine index offset based on lambda
    indexOffset = getIndexOffset(lambdaIdx);
    
    % Find valid solutions (exclude first few points)
    validIndices = and(optimalFstarVec > 0, ~isnan(optimalFstarVec));
    validIndices(1:indexOffset) = false;
    
    if sum(validIndices) <= 4
        return;  % Need at least 4 points for polynomial fit
    end
    
    % Fit polynomials
    fstarPolynomialCoeffs = polyfit(varphiGrid(validIndices), optimalFstarVec(validIndices), 3);
    fittedFunction = @(a) polyval(fstarPolynomialCoeffs, a);
    
    astarPolynomialCoeffs = polyfit(varphiGrid(validIndices), optimalAstarVec(validIndices), 3);
    astarFittedFunction = @(a) polyval(astarPolynomialCoeffs, a);
    
    % Find optimal varphi
    optimalVarphi = fminbnd(@(a) -fittedFunction(a), varphiGrid(indexOffset), varphiGrid(end));
    optimalFstar = fittedFunction(optimalVarphi);
    optimalAstar = astarFittedFunction(optimalVarphi);
    
    % Plot results
    plotVarphiResults(varphiGrid, optimalFstarVec, optimalAstarVec, ...
        optimalVarphi, optimalFstar, fittedFunction);
end

function indexOffset = getIndexOffset(lambdaIdx)
    % GETINDEXOFFSET Get index offset based on lambda index
    
    if lambdaIdx == 24
        indexOffset = 4;
    elseif lambdaIdx == 27
        indexOffset = 3;
    else
        indexOffset = 1;
    end
end

function plotVarphiResults(varphiGrid, optimalFstarVec, optimalAstarVec, ...
    optimalVarphi, optimalFstar, fittedFunction)
    % PLOTVARPHIRESULTS Plot varphi optimization results
    
    varphiVectorFine = linspace(0, 1, 100);
    
    figure(6); hold on;
    plot(varphiVectorFine, fittedFunction(varphiVectorFine))
    plot(varphiGrid, optimalFstarVec, 'ok')
    plot(optimalVarphi, optimalFstar, '*r')
    title('$f^\star(\varphi;\lambda)$', 'Interpreter', 'LateX')
    ylim([0 2])
    
    figure(7); hold on;
    plot(varphiGrid, optimalAstarVec, 'ok')
    title('$a^\star(\varphi;\lambda)$', 'Interpreter', 'LateX')
    ylim([0 4])
end

function solution = processAndPlotResults(lambdaGrid, optimalVarphiVec, optimalFstarVec, ...
    optimalAstarVec, amax, par, options, savegraph, referenceFunction, ...
    referenceFstarVec, referenceVarphiVec)
    % PROCESSANDPLOTRESULTS Process results and create final plots
    
    solution = struct();
    
    % Plot astar vs lambda
    plotAstarVsLambda(lambdaGrid, optimalAstarVec);
    
    % Create varphi star function
    solution.vstarfun = createVarphiStarFunction(lambdaGrid, optimalVarphiVec);
    plotVarphiStarFunction(lambdaGrid, optimalVarphiVec, solution.vstarfun);
    
    % Compute solution for lambda = 0
    [solution0, initialGuess] = computeLambdaZeroSolution(amax, par, options);
    
    % Compute solutions for refined lambda grid
    [lambdaGridRefined, fstarVecLambda0, astarVecLambda0, aubVecLambda0] = ...
        computeRefinedLambdaSolutions(lambdaGrid, initialGuess, par, options);
    
    % Fit functions and create final plots
    solution = fitFunctionsAndPlot(lambdaGrid, lambdaGridRefined, optimalFstarVec, ...
        fstarVecLambda0, astarVecLambda0, aubVecLambda0, solution0, solution);
    
    % Store reference function
    if exist('referenceFunction', 'var') && ~isempty(referenceFunction)
        solution.ffun0 = referenceFunction;
        solution.fstarvec0 = referenceFstarVec;
        solution.vvec0 = referenceVarphiVec;
    end
    
    % Compute varphi-dependent solutions
    solution = computeVarphiSolutions(par, options, solution);
    
    % Save graphs if requested
    if strcmp(savegraph, 'on')
        figure(1); clf(1);
        plotgraphs(par, 'graphcase', 'commit collateral numerical', ...
            'sol', solution, 'sol00', struct(), 'sol0', solution0, 'sol1', struct(), ...
            'name', 'commit_collateral_numerical', 'savegraph', 'on')
    end
    
    solution.lvec = lambdaGrid;
end

function plotAstarVsLambda(lambdaGrid, optimalAstarVec)
    % PLOTASTARVSLAMBDA Plot astar as function of lambda
    
    validIndices = find(optimalAstarVec > 0);
    
    figure(7); clf(7); hold on;
    plot(lambdaGrid(validIndices), optimalAstarVec(validIndices), 'ok')
    title('$a^\star(\lambda)$', 'Interpreter', 'LateX')
end

function varphiStarFunction = createVarphiStarFunction(lambdaGrid, optimalVarphiVec)
    % CREATEVARPHISTARFUNCTION Create interpolated function for varphi star
    
    validIndices = find(optimalVarphiVec > 0.02);
    validIndices = [validIndices(1)-1 validIndices(1:end-1)];
    
    lambdaGridTemp = linspace(lambdaGrid(validIndices(1)), 1, 10);
    
    varphiStarFunction = @(a) (a >= lambdaGrid(validIndices(1))) .* ...
        interp1(lambdaGridTemp, ...
        interp1([lambdaGrid(validIndices) 1], [optimalVarphiVec(validIndices) 1], ...
        lambdaGridTemp, 'pchip'), a, 'pchip');
end

function plotVarphiStarFunction(lambdaGrid, optimalVarphiVec, varphiStarFunction)
    % PLOTVARPHISTARFUNCTION Plot varphi star function
    
    lambdaVectorFine = linspace(0, lambdaGrid(end), 200);
    
    figure(8); clf(8); hold on;
    plot(lambdaVectorFine, varphiStarFunction(lambdaVectorFine))
    plot(lambdaGrid, optimalVarphiVec, 'ok')
    title('$\varphi^\star(\lambda)$', 'Interpreter', 'LateX')
    ylim([0 1])
end

function [solution0, initialGuess] = computeLambdaZeroSolution(alphaMax, par, options)
    % COMPUTELAMBDAZEROSOLUTION Compute solution for lambda = 0
    
    initialGuess = [alphaMax; alphaMax/2; 0];
    
    par.lambda = 0;
    [initialGuess, ~] = fmincon(@(x) -FsystCon(x, par, 'case', 'jump commit'), initialGuess, ...
        [], [], [], [], [], [], @(x) mycon(x, par, 'case', 'jump commit'), options);
    [~, solution0] = Fsyst(initialGuess, par, 'case', 'jump commit');
    disp(sum(abs(solution0.F)))
end

function [lambdaGridRefined, fstarVecLambda0, astarVecLambda0, aubVecLambda0] = ...
    computeRefinedLambdaSolutions(lambdaGrid, initialGuess, par, options)
    % COMPUTEREFINEDLAMBDASOLUTIONS Compute solutions for refined lambda grid
    
    lambdaGridRefined = unique(sort([lambdaGrid 0.2:0.025:0.30]));
    lambdaGridRefined = lambdaGridRefined(lambdaGridRefined <= 0.30);
    numLambda = length(lambdaGridRefined);
    
    fstarVecLambda0 = zeros(1, numLambda);
    astarVecLambda0 = NaN(1, numLambda);
    aubVecLambda0 = NaN(1, numLambda);
    
    for lambdaIdx = 1:numLambda
        par.lambda = lambdaGridRefined(lambdaIdx);
        
        [initialGuess, ~] = fmincon(@(x) -FsystCon(x, par, 'case', 'jump commit'), initialGuess, ...
            [], [], [], [], [], [], @(x) mycon(x, par, 'case', 'jump commit'), options);
        [~, solution0] = Fsyst(initialGuess, par, 'case', 'jump commit');
        disp(sum(abs(solution0.F)))
        
        if sum(abs(solution0.F)) < 1e-5
            aubVecLambda0(lambdaIdx) = solution0.aub;
            astarVecLambda0(lambdaIdx) = solution0.astar;
            fstarVecLambda0(lambdaIdx) = solution0.fstar;
        end
    end
end

function solution = fitFunctionsAndPlot(lambdaGrid, lambdaGridRefined, optimalFstarVec, ...
    fstarVecLambda0, astarVecLambda0, aubVecLambda0, solution0, solution)
    % FITFUNCTIONSANDPLOT Fit functions to data and create plots
    
    % Fit fstar function
    validIndices = find(~isnan(optimalFstarVec));
    polynomialCoeffs = polyfit(lambdaGrid(validIndices), log(optimalFstarVec(validIndices)), 3);
    solution.fstarfun = @(a) exp(polyval(polynomialCoeffs, a));
    
    % Fit fstar0 function
    validIndices0 = find(fstarVecLambda0 > 0);
    polynomialCoeffs0 = polyfit(lambdaGridRefined(validIndices0), fstarVecLambda0(validIndices0), 4);
    solution0.fstarfun0 = @(a) polyval(polynomialCoeffs0, a);
    solution0.lvec0 = fsolve(@(a) solution0.fstarfun0(a), 0.2);
    
    % Plot fstar vs lambda
    lambdaVectorFine = linspace(0, lambdaGrid(end), 200);
    figure(9); clf(9); hold on;
    plot(lambdaVectorFine, solution.fstarfun(lambdaVectorFine))
    plot(lambdaVectorFine, solution0.fstarfun0(lambdaVectorFine), '--b')
    plot(lambdaGridRefined, fstarVecLambda0, 'b')
    plot(lambdaGrid, optimalFstarVec, 'ok')
    title('$f^\star(\lambda)$', 'Interpreter', 'LateX')
    ylim([0 2])
    
    % Plot astar ratio
    figure(10); clf(10); hold on;
    plot(lambdaGridRefined(1:length(astarVecLambda0)), ...
        astarVecLambda0 ./ aubVecLambda0, 'b')
    title('$a^\star(\lambda)$', 'Interpreter', 'LateX')
    
    % Fit difference functions
    % Interpolate fstarVecLambda0 onto lambdaGrid to match dimensions
    validIndicesRefined = ~isnan(fstarVecLambda0) & fstarVecLambda0 > 0;
    if sum(validIndicesRefined) > 0
        % Interpolate refined grid values onto original grid
        fstarVecLambda0Interp = interp1(lambdaGridRefined(validIndicesRefined), ...
            fstarVecLambda0(validIndicesRefined), lambdaGrid, 'linear', NaN);
        
        % Find indices where both vectors have valid data
        validIndicesBoth = ~isnan(optimalFstarVec) & ~isnan(fstarVecLambda0Interp) & ...
            optimalFstarVec > 0 & fstarVecLambda0Interp > 0;
        
        if sum(validIndicesBoth) > 3
            polynomialCoeffsDiff = polyfit(lambdaGrid(validIndicesBoth), ...
                log(max(0.01, optimalFstarVec(validIndicesBoth) - fstarVecLambda0Interp(validIndicesBoth))), 3);
            solution0.dfstarfun0 = @(a) exp(polyval(polynomialCoeffsDiff, a));
        else
            % Fallback: create a simple function if not enough overlap
            solution0.dfstarfun0 = @(a) 0.01 * ones(size(a));
        end
    else
        % Fallback if no valid refined data
        solution0.dfstarfun0 = @(a) 0.01 * ones(size(a));
    end
    
    solution.fstarvec_lambda = optimalFstarVec;
    solution.sol0 = solution0;
end

function solution = computeVarphiSolutions(par, options, solution)
    % COMPUTEVARPHISOLUTIONS Compute solutions as function of varphi
    %
    % NOTE: This function requires 'jump collateral delta infinite' case in Fsyst
    % If not available, it will skip this computation
    
    % Check if the required case is available in Fsyst
    % For now, we'll skip this computation as the case is not implemented
    % Uncomment and implement the case in Fsyst.m if needed
    
    % Skip computation - case not implemented in current Fsyst
    warning('computeVarphiSolutions:CaseNotImplemented', ...
        'Skipping varphi-dependent solutions: ''jump collateral delta infinite'' case not implemented in Fsyst.m');
    
    numVarphi = 30;
    varphiVector = linspace(0, 0.98, numVarphi);
    fstarVecVarphi0 = NaN(1, numVarphi);
    
    % Alternative: Use 'jump commit' case if available
    % Uncomment below to use alternative approach
    %{
    par.lambda = 0.10;
    initialGuess = [amax; amax/2; 0];
    
    for varphiIdx = 1:numVarphi-1
        par.varphi0 = varphiVector(varphiIdx);
        
        try
            solutionVector = fsolve(@(x) Fsyst(x, par, 'case', 'jump commit'), initialGuess, options);
            [~, solution2] = Fsyst(solutionVector, par, 'case', 'jump commit');
            
            if sum(abs(solution2.F)) < 1e-5
                fstarVecVarphi0(varphiIdx) = solution2.fstar;
            end
            initialGuess = solutionVector;  % Use as next initial guess
        catch
            % Skip if computation fails
        end
    end
    %}
    
    % Plot and fit (only if we have valid data)
    if sum(~isnan(fstarVecVarphi0)) > 0
        figure(11); clf(11); hold on;
        if isfield(solution, 'ffun0')
            plot(varphiVector, solution.ffun0(varphiVector))
            plot(solution.vvec0, solution.fstarvec0, 'ok')
        end
        plot(varphiVector, fstarVecVarphi0)
        
        % Find optimal varphi
        validIndices = ~isnan(fstarVecVarphi0);
        if sum(validIndices) > 3
            polynomialCoeffs = polyfit(varphiVector(validIndices), fstarVecVarphi0(validIndices), 3);
            fittedFunction = @(a) polyval(polynomialCoeffs, a);
            solution.varphi = fminbnd(@(a) -fittedFunction(a), 0, 1);
            solution.fstar = fittedFunction(solution.varphi);
        end
    end
    
    if isfield(solution, 'ffun0')
        solution.varphi0 = fminbnd(@(a) -solution.ffun0(a), 0, 1);
        solution.fstar0 = solution.ffun0(solution.varphi0);
    end
    
    solution.vvec_ = varphiVector;
    solution.fstarvec_varphi0 = fstarVecVarphi0;
    
    % Final plot (only if we have data to plot)
    if sum(~isnan(fstarVecVarphi0)) > 0 || isfield(solution, 'ffun0')
        figure(2); clf(2);
        plotgraphs(par, 'graphcase', 'no commit collateral numerical', ...
            'sol', solution, 'name', 'nocommit_collateral_numerical', 'savegraph', 'off')
    end
end
