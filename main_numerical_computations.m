clearvars
set(0,'DefaultFigureWindowStyle','docked')
addpath('./files/')
addpath('./code/')
dbstop if error

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

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

%% END
writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20);

delete(gcp('nocreate'))

% Create parameter grids
lvec = createLambdaGrid();
nl = length(lvec);

[vvec_low, vvec_high] = createVarphiGrids();

% Create alpha grid
avec = createAlphaGrid();
na = length(avec);

% Main computation loop
for il=1:nl
    
    % Solve initial system to get starting values
    par.lambda = lvec(il);
    x1 = fsolve(@(x) Fsyst(x, par, 'case', 'jump commit'), [amax;amax/2;0], options);
    [~,sol1] = Fsyst(x1,par,'case','jump commit');
    disp(sum(abs(sol1.F)))
    
    % Set initial values for this lambda
    alb00  = 0.02;
    aub00  = sol1.aub;
    EpSa00 = sol1.EpSa;
    
    % Select varphi grid based on lambda value
    if lvec(il) <= 0.15
        vvec = vvec_low;
    else
        vvec = vvec_high;
    end
    nv = length(vvec);
    
    % Parallel loop over varphi values
    parfor iv=1:nv
        optClass = CommitCollateralClass();
        optClass.alb00 = alb00;
        optClass.aub00 = aub00;
        optClass.EpSa00 = EpSa00;
        
        % Loop over alpha values
        for ia=1:na
            filename = getSaveFileName(avec(ia), vvec(iv), lvec(il));
            
            % Check if file exists and is valid
            if exist(filename, 'file')
                f = load(filename);
                f.fstar_ = (f.estar(f.aub, f.astar_) + 1 - f.par.varphi0) / f.astar_;
                
                % Check if recomputation is needed
                needsRecompute = checkIfRecomputeNeeded(f, lvec(il), tolaub);
                
                if needsRecompute
                    optClass.optInit(lvec, il, vvec, iv, avec, ia);
                    optClass.optFunLambda(avec(ia), par, lvec, il, vvec, iv);
                end
            else
                % File doesn't exist, compute solution
                optClass.optInit(lvec, il, vvec, iv, avec, ia);
                optClass.optFunLambda(avec(ia), par, lvec, il, vvec, iv);
            end
        end
    end
end

eval('save ''./files/savefile_grid''')

%% Helper Functions

function lvec = createLambdaGrid()
    % CREATELAMBDAGRID Create refined lambda parameter grid
    %
    % Creates a grid with initial spacing, then adds midpoints twice,
    % and appends a high-resolution tail.
    
    % Base grid
    baseGrid = linspace(0.01, 0.5, 20);
    
    % Add midpoints (refine grid twice)
    refinedGrid = baseGrid;
    for refinement = 1:2
        midpoints = (refinedGrid(2:end) + refinedGrid(1:end-1)) / 2;
        refinedGrid = sort([refinedGrid, midpoints]);
    end
    
    % Combine with high-resolution tail
    lowRes = refinedGrid(1:34);
    highRes = 0.25:0.05:1;
    lvec = sort([lowRes, highRes]);
end

function [vvec_low, vvec_high] = createVarphiGrids()
    % CREATEVARPHIGRIDS Create varphi grids for different lambda ranges
    %
    % Outputs:
    %   vvec_low  - Grid for lambda <= 0.15 (sparse)
    %   vvec_high - Grid for lambda > 0.15 (dense)
    
    nv_tmp = 20;
    vvec_tmp = linspace(0.01, 0.99, nv_tmp);
    vvec_low  = vvec_tmp(1:2:nv_tmp);  % Every other point
    vvec_high = vvec_tmp(8:nv_tmp);    % Higher range
end

function avec = createAlphaGrid()
    % CREATEALPHAGRID Create refined alpha parameter grid
    %
    % Creates base grid and adds midpoints for values <= 1.5
    
    na = 10;
    avec_tmp = linspace(1, 2, na);
    avec_low = (avec_tmp(2:end) + avec_tmp(1:end-1)) / 2;
    avec_low = avec_low(avec_low <= 1.5);
    avec = sort([avec_tmp, avec_low]);
end

function needsRecompute = checkIfRecomputeNeeded(f, lambda, tolaub)
    % CHECKIFRECOMPUTENEEDED Determine if solution needs recomputation
    %
    % Inputs:
    %   f       - Loaded file structure
    %   lambda  - Current lambda value
    %   tolaub  - Tolerance for aub convergence
    %
    % Outputs:
    %   needsRecompute - Boolean indicating if recomputation needed
    
    % Check if lambda matches
    lambdaMismatch = ~strcmp(...
        num2str(round(lambda, 3)*1e3, '%.0f'), ...
        num2str(round(f.par.lambda, 3)*1e3, '%.0f'));
    
    % Check if aub converged properly
    aubConvergenceIssue = false;
    if isfield(f, 'aout') && isfield(f, 'aub0')
        aubConvergenceIssue = abs(f.aout(end) - f.aub0) / f.aub0 > tolaub;
    end
    if isfield(f, 'aub') && isfield(f, 'aub0')
        aubConvergenceIssue = aubConvergenceIssue || ...
            abs(f.aub - f.aub0) / abs(f.aub0) > tolaub;
    end
    
    needsRecompute = lambdaMismatch || aubConvergenceIssue;
end
