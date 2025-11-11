function par = parfun

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
varphi0 = 1;

% Parameter validation
if r - xmu + lambda/(xi+1) < 0
    warning('ParameterCheck:NegativePhi', 'phi is negative! This may cause numerical issues.');
end

if r + xmu - sig^2 - lambda/(xi-1) < 0
    warning('ParameterCheck:NegativeKappa', 'kappa is negative! This may cause numerical issues.');
end

% Configuration
loadsolution = 'on';
savegraph = 'off';
display = 'off';

% Parallel computation option
use_parallel = true;  % Set to false to disable parallel processing

par.names = who;
for i = 1:length(par.names)
    eval(['par.',par.names{i},' = ',par.names{i},';'])
end
