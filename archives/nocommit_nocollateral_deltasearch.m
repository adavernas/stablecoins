clearvars
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:nearlySingularMatrix')

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

% ell0 = 0.20;
% ell1 = 0.10;

% ell  = @(a) max(0,ell0 - ell1./a);
% ellp  = @(a) ell1./a.^2;
% amax  = 2*ell1/ell0;

savegraph = 'off';

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 6;
lambda = 0.1;

delt0 = (r + sig^2/2 - exp(-1)*r)*1.01;
varphi0 = 0.02;

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

%% END
writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20,'tolX',1e-20,'display','off');

x0 = amax;

x = fsolve(@(x) Fsyst(x,par,'case','stochastic no collateral'), x0, options);
[~,sol1] = Fsyst(x,par,'case','stochastic no collateral');
disp(sum(abs(sol1.FF)))

figure(1); clf(1);
sgtitle('stochastic no collateral','Interpreter','LateX')
plotfig(sol1,par)

x0 = [sol1.astar;
    sol1.aub;
    sol1.ck2;
    sol1.estar];

par.delt0 = par.xmu;

[x2,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump no collateral'), x0,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump no collateral'), options);
[~,sol2] = Fsyst(x2,par,...
    'case','jump no collateral');

figure(2); clf(2);
sgtitle(['jump no collateral $\underline{\delta}=$',num2str(sol2.delt_)],'Interpreter','LateX')
plotfig(sol2,par)

x0 = [sol2.astar;
    sol2.aub;
    sol2.ck3;
    sol2.estar;
    sol2.delt_];

[x3,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump no collateral',...
    'subcase','delta'), x0,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump no collateral',...
    'subcase','delta'), options);
[~,sol3] = Fsyst(x3,par,...
    'case','jump no collateral',...
    'subcase','delta');

figure(3); clf(3);
sgtitle(['jump no collateral $\underline{\delta}^\star=$',num2str(sol3.delt_)],'Interpreter','LateX')
plotfig(sol3,par)

x0 = [sol3.astar;sol3.aub];
[x4,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump no collateral delta infinite'), x0,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump no collateral delta infinite'), options);
[~,sol4] = Fsyst(x4,par,...
    'case','jump no collateral delta infinite');

figure(4); clf(4);
sgtitle('jump no collateral $\underline{\delta}=\infty$','Interpreter','LateX')
plotfig(sol4,par)

for ii=1:1e4
    par.r = rand;
    par.sig = rand*10;
    par.xi = rand*10;
    par.lambda = rand;
    par.xmu = -rand + rand*((par.r+par.lambda/(par.xi+1)));

    par.r = 0.29617;
    par.sig = 4.6278;
    par.xi = 9.2523;
    par.lambda = 0.21589;
    par.xmu = 0.28659;
    
    if nocommit_nocollateral_existence(par)>0
                
        [xa,~] = fmincon(@(x) -FsystCon(x,par,...
            'case','jump no collateral delta infinite'), x4,...
            [],[],[],[],[],[],@(x) mycon(x,par,...
            'case','jump no collateral delta infinite'), options);
        [~,sola] = Fsyst(xa,par,...
            'case','jump no collateral delta infinite');
        
        figure(5); clf(5);
        sgtitle('jump no collateral $\underline{\delta}=\infty$','Interpreter','LateX')
        plotfig(sola,par)

        [xb,~] = fmincon(@(x) -FsystCon(x,par,...
            'case','jump no collateral',...
            'subcase','delta'), x3,...
            [],[],[],[],[],[],@(x) mycon(x,par,...
            'case','jump no collateral',...
            'subcase','delta'), options);
        [~,solb] = Fsyst(xb,par,...
            'case','jump no collateral',...
            'subcase','delta');
        
        figure(6); clf(6);
        sgtitle(['jump no collateral $\underline{\delta}^\star=$',num2str(solb.delt_)],'Interpreter','LateX')
        plotfig(solb,par)
               
        if and(and(and(solb.fstar>sola.fstar,sola.estar>0),solb.estar>0),sum(abs(solb.F))<1e-4)
            disp(sum(abs(sola.F)))
            disp(sum(abs(solb.F)))
            disp(['delta infinite not optimal at',...
                ' r = ',num2str(par.r),...
                ' sigma = ',num2str(par.sig),...
                ' xi = ',num2str(par.xi),...
                ' lambda = ',num2str(par.lambda),...
                ' mu = ',num2str(par.xmu)])
        end
    end
end
