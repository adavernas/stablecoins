clearvars
set(0,'DefaultFigureWindowStyle','docked')

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

% ell0 = 0.20;
% ell1 = 0.10;

% ell  = @(a) max(0,ell0 - ell1./a);
% ellp  = @(a) ell1./a.^2;
% amax  = 2*ell1/ell0;

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
options = optimset('tolF',1e-20, 'tolX',1e-20);

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
disp(sum(abs(sol2.F)))

figure(2); clf(2);
sgtitle(['jump no collateral $\underline{\delta}=$',num2str(sol2.delt_)],'Interpreter','LateX') 
plotfig(sol2,par)

% figure(3); clf(3);
% avec = linspace(0,1.5*sol2.astar,1000);
% plot(avec,sol2.delta(avec),'k'); hold on;
% plot(sol2.aub,sol2.delta(sol2.aub),'*k'); hold on;


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
disp(sum(abs(sol3.F)))
disp(['delta: ',num2str(sol3.delt_)])

figure(3); clf(3);
sgtitle(['jump no collateral $\underline{\delta}^\star=$',num2str(sol3.delt_)],'Interpreter','LateX') 
plotfig(sol3,par)

[x4,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump no collateral delta infinite'), [amax;amax/2],...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump no collateral delta infinite'), options);
[~,sol4] = Fsyst(x4,par,...
    'case','jump no collateral delta infinite');
disp(sum(abs(sol4.F)))

Fmax = sol4.obj;

figure(4); clf(4);
sgtitle('jump no collateral $\underline{\delta}=\infty$','Interpreter','LateX') 
plotfig(sol4,par)

x0 = [amax;amax/2;0];

x5 = fsolve(@(x) Fsyst(x,par,'case','jump commit'), x0, options);
[~,sol5] = Fsyst(x5,par,'case','jump commit');
disp(sum(abs(sol5.F)))

figure(5); clf(5);
sgtitle('jump commit','Interpreter','LateX') 
plotfig(sol5,par)

figure(6); clf(6);
plotgraphs(par,'sol1',sol5,'sol2',sol4,'amin',0,'maxa',1.2*max(sol5.astar,sol2.astar),...
    'graphcase','uncollateralized','subcase','slides','name','uncollateralized',...
    'savegraph','off')

figure(7); clf(7);
plotgraphs(par,'sol1',sol5,'sol2',sol2,'amin',0,'maxa',1.2*max(sol5.astar,sol2.astar),...
    'graphcase','uncollateralized','subcase','slides','name','uncollateralized_delta_low',...
    'savegraph','off')