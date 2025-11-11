clearvars
set(0,'DefaultFigureWindowStyle','docked')
dbstop if error

%% PARAMETERS
r   =  0.06;
xmu  =  0.05;
sig =  0.10;

% ell0 = 0.20;
% ell1 = 0.10;

% ell  = @(a) max(0,ell0 - ell1./a);
% ellp  = @(a) ell1./a.^2;
% amax  = 2*ell1/ell0;

delt0 = 0;
ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 6;
lambda = 0.10;

varphi0 = 1;

muk = 0.05;

aubcon = NaN;
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
  
[x0,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit fullcollateral'), amax,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit fullcollateral'), options);
[~,sol0] = Fsyst(x0,par,...
    'case','jump commit fullcollateral');
disp(sum(abs(sol0.F)))

figure(1); clf(1);
sgtitle(['jump collateral commit $\varphi^\star=',num2str(sol0.varphi0),'$'],'Interpreter','LateX') 
plotfig(sol0,par)

% ni = 30;
% vvec = linspace(0,1,ni);
% svec = nan(ni,1);
% fvec = nan(ni,1);
% x0 = sol1.x;
% for ii=1:ni
%     par.varphi0 = vvec(ii);
%     
%     [x2,~] = fmincon(@(x) -FsystCon(x,par,...
%                     'case','jump collateral commit',...
%                     'subcase','constraint varphi'), x0,[],[],[],[],[],[],...
%                     @(x) mycon(x,par,...
%                     'case','jump collateral commit',...
%                     'subcase','constraint varphi'), options);
%     [~,sol2] = Fsyst(x2,par,...
%                     'case','jump collateral commit',...
%                     'subcase','constraint varphi');
%     disp(sum(abs(sol2.F)))
% 
% 
%     if and(sol2.fstar>1e-10,sum(abs(sol2.F))<1e-10)
%         fvec(ii) = sol2.fstar;
%         svec(ii) = sol2.astar;
%     end
%  end
% 
% figure(2); clf(2);
% plotgraphs(par,'xvec',vvec,'fvec',fvec','svec',svec,... 
%                'graphcase','collateralized varphi',...
%                'name','collateralized_varphi','savegraph',par.savegraph)
%            
% figure(3); clf(3);
% sgtitle(['jump collateral commit $\varphi=',num2str(sol2.varphi),'$'],'Interpreter','LateX') 
% plotfig(sol2,par,sol2.alb,max(1.1*sol2.aub,sol2.astar*1.1),'color','b')
 
x0 = [sol0.aub;
      sol0.astar];

x4 = fsolve(@(x) Fsyst(x,par,'case','jump commit no limited'), x0, options);
[~,sol4] = Fsyst(x4,par,'case','jump commit no limited');
disp(sum(abs(sol4.F)))

figure(4); clf(4);
sgtitle('jump full commit','Interpreter','LateX') 
plotfig(sol4,par)

par.varphi0 = 1;
[x5,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump collateral aub 0'), amax,[],[],[],[],[],[],...
    @(x) mycon(x,par,...
     'case','jump collateral aub 0'), options);
[~,sol5] = Fsyst(x5,par,...
    'case','jump collateral aub 0');
disp(sum(abs(sol5.F)))

figure(5); clf(5);
sgtitle('jump full commit','Interpreter','LateX') 
plotfig(sol5,par)

figure(6); clf(6);
plotgraphs(par,'sol1',sol5,'sol2',sol4,'amin',0,...
    'maxa',max(2*sol4.astar,1.2*sol5.astar),'graphcase','commit collateral',...
    'name','commit_collateral','savegraph','off')


