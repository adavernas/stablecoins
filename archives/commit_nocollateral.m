clearvars

%% PARAMETERS
r   =  0.06;
xmu =  0.05;
sig =  0.1;

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

xi = 1;
lambda = 0.01;
muk = 0.055;

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

%% END
writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20, 'display', 'on');

% astar = fsolve(@(a) ffastar(a,0,par),amax);

x0 = [amax;amax;0]*1e6;

[x0,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit'), [amax;amax/2;0],...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit'), options);
[~,sol0] = Fsyst(x0,par,...
    'case','jump commit');
disp(sum(abs(sol0.F)))
 
figure(1); clf(1);
sgtitle('jump commit','Interpreter','LateX') 
plotfig(sol0,par)

par.lambda = 0.02;

[x0,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit'), x0,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit'), options);
[~,sol0] = Fsyst(x0,par,...
    'case','jump commit');
disp(sum(abs(sol0.F)))

figure(2); clf(2);
sgtitle('jump commit','Interpreter','LateX') 
plotfig(sol0,par)

% x0 = [sol0.aub;
%       sol0.deltastar];

% fastar = fminunc( @(astar) - fvalue(astar,x0,par,options,...
%    'case','jump commit exo','delta','in stablecoins'), amax);

% fvalue(sol0.astar,x0,par,options,'case','jump commit exo','delta','in stablecoins')

% ni = 10;
% astarvec = linspace(0.8*fastar,1.2*fastar,ni);
% fstarvec = NaN(size(astarvec));
% sol1 = sol0;
% for ii=1:ni
%     x0 = [sol1.aub;
%           sol1.deltastar];
%      
%     x1 = fsolve(@(x) Fsyst(x,par,'case','jump commit exo','delta','in stablecoins',...
%         'xvar',astarvec(ii)), x0, options);
%     
%     [~,sol1] = Fsyst(x1,par,'case','jump commit exo','delta','in stablecoins',...
%         'xvar',astarvec(ii));
% 
%     fstarvec(ii) = sol1.fstar;
%     disp(sum(abs(sol1.F)))
% end
% 
% figure(2); clf(2); hold on;
% plot(astarvec,fstarvec)
% plot(sol0.astar,interp1(astarvec,fstarvec,sol0.astar),'ok')
% plot(fastar,interp1(astarvec,fstarvec,fastar),'*r')

x0 = [sol0.aub;
      sol0.astar];

x2 = fsolve(@(x) Fsyst(x,par,'case','jump commit no limited'), x0, options);
[~,sol2] = Fsyst(x2,par,'case','jump commit no limited');
disp(sum(abs(sol2.F)))

figure(3); clf(3);
sgtitle('jump full commit','Interpreter','LateX') 
plotfig(sol2,par)

figure(4); clf(4);
plotgraphs(par,'sol1',sol0,'sol2',sol2,'amin',0,...
    'maxa',max(2*sol2.astar,1.2*sol0.astar),'graphcase','commit',...
    'name','commit','savegraph','off')

% [x3,~] = fmincon(@(x) -FsystCon(x,par,...
%     'case','jump commit fullcollateral'), amax,...
%     [],[],[],[],[],[],@(x) mycon(x,par,...
%     'case','jump commit fullcollateral'), options);
% [~,sol3] = Fsyst(x3,par,...
%     'case','jump commit fullcollateral');
% disp(sum(abs(sol3.F)))

% ff = @(a) ellp(a).*a - (muk-r+ell(a));
% 
% avec = linspace(0,3,100);
% figure(5); clf(5);
% plot(avec,ff(avec));
% 
% [x4,F] = fsolve(@(x) ellp(x)*x - (muk-r+ell(x)), amax, options);





