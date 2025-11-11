clearvars

%% PARAMETERS
r   =  0.06;
xmu =  0.05;
sig =  0.1;

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 1;
lambda = 0.02;
muk = 0.055;

varphi0 = 1;

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

[x1,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit fullcollateral'), amax,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit fullcollateral'), options);
[~,sol1] = Fsyst(x1,par,...
    'case','jump commit fullcollateral');
disp(sum(abs(sol1.F)))

figure(2); clf(2);
sgtitle(['jump collateral commit $\varphi^\star=',num2str(sol1.varphi0),'$'],'Interpreter','LateX') 
plotfig(sol1,par)

figure(3); clf(3);
plotgraphs(par,'sol1',sol0,'sol2',sol1,'amin',0,...
    'maxa',max(2*sol1.astar,1.2*sol1.astar),'graphcase','vincent',...
    'name','commit_collateral','savegraph','off')


