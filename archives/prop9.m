clearvars
set(0,'DefaultFigureWindowStyle','docked')

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;

xi = 1;
lambda = 0.010;

gk3 = NaN;

%% END

writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20, 'display','on');

ng = 40;
gk3vec = linspace(-1.1,-10,ng);

vars = {'dLdg','FF','phi','Chat','Theta'};
nv = length(vars);
for iv=1:nv
    eval([vars{iv},'vec = NaN(1,ng);'])
end

x0 = [0.5;0.5;1];
for ig=1:ng
    par.gk3 = gk3vec(ig);
    [x,~] = fsolve(@(x) systEqn(x,par,'case','phi>0'), x0, options);
    [F,sol] = systEqn(x,par,'case','phi>0'); %#ok<ASGLU>
    
    if sum(abs(F))<1e-8
        for iv=1:nv
            eval([vars{iv},'vec(ig) = sol.',vars{iv},';'])
        end
        x0 = x;
    end
end

x0 = [0.5;0.5;1;-2];
[x,~] = fsolve(@(x) systEqn(x,par,'case','phi>0','subcase','phi>0,gk3'), x0, options);
[F,sol] = systEqn(x,par,'case','phi>0','subcase','phi>0,gk3');
    
nl = 2;
nc = 3;
figure(1); clf(1);

hf = 1;
subplot(nl,nc,hf); hold on;
plot(gk3vec,dLdgvec)
plot(sol.gk3,sol.dLdg,'*r')
xlabel('$\gamma$','Interpreter','LateX')
title('$\partial L/\partial \gamma$','Interpreter','LateX')

hf = hf+1;
subplot(nl,nc,hf); hold on;
plot(gk3vec,FFvec)
plot(sol.gk3,sol.FF,'*r')
xlabel('$\gamma$','Interpreter','LateX')
title('$FF$','Interpreter','LateX')

hf = hf+1;
subplot(nl,nc,hf); hold on;
plot(gk3vec,phivec)
plot(sol.gk3,sol.phi,'*r')
xlabel('$\gamma$','Interpreter','LateX')
title('$\phi$','Interpreter','LateX')

hf = hf+1;
subplot(nl,nc,hf); hold on;
plot(gk3vec,Chatvec)
plot(sol.gk3,sol.Chat,'*r')
xlabel('$\gamma$','Interpreter','LateX')
title('$\widehat{C}$','Interpreter','LateX')

hf = hf+1;
subplot(nl,nc,hf); hold on;
plot(gk3vec,Thetavec)
plot(sol.gk3,sol.Theta,'*r')
xlabel('$\gamma$','Interpreter','LateX')
title('$\Theta$','Interpreter','LateX')

disp(['phi = ',num2str(sol.phi),...
    ' Theta = ',num2str(sol.Theta),...
    ' Ctilde = ',num2str(sol.Ctilde),...
    ' FF = ',num2str(sol.FF)])