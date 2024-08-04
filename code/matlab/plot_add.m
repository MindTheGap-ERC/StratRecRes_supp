M = 10;
S = 0.05;
L = 10;
t = linspace(0, 500, 100);
d = linspace(0, 1.1 * L, 100);

epsilon = 0.1;
init_distr = @(x) 2/(epsilon*sqrt(2*pi))*exp(-0.5*(x/epsilon).^2);
odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-6); 

%%
plot(d, init_distr(d))

%%
ADD = get_ADD(t,d, S, M, L, init_distr, odeoptions);

%%
contourf( t,d,  ADD', -0.01:0.01:1, 'ShowText','on')
set(gca, 'YDir','reverse')
set(gca, "XAxisLocation","top")
set(gca,'ColorScale','log')
xlabel("Age (years)")
ylabel("Depth (cm)")



%%
function ADD=get_ADD(t, d, S, M, L, init_distr, odeoptions)
% see https://www.mathworks.com/help/matlab/ref/pdepe.html

function [pL,qL,pR,qR] = boundaryconditions(~,~,~,uR,~) % left = surface, right = bottom
% surface: no flux
pL = 0;
qL = 1;
% bottom : no diffusive flux
pR = uR;
qR = 1;
end

function [c,f,s] = pdefun(x,~,u,dudx)
c = 1; % symmetry constant
f =  M * double(x <= L) *dudx -S * u; % flux term
s = 0; % source term
end

ADD = pdepe(0,@pdefun,init_distr,@boundaryconditions,d,t,odeoptions);
end