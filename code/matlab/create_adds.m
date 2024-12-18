%% define filepath for data storage
path_parts = split(mfilename('fullpath'), filesep);
relevant_path_parts = path_parts(1:end-3);
new_path = join(relevant_path_parts, filesep);
data_path = join([new_path{1},"data","matlab_outputs"], filesep);
filename = "ADDs_from_matlab.mat";
file_path = strcat(data_path, filesep, filename);


% set dimensionless times and depths
t_dimless = sort(unique([0:10^-6:0.001,1.6:10^-3:9,0.001:10^-4:1.6])); %times of evaluation
depth_subdivisions = 200;
d_dimless = linspace(0,1.1,depth_subdivisions); %depths of evaluation
odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-9); %settings for solution of pde
peclet_numbers = 10.^(-5:0.25:3); % relation mixing /transport

% initialize storage
tavg_below_sml = zeros(length(t_dimless), length(peclet_numbers));

% initial condition
epsilon = 0.005;
init_distr = @(x) 2/(epsilon*sqrt(2*pi))*exp(-0.5*(x/epsilon).^2);


%% calculate ADDs
disp("Starting calculation...")
for i=1:length(peclet_numbers)
    peclet_no = peclet_numbers(i);
    disp(strcat("Calculating ADD for Peclet number ", num2str(peclet_no)))
    ADD = get_ADD_dimless(t_dimless, d_dimless, peclet_no, init_distr, odeoptions);
    tavg_below_sml(:, i) = ADD(:, find(d_dimless>1, 1));
    %ADD_res(:,:,i) = ADD;
end




%% Simulate ADD for Po4 core
ages = 0:200; % age in years
depths = 0:125; % depth in cm

m_depths = [0,17, 18,120,125]; % depths where mixing values change
m_vals = [780, 780, 20, 20, 0]; % mixing values at depths in m_depths

f = @(x) interp1(m_depths, m_vals, x);
S = 2.2; % sedimentation rate (cm/y) for Po4 core based on tomasovych et al 2018

%%
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-6); %settings for solution of pde

epsilon = 1;
init_distr = @(x) 2/(epsilon*sqrt(2*pi))*exp(-0.5*(x/epsilon).^2);

%%
S_fun = @(d) S; %sedimentation rate
M_fun = @(d) interp1(m_depths, m_vals, d); %mixing
L_fun = @(d) 0 ; % no disintegration
fZero = init_distr;
surfaceInflux =  0;
surfaceLossRate =  0;

%%

u=PartiMoDe_StEn(ages, depths, S_fun, M_fun, L_fun, fZero, surfaceInflux, surfaceLossRate, odeoptions);

%%
%plot(depths, u(:,60))

%%
%contourf(ages, depths, u, "LevelStep",0.001)

%%
save(file_path, ...
    "tavg_below_sml", ...
    "t_dimless", ...
    "epsilon", ...
    "odeoptions", ...
    "peclet_numbers", ...
    "ages", ...
    "depths", ...
    "u")

disp(strcat("Done. Output saved into data/matlab_outputs/", filename))
%%
function u=PartiMoDe_StEn(times, depths, S, M, L, fZero, surfaceInflux, surfaceLossRate, odeOptions)
%% Description
% Implementation of the PartiMode model where particle movement and 
% destruction depends on location, but not on time 
% Author: Niklas Hohmann
% email: N.H.Hohmann (at) uu.nl , ORCID: https://orcid.org/0000-0003-1559-1838
%% Inputs:
% times: Vector with strictly increasing positive numbers. The points in
%       time where the model outputs are determined
% depths: vector of strictly increasing numbers with depths(1)=0. The
%       depths/locations where the model outputs are determined
% S:    Function handle. S=S(x) is a function of one variable that  determines
%       the advective flux at depth x
% M:    Function handle. M=M(x) is a function of one variable that determines
%       the diffusive flux at depth x
% L:    Function handle. L=L(x) is a function of one variable that determines the disintegration rate of
%       particles at depth x
% fZero: Function handle. fZero=fZero(x) is a function of one variable that describes the initial
%       condition, e.g. the state of the system at t=0
% surfaceInflux: a scalar >= 0 describing the influx of new 
%       particles through the sedimen surface per time unit 
% surfaceLossRate: a scalar >= 0 describing the rate with which with which
%       particles are lost through the sediment surface
% odeOption: options structure to be handed over to the ODE solver ode15s.
%       For details type
%       doc odeset
%       into the command line and/or see
%       https://www.mathworks.com/help/matlab/math/summary-of-ode-options.html

%% Outputs
% u: Matrix with length(depths) rows and length(times) colums. u(i,j) is
%       the value of u(x,t) at depth x=depths(i) and time t=times(j)

%% See also
% * File "Example_PartiMode_StEn.m" for an example
% * Function PartiMoDe_LoFi: particle input changes with time
% * Function PartiMoDe: effects are depth and time-dependent
% * Function PartiMoDe_Interpol: interpolation of results
% * Function StratProfile: determines particle density as a function of
%       depth

%% Define Flux and Source/Sink Terms
% See Matlab help page of the pdepe function for details on how to code
%       partial differential equations:
%       https://www.mathworks.com/help/matlab/ref/pdepe.html
function [c,f,s] = pdefun(x,~,u,dudx)
% coefficient on left side (in front of d/dt u(x,t)) is always 1
c = 1;
% flux = diffusive flux + advective flux  (negative sign because of
%       downcore advection)
f =  M(x) * dudx - S(x) * u;
% Sink term: loss of particles because of the disintegration rate
s = - L(x)*u;
end

%% Define Boundary Conditions
function [pL,qL,pR,qR] = boundaryconditions(~,uL,xR,uR,~)
% See Matlab help page of the pdepe function for details on how to code
%       partial differential equations
%       https://www.mathworks.com/help/matlab/ref/pdepe.html
%       Here left ("L") is the sediment surface, and R is the bottom of the
%       observed interval

% Standard form for boundary conditions at the surface is
%       pL + qL*flux=0 
%       Particle flux through the sediment surface is determined by the
%       surface influx and the surface loss rate
%       flux = surface influx - surface loss rate * u
pL = surfaceInflux - surfaceLossRate*uL;
qL = 1;

% Standard form for boundary conditions at bottom is
%       pR + qR*flux=0 
% No diffusive flux of particles at the bottom. Since
%       flux = flux_adv + flux_diff
%       we get
pR = S(xR)*uR;
qR = 1;
end

%% Solve PDE
% Hand functions over to pdepe to determine u(x,t) at the specified
%       times and depths
u = pdepe(0,@pdefun,fZero,@boundaryconditions,depths,times,odeOptions);
u=u';
end

%% Function to determine ADD
function ADD=get_ADD_dimless(t_dimless, d_dimless, peclet_no ,init_distr, odeoptions)
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
f =  (1/peclet_no) * double(x <= 1) *dudx -u; % flux term
s = 0; % source term
end

ADD = pdepe(0,@pdefun,init_distr,@boundaryconditions,d_dimless,t_dimless,odeoptions);
end


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

function [v, t] = cumul_int(d, t, normalize)
% determines cumulative integral over a density using the trapezoidal rule
% inputs
% d densities
% t times at which densities are observed. vector of the same length as 
t_diff = diff(t);
d_mean = 0.5 * (d(1:end - 1) + d(2:end) );
vals = cumsum([0, t_diff .* d_mean']);
if normalize
    vals = vals/max(vals);
end
v = vals;
end

