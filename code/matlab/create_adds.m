%% define filepath for data storage
path_parts = split(mfilename('fullpath'), filesep);
relevant_path_parts = path_parts(1:end-3);
new_path = join(relevant_path_parts, filesep);
data_path = join([new_path{1},"data","matlab_outputs"], filesep);
filename = "tavg_from_matlab.mat";
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

save(file_path, ...
    "tavg_below_sml", ...
    "t_dimless", ...
    "epsilon", ...
    "odeoptions", ...
    "peclet_numbers")
disp(strcat("Done. Output saved into data/matlab_outputs/", filename))

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

