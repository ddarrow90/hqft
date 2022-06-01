clear all
clc
close all

mkdir('results')

add_to_path

%% Set up ICs
% Fill-in uniform distribution of impact parameters (positive values)
n=36;       % number of trajectories (multiples of 8 best)
maxIP = 13;
deltaIP = 0.25;

xIP = 0:deltaIP:maxIP     %Array of impact parameters

% Distribute these on a ciruclar domain
Rdomain = 22;
yIP = -sqrt(Rdomain.^2-xIP.^2);

% Initial velocity
ui = 0;
vi = 0.039;

% Convert accoriding to \lamda_shallow
% If nothing, the IP is referred to \lambda_deep
convert_ICs_to_lshallow = true;

%% Run
% Get the current pool and shut it down (in case it was open)
poolobj = gcp('nocreate');delete(poolobj)

% Start new parallel pool
parpool(5) % Set number of cores

tic
parfor i=1:size(xIP,2) 
    xi = xIP(i);
    yi = yIP(i);
    hole(xi,yi,ui,vi,convert_ICs_to_lshallow);
end

%reduce_filesize

poolobj = gcp('nocreate');delete(poolobj)
toc