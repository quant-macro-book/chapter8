% Purpose:
%  Solve Krusell and Smith (1998) model.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 09/11/2019  T. Yamada   Original code

clear;
clear global;
close all;
format short;

load calibration.mat;
global Params

% ***** MAIN VARIABLES IN THE LOOP ******
policy = zeros(Params.nd, Params.ne, Params.nk, Params.nz);
k_path = zeros(Params.nums, 1);
z_path = zeros(Params.nums, 1);
%========================================

% ***** INITIAL GUESS *****
icept  = zeros(Params.nz, 1);
icept1 = zeros(Params.nz, 1);
icept2 = zeros(Params.nz, 1);
slope  = ones(Params.nz, 1);
slope1 = ones(Params.nz, 1);
slope2 = ones(Params.nz, 1);
R2     = zeros(Params.nz, 1);
%==========================

% good initial guess
icept(1) = 0.1391;
icept(2) = 0.1297;
slope(1) = 0.9617;
slope(2) = 0.9630;

disp(' ');
disp('-+-+-+- Solving Krusell and Smith model -+-+-+-');

tic

%% MAIN LOOP TO FIND AN FIXED POINT
while (it<=maxit && metric>toler)

    fprintf('----- main iteration: %d ----- \n' , it);
    disp(' ')

    % compute policy function using endogenous grid method
    [policy] = end_grid_method(icept, slope);
toc

% for debug use
%save tmp_policy.mat;
%load tmp_policy.mat;

    % law of motion by simulation
    [k_path, z_path] = law_of_motion_sim(policy);
toc

% for debug use
%save tmp_simulation.mat;
%load tmp_simulation.mat;

    % new coefficients from regression
    [icept1, slope1, R2] = regress_KS(k_path, z_path);

    % iteration error
    metric1 = max(abs((icept-icept1)./icept));
	metric2 = max(abs((slope-slope1)./slope));
    metric = max(metric1, metric2);
    it = it + 1;
    
    % update coefficients
    icept2 = adj*icept + (1.0-adj)*icept1;
    slope2 = adj*slope + (1.0-adj)*slope1;
    icept  = icept2;
    slope  = slope2;

    disp('error (%):')
    disp(metric*100.0)
    disp('intercept:')
    disp(icept)
    disp('slope:')
    disp(slope)
    disp('R^2:')
    disp(R2)
    
end

toc

%% ===================================%;
%            SAVE RESULTS             %;
%=====================================%;

save final_results.mat;

return
