% Collects all calibration parameters.
% All parameters are taken from Krusell and Smith (1998,JPE).
% All the variables are packed in "Params".
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 09/11/2019  T. Yamada   Original code

clear;
clear global;
close all;
format short;

global Params

%% ***** CALIBRATION: PREFERENCE *****
Params.beta = 0.99; % discount factor
Params.gamma = 1.0; % relative risk aversion
%=====================================

% ***** CALIBRATION: PRODUCTION *****
Params.alpha = 0.36;  % capital share
Params.delta = 0.025; % depreciation rate
%====================================

% ***** CALIBRATION: IDIOSYNCRATIC UNEMPLOYMENT SHOCK *****
Params.ne = 2;              % employed, unemployed
Params.endow = [1.0, 0.05]; % labor endowment: UI = 0.05
%==========================================================

% ***** AGGREGATE PRODUCTIVITY SHOCK AND UNEMPLOYMENT DYNAMICS *****
Params.nz = 2;             % #aggregate state = 2
Params.unempg = 0.04;      % unemployment rate when good
Params.unempb = 0.1;       % unemployment rate when bad
Params.durug = 1.5;        % duration of unemployment when good
Params.durub = 2.5;        % duration of unemployment when bad
Params.durgd = 8.0;        % duration of boom (8 quarters)
Params.durbd = 8.0;        % duration of recession (8 quarters)
Params.tfp = [1.01, 0.99]; % TFP shock
%===================================================================

% ***** GRIDS AND SUMULATION *****
Params.na = 21;      % #grid for VFI/EGM
Params.nd = 501;     % #grid for policy function
Params.nk = 6;       % #grid for aggregate capital
Params.amax = 300.0; % max of asset holding
Params.amin = 0.0;   % min of asset holding
Params.kmax = 40.0;  % max of aggregate capital grid
Params.kmin = 30.0;  % min of aggregate capital grid
Params.nums = 6000;  % simulation period
Params.numi = 1000;  % #individuals in the simulation
%=================================

% %% ***** VARIABLES IN THE MAIN LOOP *****
% Params.policy = zeros(Params.na, Params.ne, Params.nk, Params.nz);
% Params.k_path = zeros(Params.nums);
% Params.z_path = zeros(Params.nums);

% % ***** INITIAL GUESS *****
% Params.icept  = zeros(Params.nz, 1);
% Params.icept1 = zeros(Params.nz, 1);
% Params.icept2 = zeros(Params.nz, 1);
% Params.slope  = ones(Params.nz, 1);
% Params.slope1 = ones(Params.nz, 1);
% Params.slope2 = ones(Params.nz, 1);
% Params.R2     = zeros(Params.nz, 1);

% ***** TOLERANCE OF ERROR *****
it = 1;
maxit = 1000;
metric1 = 1.0;
metric2 = 1.0;
metric = 1.0;
adj = 0.5;
toler = 1.0e-004;
%===============================

%% ***** GRID *****
% Params.aprime = linspace(Params.amin, Params.amax, Params.na);
% Params.grid = linspace(Params.amin, Params.amax, Params.nd);
Params.aprime = grid_exp3(Params.amin, Params.amax, Params.na);
Params.grid = grid_exp3(Params.amin, Params.amax, Params.nd);
Params.kgrid = linspace(Params.kmin, Params.kmax, Params.nk);
%==================

%% ***** IDIOSYNCRATIC UNEMPLOYMENT RISK *****
[Params.prob, Params.tran_zz] = trans_prob(Params);
%=======================================

%% ***** AGGREGATE LABOR *****
Params.pi_dist(:, 1) = [1.0-Params.unempg, Params.unempg];
Params.pi_dist(:, 2) = [1.0-Params.unempb, Params.unempb];
Params.L_agg = zeros(2,1);
for z = 1: Params.nz
    Params.L_agg(z)  = Params.endow*Params.pi_dist(:, z);
end
%=============================

%% SAVE VARIABLES AS .mat FORMAT
save calibration.mat;

return;
