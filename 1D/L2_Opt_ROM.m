%% 2D steady Stokes equation
% min \lambda^2 ||yu - \hat yu|| + ||yp - y p||
%
%% Main Code
% parameters
format short e

clc;clear all;

addpath(genpath(strcat(pwd,'/Example')));

load('train_data.mat');
load('Model_FOM.mat');
load('Opts.mat');

N = 1:1:5;
Opts.para.lambda = 1;

tic;

train_data = Modified_data(train_data, Opts);
Opts.train_data = train_data;
Opts.Int.point = Opts.Xi.xi_train;
Opts.blk = Model_FOM.blk;

Opts.L2_Opts.maxit = 10000;
Opts.L2_Opts.tol = 1e-6;
Opts.L2_Opts.Armijo.rho = 0.5;
Opts.L2_Opts.Armijo.sigma = 1e-5;
Opts.L2_Opts.init_step_size = 1;
Opts.L2_Opts.m = 3;

fprintf('\n--------------------------------------------------------------');
fprintf('\nL2-Opt-ROM');
fprintf('\n--------------------------------------------------------------');

for i_N = 1:length(N)
    fprintf('\nStart: N = %d', N(i_N));
    load(strcat('POD_result_', num2str(N(i_N)), '.mat'));

    Opts.para.N = N(i_N);
    Opts.para.rp = N(i_N);
    Opts.para.ru = N(i_N);
    Opts.para.rus = 2 * Opts.para.ru;
    Opts.initial_point = Modified_Model(result.Model, Opts);

    Model_L2_Opt_ROM_POD = L2_Opt_ROM_Solver(Opts);

    Model_L2_Opt_ROM_POD.base = result.Model.base;

    data_L2_Opt_ROM_POD = compute_xy_from_ROM_Model_ROM(Model_L2_Opt_ROM_POD, Opts, Opts.Xi.xi_true);
end
