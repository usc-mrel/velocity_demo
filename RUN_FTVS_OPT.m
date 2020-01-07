%RUN_VS_OPT designs FT-VS pulse with optimized at a user defined range
% of B0 and B1.
%
%   Author: Vanessa, adapted from Terrence Jao
%   Date: ...
%   References:

%% Pulse Timing Calculation
clear all; close all; clc;
addpath('useful'); addpath('rf_tools'); addpath(genpath('bloch_mex'));
addpath('rf_tools/mex');

% -------------------------------------------------------------------------
% Do we want to plot?
ploton = 0;
% Do we want to display messages?
verbose = 1;        %Turn off verbose for slightly faster runtime

B1_val_hp = 1.4; %1.72; %2; %10.20; % uT
B1_val_inv = 10.2; %12; % uT
Grad_val = 1.45; % mT/cm

% durations (mS)
grad_ramp = 0.302;
Tgap = 0; %0.1; % primarily for eddy current artifacts

% sampling periods
dtGz = 0.002;

comp_hp = 0;
comp_inv = 0; % 90-180-90, 90-240-90, 90-360-90
single_refocus = 0; %0,1
grad_var = 3;
sinc_weight = 1;
num_sp = 9; % 5,9

%rf_weight = dzrf(num_sp,8,'inv','max',0.1,0.001);
rf_weight = dzrf(num_sp,6,'inv','pm',0.1,0.01);
hpscale = real(num_sp/sum(rf_weight)*rf_weight);
init_vals = real(hpscale);

%% Optimizaiton
lambda = [0.001 0.01 0.1 1 10 100 1000]; %10.^(-3:2); %[0 logspace(0,6)];
%lambda = [1 1]; %[10.^(-3:3) 0.5 0.8 0.3]; %1/4:3);
est_vec = zeros(length(lambda),length(hpscale)); fval_vec = est_vec;

for itr=3:length(lambda)-1

model = @(x) opt_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp, itr,lambda(itr), x);

options = optimoptions(@fmincon,'MaxIterations',1000, ...
    'MaxFunctionEvaluations',10000, ...
    'OptimalityTolerance',1e-9, ...
    'ConstraintTolerance',1e-9,...
    'Display','iter-detailed',...
    'PlotFcn','optimplotfval',...
    'UseParallel',1);

lb = -10*ones(1,num_sp); %1.2*max(hpscale)*ones(1,num_sp);
ub = 10*ones(1,num_sp); %1.2*max(hpscale)*ones(1,num_sp);

[est, fval]= fmincon(model,init_vals,[],[],[],[],lb,ub,[],options);

est_vec(itr,:) = est;
fval_vec(itr) = fval;
end
%%
disp(['init_vals = ' num2str(hpscale)])
disp(['est = ' num2str(est)])
disp(['fval = ' num2str(fval)])
%%
