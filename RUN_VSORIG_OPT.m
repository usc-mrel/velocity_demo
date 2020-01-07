%RUN_VS_OPT designs VS BIR4 or BIR8 with optimized at a user defined range
% of B0 and B1.
%
%   Author: Terrence Jao
%   Date: ...
%   References:

%% Pulse Timing Calculation
clear all; close all;
% -------------------------------------------------------------------------
% Do we want to plot?
ploton = 0;
% Do we want to display messages?
verbose = 1;        %Turn off verbose for slightly faster runtime

%--------------------------------------------------------------------------
%Parameters of BGS Optimization Script - Tweak Here!
%--------------------------------------------------------------------------


% ---- Tissue Paramters ---- %
T1 = 1.932;         % T1 Blood at 3T (x)
T2 = 0.275;         % T2 blood 3T (s)

% ---- Pulse Paramters ---- %
Vc = 10;            % Cutoff Velocity (cm/s)
b1nom = .15;        % b1 scale
ramp_time = 500;    % length of ramps (us)
% rf_width = 1500;    % length of RF subpulse (us)
grad_wait = 0;      % Delay between gradient and RF (us)
dtRF = 4;           % time step RF (us)
dtG = 4;            % time step Gradient (us)
type = 8;           % 4 - BIR4 or 8 - BIR8


% ---- VS Optimization Paramters ---- %
z0 = 0;
Nv = 51;            % Number of velocities to simulate
Nf = 51;           % Number of frequencies to simulate
Nb1 = 1;           % Number of b1 scales to simulate

delf = 500;         % off-resonance range to simulate (-delf, + delf) (Hz);

v = linspace(-2*Vc, 2*Vc, Nv);
df = linspace(-delf,delf,Nf);
% b1scale = linspace(.6,1.2,Nb1);
b1scale = 1;

maxiter = 20;
outIter = 10;
opttype = 'sqp';        %Optimiation Algorithm

if verbose == 1
    %     fprintf(['Simulating BIR%d with:\nVc: %0.1f, b1nom %0.1f, ramptime: %d\n'...
    %         'RF_width: %d, grad_wait = %d, dtRF = %d, dtG = %d\n'], ...
    %        type, Vc, b1nom, ramp_time, rf_width, grad_wait, dtRF, dtG);
    fprintf(['Simulating BIR%d with:\nVc: %0.1f, b1nom %0.1f, ramptime: %d\n'...
        'grad_wait = %d, dtRF = %d, dtG = %d\n'], ...
        type, Vc, b1nom, ramp_time, grad_wait, dtRF, dtG);
end

%--------------------------------------------------------------------------
%Running the Optimization
%--------------------------------------------------------------------------

if verbose; fprintf('Performing Optimization using %s\n', opttype); end;

options = optimset('Algorithm', opttype, 'display', 'iter', 'MaxIter', maxiter);

kappaParam = 60;
zeta = 15;
omega = 39.8; %kHz
rf_width = 1500;

%Initial Guess for Tm (remember you may get stuck at local minima)
tm = [kappaParam, zeta, omega, rf_width];
fvalprev = inf;

for jj = 1:outIter
    fprintf('-------------------------------\n');
    fprintf('Outer Iteration %d\n', jj);
    fprintf('-------------------------------\n');
    %Optimizing metric function
    %     f = @(x) vsimFun(x, b1nom, Vc, ramp_time, rf_width, grad_wait, ...
    %         dtRF, dtG, type, df, b1scale, v, T1, T2);
    
    f = @(x) vsimFun(x, b1nom, Vc, ramp_time, grad_wait, ...
        dtRF, dtG, type, df, b1scale, v, T1, T2);
    
    %Constraint function
    %c = @(x) confun(x,TS(k),flipSchd, TR);
    
    %Optimization With Constraints
    %     [tm_opt,fval] = fminsearch(f,tm, [], [], [], [], [],[], ...
    %         c, options);
    [tm_opt,fval] = fminsearch(f,tm, options);
    tm = tm_opt;
    
    %Stopping Criteria
    if abs(fvalprev-fval) < 1e-5
        break;
    else
        fvalprev = fval;
    end
end
if verbose; fprintf('\nOptimization Complete\n'); end;


%% ------------------------------------------------------------------------
%Running Simulation with Optimal Parameters
%--------------------------------------------------------------------------
kappa = atan2(tm_opt(1),1);
zeta = tm_opt(2);
omega = tm_opt(3); %kHz
rf_width = tm_opt(4);

rf_width = ceil(rf_width);
rf_width = rf_width - mod(rf_width,4);

fprintf('Optimal kappaParam = %0.2f kappa = %0.2f, zeta = %0.2f omega = %0.3f rf_width = %0.2f\n', tm(1), kappa, zeta, omega, rf_width);


% Generate Pulse
[b1, gz, ~] = genVSBIR(zeta, kappa, omega, Vc, 'type', type, ...
    'Smax', 7.5, 'dt', [dtRF dtG], 'b1max', b1nom, 'sym', 1, ...
    'ramp_time', ramp_time, 'grad_wait', grad_wait, 'RF_width', rf_width);

Nb1 = 11;
b1scale_opt = linspace(.5, 1.5, Nb1);
% Initial Magnetization
mx0 = zeros(Nv, Nf, Nb1); my0 = zeros(Nv, Nf, Nb1); mz0 = ones(Nv, Nf, Nb1);

% Bloch Simulation over range of B1 and Off-resonance
[~,~,mz] = blochvb1(b1,gz,dtRF*1e-6,T1,T2,df,0,v,b1scale_opt,0,mx0,my0,mz0);

figure(2);
subplot(1,3,1); imagesc(b1scale_opt*b1nom, v, squeeze(mz(:,df==0,:))); axis square; colorbar;
xlabel('B1'); ylabel('Velocity');
subplot(1,3,2); imagesc(df, b1scale_opt*b1nom, squeeze(mz(v==0,:,:)).'); axis square; colorbar;
xlabel('Frequency'); ylabel('B1');
subplot(1,3,3); imagesc(df, v, squeeze(mz(:,:,b1scale_opt==1))); axis square; colorbar;
xlabel('Frequency'); ylabel('Velocity');
