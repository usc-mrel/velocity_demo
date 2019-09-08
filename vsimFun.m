% function J = vsimFun(tm, b1nom, Vc, ramp_time, rf_width, grad_wait, ...
%     dtRF, dtG, type, df, b1scale, v, T1, T2)
function J = vsimFun(tm, b1nom, Vc, ramp_time, grad_wait, ...
    dtRF, dtG, type, df, b1scale, v, T1, T2)
%VSIMFUM passes velocity Selective BIR4 or BIR8 pulse parameters to a bloch
%   simulator and returns the objective function J, the squared residual
%   magnetization across a range of B1 scales and Off-resonances subtracted
%   by it's ideal profile with no B1 inhomogeneity or off-resonance for 
%   use in a non-linear optimization scheme (fmincon). 
%   
%   J = VSIMFUM(tm, 1nom, Vc, ramp_time, rf_width, grad_wait, ...
%    dtRF, dtG, type, df, b1scale, v, T1, T2)    
%   returns the objective function, J, that a non-linear optimization 
%   scheme is trying to minimize by finding optimal pulse paramters, tm,  
%   which correspond to kappa, zeeta, and omega, for a VS-BIR4 or BIR8
%   pulse
%   
%   Imaging parameters: 
%
%   Author: Terrence Jao
%   Date: ...

kappa = atan2(tm(1),1);
zeta = tm(2);
omega = tm(3); %kHz
rf_width = tm(4);

Nv = length(v);
Ndf = length(df);
Nb1 = length(b1scale);

% Generate Pulse
[b1, gz, ~] = genVSBIR(zeta, kappa, omega, Vc, 'type', type, ...
    'Smax', 7.5, 'dt', [dtRF dtG], 'b1max', b1nom, 'sym', 1, ...
    'ramp_time', ramp_time, 'grad_wait', grad_wait, 'RF_width', rf_width); 

NRF = length(b1);
NGz = length(gz);

rho = abs(b1);
phi = angle(b1);

%RF and Gradient are at different sampling frequencies - convert to ms
tRF = 1e-3*(0:dtRF:dtRF*(NRF-1));       %In milliseconds
tGz = 1e-3*(0:dtG:dtG*(NGz-1));       %In milliseconds

h = figure(1);
set(h, 'color', [1 1 1]);
subplot(3,1,1);
plot(tRF, rho, 'LineWidth', 1.5, 'color', 'blue');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1, 'XTick', 0:2:tRF(end), ...
    'XLim', [0 tRF(end)]);
xlabel('time (ms)');
ylabel('B1 Magnitude');

subplot(3,1,2);
plot(tRF, phi, 'LineWidth', 1.5, 'color', 'red');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1, 'XTick', 0:2:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-pi pi], 'YTick', -pi:pi/2:pi);
xlabel('time (ms)');
ylabel('B1 Phase');

subplot(3,1,3);
plot(tGz,gz, 'LineWidth', 1.5, 'color', 'green');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1, 'XTick', 0:2:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-6 6]);
xlabel('time (ms)');
ylabel('Gradient (G/cm)')

drawnow;

% Initial Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv, Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);

% Bloch Simulation over range of B1 and Off-resonance

% Reference Profile
% [~,~,mz_ref] = blochvb1(b1,gz,dtRF*1e-6,T1,T2,0,0,v,1,0,mx0,my0,mz0);
mz_ref = cos(pi*v/Vc).';
mz_ref = repmat(mz_ref, [1 Ndf Nb1]);

% Pulse Profile
[~,~,mz] = blochvb1(b1,gz,dtRF*1e-6,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
% mz = reshape(mz, [Nv, Ndf*Nb1]);

figure(2);
subplot(1,2,1); imagesc(mz_ref); caxis([-1 1]); axis square; 
subplot(1,2,2); imagesc(mz); caxis([-1 1]); axis square;
drawnow;

%Weight
% weight static tissue more heavily - can't have spurious labeling

stdGauss = 2*Vc;
gauss = @(diffX) exp(-((diffX).^2./(2*stdGauss^2)));
w = gauss(v).';
w = repmat(w, [1 Ndf]);
% w = repmat(w, [1 Ndf*Nb1]);
w = w.*(w.');
w = repmat(w, [1 1 Nb1]);

J = sum(mean(w.*(mz-mz_ref).^2,2))/Nv;
