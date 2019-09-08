% function J = vsimFun(tm, b1nom, Vc, ramp_time, rf_width, grad_wait, ...
%     dtRF, dtG, type, df, b1scale, v, T1, T2)
function J = vsvlimFun(tm, b1nom, Vc, ramp_time, grad_wait, ...
    dtRF, dtG, df, b1scale, v, T1, T2)
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
%   Author: Vanessa, adapted from Terrence Jao
%   Date: ...

phi_sp_nom = zeros(1,9);
phi_sp = tm(1:9); 
phi_inv_nom = [0 0 pi pi pi 0 0 pi pi pi 0 0 0 pi pi 0];
phi_inv = tm(10:25); 

delay_nom = 0; 
delay = delay_nom; %tm(end);

inv_nom = ones(1,16);
inv_amp = inv_nom; %tm(1:16); 

Nv = length(v);
Ndf = length(df);
Nb1 = length(b1scale);

[b1, gz] = genVSVLpulse(phi_sp, phi_inv, delay, Vc, inv_amp, ...
    'Smax', 7.5, 'dt', [dtRF dtG], 'b1max', b1nom, 'sym', 1, ...
    'ramp_time', ramp_time, 'grad_wait', grad_wait); 

[b10, ~] = genVSVLpulse(phi_sp_nom, phi_inv_nom, delay_nom, Vc, inv_nom,  ...
    'Smax', 7.5, 'dt', [dtRF dtG], 'b1max', b1nom, 'sym', 1, ...
    'ramp_time', ramp_time, 'grad_wait', grad_wait); 

NRF = length(b1);
NGz = length(gz);

rho = abs(b1);
phi = angle(b1);

%RF and Gradient are at different sampling frequencies - convert to ms
tRF = 1e-3*(0:dtRF:dtRF*(NRF-1));       %In milliseconds
tGz = 1e-3*(0:dtG:dtG*(NGz-1));       %In milliseconds
 
% Initial Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv, Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);

% Reference Profile
[~,~,mz_ref] = blochvb1(b10,gz,dtRF*1e-6,T1,T2,0,0,v,1,0,mx0,my0,mz0);
%%
mz_ref = repmat(mz_ref, [1 Ndf Nb1]);
for itr = 1:Nb1 
    mz_ref(18:32,:,itr) = mz_ref(18:32,:,itr)*b1scale(itr);
end

 
% Pulse Profile 
[~,~,mz] = blochvb1(b1,gz,dtRF*1e-6,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);

%Weight 
% weight static tissue more heavily - can't have spurious labeling
%w = [0.2*ones(1,18) ones(1,15) 0.2*ones(1,18)];
w = [0.2*ones(1,17) ones(1,6) 5*ones(1,2) 10 5*ones(1,2) ones(1,6) 0.2*ones(1,17)];
w = repmat(w.', [1 Ndf Nb1]); 

J = mean(mean(mean(w.*(mz-mz_ref).^2)));
%J = max(max(max((mz(18:32)-mz_ref(18:32)).^2))); 
