function [b1, gz] = genVSVLpulse(phi_sp, phi_inv, delay, Vmax, inv_amp, ...
    varargin)
%   generate Qin Qin's pulse, with adapted sub-pulse envelope!
%
%   Author: Vanessa, adapted from Terrence Jao
%   Date:   1-29-2014
%   Ref:
%   [1] Qin Qin's pulse

%--------------------------------------------------------------------------
%Input Parser
%--------------------------------------------------------------------------
p = inputParser;
p.addRequired('phi_sp');
p.addRequired('phi_inv');
p.addRequired('delay', @isscalar);
p.addRequired('Vmax', @isscalar);
p.addParamValue('sym', 1, @(x) x == 1 || x==0);
p.addParamValue('dt', [2 4], @(x) length(x) == 2);
p.addParamValue('Gmax', 2, @isscalar);          %Max Grad 4 G/cm
p.addParamValue('Smax', 15, @isscalar);         %Max Slew 15 G/cm/ms
p.addParamValue('B1max', 0.16, @isscalar);      %Max B1 ampliutde 0.16G
p.addParamValue('gamma',.026752, @isscalar);    %Gamma radians/us/G
p.addParamValue('ramp_time', 300, @isscalar);   %Ramp Time us
p.addParamValue('grad_wait', 0, @isscalar);     %Gradient Wait
p.parse(phi_sp, phi_inv, delay, Vmax, varargin{:});
dt = p.Results.dt;
gmax = p.Results.Gmax;
smax = p.Results.Smax;
b1max = p.Results.B1max; 
gamma = p.Results.gamma;
r = p.Results.ramp_time;
grad_wait = p.Results.grad_wait;

%dtRF = dt(1);   dtGz = dt(2);
smax = 1e-3*smax;                %conversion from G/cm/ms to G/cm/us
%--------------------------------------------------------------------------
%Generate RF and Gradient
%--------------------------------------------------------------------------
gam = 42.58; % MHz/T

% amplitudes
B1_val_hp = 3.13; %%3.88; % %5.2; % uT
B1_val_inv = 10; % 11.5?? % uT
Grad_val = 2.2; %1.45; % mT/cm

% durations (mS)
hp_dur = 20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;
grad_ramp = 0.302;
Tgap = delay; %.34;

% sampling periods
dtRF = 0.002*51;
dtGz = 0.002*51;
 
% res
hp_res = round(hp_dur/dtRF);
inv_res = round(inv_dur/dtRF);
grad_ramp_res = round(grad_ramp/dtGz);
Tgap_res = round(Tgap/dtRF); 

% sub-waveforms
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
grad_wait = zeros(1, length(grad));

inv = ones(1,inv_res)*B1_val_inv;
inv_wait = zeros(1,inv_res); 
gap = zeros(1,Tgap_res);

phi = phi_inv; %[0 0 pi pi pi 0 0 pi  pi pi 0 0  0 pi pi 0];

rf_sub = ones(1,hp_res)*B1_val_hp;
hpscale = 9/sum(dzrf(9,2,'inv','max'))*dzrf(9,2,'inv','max'); 

  
b1 = []; gz = []; 
for step = 1:8 
        b1 = [b1      hpscale(step)*rf_sub.*exp(1i*phi_sp(step))       grad_wait  gap  ...
            inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait  zeros(1,length(rf_sub))    grad_wait  gap  ...
            inv.*exp(1i*phi(2*step)).*inv_amp(2*step)          gap grad_wait ];
        gz = [gz       zeros(1,length(rf_sub))     -grad      gap  ...
            inv_wait  gap -grad          zeros(1,length(rf_sub)) -grad      gap  ...
            inv_wait           gap -grad ]; 
end
b1 = [b1 hpscale(end)*rf_sub.*exp(1i*phi_sp(end))  ]*1e-2; 
gz = [gz zeros(1,length(rf_sub))  ];
 