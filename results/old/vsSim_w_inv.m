%% ------------------------------------------------------------------------
% Demo Pulse Generation
%--------------------------------------------------------------------------

% close all; clear all;
clear all;


kappa = atan2(62.96,1);
zeta = 20.50;
omega = 21.60; %kHz
RF_width = 2236;

%Tried Optimizattion
Vmax = 10;       
b1max = .16;
Gmax = 10;
ramp_time = 500;

grad_wait = 748;
type = 8;

%% TJAOs TRY ---Pulse Generation---%
% Vmax = 10; 
% [b1, gz, dt] = genVSBIR(zeta, kappa, omega, Vmax, 'type', type, ...
%     'Smax', 7.5, 'dt', [4 4], 'b1max', b1max, 'sym', 1, 'ramp_time', ramp_time, 'grad_wait', grad_wait, 'RF_width', RF_width, 'Gmax', Gmax); %'ramp_time', 1000,
% NRF = length(b1);
% NGz = length(gz);
% dtRF = dt(1); dtGz = dt(2);
% 
% rho = abs(b1);
% phi = angle(b1);
% 
% b10 = b1;

%% VL TRIES!

opt = 4; % hard (1) or sinc (2) pulse envelope

gam = 42.58; % MHz/T

% amplitudes
B1_val_hp = 3.13; %%3.88; % %5.2; % uT
B1_val_inv = 10.02;  % 11.5?? % uT
Grad_val = 1.45; % mT/cm

% durations (mS)
hp_dur = 20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;
grad_ramp = 0.302; 
Tgap = 0; %.34;

% sampling periods
dtRF = 0.002; %*51;
dtGz = 0.002; %*51;

% res
hp_res = round(hp_dur/dtRF);
inv_res = round(inv_dur/dtRF);
grad_ramp_res = round(grad_ramp/dtGz);
if opt == 1
    %grad_ramp_res = round(.45*grad_ramp_res);
    grad_ramp_res = round(grad_ramp_res);
    Grad_val = Grad_val/1.40;
elseif opt == 2
    grad_ramp_res = round(9/8*grad_ramp_res);
elseif opt == 3
    grad_ramp_res = 1/2*round(9/8*grad_ramp_res);
elseif opt == 4
    grad_ramp_res = round(grad_ramp_res); %round(9/8*grad_ramp_res);
    
end
Tgap_res = round(Tgap/dtRF);

% sub-waveforms
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
grad2 = [1:grad_ramp_res grad_ramp_res grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
grad_wait = zeros(1, length(grad));

hp = ones(1,hp_res)*B1_val_hp;
if (opt==4)
    inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
else
    inv = [ones(1,inv_res)]*B1_val_inv;
end
hp_wait = zeros(1,hp_res);
inv_wait = zeros(1,length(inv));
gap = zeros(1,Tgap_res);
inv_amp = ones(1,18);
phi = [0 0 pi pi pi 0 0 pi] + 0;
phi = [phi pi pi 0 0  0 pi pi 0];
phi_hp = zeros(1,9); 


    
 %  phi_hp = est_opts(itr,1:9);
 %  phi = est_opts(itr,10:25);

%% 
if opt == 1;
    hpscale = ones(1,9);
elseif opt == 2;
    hp = ones(1,2*hp_res)*1/2*B1_val_hp;
    hp_wait = zeros(1,2*hp_res);
    hpscale = 9/sum(dzrf(9,4,'st'))*dzrf(9,4,'st'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
elseif opt == 3;
    hp = ones(1,2*hp_res)*1/2*B1_val_hp;
    hp_wait = zeros(1,2*hp_res);
    hpscale = 9/sum(dzrf(9,2,'st'))*dzrf(9,2,'st'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
elseif opt==4;
    hp = ones(1,hp_res)*B1_val_hp; %ones(1,2*hp_res)*1/2*B1_val_hp;
    hp_wait = zeros(1,hp_res);
    hpscale = 9/sum(dzrf(9,2,'inv','max'))*dzrf(9,2,'inv','max'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);

end
b1 = []; gz = []; gz_flip = []; gz_off = [];

%% double refocusing
rf_sub = hp;
for step = 1:8
        b1 = [b1       hpscale(step)*rf_sub.*exp(1i*phi_hp(step))          grad_wait  gap  ...
            inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait  zeros(1,length(rf_sub))    grad_wait  gap  ...
            inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait ];

        gz = [gz       zeros(1,length(rf_sub))       -grad     gap  ...
            inv_wait  gap grad         zeros(1,length(rf_sub)) -grad      gap  ...
            inv_wait           gap grad ]; 
        
        gz_flip = [gz_flip       zeros(1,length(rf_sub))       -grad     gap  ...
            inv_wait  gap -grad         zeros(1,length(rf_sub)) -grad      gap  ...
            inv_wait           gap -grad ]; 
        
        gz_off = [gz_off       zeros(1,length(rf_sub))       -grad*0     gap  ...
            inv_wait  gap 0*grad         zeros(1,length(rf_sub)) -grad*0      gap  ...
            inv_wait           gap 0*grad ]; 
        
     
end
b1 = [b1 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; % add zeros(1,9000) to make add 1 sec after pulse!
sc2 = 1; %0.5; %0.4;
gz = sc2*[gz zeros(1,length(rf_sub)) ]; 
gz_flip = sc2*[gz_flip zeros(1,length(rf_sub)) ]; 
gz_off = sc2*[gz_off zeros(1,length(rf_sub))]; 

%%
NRF = length(b1);
NGz = length(gz);

rho = abs(b1);
phi = angle(b1);

%% end VL tries!

%RF and Gradient are at different sampling frequencies - convert to ms
tRF = (0:dtRF:dtRF*(NRF-1));       %In milliseconds
tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds

%---Plotting---%
h = figure(1);
set(h, 'color', [1 1 1]);
subplot(3,1,1);
plot(tRF, rho, 'LineWidth', 1.5, 'color', 'blue');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)]);
xlabel('time (ms)');
ylabel('B1 Magnitude');
colormap jet;

subplot(3,1,2);
plot(tRF, phi, 'LineWidth', 1.5, 'color', 'red');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-pi pi], 'YTick', -pi:pi/2:pi);
xlabel('time (ms)');
ylabel('B1 Phase');
colormap jet;

subplot(3,1,3);
plot(tGz,gz, 'LineWidth', 1.5, 'color', 'green');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-6 6]);
xlabel('time (ms)');
ylabel('Gradient (G/cm)')
colormap jet;

drawnow;
%% ------------------------------------------------------------------------
% Demo Pulse Simulation
%--------------------------------------------------------------------------

z0 = 0;
v = linspace(-Vmax,Vmax,5);
%df = 0;
df = linspace(-10,10,3);
%b1scale = 1; 
b1scale = linspace(.9,1.1,3);

Nb1 = length(b1scale);
Nv = length(v);
Ndf = length(df);

HR = 80; %:120;
for itr2 = 1:length(HR)
% Intialize Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv, Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);
% Arterial blood T1 and T2
%T1 = 1.664; 
%T2 = 0.2; %.150;

T1 = 1.25; %1.932; % 1.420; % calf,  1.932; % myocardium
T2 = 0.05; %0.275; %0.031; % calf, 0.275; % myocardium


[mx,my,mz] = blochvb1(b1,gz,dtRF*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz_flip] = blochvb1(b1,gz_flip,dtRF*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz_off] = blochvb1(b1,gz_off,dtRF*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);





TR = 60/HR(itr2)/dtRF;
[mx_ref,my_ref,mz_ref] = blochvb1(zeros(1,round(TR)),zeros(1,round(TR)),dtRF,T1,T2,df,z0,v,b1scale,0,mx,my,mz);
mzref_0 = mz_ref(3,2,2);
npts = 2001; 
mzvec = zeros(1,npts); mzvec(1) = mz(3,2,2);
mx_pt = mx; my_pt = my; mz_pt = mz;
for itr = 1:TR
    [mx_u,my_u,mz_u] = blochvb1(0,0,dtRF,T1,T2,df,z0,v,b1scale,0,mx_pt,my_pt,mz_pt);
    mzvec_pt(itr) = mz_u(3,2,2);   
    mx_pt = mx_u;
    my_pt = my_u;
    mz_pt = mz_u;
end
[~, inv_pt] = min(abs(mzref_0 - mzvec_pt));
% contstrained optimization: inv_pt such that mz_ref_0 = 0 !
for itr = 1:npts
    [mx_u,my_u,mz_u] = blochvb1(0,0,dtRF,T1,T2,df,z0,v,b1scale,0,mx,my,mz);
    mzvec(itr) = mz_u(3,2,2);
    if itr == inv_pt
        mx = -mx_u;
        my = -my_u;
        mz = -mz_u;
    else       
        mx = mx_u;
        my = my_u;
        mz = mz_u;
    end
end
    inv_time(itr2) = inv_pt*dtRF;
end

%%
figure(1010110); clf; plot((1:npts),mzvec); xlabel('sec'); ylabel('Mz');
hold on; plot(inv_pt,mzvec(inv_pt),'r*'); plot(2*TR,mzvec(2*TR),'b*');
figure(500); plot(HR, inv_time); xlabel('HR'); ylabel('ms'); title('inversion time');