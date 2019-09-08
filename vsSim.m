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
Vmax = 20;       
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

opt = 4; % hard (1) or max-phase (4) pulse envelope

gam = 42.58; % MHz/T

% amplitudes
B1_val_hp = 3.13; %3.13; %%3.88; % %5.2; % uT
B1_val_inv = 9.98; %%10.02;  % 11.5?? % uT
Grad_val = 1.45; % mT/cm

% durations (mS)
hp_dur = 20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;
grad_ramp = 0.302; 
Tgap = 0.1; %.34;

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
    Grad_val = Grad_val*.621*5/10; % Vc=10, hard inversions
    %Grad_val = 0.3965*Grad_val*6/10; %Vc=10, composite inversions
elseif opt == 2
    grad_ramp_res = round(9/8*grad_ramp_res);
elseif opt == 3
    grad_ramp_res = 1/2*round(9/8*grad_ramp_res);
elseif opt == 4
    %grad_ramp_res = 0.67*round(grad_ramp_res); %round(9/8*grad_ramp_res);
    %grad_ramp_res = 0.96*round(grad_ramp_res); % Vc=6
    Grad_val = 0.96*Grad_val*5/10; %*8/10;
    %grad_ramp_ref = round(grad_ramp_res); 
    %Grad_val = Grad_val*0.621; 
    
end
Tgap_res = round(Tgap/dtRF);

% sub-waveforms
gap_bp = zeros(1,100);
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
%grad_bp = 0*[decimate(grad,2) gap_bp gap_bp -decimate(grad,2)];
%grad = 2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;

grad_wait = zeros(1, length(grad));

hp = ones(1,hp_res)*B1_val_hp;
if (opt==4)
  %  inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
     inv = [ones(1,1/2*inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,inv_res) zeros(1,162) ones(1,1/2*inv_res).*exp(1i*pi/2)]*B1_val_inv;
  %  inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*2).*exp(1i*2*pi/3) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
  %  inv = [ones(1,inv_res)]*B1_val_inv;
  %  inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*4/3).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;


else
 %   inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
    inv = [ones(1,inv_res)]*B1_val_inv;
 %   inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*2).*exp(1i*2*pi/3) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
 %   inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*4/3).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
 %   inv = [ones(1,1/2*inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,inv_res*4/3) zeros(1,162) ones(1,1/2*inv_res).*exp(1i*pi/2)]*B1_val_inv;

end
hp_wait = zeros(1,hp_res);
inv_wait = zeros(1,length(inv));
gap = zeros(1,Tgap_res);
inv_amp = ones(1,18);
phi = [0 0 pi pi pi 0 0 pi] + 0;
phi = [phi pi pi 0 0  0 pi pi 0];
phi_hp = zeros(1,9); 
%%
for itr = 5; %1:size(est_opts,1)
    
 %  phi_hp = est_opts(itr,1:9);
 %  phi = est_opts(itr,10:25);


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
    %rf = dzrf(9,2,'inv','max');
    %hpscale = 9/sum(rf)*rf; 
    hpscale = ones(1,9);

end
b1 = []; gz = []; gz_flip = []; gz_off = [];
b12 = [];
%% double refocusing
rf_sub = hp;
       
for step = 1:8
%        if step==1
%            b1_ud = [    hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
%                inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  zeros(1,length(rf_sub))  gap   grad_wait  gap  ...
%                inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
%            b1 = [b1   b1_ud];
%            
%            gz_ud = [zeros(1,length(rf_sub))    gap   grad_wait     gap  ...
%                inv_wait  gap -grad    gap     zeros(1,length(rf_sub)) gap grad      gap  ...
%                inv_wait      gap    grad_bp gap ];
%            gz = [gz gz_ud];
%        elseif step==8
%            b1_ud = [    hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
%                inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  zeros(1,length(rf_sub))  gap   grad_wait  gap  ...
%                inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
%            b1 = [b1   b1_ud];
%            
%            gz_ud = [zeros(1,length(rf_sub))    gap   -grad_bp     gap  ...
%                inv_wait  gap -grad    gap     zeros(1,length(rf_sub)) gap grad      gap  ...
%                inv_wait      gap    grad_wait gap ];
%            gz = [gz gz_ud];
%        else
           b1_ud = [    hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
               inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  zeros(1,length(rf_sub))  gap   grad_wait  gap  ...
               inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
%                    b1_ud = [hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
%                         inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap    grad_wait gap  zeros(1,length(rf_sub))  gap grad_wait  gap  ...
%                         inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
           
           b1 = [b1   b1_ud];
           
           gz_ud = [zeros(1,length(rf_sub))    gap   -grad     gap  ...
            inv_wait  gap grad    gap     zeros(1,length(rf_sub)) gap -grad      gap  ...
            inv_wait      gap    grad gap ];
%         gz_ud = [zeros(1,length(rf_sub))    gap   -grad_bp    gap  ...
%              inv_wait  gap -grad  gap     zeros(1,length(rf_sub)) gap grad gap  ...
%              inv_wait      gap    grad_bp gap ];
         gz = [gz gz_ud]; 
 %      end
                  
end

b1 = [b1 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; % add zeros(1,9000) to make add 1 sec after pulse!
b12 = b1; % add zeros(1,9000) to make add 1 sec after pulse!


gz = [gz zeros(1,length(rf_sub)) ]; 
gz_flip = -abs(gz);
gz_off = zeros(size(gz));



%%
NRF = length(b1);
NGz = length(gz);

rho = abs(b1);
phi = angle(b1);


%RF and Gradient are at different sampling frequencies - convert to ms
tRF = (0:dtRF:dtRF*(NRF-1));       %In milliseconds
tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds

%---Plotting---%
h = figure(1);
set(h, 'color', [1 1 1]);
subplot(3,1,1);
plot(tRF, rho, 'LineWidth', 1.5, 'color', 'blue');
grid on;
set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)] , 'YLim', [0 0.105]);
xlabel('time (ms)');
ylabel('|B1| (G)');
colormap jet;

subplot(3,1,2);
plot(tRF, phi, 'LineWidth', 1.5, 'color', 'red');
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-pi pi], 'YTick', -pi:pi:pi,'YTickLabel', {'-\pi  ','0','\pi '});
xlabel('time (ms)');
ylabel('\phiB1 (rad)');
colormap jet;

subplot(3,1,3);
plot(tGz,gz, 'LineWidth', 1.5, 'color', 'green');
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tRF(end), ...
    'XLim', [0 tRF(end)], 'YLim', [-6 6]);
xlabel('time (ms)');
ylabel('Gv (G/cm)')
colormap jet;

drawnow;
%%
% id1 = 1267, id2 = 1854 mid=1561
% id1 = 4179, id2 = 4766 mid=4473
% dist = 2912
%gz_old = gz;
%gz = gz_flip;
mid_inv = [hp_wait gap grad_wait gap zeros(1,length(inv_wait)/2)];
%inv_dist= [hp_wait gap grad_wait gap 
inv_dist= [inv_wait gap grad_wait gap ...
           hp_wait gap grad_wait gap]; % inv_wait gap grad_wait gap];
mid_inv = length(mid_inv);
inv_dist= length(inv_dist);
inv_flag=1;
m0_val = gz(1)*dtGz/1000/10000;
m1_val = gz(1)*tGz(1)/1000/10000*dtGz/1000;
df = 100;
for itr=1:length(gz)
    m0(itr) = (m0_val+gz(itr)*dtGz/1000/10000); %G/cm*s*1T/10000G
    m1(itr) = (m1_val + gz(itr)*tGz(itr)/1000/10000*dtGz/1000); % G/cm*ms*1s/1000ms*1T/10000G= T/cm*s*s
    if (mod(itr,inv_dist)-mid_inv)==0
        m0(itr) = -0.6*m0(itr);
        m1(itr) = -0.6*m1(itr);
    end
    m0_val = m0(itr);
    m1_val = m1(itr);
end
gamma = 42.58e6; % Hz/T
kv = gamma*m1; % Hz/T*T/cm*s*s = s/cm
m0 = gamma*m0; % T/cm*s*Hz/T = 1/cm
figure(200);
subplot(3,1,1)
plot(tGz,gz,tGz,abs(b1)*10,tGz(2*inv_dist + length(hp_wait)/2),m0(2*inv_dist + length(hp_wait)/2)/100,'r*'); xlabel('time (ms)'); 
legend('gz','b1');

kv_locs = [length(hp_wait)/2 + [0:8]*(2*inv_dist)];

subplot(3,1,2)
plot(tGz,m0,tGz(kv_locs),m0(kv_locs),'r*'); title('0th moment'); xlabel('time (ms)'); ylabel('1/cm');
legend('m0');

subplot(3,1,3)
plot(tGz,kv,tGz(kv_locs),kv(kv_locs),'r*',tGz,zeros(1,length(kv))); title('kv'); xlabel('time (ms)'); ylabel('kv (s/cm)');

drawnow;
%% ------------------------------------------------------------------------
% Demo Pulse Simulation
%--------------------------------------------------------------------------

%z0 = 0; 
z0 = linspace(-5,5,61);
v = linspace(-2*Vmax,2*Vmax,51);
%v=0;
%df = 0;
df = linspace(-200,200,11);
%b1scale = 1;  
b1scale = linspace(.6,1.4,15);

Nb1 = length(b1scale);
Nv = length(v);
Ndf = length(df);
Nz = length(z0);

% Intialize Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv,Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);

% Arterial blood T1 and T2
%T1 = 1.664; 
%T2 = 0.2; %.150;

% no relaxation
T1 = 100;
T2 = 100;

% silicone oil
%T1 = 1.111;
%T2 = 227;

%T1 = 1.25; %1.932; % 1.420; % calf,  1.932; % myocardium
%T2 = 0.05; %0.275; %0.031; % calf, 0.275; % myocardium
%% test

tic;
mx = zeros(Nz,Ndf,Nb1); my=mx; mz=mx;
for itr = 1:Nb1
    itr
    [mx(:,:,itr),my(:,:,itr),mz(:,:,itr)] = bloch(b1*b1scale(itr),gz,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_flip(:,:,itr),my_flip(:,:,itr),mz_flip(:,:,itr)] = bloch(b12*b1scale(itr),gz_flip,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_off(:,:,itr),my_off(:,:,itr),mz_off(:,:,itr)] = bloch(b1*b1scale(itr),gz_off,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
end


[mx2,my2,mz2] = blochvb1(b1,gz,dtRF*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
[~,~,mz2_flip] = blochvb1(b1,gz_flip,dtRF*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
[~,~,mz2_off]  = blochvb1(b1,gz_off,dtRF*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
toc;
%%
%mz = mz_flip;
% Separately show simulations to check what you expect to see with each of the changed features?i.e.
% 	1) sharpness of the *V* selectivity (at perfect B0 and B1+) achieved using max-phase sinc envelope
% 	2) B1+ insensitivity using the MLEV refocusing pulses
% 	3) B1 scale pre-emptively increased?	
%   4) Gradients flipped vs. off.

% %% sharpness of the Vcut selectivity 
% figure(2); clf;
% %subplot(1,4,1);
% plot(v,squeeze(mz(:,df==0,b1scale==1))); axis square; 
% hold on;
% h1 = plot(-10,0,'b*');
% h2 = plot(3,squeeze(mz(v==2,df==0,b1scale==1)),'r*');
% legend([h1 h2],'Vc=-6' ,'Vc=2');
% xlabel('Vc (cm/s)'); ylabel('Mz/M0'); ylim([-1.02 1.02]); xlim([-40 40]); grid on;
% set(gca,'FontSize',14);
% title('\Delta f=0, B1scale=1'); colorbar;

%% Vs stripe artifact
figure(2); clf;
df = linspace(-200,200,11);

% subplot(1,3,1);
% imagesc(df,z0,squeeze(mz(:,:,b1scale==1))); axis square; colorbar; caxis([-1 0.4]);
% xlabel('\Delta f (Hz)'); ylabel('z (cm)');

toplot = squeeze(mz(:,df==0,:));
subplot(3,2,1);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.4]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Control');

toplot = squeeze(mz(:,df==0,:));
subplot(3,2,2);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.4]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Control');

toplot = squeeze(mz_flip(:,df==0,:));
subplot(3,2,3);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.4]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label')

toplot = squeeze(mz_off(:,df==0,:));
subplot(3,2,4);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.4]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label')

toplot = mz(:,df==0,:)-mz_flip(:,df==0,:);
toplot = squeeze(toplot);
subplot(3,2,5);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-0.2 0.2]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label - Control');

toplot = mz(:,df==0,:)-mz_off(:,df==0,:);
toplot = squeeze(toplot);
subplot(3,2,6);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-0.2 0.2]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label - Control');

% complex sum!? 
%% Vs profile

figure(3); clf;
subplot(1,2,1);
imagesc(df,v,squeeze(mz2(:,:,b1scale==1))); axis square; colorbar; caxis([-1 1]);
xlabel('\Delta f (Hz)'); ylabel('v (cm/s)');
subplot(1,2,2);
imagesc(b1scale,v,squeeze(mz2(:,df==0,:))); axis square; colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('v (cm/s)');
% hold on;
% h1 = plot(-10,0,'b*');
% h2 = plot(3,squeeze(mz(v==2,df==0,b1scale==1)),'r*');
% legend([h1 h2],'Vc=-6' ,'Vc=2');
% xlabel('Vc (cm/s)'); ylabel('Mz/M0'); ylim([-1.02 1.02]); xlim([-40 40]); grid on;
% set(gca,'FontSize',14);
% title('\Delta f=0, B1scale=1'); colorbar;

% subplot(1,4,2); imagesc(df, b1scale, squeeze(mz(v==-6,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1scale');  caxis([-1 1]); title('Vc=-6 cm/s');
% set(gca,'FontSize',14);
% subplot(1,4,3); imagesc(df, b1scale, squeeze(mz(v==0,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1scale');  caxis([-1 1]); title('Vc=0 cm/s');
% set(gca,'FontSize',14);
% subplot(1,4,4); imagesc(df, b1scale, squeeze(mz(v==3,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1scale');  caxis([-1 1]); title('Vc=3 cm/s');
% set(gca,'FontSize',14);
% 
% %% B1 scale pre-emptively increased
% figure(3); clf;
% subplot(1,3,1); imagesc(df, v, squeeze(mz(:,:,b1scale==0.7))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title('B1scale=0.7'); %title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);
% subplot(1,3,2); imagesc(df, v, squeeze(mz(:,:,b1scale==1))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title('B1scale=1'); %title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);
% subplot(1,3,3); imagesc(df, v, squeeze(mz(:,:,round(b1scale*100)/100==1.3))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title('B1scale=1.3'); %title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);

%subplot(3,3,1); imagesc(b1scale, v, squeeze(mz(:,df==-150,:))); axis square; colorbar; colormap jet; 
%xlabel('B1scale'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=-150 Hz');
%set(gca,'FontSize',14);
%subplot(3,3,4); imagesc(b1scale, v, squeeze(mz(:,df==0,:))); axis square; colorbar; colormap jet; 
%xlabel('B1scale'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=0 Hz');
%set(gca,'FontSize',14);
%subplot(3,3,7); imagesc(b1scale, v, squeeze(mz(:,df==100,:))); axis square; colorbar; colormap jet; 
%xlabel('B1scale'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=100 Hz');
%set(gca,'FontSize',14);

% figure(101); clf; 
% figure(102); clf;
% sc = 0.4:0.05:1.05;
% for itr2 = 1:14
% [~, b1idx] = min(abs(b1scale-sc(itr2)));
% % figure(101);
% % subplot(1,9,itr2)
% % %plot(v,squeeze(mz(:,df==-100,b1idx) - mz_flip(:,df==-100,b1idx)),'b', ...
% % %     v,squeeze(mz(:,df==-50,b1idx)  - mz_flip(:,df==-50,b1idx)),'c', ...
% % plot(v,squeeze(mz(:,df==0,b1idx)    - mz_flip(:,df==0,b1idx)),'k','LineWidth',2); %, ...
% %      %v,squeeze(mz(:,df==50,b1idx)   - mz_flip(:,df==50,b1idx)),'m', ...
% %      %,squeeze(mz(:,df==100,b1idx)  - mz_flip(:,df==100,b1idx)),'y','LineWidth',2);
% % xlabel('velocity (cm/s)'); ylabel('Mz'); ylim([-0.5 2]); xlim([-4*Vmax 4*Vmax]);
% % hold on;
% % %vidx = v(v>-5&v<5);
% % vidx = v(v<-10);
% % h = area(vidx,2*ones(size(vidx)),-0.5);
% % % h.FaceColor = [0 0.35 0.95];
% % % h.EdgeColor = [0 0.55 0.95];
% % h.FaceColor = [0 0.95 0.35];
% % h.EdgeColor = [0 0.95 0.75];
% % h.FaceAlpha = [0.1];
% % 
% % vidx2 = v(v>10);
% % h2 = area(vidx2,2*ones(size(vidx2)),-0.5);
% % h2.FaceColor = [0 0.95 0.35];
% % h2.EdgeColor = [0 0.95 0.75];
% % h2.FaceAlpha = [0.1];
% % 
% % set(gca,'FontSize',14,'FontWeight','bold');
% % % if (itr2==7) 
% % %     legend('\Deltaf=-100', ...
% % %        '\Deltaf=-50', ...
% % %        '\Deltaf=0', ...
% % %        '\Deltaf=50', ...
% % %        '\Deltaf=100','Location','best');
% % % end
% % title(['b1scale=' num2str(sc(itr2)) ' ']);
% 
% figure(102);
% % subplot(3,6,itr2)
% % plot(v,squeeze(mz(:,df==-100,b1idx)),'b', ...
% %      v,squeeze(mz(:,df==-50,b1idx)),'c', ...
% %      v,squeeze(mz(:,df==0,b1idx)),'k', ...
% %      'LineWidth',2);
% % %     v,squeeze(mz(:,df==50,b1idx)),'m', ...
% % %     v,squeeze(mz(:,df==100,b1idx)),'y',
% % title(['b1scale=' num2str(sc(itr2)) ' ']); set(gca,'FontSize',14,'FontWeight','bold');
% % xlabel('velocity (cm/s)'); ylabel('Tag'); ylim([-1 1]); xlim([-4*Vmax 4*Vmax]);
% % hold on;
% % 
% % subplot(3,6,itr2 + 6)
% % plot(v,squeeze(mz_off(:,df==-100,b1idx)),'b', ...
% %      v,squeeze(mz_off(:,df==-50,b1idx)),'c', ...
% %      v,squeeze(mz_off(:,df==0,b1idx)),'k', ...
% %      'LineWidth',2);
% % %     v,squeeze(mz_off(:,df==50,b1idx)),'m', ...
% % %     v,squeeze(mz_off(:,df==100,b1idx)),'y',
% % set(gca,'FontSize',14,'FontWeight','bold');
% % xlabel('velocity (cm/s)'); ylabel('Control'); ylim([-1 1]); xlim([-4*Vmax 4*Vmax]);
% % hold on;
% toplot = squeeze(mz(:,df==0,b1idx)    - mz_flip(:,df==0,b1idx));
% subplot(2,7,itr2)
% %plot(v,squeeze(mz(:,df==-100,b1idx) - mz_flip(:,df==-100,b1idx)),'b', ...
% %     v,squeeze(mz(:,df==-50,b1idx)  - mz_flip(:,df==-50,b1idx)),'c', ...
% plot(v,toplot,'k','LineWidth',2); %, ...
% %   v,squeeze(mz(:,df==50,b1idx)   - mz_off(:,df==50,b1idx)),'m', ...
% %   v,squeeze(mz(:,df==100,b1idx)  - mz_off(:,df==100,b1idx)),'y',
% 
% xlabel('velocity (cm/s)'); ylabel(['Tag - Control']);
% ylim([-0.5 2]); xlim([-4*Vmax 4*Vmax]);
% hold on;
% %vidx = v(v>-5&v<5);
% vidx = v(v<-10);
% h = area(vidx,2*ones(size(vidx)),-0.5);
% % h.FaceColor = [0 0.35 0.95];
% % h.EdgeColor = [0 0.55 0.95];
% h.FaceColor = [0 0.95 0.35];
% h.EdgeColor = [0 0.95 0.75];
% h.FaceAlpha = [0.1];
% 
% vidx2 = v(v>10);
% h2 = area(vidx2,2*ones(size(vidx2)),-0.5);
% h2.FaceColor = [0 0.95 0.35];
% h2.EdgeColor = [0 0.95 0.75];
% h2.FaceAlpha = [0.1];
% legend(['relative signal= ' num2str(mean(toplot(v>10))) ' ']); 
% %set(gca,'FontSize',14,'FontWeight','bold');
% if (itr2==6) 
% %    legend('\Deltaf=-100', ...
% %       '\Deltaf=-50', ...
% %       '\Deltaf=0', ...
% %       'Location','best');
% %       '\Deltaf=50', ...
% %      '\Deltaf=100','Location','best');
% end
% 
% title(['b1scale=' num2str(sc(itr2)) ' ']);
% 
% end


%%
% figure(3);
% subplot(1,3,1); imagesc(b1scale*max(abs(b1)), v, squeeze(mz(:,df==0,:))); axis square; colorbar; colormap jet; 
% xlabel('B1 (G)'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=0 Hz');
% set(gca,'FontSize',14);
% subplot(1,3,2); imagesc(df, b1scale*max(abs(b1)), squeeze(mz(v==0,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1 (G)');  caxis([-1 1]); title('Vc=0 cm/s');
% set(gca,'FontSize',14);
% subplot(1,3,3); imagesc(df, v, squeeze(mz(:,:,round(b1scale*1000)/1000==1))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);
% figure(101); clf; 
% sc = 0.7:0.1:1.3;
% for itr2 = 1:7
% figure(101);
% subplot(1,7,itr2)
% [~, b1idx] = min(abs(b1scale-sc(itr2)));
% plot(v,squeeze(mz(:,df==-100,b1idx)),'b', ...
%      v,squeeze(mz(:,df==-50,b1idx)),'c', ...
%      v,squeeze(mz(:,df==0,b1idx)),'k', ...
%      v,squeeze(mz(:,df==50,b1idx)),'m', ...
%      v,squeeze(mz(:,df==100,b1idx)),'y','LineWidth',2);
% xlabel('velocity (cm/s)'); ylabel('Mz'); ylim([-1 1]); xlim([-4*Vmax 4*Vmax]);
% hold on;
% %vidx = v(v>-5&v<5);
% vidx = v(v<-10);
% h = area(vidx,2*ones(size(vidx)),-1);
% %h.FaceColor = [0 0.35 0.95];
% %h.EdgeColor = [0 0.55 0.95];
% h.FaceColor = [0 0.95 0.35];
% h.EdgeColor = [0 0.95 0.75];
% h.FaceAlpha = [0.1];
% 
% vidx2 = v(v>10);
% h2 = area(vidx2,2*ones(size(vidx2)),-1);
% h2.FaceColor = [0 0.95 0.35];
% h2.EdgeColor = [0 0.95 0.75];
% h2.FaceAlpha = [0.1];
% 
% set(gca,'FontSize',14,'FontWeight','bold');
% if (itr2==7) 
%     legend('\Deltaf=-100', ...
%        '\Deltaf=-50', ...
%        '\Deltaf=0', ...
%        '\Deltaf=50', ...
%        '\Deltaf=100','Location','best');
% end
% title(['b1scale=' num2str(sc(itr2)) ' ']);
% 
% end


%%
% idx = find((squeeze(mz(:,df==0,b1scale==1))<0.8)&(squeeze(mz(:,df==0,b1scale==1))>-0.8));
% vals = v(idx);
% vals = vals(vals<0)
% v_prof_width = max(vals) - min(vals)
% v_cut_idx = find(abs(squeeze(mz(50:150,df==0,b1scale==1)))==min(abs(squeeze(mz(50:150,df==0,b1scale==1)))));
% v_cut = v(v_cut_idx + 49)
end

