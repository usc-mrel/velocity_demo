%% ------------------------------------------------------------------------
% Demo Pulse Generation
%--------------------------------------------------------------------------

% close all; clear all;
clear all;
  

%% VL TRIES!

opt = 4; % hard (1) or max-phase (4) pulse envelope

gam = 42.58; % MHz/T

% amplitudes
B1_val_hp = 3.13; %%3.88; % %5.2; % uT
B1_val_inv = 10.02;  % 11.5?? % uT
Grad_val = 1.45; % mT/cm
if opt==1; Grad_val = Grad_val/1.4; end;


% durations (mS)
hp_dur = 20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;
grad_ramp = 0.302; 
Tgap = 0.1; %.34;

% sampling periods
dtRF = 0.002; 
dtGz = 0.002; 

% res
hp_res = round(hp_dur/dtRF);
inv_res = round(inv_dur/dtRF);
grad_ramp_res = round(grad_ramp/dtGz);
Tgap_res = round(Tgap/dtRF);

% sub-waveforms
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
grad2 = [1:grad_ramp_res grad_ramp_res grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
grad_wait = zeros(1, length(grad));

hp = ones(1,hp_res)*B1_val_hp;
if opt==4; 
    hpscale = 9/sum(dzrf(9,2,'inv','max'))*dzrf(9,2,'inv','max'); 
    inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
else
    hpscale = ones(1,9);
    inv = ones(1,inv_res)*B1_val_inv;
end


hp_wait = zeros(1,hp_res);
inv_wait = zeros(1,length(inv));
gap = zeros(1,Tgap_res);
inv_amp = ones(1,18);
phi_hp = zeros(1,9); 
phi = [0 0 pi pi pi 0 0 pi] + 0;
phi = [phi pi pi 0 0  0 pi pi 0];

b1 = []; gz = []; gz_flip = []; gz_off = [];
b12 = []; b13 = []; b14 = [];

% b1 = dzrf(500,4,'se');
% gz = 2*ones(size(b1));
rf_sub = hp;
for step = 1:8
        b1_ud = [    hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
            inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  zeros(1,length(rf_sub))  gap   grad_wait  gap  ...
            inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
        b1 = [b1   b1_ud];

        
        gz_ud = [zeros(1,length(rf_sub))    gap   -grad     gap  ...
            inv_wait  gap grad    gap     zeros(1,length(rf_sub)) gap -grad      gap  ...
            inv_wait      gap    grad gap ];
        gz = [gz gz_ud]; 
        
        %phi2 = gam*cumsum(pi/4*gz_ud)/length(gz_ud);
        phi2 = [    0*hpscale(step)*rf_sub.*exp(1i*phi_hp(step))      gap    grad_wait  gap  ...
            pi/2*ones(1,length(inv)) gap     grad_wait gap  zeros(1,length(rf_sub))  gap   grad_wait  gap  ...
            pi/2*ones(1,length(inv))         gap grad_wait gap ];
        b12 = [b12 b1_ud.*exp(1i*phi2)];
        
        %phi3 = gam*cumsum(pi/2*gz_ud)/length(gz_ud);
        phi3 = 2*phi2;
        b13 = [b13 b1_ud.*exp(1i*phi3)];
        
        %phi4 = gam*cumsum(3*pi/4*gz_ud)/length(gz_ud);
        phi4 = 3*phi2;
        b14 = [b14 b1_ud.*exp(1i*phi4)];

end
b1 = [b1 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; % add zeros(1,9000) to make add 1 sec after pulse!
b12 = [b12 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; 
b13 = [b13 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2;
b14 = [b14 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; 
gz = [gz zeros(1,length(rf_sub)) ]; 



NRF = length(b1);
NGz = length(gz);
%%
rho = abs(b12);
phi = angle(b12);


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
    'XLim', [0 tRF(end)] , 'YLim', [0 B1_val_inv/100]);
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
%% ------------------------------------------------------------------------
% Demo Pulse Simulation
%--------------------------------------------------------------------------

%z0 = 0; 
z0 = linspace(-2.5,2.5,201);
%df = 0;
df = linspace(-100,100,51);
b1scale = linspace(0.7,1.3,51);  

Ndf = length(df);
Nz = length(z0);

% Intialize Magnetization
mx0 = zeros(Nz, Ndf); my0 = zeros(Nz,Ndf); mz0 = ones(Nz, Ndf);

% Arterial blood T1 and T2
%T1 = 1000; %.664; 
%T2 = 1000; %0.2; %.150;

T1 = 1.25; %1.932; % 1.420; % calf,  1.932; % myocardium
T2 = 0.05; %0.275; %0.031; % calf, 0.275; % myocardium
%% test

tic;
[mx,my,mz] =     bloch(b1,gz,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
[mx,my,mz2] =     bloch(b12,gz,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
[mx,my,mz3] =     bloch(b13,gz,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
[mx,my,mz4] =     bloch(b14,gz,dtRF*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
toc;

tic;
for itr = 1:length(b1scale);
[mx,my,Mz] =     bloch(b1.*b1scale(itr),gz,dtRF*1e-3,T1,T2,0 ,z0,0,mx0,my0,mz0);
MZ(:,itr) = Mz;
[mx,my,Mz2] =     bloch(b12.*b1scale(itr),gz,dtRF*1e-3,T1,T2,0 ,z0,0,mx0,my0,mz0);
MZ2(:,itr) = Mz2;
[mx,my,Mz3] =     bloch(b13.*b1scale(itr),gz,dtRF*1e-3,T1,T2,0 ,z0,0,mx0,my0,mz0);
MZ3(:,itr) = Mz3;
[mx,my,Mz4] =     bloch(b14.*b1scale(itr),gz,dtRF*1e-3,T1,T2,0 ,z0,0,mx0,my0,mz0);
MZ4(:,itr) = Mz4;
end
toc;
%%
figure(1010101)
subplot(5,1,1)
imagesc(z0,df,mz.'); %caxis([-0.55 -0.45]); %caxis([-1 -0.9]); %
caxis([-0.45 -0.35]);

set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,2)
imagesc(z0,df,mz2.'); %caxis([-0.55 -0.45]); %caxis([-1 -0.9]);%
caxis([-0.45 -0.35]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,3)
imagesc(z0,df,mz3.'); %caxis([-0.55 -0.45]); %caxis([-1 -0.9]);%
caxis([-0.45 -0.35]); 
ylabel('Off-resonance (Hz)');
set(gca,'FontSize',12,'FontWeight','bold');
colorbar;

subplot(5,1,4)
imagesc(z0,df,mz4.'); %caxis([-0.55 -0.45]);%caxis([-1 -0.9]);%
caxis([-0.45 -0.35]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,5)
imagesc(z0,df,(mz + mz2 + mz3 + mz4).'/4); %caxis([-0.55 -0.45]);%caxis([-1 -0.9]);%
caxis([-0.45 -0.35]);
xlabel('Position [cm]'); 
%title('mean');
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 


figure(1010102)
subplot(5,1,1)
imagesc(z0,b1scale,MZ.'); %caxis([-0.5 0.1]);%caxis([-1 -0.5]);%caxis([-0.4 -0.2]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,2)
imagesc(z0,b1scale,MZ2.'); caxis([-0.5 0.1]);%caxis([-1 -0.5]); % 
caxis([-0.4 -0.2]);
%xlabel('Position [cm]'); ylabel('B1 scale');
%title('RF2');
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,3)
imagesc(z0,b1scale,MZ3.'); caxis([-0.5 0.1]);%caxis([-1 -0.5]);%
caxis([-0.4 -0.2]); 
ylabel('B1 scale');
set(gca,'FontSize',12,'FontWeight','bold');
colorbar;

subplot(5,1,4)
imagesc(z0,b1scale,MZ4.'); caxis([-0.5 0.1]);%caxis([-1 -0.5]);%
caxis([-0.4 -0.2]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 

subplot(5,1,5)
imagesc(z0,b1scale,(MZ + MZ2 + MZ3 + MZ4 ).'/4); caxis([-0.5 0.1]);%caxis([-1 -0.5]);%
caxis([-0.4 -0.2]);
xlabel('Position [cm]'); set(gca,'FontSize',12,'FontWeight','bold');
colorbar; 



%%
% figure(3);
% subplot(1,3,1); imagesc(b1scale*max(abs(b1)), v, squeeze(mz(:,df==0,:))); axis square; colorbar; colormap jet; 
% xlabel('B1 (G)'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=0 Hz');
% set(gca,'FontSize',14);
% subplot(1,3,2); imagesc(df, b1scale*max(abs(b1)), squeeze(mz(v==0,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1 (G)');  caxis([-1 1]); title('Vc=0 cm/s');
% set(gca,'FontSize',14);
% subplot(1,3,3); imagesc(df, v, squeeze(mz(:,:,b1scale==1))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);
