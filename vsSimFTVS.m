%% ------------------------------------------------------------------------
% Demo FT-VS Pulse Generation and Bloch Simulations
%--------------------------------------------------------------------------
close all; clc; clear all;
addpath('useful'); addpath('rf_tools'); addpath('bloch_mex');

B1_val_hp = 1.4; %2; %10.20; % uT
B1_val_inv = 10.2; % uT
Grad_val = 1.45; % mT/cm

% durations (mS)
grad_ramp = 0.302; 
Tgap = 0.25; %0.1; % primarily for eddy current artifacts

% sampling periods
dtGz = 0.002;

% pulse design options: 
% note 1: other options not thoroughly tested, but may provide inspiration
% for other potential pulse designs
% note 2: sinc_weight=1 uses sub-pulse amplitudes that were optimized and 
% determined from previously running RUN_FTVS_OPT.m. See variables 
% 'hpscale' and 'rf_weight' in useful/gen_FVEVS to play with other subpulse
% amplitudes.
comp_hp = 0; % composite hard pulses: 0=off, 1=on
comp_inv = 0; % composite inversion: 0=off, 1=90-180-90, 2=90-240-90, 3=90-360-90
single_refocus = 0; %0=double-refocusing,1=single-refocusing
grad_var = 3;
sinc_weight = 1; 
num_sp = 9; % 5,9
% grad_var=1: 4 unipolars
% grad_var=3: 2 unipolars (between 180s)
% grad_var=4: 2 unipolars, (outer pairs) (similar label-control without gradient compensation)
% grad_var=5: 2 bipolars (outer pairs) (does nothing)
% grad_var=6: missing the dephasing unipolars
% grad_var=7: 4 bipolars  

if(single_refocus==1)
    [b1, gz, gz_flip, gz_off,inv_start,inv_dist,kv_locs] = gen_FVEVS_singleinv(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap,comp_inv,sinc_weight,num_sp);
else
    [b1, gz, gz_flip, gz_off,inv_start,inv_dist,kv_locs,hpscale] = gen_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap,comp_inv,grad_var,comp_hp,single_refocus,sinc_weight,num_sp);
end

% 0th and 1st moment calculations
t = [0:length(b1)-1]*dtGz;
[m0,t] = grad_moment(gz_flip,dtGz,0,inv_start,inv_dist,1);
[m1,~] = grad_moment(gz,dtGz,1,inv_start,inv_dist,1);
deltam0 = round((m0(kv_locs(2)) - m0(kv_locs(1)))*1000)/1000;
deltam1 = round((m1(kv_locs(2)) - m1(kv_locs(1)))*1000)/1000;
FOVm1 = round(1/deltam1/10)*10;
FOVm0 = round(1/deltam0*10)/10;

% plot kv-encoding and moments
figure(2);
subplot(3,1,3)
stem(m1(kv_locs),real(hpscale),'*r')
ylabel('B1 (G)'); xlabel('kv (s/cm)');
legend('kv-encoding');

subplot(3,1,2)
plot(t,m0,t(kv_locs),m0(kv_locs),'r*'); 
title(['0th moment (Grad flip): \Deltakz=' num2str(deltam0) ' 1/cm, \DeltaFOVz=' num2str(FOVm0) ' cm']); xlabel('time (ms)'); ylabel('kz (1/cm)');
legend('kz','kz-encoding');

subplot(3,1,1)
plot(t,m1,t(kv_locs),m1(kv_locs),'r*'); 
title(['1st moment: \Deltakv=' num2str(deltam1) ' s/cm , \DeltaFOVv=' num2str(FOVm1) ' cm/s']); xlabel('time (ms)'); ylabel('kv (s/cm)');
legend('kv','kv-encoding');

drawnow;
%% ------------------------------------------------------------------------
% Demo Bloch Simulation
%--------------------------------------------------------------------------
 
z0 = linspace(-8,8,7); % adjust for Bloch Sim of different z-values (cm)
v = linspace(-FOVm1,FOVm1,261); % adjust for different velocities of blood (cm/s)
v2 = linspace(-3,3,11); % adjust for differnt velocities of myocardium (cm/s)
df = linspace(-125,125,11); % adjust for differetn off-resonance values (Hz)
b1scale = linspace(.7,1.2,11); % adjust for different relative b1 scale

Nb1 = length(b1scale); Nv = length(v); Ndf = length(df); 
Nz = length(z0); Nv2 = length(v2);

% Intialize Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv,Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);
mx02 = zeros(Nv2, Ndf, Nb1); my02 = zeros(Nv2,Ndf, Nb1); mz02 = ones(Nv2, Ndf, Nb1);

% Arterial blood T1 and T2
%T1 = 1.664; 
%T2 = 0.2; %.150;

% no relaxation
T1 = 10; % 1.664;
T2 = 10; % 0.2;

% silicone oil
%T1 = 1.111;
%T2 = 227;

%T1 = 1.25; % calf,  1.932; % myocardium
%T2 = 0.05; % calf, 0.275; % myocardium
%% bloch simulations (space)
tic;
 for itr = 1:Nb1
     [~,~,mz_off(:,:,itr)] = bloch(b1*b1scale(itr),gz_off,dtGz*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
     [~,~,mz_flip(:,:,itr)] = bloch(b1*b1scale(itr),gz_flip,dtGz*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
 end

%% bloch simulations (velocity)

[~,~,mz2] = blochvb1(b1,gz,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
[~,~,mz22] = blochvb1(b1,gz,dtGz*1e-3,T1,T2,df,0,v2,b1scale,0,mx0,my0,mz0);

[~,~,mz2_flip] = blochvb1(b1,gz_flip,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
[~,~,mz2_off]  = blochvb1(b1,gz_off,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
toc;

%% Uncomment to check for VS-stripe artifact.
% relevant only if using a control pulse with non-zero gradients with same 
% amplitude (called control-uni, below). Note: make sure z0 is sampled
% finely!
%
% figure(3); clf;
% 
% toplot = squeeze(mz_flip(:,:,b1scale==1));
% subplot(3,2,1);
% imagesc(df,z0,toplot); 
%  colorbar; caxis([-1 1]);
% xlabel('\Delta f'); ylabel('z (cm)');
% title('Control-uni');
% 
% toplot = squeeze(mz2(:,:,b1scale==1));
% subplot(3,2,5);
% imagesc(df,v,toplot); 
% colorbar; caxis([-1 1]);
% xlabel('\Delta f'); ylabel('v (cm/s)');
% title('Label');
% 
% toplot = squeeze(mz_flip(:,df==0,:));
% subplot(3,2,3);
% imagesc(b1scale,z0,toplot); 
% colorbar; caxis([-1 1]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Control-uni');
% 
% toplot = squeeze(mz2_off(:,:,b1scale==1));
% subplot(3,2,2);
% imagesc(df,v,toplot); 
% colorbar; caxis([-1 1]);
% xlabel('\Delta f'); ylabel('v (cm/s)');
% title('Control-off');
% 
% toplot = squeeze(mz2(:,df==0,:));
% subplot(3,2,6);
% imagesc(b1scale,v,toplot); 
% colorbar; caxis([-1 1]);
% xlabel('B1 scale'); ylabel('v (cm/s)');
% title('Label');
% 
% toplot = squeeze(mz2_off(:,df==0,:));
% subplot(3,2,4);
% imagesc(b1scale,v,toplot); 
% colorbar; caxis([-1 1]);
% xlabel('B1 scale'); ylabel('v (cm/s)');
% title('Control-off');


%% ASL signal over myocardium and blood for different velocitites
b1scale2 = b1scale - 0.2; % RF pulse was designed with lower B1 than peak
% B1 of scanner and then scaled, to account for anticipated lower B1 in the
% heart

mz_0v = permute(repmat(squeeze(mz_off(z0==0,df==0,:)),[1 length(v)]),[2 1]);
mz_0v_df = repmat(squeeze(mz_off(z0==0,:,b1scale2==1)),[length(v) 1]);
mz_0v2 = permute(repmat(squeeze(mz2_off(z0==0,df==0,:)),[1 length(v2)]),[2 1]);
mz_0v_df2 = repmat(squeeze(mz2_off(z0==0,:,b1scale2==1)),[length(v2) 1]);

for opt = 1:2;
    figure(3+opt); clf;

x = b1scale2;   y = v;               z=mz_0v - squeeze(mz2(:,df==0,:));
x2 = b1scale2;  y2 = v2;             z2=mz_0v2 - squeeze(mz22(:,df==0,:));
x3 = b1scale2;  y3 = v(v>-130&v<-10); z3 = z(v>-130&v<-10,:);
x4 = b1scale2;  y4 = v(v>10&v<130);   z4 = z(v>10&v<130,:);
if (opt==2)
    x = df;  z=mz_0v_df - squeeze(mz2(:,:,b1scale2==1));
    x2 = df; z2=mz_0v_df2 - squeeze(mz22(:,:,b1scale2==1));
    x3 = df; z3 = z(v>-130&v<-10,:);
    x4 = df; z4 = z(v>10&v<130,:);
end
% create two axes
ax1 = axes;
imagesc(ax1,x,y,z)
ax2 = axes;
imagesc(ax2,x2,y2,z2)
ax3 = axes;
imagesc(ax3,x3,y3,z3)
ax4 = axes;
imagesc(ax4,x4,y4,z4)
% link them together
linkaxes([ax1,ax2,ax3,ax4])
% hide the top axes
xlabel(ax1,'B1 scale'); ylabel(ax1,'v (cm/s)')
set(ax1,'Ytick',-140:40:140);
title(ax1,'ASL Signal (\Deltaf=0)');
if (opt==2) 
    title(ax1,'ASL Signal (b1=1)');
    xlabel(ax1,'\Deltaf');
end

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];
% give each one its own colormap 
caxis(ax1,[-0.3 2]); colormap(ax1,'gray'); 
 caxis(ax2,[-0.02 0.02]); colormap(ax2,'cool');
caxis(ax3,[0 2]); colormap(ax3,'winter'); 
caxis(ax4,[0 2]); colormap(ax4,'winter'); 

% add colorbars and get everything lined up
set([ax1,ax2,ax3,ax4],'Position',[.18 .11 .685 .815]);
cb1 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 
cb2 = colorbar(ax3,'Position',[.04 .11 .0675 .815]); 
cbtitle1 = get(cb1,'Title');
set(cbtitle1,'String','myocardium');
set(cb1,'YTick',-0.02:0.01:0.02)
cbtitle2 = get(cb2,'Title');
set(cbtitle2,'String','blood');
set(cb2,'YTick',0:0.5:2)
% make font size bigger
set(ax1,'FontSize',16);
set(ax2,'FontSize',16);
set(ax3,'FontSize',16);
set(ax4,'FontSize',16);
end
vs_lab_blood_idx = (v>10&v<130)|(v>-130&v<-10);
blood_signal = repmat(mz2(v==0,:,:),[sum(vs_lab_blood_idx) 1 1])-mz2(vs_lab_blood_idx,:,:); % ideal to have  Mz=-1 (inversion) of blood

vs_lab_myo_idx = v>-3&v<3;
myo_signal = repmat(mz2(v==0,:,:),[sum(vs_lab_myo_idx) 1 1]) - mz2(vs_lab_myo_idx,:,:); % ideal to have no signal over myocardium

max_m = max(abs(myo_signal(:)));
min_m = min(abs(myo_signal(:)));
mean_m = mean(abs(myo_signal(:)));
std_m = std(abs(myo_signal(:)));
max_b = max(abs(blood_signal(:)));
min_b = min(abs(blood_signal(:)));
mean_b = mean(abs(blood_signal(:)));
std_b = std(abs(blood_signal(:)));
disp(['ASL signal (myo): max=' num2str(max_m) ' min=' num2str(min_m) ' mean=' num2str(mean_m) ' std=' num2str(std_m)]);
disp(['ASL signal (blood): max=' num2str(max_b) ' min=' num2str(min_b) ' mean=' num2str(mean_b) ' std=' num2str(std_b)]);

%% ASL signal over the myocardium
    
figure(6)
suptitle('ASL signal (myocardium)');
subplot(2,1,1)
imagesc(v2,df,mz_0v_df2 - squeeze(mz22(:,:,b1scale2==1))); axis square; axis off;
caxis([-0.02 0.02]); colorbar; colormap('cool');
ylabel('v'); xlabel('\Deltaf'); title('b1=1');

subplot(2,1,2)
imagesc(v2,b1,mz_0v2 - squeeze(mz22(:,df==0,:))); axis square; axis off;
caxis([-0.02 0.02]); colorbar; colormap('cool');
ylabel('v'); xlabel('b1'); title('\Deltaf=0');
%% velocity profile and k-v encoding before/after subpulse optimization
deltakv = m1(kv_locs(2))-m1(kv_locs(1));

% bloch simulations
v_blc = linspace(-FOVm1,FOVm1,211);

load('useful/mz2_before_mod_140cmsvfov.mat'); load('useful/mz2_after_mod_140cmsvfov.mat'); 
load('useful/hpscale_before_mod_140cmsvfov.mat'); load('useful/hpscale_after_mod_140cmsvfov.mat');

figure(7) 
subplot(2,1,1)
plot(v_blc,v_prof_blc_after,'k', v_blc,v_prof_blc_before,'b', 'LineWidth',2)
ylabel('Mz/M0'); xlabel('v (cm/s)');
xlim([-70 70]); legend('after','before');
title('Velocity profile (\Deltaf=0,b1=1)')
set(gca,'FontSize',16);
subplot(2,1,2)
h = stem([m1(kv_locs); m1(kv_locs)].',[real(hpscale_after); real(hpscale_before)].','*','LineWidth',2)
h(1).Color = 'black'; h(2).Color = 'blue';
ylabel('b1'); xlabel('kv (s/cm)');
title(['Kv-encoding (\Deltakv = ' num2str(round(deltakv*1000)/1000) ', FOVv= ' num2str(round(1./deltakv)) ')']);
set(gca,'FontSize',16);
xlim([-0.06 0.002]);