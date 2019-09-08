%% ------------------------------------------------------------------------
% Demo Pulse Generation
%--------------------------------------------------------------------------

% close all; clear all;
clear all;

% amplitudes
B1_val_hp = 3.2; %10.20; % uT
B1_val_inv = 10.20; % uT
Grad_val = 1.45; % mT/cm

% durations (mS)
grad_ramp = 0.302; 
Tgap = 0; %0.1; % primarily for eddy current artifacts

% sampling periods
dtGz = 0.002;

comp_hp = 0;
comp_inv = 0; % 90-180-90, 90-240-90, 90-360-90
single_refocus = 0;
grad_var = 3;
sinc_weight = 1;
% variation 1: 4 unipolars
% variation 2: 2 unipolars (middle 180s)
% variation 3: 2 unipolars, 2 bipolars
% variation 4: 2 unipolars, (outer pairs) (similar label-control without gradient compensation)
% variation 5: 2 bipolars (outer pairs) (does nothing)
% variation 6: missing the dephasing unipolars
% variation 7: 4 bipolars  

[b1, gz, gz_flip, gz_off,inv_start,inv_dist,kv_locs] = gen_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap,comp_inv,grad_var,comp_hp,single_refocus,sinc_weight);
%[b1, gz, gz_flip, gz_off,inv_start,inv_dist,kv_locs] = gen_FVEVS_singleinv(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap,comp_inv,sinc_weight);
t = [0:length(b1)-1]*dtGz;
% kv_locs(1)
% 0th and 1st moment calculations

% [moment_vec,t] = grad_moment(grad,dt,order,inv_start,inv_dist,b1scale)
[m0,t] = grad_moment(gz,dtGz,0,inv_start,inv_dist,1);
[m1,~] = grad_moment(gz,dtGz,1,inv_start,inv_dist,1);

figure(200);
subplot(3,1,1)
plot(t,gz,t,abs(b1)*10); xlabel('time (ms)'); 
legend('gz','b1');

subplot(3,1,2)
plot(t,m0,t(kv_locs),m0(kv_locs),'r*'); 
title('0th moment'); xlabel('time (ms)'); ylabel('kz (1/cm)');

subplot(3,1,3)
plot(t,m1,t(kv_locs),m1(kv_locs),'r*'); 
title('1st moment'); xlabel('time (ms)'); ylabel('kv (s/cm)');

drawnow;
%% ------------------------------------------------------------------------
% Demo Pulse Simulation
%--------------------------------------------------------------------------
 
z0 = linspace(-2.5,2.5,15);
v = linspace(-100,100,111);
df = linspace(-200,200,15); 
b1scale = linspace(.6,1.4,11);

Nb1 = length(b1scale); Nv = length(v); Ndf = length(df); Nz = length(z0);

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

%T1 = 1.25; % calf,  1.932; % myocardium
%T2 = 0.05; % calf, 0.275; % myocardium
%% bloch simulations (space)
tic;
mx = zeros(Nz,Ndf,Nb1); my=mx; mz=mx;
for itr = 1:Nb1
    itr
    [mx(:,:,itr),my(:,:,itr),mz(:,:,itr)] = bloch(b1*b1scale(itr),gz,dtGz*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_flip(:,:,itr),my_flip(:,:,itr),mz_flip(:,:,itr)] = bloch(b1*b1scale(itr),gz_flip,dtGz*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_off(:,:,itr),my_off(:,:,itr),mz_off(:,:,itr)] = bloch(b1*b1scale(itr),gz_off,dtGz*1e-3,T1,T2,df,z0,0,mx0,my0,mz0);
end

%% bloch simulations (velocity)

[mx2,my2,mz2] = blochvb1(b1,gz,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_flip] = blochvb1(b1,gz_flip,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_off]  = blochvb1(b1,gz_off,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
toc;

%% Vs stripe artifact
figure(2); clf;

% subplot(1,3,1);
% imagesc(df,z0,squeeze(mz(:,:,b1scale==1))); axis square; colorbar; caxis([-1 0.4]);
% xlabel('\Delta f (Hz)'); ylabel('z (cm)');

subplot(5,1,1)
plot(t,gz,t,abs(b1)*10); xlabel('time (ms)'); 
legend('gz','b1');

% toplot = squeeze(mz(:,df==0,:));
% subplot(5,2,2);
% imagesc(b1scale,z0,toplot); 
%  colorbar; caxis([-1 0.2]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Label');

toplot = squeeze(mz(:,df==0,:));
subplot(5,1,2);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.2]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label');

% toplot = squeeze(mz_flip(:,df==0,:));
% subplot(5,2,5);
% imagesc(b1scale,z0,toplot); 
%  colorbar; caxis([-1 0.2]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Control')

toplot = squeeze(mz_off(:,df==0,:));
subplot(5,1,3);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 0.2]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Control')

toplot = mz(:,df==0,:)-mz_off(:,df==0,:);
toplot = squeeze(toplot);
subplot(5,1,4);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-0.05 0.05]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label - Control');

% toplot = mz(:,df==0,:)-mz_off(:,df==0,:);
% toplot = squeeze(toplot);
% subplot(5,2,8);
% imagesc(b1scale,z0,toplot); 
%  colorbar; caxis([-0.2 0.2]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Label - Control');

% complex sum!? 
% Vs profile

%figure(3); clf;
subplot(5,2,9);
imagesc(df,v,squeeze(mz2(:,:,b1scale==1))); colorbar; caxis([-1 1]);
xlabel('\Delta f (Hz)'); ylabel('v (cm/s)');
subplot(5,2,10);
imagesc(b1scale,v,squeeze(mz2(:,df==0,:))); colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('v (cm/s)');
%%
idx=2:7; idx2=10:15;
mz_0v = permute(repmat(squeeze(mz_off(z0==0,df==0,idx)),[1 length(v)]),[2 1]);
mz_0v_df = repmat(squeeze(mz_off(z0==0,idx2,b1scale==1)),[length(v) 1]);
figure(3)
subplot(1,3,1)
plot(v,squeeze(mz2(:,df==0,b1scale==1)),[-15 15],[0 0],'r*'); axis square; 
xlabel('v (cm/s)'); ylabel('Mz'); axis tight; legend('Mz','v=15');
subplot(1,3,2)
plot(v,squeeze(mz2(:,df==0,idx)) - mz_0v); axis square; 
xlim([-2 2]); ylim([-0.02 0.02]);
xlabel('v (cm/s)'); title('Label-Control')
legend(['B1sc = ' num2str(b1scale(2))],['B1sc = ' num2str(b1scale(3))], ...
       ['B1sc = ' num2str(b1scale(4))],['B1sc = ' num2str(b1scale(5))], ...
       ['B1sc = ' num2str(b1scale(6))],['B1sc = ' num2str(b1scale(7))]);
subplot(1,3,3)
plot(v,squeeze(mz2(:,idx2,b1scale==1)) - mz_0v_df); axis square; 
xlim([-2 2]); ylim([-0.02 0.02]);
xlabel('v (cm/s)'); title('Label-Control')
legend(['\Deltaf = ' num2str(df(10))],['\Deltaf = ' num2str(df(11))], ...
       ['\Deltaf = ' num2str(df(12))],['\Deltaf = ' num2str(df(13))], ...
       ['\Deltaf = ' num2str(df(14))],['\Deltaf = ' num2str(df(15))]);


   

