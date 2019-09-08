%% ------------------------------------------------------------------------
% Demo Pulse Generation
%--------------------------------------------------------------------------

% close all; clear all;
clear all;

%---BIR-4 Parameters
% omega = 18.430;		
% zeta = 10.46; 		
% kappa = atan2(87.24,1); 		
% RF_width = 1816;			

%---BIR-8 Parameters
%kappa = atan2(60,1);
%zeta = 15;
%omega = 39.8; %kHz
%RF_width = 1500;

kappa = atan2(62.96,1);
zeta = 20.50;
omega = 21.60; %kHz
RF_width = 2236;

%Tried Optimizattion
Vmax = 30;       
b1max = .16;
Gmax = 10;
ramp_time = 500;

grad_wait = 748;
type = 8;

%% TJAOs TRY ---Pulse Generation---%

[b1, gz, dt] = genVSBIR(zeta, kappa, omega, Vmax, 'type', type, ...
    'Smax', 7.5, 'dt', [4 4], 'b1max', b1max, 'sym', 1, 'ramp_time', ramp_time, 'grad_wait', grad_wait, 'RF_width', RF_width, 'Gmax', Gmax); %'ramp_time', 1000,
t = [0:length(b1)-1]*dt(1);

gz_flip = -abs(gz);
gz_off = zeros(size(gz));

NRF = length(b1);
NGz = length(gz);
dtRF = dt(1); dtGz = dt(2);

rho = abs(b1);
phi = angle(b1);

b10 = b1;

%% VL TRIES!

% opt = 4; % hard (1) or sinc (2) pulse envelope
% 
% gam = 42.58; % MHz/T
% 
% % amplitudes
% B1_val_hp = 3.13; %%3.88; % %5.2; % uT
% B1_val_inv = 11.5; %11.5; % 11.5?? % uT
% Grad_val = 2.2; %1.45; % mT/cm
% 
% % durations (mS)
% hp_dur = 20/360/gam/B1_val_hp*10^3;
% inv_dur = 180/360/gam/B1_val_inv*10^3;
% grad_ramp = 0.302; 
% Tgap = 0; %.34;
% 
% % sampling periods
% dtRF = 0.002*51;
% dtGz = 0.002*51;
% 
% % res
% hp_res = round(hp_dur/dtRF);
% inv_res = round(inv_dur/dtRF);
% grad_ramp_res = round(grad_ramp/dtGz);
% if opt == 1
%     %grad_ramp_res = round(.45*grad_ramp_res);
%     grad_ramp_res = round(grad_ramp_res);
% elseif opt == 2
%     grad_ramp_res = round(9/8*grad_ramp_res);
% elseif opt == 3
%     grad_ramp_res = 1/2*round(9/8*grad_ramp_res);
% elseif opt == 4
%     grad_ramp_res = round(grad_ramp_res); %round(9/8*grad_ramp_res);
%     
% end
% Tgap_res = round(Tgap/dtRF);
% 
% % sub-waveforms
% grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
% grad2 = [1:grad_ramp_res grad_ramp_res grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
% grad_wait = zeros(1, length(grad));
% 
% hp = ones(1,hp_res)*B1_val_hp;
% inv = [ones(1,1/2*inv_res) ones(1,inv_res).*exp(1i*pi/2) ones(1,1/2*inv_res)]*B1_val_inv;
% %inv = [ones(1,inv_res)]*B1_val_inv;
% hp_wait = zeros(1,hp_res);
% inv_wait = zeros(1,2*inv_res);
% %inv_wait = zeros(1,inv_res);
% gap = zeros(1,Tgap_res);
% inv_amp = ones(1,18);
% phi = [0 0 pi pi pi 0 0 pi] + 0;
% phi = [phi pi pi 0 0  0 pi pi 0];
% phi_hp = zeros(1,9); 
% %%
% %load('est');
% %phi_hp = est;
% %phi = phi + (pi - [est(1) est(1) est(2) est(2) est(3) est(3) est(4) est(4) ...
% %    est(5) est(5) est(6) est(6) est(7) est(7) est(8) est(8)]) ;
%     
%                     
% 
% %load('est_opts') % <- mean Generate initial guesses!?
% %load('est_opts2') % <- mean w/0 b1 scaling...
% %load('est_opts3') % <- max w/0 b1 scaling...
% %load('est_opts4') % <- max, over inv w/0 b1 scaling...
% %load('est_opts5') % <- mean, with inv phase equal to sub-pulse phase...
% 
% 
% for itr = 5; %1:size(est_opts,1)
%     
%  %  phi_hp = est_opts(itr,1:9);
%  %  phi = est_opts(itr,10:25);
%     
% 
% %% optmizing over B1 only (no df)
% % phi_hp =  [-0.2160   -0.4066   -1.9293 0.0047    1.8630   -2.3865 ...
% %             2.4916    0.7299   -0.8933];
% % phi =   [ -0.4273   -0.5228    2.0189 0.9676    0.5584    0.3689 ...
% %       -0.8347    1.6182    0.8704 -1.4080   -1.9820    0.5488 ...
% %       0.0155    2.0742    1.5645 -2.0424]; % best!
%     
% %% optimizing over df only (using only B1 as initial guess) - terrible
% 
% %% 
% if opt == 1;
%     hpscale = ones(1,9);
% elseif opt == 2;
%     hp = ones(1,2*hp_res)*1/2*B1_val_hp;
%     hp_wait = zeros(1,2*hp_res);
%     hpscale = 9/sum(dzrf(9,4,'st'))*dzrf(9,4,'st'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
% elseif opt == 3;
%     hp = ones(1,2*hp_res)*1/2*B1_val_hp;
%     hp_wait = zeros(1,2*hp_res);
%     hpscale = 9/sum(dzrf(9,2,'st'))*dzrf(9,2,'st'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
% elseif opt==4;
%     hp = ones(1,hp_res)*B1_val_hp; %ones(1,2*hp_res)*1/2*B1_val_hp;
%     %hp = 10*hp*sech(0.1*hp_res).^(1 + 1i*0.1);
%     hp_wait = zeros(1,hp_res);
%     hpscale = 9/sum(dzrf(9,2,'inv','max'))*dzrf(9,2,'inv','max'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
% 
% end
% b1 = []; gz = [];
% 
% %%
% %rf_sub = hp;
% %hpscale = ones(1,9); 
% % beta = 300;
% % mu = 2;
% % pwsech = 25000;
% % dtus = 116.5;
% % cexp = 0.01;
% % ssinvthick = 4;
% % b1_scale = 1;
% % [rf_im,grf,bw] = gensech(beta,  mu,  pwsech,  dtus,  cexp, ssinvthick, b1_scale);
% % rf_sub = rf_im*100; 
% % hp_wait = zeros(1,length(rf_sub));
% 
% % hpscale = ones(1,9); 
% % beta = 300;
% % mu = 2;
% % pwsech = 25000;
% % dtus = 116.5/9;
% % cexp = 0.01;
% % ssinvthick = 4;
% % b1_scale = 1;
% % [rf_im,grf,bw] = gensech(beta,  mu,  pwsech,  dtus,  cexp, ssinvthick, b1_scale);
% % [rf_im2,~,~] = gensech(beta,  mu,  22718*2,  2, cexp, ssinvthick, b1_scale);
% % rf_im = ones(1,length(rf_im))*B1_val_inv*1e-2;
% % %rf_sub = rf_im*100; 
% % for step = 1:8
% %         rf_sub = abs(rf_im(1 + (step-1)*hp_res:step*hp_res))*100;
% %         rf_wait = zeros(1,length(rf_sub));
% %         b1 = [b1       rf_sub          grad_wait  gap  ...
% %             inv gap     grad_wait  rf_wait    grad_wait  gap  ...
% %             inv         gap grad_wait ];
% %         gz = [gz       rf_wait       -grad      gap  ...
% %             inv_wait  gap grad          rf_wait -grad      gap  ...
% %             inv_wait           gap grad ];
% % end
% % rf_sub = rf_im(1+8*hp_res:end)*100;
% % b1 = [b1 rf_sub ]*1e-2; 
% % b1 = b1.*rf_im2(1:length(b1)); 
% 
% %% double refocusing
% rf_sub = hp;
% for step = 1:8
%         b1 = [b1       hpscale(step)*rf_sub.*exp(1i*phi_hp(step))          grad_wait  gap  ...
%             inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait  zeros(1,length(rf_sub))    grad_wait  gap  ...
%             inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait ];
% 
%         gz = [gz       zeros(1,length(rf_sub))       -grad     gap  ...
%             inv_wait  gap grad         zeros(1,length(rf_sub)) -grad      gap  ...
%             inv_wait           gap grad ]; 
%         
%      
% end
% b1 = [b1 hpscale(end)*rf_sub.*exp(1i*phi_hp(end)) ]*1e-2; % add zeros(1,9000) to make add 1 sec after pulse!
% sc2 = 1; %0.5; %0.4;
% gz = sc2*[gz zeros(1,length(rf_sub)) ]; 
% %% single refocusing
% % B1_val_inv = 12.04; % 11.5?? % uT
% % inv_dur = 180/360/gam/B1_val_inv*10^3;
% % inv_res = round(inv_dur/dtRF);
% % inv = ones(1,inv_res).*B1_val_inv;
% % 
% % grad_ramp_res = 2*grad_ramp_res;
% % % sub-waveforms
% % grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
% % grad_wait = zeros(1, length(grad));
% % hpscale = 6/sum(dzrf(6,2,'inv','max'))*dzrf(6,2,'inv','max'); %9/sum(sinc((-4:4)/3))*sinc((-4:4)/3);
% % phi = [0 0 pi pi 0];
% % rf_sub = hp;
% % %gap = zeros(1,10);
% % for step = 1:5
% %         b1 = [b1       hpscale(step)*rf_sub          grad_wait  gap  ...
% %             [ones(1,1/2*inv_res)*B1_val_inv inv.*exp(1i*pi/2) ones(1,1/2*inv_res)*B1_val_inv].*exp(1i*phi(step)) gap     grad_wait];
% %         gz = [gz       zeros(1,length(rf_sub))       grad*(-1)^step     gap  ...
% %             inv_wait inv_wait  gap grad*(-1)^step]; 
% % end
% % b1 = [b1 hpscale(end)*rf_sub ]*1e-2; % add zeros(1,9000) to make add 1 sec after pulse!
% % sc2 = 3; %0.5; %0.4;
% % gz = sc2*[gz zeros(1,length(rf_sub)) ]; 
% %%
% NRF = length(b1);
% NGz = length(gz);
% 
% rho = abs(b1);
% phi = angle(b1);

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

z0 = linspace(-2.5,2.5,5);
v = linspace(-80,80,311);
df = linspace(-200,200,13); 
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

%%
% figure(2);
% subplot(1,3,1); imagesc(b1scale*max(abs(b1)), v, squeeze(mz(:,df==0,:))); axis square; colorbar; colormap jet; 
% xlabel('B1 (G)'); ylabel('v (cm/sec)');  caxis([-1 1]); title('\Deltaf=0 Hz');
% set(gca,'FontSize',14);
% subplot(1,3,2); imagesc(df, b1scale*max(abs(b1)), squeeze(mz(v==0,:,:)).'); axis square; colorbar; colormap jet;
% xlabel('\Delta f (Hz)'); ylabel('B1 (G)');  caxis([-1 1]); title('Vc=0 cm/s');
% set(gca,'FontSize',14);
% subplot(1,3,3); imagesc(df, v, squeeze(mz(:,:,b1scale==1))); axis square; colorbar;
% xlabel('\Delta f (Hz)'); ylabel('v (cm/sec)');  caxis([-1 1]);colormap jet; title(['B1max=' num2str(max(abs(b1))) ' G']);
% set(gca,'FontSize',14);
% 
% sc = 0.8:0.1:1.3;
% figure(101);
% for itr2 = 1:6
% subplot(1,6,itr2)
% plot(v,squeeze(mz(:,df==-50,b1scale==sc(itr2))),'b', ...
%      v,squeeze(mz(:,df==-25,b1scale==sc(itr2))),'c', ...
%      v,squeeze(mz(:,df==0,b1scale==sc(itr2))),'k', ...
%      v,squeeze(mz(:,df==25,b1scale==sc(itr2))),'m', ...
%      v,squeeze(mz(:,df==50,b1scale==sc(itr2))),'y','LineWidth',2);
% xlabel('velocity (cm/s)'); ylabel('Mz');
% set(gca,'FontSize',14,'FontWeight','bold');
% if (itr2==6) 
%     legend('\Deltaf=-50', ...
%        '\Deltaf=-25', ...
%        '\Deltaf=0', ...
%        '\Deltaf=25', ...
%        '\Deltaf=50','Location','best');
% end
% title(['b1scale=' num2str(sc(itr2)) ' ']);
% end
% 
% 
% % figure(10101);
% % plot(v,squeeze(mz(:,df==0,b1scale==1)),'k',...
% %      v,squeeze(mz(:,df==0,b1scale==1.2)),'g', ...
% %      v,squeeze(mz(:,df==0,b1scale==0.8)),'b','LineWidth',2);
% % xlabel('velocity (cm/s)'); ylabel('Mz');
% % set(gca,'FontSize',14,'FontWeight','bold');
% % legend(['B1=' num2str(1*max(abs(b1))) 'G'], ...
% %        ['B1=' num2str(1.2*max(abs(b1))) 'G'], ...
% %        ['B1=' num2str(0.8*max(abs(b1))) 'G']); 
% % title('\Deltaf=0 Hz');
% 
% % filename_rho = 'VSVL.rho';
% % filename_theta = 'VSVL.theta';
% % filename_gz = 'VSVL.gz';
% % signa(abs(b1),filename_rho);
% % signa(angle(b1),filename_theta);
% % signa(gz,filename_gz);
% %     
% % idx = find((squeeze(mz(:,df==0,b1scale==1))<0.8)&(squeeze(mz(:,df==0,b1scale==1))>-0.8));
% % vals = v(idx);
% % vals = vals(vals<0)
% % v_prof_width = max(vals) - min(vals)
% % v_cut_idx = find(abs(squeeze(mz(:,df==0,b1scale==1)))==min(abs(squeeze(mz(:,df==0,b1scale==1)))));
% % v_cut = v(v_cut_idx)
% %end
% %% hyperbolic secant
% % beta = 300;
% % mu = 2;
% % pwsech = 25000;
% % dtus = 116.5;
% % cexp = 0.01;
% % ssinvthick = 4;
% % b1_scale = 1;
% % [rf_im,grf,bw] = gensech(beta,  mu,  pwsech,  dtus,  cexp, ssinvthick, b1_scale);
% % figure(10110)
% % subplot(3,1,1)
% % plot(abs(rf_im))
% % subplot(3,1,2)
% % plot(angle(rf_im))
% % subplot(3,1,3)
% % plot(grf)

tic;
mx = zeros(Nz,Ndf,Nb1); my=mx; mz=mx;
for itr = 1:Nb1
    itr
    [mx(:,:,itr),my(:,:,itr),mz(:,:,itr)] = bloch(b1*b1scale(itr),gz,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_flip(:,:,itr),my_flip(:,:,itr),mz_flip(:,:,itr)] = bloch(b1*b1scale(itr),gz_flip,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_off(:,:,itr),my_off(:,:,itr),mz_off(:,:,itr)] = bloch(b1*b1scale(itr),gz_off,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
end

%% bloch simulations (velocity)

[mx2,my2,mz2] = blochvb1(b1,gz,dtGz*1e-6,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_flip] = blochvb1(b1,gz_flip,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_off]  = blochvb1(b1,gz_off,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
toc;

%% Vs stripe artifact
figure(2); clf;

% subplot(1,3,1);
% imagesc(df,z0,squeeze(mz(:,:,b1scale==1))); axis square; colorbar; caxis([-1 0.4]);
% xlabel('\Delta f (Hz)'); ylabel('z (cm)');

subplot(3,1,1)
plot(t,gz,t,abs(b1)*10); xlabel('time (ms)'); 
legend('gz','b1');

% toplot = squeeze(mz(:,df==0,:));
% subplot(5,2,2);
% imagesc(b1scale,z0,toplot); 
%  colorbar; caxis([-1 0.2]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Label');

toplot = squeeze(mz(:,df==0,:));
subplot(3,3,4);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label');

% toplot = squeeze(mz_flip(:,df==0,:));
% subplot(5,2,5);
% imagesc(b1scale,z0,toplot); 
%  colorbar; caxis([-1 0.2]);
% xlabel('B1 scale'); ylabel('z (cm)');
% title('Control')

toplot = squeeze(mz_off(:,df==0,:));
subplot(3,3,5);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Control')

toplot = mz(:,df==0,:)-mz_off(:,df==0,:);
toplot = squeeze(toplot);
subplot(3,3,6);
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
subplot(3,2,5);
imagesc(df,v,squeeze(mz2(:,:,b1scale==1))); colorbar; caxis([-1 1]);
xlabel('\Delta f (Hz)'); ylabel('v (cm/s)');
subplot(3,2,6);
imagesc(b1scale,v,squeeze(mz2(:,df==0,:))); colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('v (cm/s)');
%%
idx=1:2:11; idx2=3:2:13;
mz_0v = permute(repmat(squeeze(mz_off(z0==0,df==0,idx)),[1 length(v)]),[2 1]);
mz_0v_df = repmat(squeeze(mz_off(z0==0,idx2,b1scale==1)),[length(v) 1]);
figure(3); clf;
subplot(1,3,1)
plot(v,squeeze(mz2(:,df==0,b1scale==1)),[-15 15],[0 0],'r*'); axis square; 
xlabel('v (cm/s)'); ylabel('Mz'); axis tight; legend('Mz','v=15');
subplot(1,3,2)
plot(v,squeeze(mz2(:,df==0,idx)) - mz_0v); axis square; 
xlim([-6 6]); ylim([-0.04 0.04]);
xlabel('v (cm/s)'); title('Label-Control')
legend(['B1sc = ' num2str(b1scale(1))],['B1sc = ' num2str(b1scale(3))], ...
       ['B1sc = ' num2str(b1scale(5))],['B1sc = ' num2str(b1scale(7))], ...
       ['B1sc = ' num2str(b1scale(9))],['B1sc = ' num2str(b1scale(11))],'AutoUpdate','off');
hold on; 
plot([-2 2 2 -2 -2],[-0.02 -0.02 0.02 0.02 -0.02],'k--');

subplot(1,3,3)
plot(v,squeeze(mz2(:,idx2,b1scale==1)) - mz_0v_df); axis square; 
xlim([-6 6]); ylim([-0.04 0.04]);
xlabel('v (cm/s)'); title('Label-Control')
legend(['\Deltaf = ' num2str(df(3))],['\Deltaf = ' num2str(df(5))], ...
       ['\Deltaf = ' num2str(df(7))],['\Deltaf = ' num2str(df(9))], ...
       ['\Deltaf = ' num2str(df(11))],['\Deltaf = ' num2str(df(13))],'AutoUpdate','off');
hold on; 
plot([-2 2 2 -2 -2],[-0.02 -0.02 0.02 0.02 -0.02],'k--');

   

