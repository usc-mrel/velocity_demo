%% ------------------------------------------------------------------------
% Demo VS Pulse Generation and Bloch Simulations using adiabatic BIR pulses
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

%% ---Pulse Generation---%

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
% Demo Pulse Bloch Simulation
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

tic;
mx = zeros(Nz,Ndf,Nb1); my=mx; mz=mx;
for itr = 1:Nb1
    itr
    [mx(:,:,itr),my(:,:,itr),mz(:,:,itr)] = bloch(b1*b1scale(itr),gz,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_flip(:,:,itr),my_flip(:,:,itr),mz_flip(:,:,itr)] = bloch(b1*b1scale(itr),gz_flip,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
    [mx_off(:,:,itr),my_off(:,:,itr),mz_off(:,:,itr)] = bloch(b1*b1scale(itr),gz_off,dtGz*1e-6,T1,T2,df,z0,0,mx0,my0,mz0);
end

% bloch simulations (velocity)

[mx2,my2,mz2] = blochvb1(b1,gz,dtGz*1e-6,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_flip] = blochvb1(b1,gz_flip,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
%[~,~,mz2_off]  = blochvb1(b1,gz_off,dtGz*1e-3,T1,T2,df,z0,v,b1scale,0,mx0,my0,mz0);
toc;

%% plot pulse and Bloch simulations
figure(2); clf;

subplot(3,1,1)
plot(t,gz,t,abs(b1)*10); xlabel('time (ms)'); 
legend('gz','b1');

toplot = squeeze(mz(:,df==0,:));
subplot(3,3,4);
imagesc(b1scale,z0,toplot); 
 colorbar; caxis([-1 1]);
xlabel('B1 scale'); ylabel('z (cm)');
title('Label');

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

   

