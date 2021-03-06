function J = opt_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp, itr, lambda, hpscale_0)
% Generate double-inversion FVE-VS pulse
%  hpscale = opt_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp, hpscale_0)


%num_sp=9;
gam = 42.58; % MHz/T
dtGz = 0.02; % ms

hp_dur = (360/num_sp)/2/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;

% res
hp_res = round(hp_dur/dtGz);
inv_res = round(inv_dur/dtGz);
grad_ramp_res = round(grad_ramp/dtGz);
  
Tgap_res = round(Tgap/dtGz);

% sub-waveforms
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
if grad_var==2 || grad_var==3 || grad_var==5 || grad_var==7
    gap_bp = zeros(1,round(100*(0.002/dtGz)));
    grad_bp = 0*[decimate(grad,2) gap_bp gap_bp -decimate(grad,2)];
    %grad = 1.67*2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
    if (sinc_weight==0)
        grad = 1/3*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad)
    else
        %grad = 1.6; %1.1*[1:grad_ramp_res grad_ramp_res*ones(1,round(200*(0.002/dtGz))) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        %grad = 0.915*
        grad = 140/80*0.84*[1:grad_ramp_res grad_ramp_res*ones(1,round(200*(0.002/dtGz))) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad); 
    end
end
if grad_var==2 
    grad_bp = 0*grad_bp;
end

grad_wait = zeros(1, length(grad));

if (comp_hp)
    rf_sub = [ones(1,1/2*hp_res) ones(1,hp_res).*exp(1i*pi/2) ones(1,1/2*hp_res)]*B1_val_hp;
else
    rf_sub = ones(1,hp_res)*B1_val_hp;
end

if (comp_inv==1)
    inv = [ones(1,1/2*inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,inv_res) zeros(1,162) ones(1,1/2*inv_res).*exp(1i*pi/2)]*B1_val_inv;
elseif (comp_inv==2)
    inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*4/3).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
elseif (comp_inv==3)
    inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*2).*exp(1i*4*pi/6) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
else
    inv = [ones(1,inv_res)]*B1_val_inv;
end

rfsub_wait = zeros(1,length(rf_sub));
inv_wait = zeros(1,length(inv));
gap = zeros(1,Tgap_res);
inv_amp = ones(1,18);
phi = [0 0 pi pi pi 0 0 pi pi pi 0 0  0 pi pi 0];

% if (sinc_weight)
%     if (num_sp>8)
%     rf_weight = dzrf(num_sp,8,'inv','max',0.1,0.001);
%     else
%     rf_weight = dzrf(9,4,'inv','max',0.1,0.01); 
%     end
%     hpscale = num_sp/sum(rf_weight)*rf_weight;
%     disp(hpscale)
% else
%     hpscale = ones(1,num_sp);
% end
hpscale = real(hpscale_0);
phi_sp = -1j*[0:pi:(num_sp-1)*pi];

% b1 calc
b1 = []; 
for step=1:num_sp-1   
    b1_ud = [hpscale(step)*rf_sub.*exp(phi_sp(step))      gap    grad_wait  gap  ...
        inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  rfsub_wait  gap   grad_wait  gap  ...
        inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
    b1 = [b1   b1_ud];
end
b1 = [b1 hpscale(end)*rf_sub.*exp(phi_sp(step))]*1e-2;
%%% scale up B1 %%%
    %b1 = 1.2*b1;
%%%%%%%%%%%%%%%%%%%
rho = abs(b1);
phi = angle(b1);

% gv calc
gz = [];
if grad_var==1
    for step = 1:num_sp-1    
        gz_ud = [rfsub_wait    gap   -grad     gap  ...
            inv_wait  gap grad    gap     rfsub_wait gap -grad      gap  ...
            inv_wait      gap    grad gap ];
        gz = [gz gz_ud];
        
    end
elseif grad_var==2 || grad_var==3
    for step = 1:num_sp-1
         gz_ud = [rfsub_wait    gap   -grad_bp     gap  ...
                inv_wait  gap -grad    gap     rfsub_wait gap grad      gap  ...
                inv_wait      gap    grad_bp gap ];
         gz = [gz gz_ud];
    end    
elseif grad_var==4
    for step = 1:num_sp-1
        gz_ud = [rfsub_wait gap   -grad      gap  ...
            inv_wait  gap grad_wait  gap     rfsub_wait gap grad_wait gap ...
            inv_wait  gap grad           gap ];
        gz = [gz gz_ud];
    end
elseif grad_var==5
    for step = 1:num_sp-1
        gz_ud = [rfsub_wait gap   -grad_bp      gap  ...
            inv_wait  gap grad_wait  gap     rfsub_wait gap grad_wait gap ...
            inv_wait  gap grad_bp           gap ];
        gz = [gz gz_ud];
    end
elseif grad_var==6
    gz_area = 1/2*length(grad);
    gz_scale = gz_area/length(grad) + 1
    for step = 1:num_sp-1
        if step==1
            gz_ud = [rfsub_wait gap   grad_wait      gap  ...
                inv_wait  gap grad*gz_scale  gap   rfsub_wait gap -grad*gz_scale gap ...
                inv_wait  gap grad           gap ];
        elseif step==8            
            gz_ud = [rfsub_wait gap -grad        gap  ...
                inv_wait  gap grad*gz_scale  gap  rfsub_wait gap -grad*gz_scale gap ...
                inv_wait  gap grad_wait           gap ];
        else            
            gz_ud = [rfsub_wait gap -grad gap ...
                inv_wait  gap grad gap rfsub_wait gap -grad gap ...
                inv_wait  gap grad gap ];
        end
        gz = [gz gz_ud];
    end
elseif grad_var==7
    for step = 1:num_sp-1
         gz_ud = [rfsub_wait    gap  grad_bp  gap  ...
                inv_wait  gap -grad_bp  gap rfsub_wait gap -grad_bp gap  ...
                inv_wait      gap    grad_bp gap ];
         gz = [gz gz_ud];
    end  
end
gz = [gz rfsub_wait];

if (single_refocus)
    b1 = [];
    for step=1:num_sp-1
        b1_ud = [hpscale(step)*rf_sub      gap    grad_wait  gap  ...
            inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1)    gap   grad_wait  gap];
        b1 = [b1   b1_ud];
    end
    b1 = [b1 hpscale(end)*rf_sub]*1e-2;
    rho = abs(b1);
    phi = angle(b1);
    
    gz = [];
        for step = 1:num_sp-1
            gz_ud = [rfsub_wait    gap   -grad     gap  ...
                inv_wait  gap    grad gap ];
            gz = [gz gz_ud];
            
        end
        gz = [gz rfsub_wait];        
end


%Gradient samples and time fector
NGz = length(gz);
tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds

%---Plotting---%
invs = b1.*(b1==max(abs(b1))|b1==max(b1));
%---Plotting---%
h = figure(1);
set(h, 'color', [1 1 1]);
subplot(2,1,1);
plot(tGz, real(b1),'b', tGz,real(invs),'g','LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)] , 'YLim', [-0.125 0.125]);
xlabel('time (ms)'); ylabel('B1 (G)'); colormap jet;
legend('velocity selection','MLEV refocusing');

% subplot(3,1,2);
% plot(tGz, phi, 'LineWidth', 1.5, 'color', 'red');
% grid on;
% set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tGz(end), ...
%     'XLim', [0 tGz(end)], 'YLim', [-pi pi], 'YTick', -pi:pi:pi,'YTickLabel', {'-\pi  ','0','\pi '});
% xlabel('time (ms)'); ylabel('\phiB1 (rad)'); colormap jet;

subplot(2,1,2);
plot(tGz,gz,'r','LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)], 'YLim', [-6 6]);
xlabel('time (ms)'); ylabel('Gz (G/cm)'); colormap jet;

%%
%%%%%%%
% bloch simulations
%%%%%%%

v = linspace(-90,90,311);
df = linspace(-125,125,7); 
%b1scale = linspace(0.6,1,6);
b1scale = linspace(0.7,1.2,6);

Nb1 = length(b1scale); Nv = length(v); Ndf = length(df);

% Intialize Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv,Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);

% no relaxation
T1 = 10; %1.664;
T2 = 10; %0.2;

tic;
[~,~,mz2] = blochvb1(b1,gz,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
toc;
%%
b1scale2 = b1scale-0.2;
figure(2); %'units','normalized','outerposition',[0 0 0.5 0.5])
Figgy2 = gcf;
Figgy2.Position = [0, 0, 1500, 300];
subplot(1,5,1)
plot(v,squeeze(mz2(:,df==0,b1scale2==1)),[-15 15],[0 0],'r*'); 
xlabel('v (cm/s)'); ylabel('Mz'); axis tight; legend('Mz','v=15');
mz_0v = permute(repmat(squeeze(mz2(v==0,df==0,:)),[1 length(v)]),[2 1]);
mz_0v_df = repmat(squeeze(mz2(v==0,:,b1scale==1)),[length(v) 1]);
subplot(1,5,2)
imagesc(b1scale2,v, mz_0v - squeeze(mz2(:,df==0,:))); 
ylim([-3 3]); % ylim([-0.02 0.02]); 
caxis([-0.02 0.02]);
colorbar; title('Control-Label (\Deltaf=0)');
xlabel('B1scale'); ylabel('v (cm/s)');

subplot(1,5,3)
imagesc(df,v, mz_0v_df - squeeze(mz2(:,:,b1scale==1))); 
ylim([-3 3]); % ylim([-0.02 0.02]); 
caxis([-0.02 0.02]); 
colorbar; title('Control-Label (B1scale=1)');
xlabel('\Delta f (Hz)'); ylabel('v (cm/s)');

subplot(1,5,4);
imagesc(b1scale2,v,mz_0v - squeeze(mz2(:,df==0,:))); colorbar; caxis([0 2]);
xlabel('B1scale'); ylabel('v (cm/s)');
title('Control-Label (\Deltaf=0)');

subplot(1,5,5);
imagesc(df,v,mz_0v_df - squeeze(mz2(:,:,b1scale2==1))); colorbar; caxis([0 2]);
xlabel('\Delta f (Hz)'); ylabel('v (cm/s)');
title('Control-Label (B1scale=1)');
%% Vs profiles

% control-label signal (blood)
%blood_idx = (v>10&v<104)|(v>-104&v<-10);
%blood_signal = repmat(mz2(v==0,:,:),[sum(blood_idx) 1 1])-mz2(blood_idx,:,:); % ideal to have Mz=-1 (inversion) of blood
vs_lab_blood_idx = (v>10&v<70)|(v>-70&v<-10);
blood_signal = repmat(mz2(v==0,:,:),[sum(vs_lab_blood_idx) 1 1])-mz2(vs_lab_blood_idx,:,:); % ideal to have  Mz=-1 (inversion) of blood

vs_lab_myo_idx = v>-2&v<2;
myo_signal = repmat(mz2(v==0,:,:),[sum(vs_lab_myo_idx) 1 1]) - mz2(vs_lab_myo_idx,:,:); % ideal to have no signal over myocardium

penalty = max(abs(myo_signal(:)));

%lambda = 0.5;  %0.5;
%J = (lambda)*mean(abs(myo_signal(:))) - (mean(abs(blood_signal(:))));
%J = (1)*mean(abs(myo_signal(:))) - (mean(abs(blood_signal(:)))) + ...
%    (lambda)*penalty;
J =  - (mean(abs(blood_signal(:)))) + ...
    (lambda)*penalty;
% look for weighting with mean???
% % weight static tissue more heavily - can't have spurious labeling
% 
% stdGauss = 2*Vc;
% gauss = @(diffX) exp(-((diffX).^2./(2*stdGauss^2)));
% w = gauss(v).';
% w = repmat(w, [1 Ndf]);
% % w = repmat(w, [1 Ndf*Nb1]);
% w = w.*(w.');
% w = repmat(w, [1 1 Nb1]);
% 
% J = sum(mean(w.*(mz-mz_ref).^2,2))/Nv;

results.max_myo_signal = max(myo_signal(:));
results.mean_myo_signal = mean(myo_signal(:));
results.std_myo_signal = std(myo_signal(:));
results.max_blood_signal = max(blood_signal(:));
results.mean_blood_signal = mean(blood_signal(:));
results.std_blood_signal = std(blood_signal(:));
results.lambda = lambda; 
results.rfweight = hpscale*B1_val_hp; 
if lambda < 1
    filename = ['./results/resultsmat_itr_' num2str(itr) '_lambda_'];
else
    filename = ['./results/resultsmat_itr_' num2str(itr) '_lambda_' num2str(round(lambda))];
end
save(filename,'results');

end

