function J = opt_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp, hpscale_0)
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
    gap_bp = zeros(1,round(200*(0.002/dtGz)));
    grad_bp = 0*[decimate(grad,2) gap_bp gap_bp -decimate(grad,2)];
    %grad = 1.67*2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
    if (sinc_weight==0)
        grad = 1/3*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad)
    else
        grad = 1.4*[1:grad_ramp_res grad_ramp_res*ones(1,round(400*(0.002/dtGz))) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
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
hpscale = hpscale_0;
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
h = figure(1);
set(h, 'color', [1 1 1]);
subplot(3,1,1);
plot(tGz, rho, 'LineWidth', 1.5, 'color', 'blue');
grid on;
set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)] , 'YLim', [0 0.105]);
xlabel('time (ms)'); ylabel('|B1| (G)'); colormap jet;

subplot(3,1,2);
plot(tGz, phi, 'LineWidth', 1.5, 'color', 'red');
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)], 'YLim', [-pi pi], 'YTick', -pi:pi:pi,'YTickLabel', {'-\pi  ','0','\pi '});
xlabel('time (ms)'); ylabel('\phiB1 (rad)'); colormap jet;

subplot(3,1,3);
plot(tGz,gz, 'LineWidth', 1.5, 'color', 'green');
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)], 'YLim', [-6 6]);
xlabel('time (ms)'); ylabel('Gv (G/cm)'); colormap jet;

%%
%%%%%%%
% bloch simulations
%%%%%%%

v = linspace(-120,120,211);
df = linspace(-200,200,15); 
b1scale = linspace(.6,1.4,11);

Nb1 = length(b1scale); Nv = length(v); Ndf = length(df);

% Intialize Magnetization
mx0 = zeros(Nv, Ndf, Nb1); my0 = zeros(Nv,Ndf, Nb1); mz0 = ones(Nv, Ndf, Nb1);

% no relaxation
T1 = 10; %1.664;
T2 = 10; %0.2;

tic;
[~,~,mz2] = blochvb1(b1,gz,dtGz*1e-3,T1,T2,df,0,v,b1scale,0,mx0,my0,mz0);
toc;

%% Vs profiles

% control-label signal (blood)
vs_lab_blood_idx = (v>20&v<104)|(v>-104&v<-20);
blood_signal = mz2(vs_lab_blood_idx,:,:) - 1; % ideal to have Mz=-1 (inversion) of blood

% contorl-label signal (myocardium)
vs_lab_myo_idx = v>-2&v<2;
myo_signal = mz2(vs_lab_myo_idx,:,:) - mz2(v==0,:,:); % ideal to have no labeled signal over myo
   
lambda = 1;  %0.5;
J = (lambda)*mean(myo_signal(:).^2) % + (1-lambda)*1/mean(blood_signal(:).^2)

end

