function [b1, gz, gz_flip, gz_off, inv_start, inv_dist, kv_locs] = gen_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

gam = 42.58; % MHz/T
dtGz = 0.002; % ms

hp_dur = 20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;

% res
hp_res = round(hp_dur/dtGz);
inv_res = round(inv_dur/dtGz);
grad_ramp_res = round(grad_ramp/dtGz);
  
Tgap_res = round(Tgap/dtGz);

% sub-waveforms
grad = [1:grad_ramp_res grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
if grad_var==2 || grad_var==3 || grad_var==5 || grad_var==7
    gap_bp = zeros(1,100);
    grad_bp = 1*[decimate(grad,2) gap_bp gap_bp -decimate(grad,2)];
    %grad = 1.67*2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
    if (sinc_weight==0)
        grad = 1/3*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad)
    else
        grad = 1/1.67*2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad)
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

if (sinc_weight)
    rf_weight = dzrf(9,3,'inv','max');
    %rf_weight = dzrf(9,4,'inv','pm',0.01,0.01); 
    hpscale = 9/sum(rf_weight)*rf_weight;
    disp(hpscale)
else
    hpscale = ones(1,9);
end

% b1 calc
b1 = []; 
for step=1:8    
    b1_ud = [hpscale(step)*rf_sub      gap    grad_wait  gap  ...
        inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1) gap     grad_wait gap  rfsub_wait  gap   grad_wait  gap  ...
        inv.*exp(1i*phi(2*step)).*inv_amp(2*step)         gap grad_wait gap ];
    b1 = [b1   b1_ud];
end
b1 = [b1 hpscale(end)*rf_sub]*1e-2;
rho = abs(b1);
phi = angle(b1);

% gv calc
gz = [];
if grad_var==1
    for step = 1:8       
        gz_ud = [rfsub_wait    gap   -grad     gap  ...
            inv_wait  gap grad    gap     rfsub_wait gap -grad      gap  ...
            inv_wait      gap    grad gap ];
        gz = [gz gz_ud];
        
    end
elseif grad_var==2 || grad_var==3
    for step = 1:8
         gz_ud = [rfsub_wait    gap   -grad_bp     gap  ...
                inv_wait  gap -grad    gap     rfsub_wait gap grad      gap  ...
                inv_wait      gap    grad_bp gap ];
         gz = [gz gz_ud];
    end    
elseif grad_var==4
    for step = 1:8
        gz_ud = [rfsub_wait gap   -grad      gap  ...
            inv_wait  gap grad_wait  gap     rfsub_wait gap grad_wait gap ...
            inv_wait  gap grad           gap ];
        gz = [gz gz_ud];
    end
elseif grad_var==5
    for step = 1:8
        gz_ud = [rfsub_wait gap   -grad_bp      gap  ...
            inv_wait  gap grad_wait  gap     rfsub_wait gap grad_wait gap ...
            inv_wait  gap grad_bp           gap ];
        gz = [gz gz_ud];
    end
elseif grad_var==6
    gz_area = 1/2*length(grad);
    gz_scale = gz_area/length(grad) + 1
    for step = 1:8
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
    for step = 1:8
         gz_ud = 8*[rfsub_wait    gap  grad_bp  gap  ...
                inv_wait  gap -grad_bp  gap rfsub_wait gap -grad_bp gap  ...
                inv_wait      gap    grad_bp gap ];
         gz = [gz gz_ud];
    end  
end
gz = [gz rfsub_wait];

if (single_refocus)
    b1 = [];
    for step=1:8
        b1_ud = [hpscale(step)*rf_sub      gap    grad_wait  gap  ...
            inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1)    gap   grad_wait  gap];
        b1 = [b1   b1_ud];
    end
    b1 = [b1 hpscale(end)*rf_sub]*1e-2;
    rho = abs(b1);
    phi = angle(b1);
    
    gz = [];
        for step = 1:8
            gz_ud = [rfsub_wait    gap   -grad     gap  ...
                inv_wait  gap    grad gap ];
            gz = [gz gz_ud];
            
        end
        gz = [gz rfsub_wait];        
end


gz_flip = -abs(gz);
gz_off = zeros(size(gz));

%Gradient samples and time fector
NGz = length(gz);
tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds

inv_start = length([rfsub_wait gap grad_wait gap zeros(1,length(inv_wait)/2)]);
inv_dist  = length([inv_wait gap grad_wait gap rfsub_wait gap grad_wait gap]);

if (single_refocus)
    kv_locs = length(rfsub_wait)/2 + [0:8]*(inv_dist);
else
    kv_locs = length(rfsub_wait)/2 + [0:8]*(2*inv_dist);
end

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

end

