function [b1, gz, gz_flip, gz_off, inv_start, inv_dist, kv_locs, rf_weight] = gen_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp)
% Generate double-inversion FVE-VS pulse
%  [b1, gz, gz_flip, gz_off, inv_start, inv_dist, kv_locs] = gen_FVEVS(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv, grad_var, comp_hp,single_refocus,sinc_weight,num_sp)


%num_sp=9;
gam = 42.58; % MHz/T
dtGz = 0.002; % ms

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
    gap_bp = zeros(1,100);
    grad_bp = 0*[decimate(grad,2) gap_bp gap_bp -decimate(grad,2)];
    %grad = 1.67*2*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
    if (sinc_weight==0)
        grad = 1/3*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        max(grad)
    else
         % 1.1 scale before
        grad = 0.85*[1:grad_ramp_res grad_ramp_res*ones(1,200) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
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
    %rf_weight = dzrf(9,4,'inv','max',0.1,0.01); 
    %hpscale = num_sp/sum(rf_weight)*rf_weight;
    %hpscale(2:2:end) = -hpscale(2:2:end);
    %hpscale = [-0.715172173301665,3.211163965528062,2.638765485759173,-0.576309590203551,6.723038271406336,3.915309605235483,-1.370090730629797,-0.478035545270557,1.124736670821164];
    hpscale = [0.0043015      1.9544     -1.2015     0.56621      6.8067      2.2401     0.34743       -1.05       2.173];
    rf_weight = hpscale*B1_val_hp*10e-2;  
    %hpscale = [0.0164 2.6874 -2.0084 0.7862 9.7065 3.3431 0.5917 -1.5175 3.0098]*2/3;

else
    hpscale = ones(1,num_sp);
end
phi_sp = -1j*[0:pi:(num_sp-1)*pi];

% b1 calc 
b1 = []; 
for step=1:num_sp-1   
    b1_ud = [hpscale(step)*rf_sub.*exp(phi_sp(step))          grad_wait  gap  ...
        inv.*exp(1i*phi(2*step-1)).*inv_amp(2*step-1)      grad_wait gap  rfsub_wait     grad_wait  gap  ...
        inv.*exp(1i*phi(2*step)).*inv_amp(2*step)          grad_wait gap ];
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
         gz_ud = [rfsub_wait      -grad_bp     gap  ...
                inv_wait  gap -grad        rfsub_wait grad      gap  ...
                inv_wait         grad_bp gap ];
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


gz_flip = -abs(gz);
gz_off = zeros(size(gz));

%Gradient samples and time fector
NGz = length(gz);
tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds

inv_start = length([rfsub_wait  grad_wait gap zeros(1,round(length(inv_wait)/2))]);
inv_dist  = length([inv_wait gap grad_wait rfsub_wait grad_wait gap]);

if (single_refocus)
    kv_locs = ceil(length(rfsub_wait)/2 + [0:num_sp-1]*(inv_dist));
else
    kv_locs = ceil(length(rfsub_wait)/2 + [0:num_sp-1]*(2*inv_dist));
end
invs = b1.*(real(b1)==0.102|real(b1)==-0.102);
subpulses = real(b1).*(real(b1)~=0.102&real(b1)~=-0.102);
subpulses(real(b1)==0.102|real(b1)==-0.102) = abs(invs(real(b1)==0.102|real(b1)==-0.102));
%---Plotting---%
h = figure(1); clf;
set(h, 'color', [1 1 1]);
% x = linspace(0,10);
% y = sin(3*x);
% yyaxis left
% plot(x,y)
% 
% z = sin(3*x).*exp(0.5*x);
% yyaxis right
% plot(x,z)
% ylim([-150 150])

left_color = [0 0 1];
right_color = [1 0 0];
set(h,'defaultAxesColorOrder',[left_color; right_color]);

%subplot(2,1,1); 
yyaxis left
plot(tGz, subpulses,'b', tGz,abs(invs),'g-','LineWidth', 1.5);
txt = ['+180x' ;'+180x' ;'-180x' ;'-180x'; '-180x' ;'+180x'; '+180x' ;'-180x'; '-180x' ;'-180x'; '+180x'; '+180x'; '+180x'; '-180x'; '-180x' ;'+180x'];
invs_diff = diff(abs(invs));
invs_vals = invs_diff==max(invs_diff);
invs_locs = find(invs_vals==1);
for itr=1:16
text(tGz(invs_locs(itr)) - 0.5,abs(invs(invs_locs(itr))) + 0.107,txt(itr,:)); hold on;
end
grid on;
set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)] , 'YLim', [-0.05 0.125]);
xlabel('time (ms)'); ylabel('B1 (G)'); colormap jet;

% subplot(3,1,1); 
% plot(tGz, subpulses,'b','LineWidth', 1.5);
% set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tGz(end), ...
%     'XLim', [0 tGz(end)] , 'YLim', [-0.05 0.125]);
% xlabel('time (ms)'); ylabel('B1 (G)'); colormap jet;
% legend('subpulse');
% 
% subplot(3,1,2); 
% plot(tGz,abs(invs),'g','LineWidth', 1.5);
% txt = ['+180x' ;'+180x' ;'-180x' ;'-180x'; '-180x' ;'+180x'; '+180x' ;'-180x'; '-180x' ;'-180x'; '+180x'; '+180x'; '+180x'; '-180x'; '-180x' ;'+180x'];
% invs_diff = diff(abs(invs));
% invs_vals = invs_diff==max(invs_diff);
% invs_locs = find(invs_vals==1);
% for itr=1:16
% text(tGz(invs_locs(itr)) - 0.5,abs(invs(invs_locs(itr))) + 0.11,txt(itr,:)); hold on;
% end
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight','bold','LineWidth', 2, 'XTick', 0:4:tGz(end), ...
%     'XLim', [0 tGz(end)] , 'YLim', [-0.05 0.125]);
% xlabel('time (ms)'); ylabel('B1 (G)'); colormap jet;
%legend('inversion');

%subplot(2,1,2);
yyaxis right
plot(tGz,gz,'r', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14,  'FontWeight','bold', 'LineWidth', 2, 'XTick', 0:4:tGz(end), ...
    'XLim', [0 tGz(end)], 'YLim', [-1.5 3.75]);
xlabel('time (ms)'); ylabel('Gz (G/cm)'); colormap jet;
legend('velocity selection','MLEV refocusing','Gz');


end

