function [b1, gz, gz_flip, gz_off, inv_start, inv_dist, kv_locs] = gen_FVEVS_singleinv(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv,sinc_weight,num_sp)
% Calculate FVE-VS pulse with single inversions!
%
% [b1, gz, gz_flip, gz_off, inv_start, inv_dist, kv_locs] = gen_FVEVS_singleinv(grad_ramp,Grad_val,B1_val_hp,B1_val_inv,Tgap, comp_inv,sinc_weight)
%
%
ftype='max';
if (num_sp==9)
   tbw=7; 
else
    tbw=4;
end
d1=0.01;
vfov_pres=100;
d2=0.001;
gam = 42.58; % MHz/T
dtGz = 0.002; % ms

grad_ramp_res = round(grad_ramp/dtGz);

Tgap_res = round(Tgap/dtGz);
gap = zeros(1,Tgap_res);

% durations
hp_dur = (9/num_sp)*20/360/gam/B1_val_hp*10^3;
inv_dur = 180/360/gam/B1_val_inv*10^3;
% B1_val_hp=hp_dur/inv_dur*B1_val_hp;
% hp_dur=inv_dur;

% res
hp_res = round(hp_dur/dtGz);
inv_res = round(inv_dur/dtGz);
vfov=0;
grad_bp_scale=1;
%% define pulse sub-components
while(vfov>vfov_pres+5 || vfov<vfov_pres-5)
    %bipolar gradient pair
    if(vfov>vfov_pres+5)
        grad_bp_scale=grad_bp_scale-0.01;
    else
        grad_bp_scale=grad_bp_scale+0.01;
    end
    if (sinc_weight==0)
        grad = [1:grad_ramp_res grad_ramp_res*ones(1,400) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        grad_bp = 1.06*[decimate(grad,2) zeros(1,inv_res-hp_res) -decimate(grad,2)];
    else
        grad = [1:grad_ramp_res grad_ramp_res*ones(1,436) grad_ramp_res:-1:1]/grad_ramp_res*Grad_val;
        grad_bp = grad_bp_scale*[decimate(grad,2)  -decimate(grad,2)];
        % zeros(1,inv_res-hp_res)
    end
    grad_wait = zeros(1, length(grad_bp));
    
    % sub-pulse envelope
    if (sinc_weight)
        % rf_weight = dzrf(9,4,'inv','max',0.1,0.0001);
        %rf_weight = dzrf(9,4,'inv','pm',0.1,0.001);
        rf_weight = dzrf(num_sp,tbw,'inv',ftype,d1,d2);
        hpscale = num_sp/sum(rf_weight)*rf_weight;
        %disp(hpscale)
    else
        hpscale = ones(1,num_sp);
    end
    
    % inversions
    if (comp_inv==1)
        inv = [ones(1,1/2*inv_res).*exp(1i*pi/2) zeros(1,162) ones(1,inv_res) zeros(1,162) ones(1,1/2*inv_res).*exp(1i*pi/2)]*B1_val_inv;
    elseif (comp_inv==2)
        inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*4/3).*exp(1i*pi/2) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
    elseif (comp_inv==3)
        inv = [ones(1,1/2*inv_res) zeros(1,162) ones(1,inv_res*2).*exp(1i*4*pi/6) zeros(1,162) ones(1,1/2*inv_res)]*B1_val_inv;
    else
        inv = [ones(1,inv_res)]*B1_val_inv;
    end
    inv_wait = zeros(1,length(inv));
    
    %% b1 calc
    phi = [0 0 pi pi pi 0 0 pi pi pi 0 0 0 pi pi 0]; % MLEV phase cycling pattern
    
    rf_sub = ones(1,hp_res)*B1_val_hp;
    rfsub_wait = zeros(1,length(rf_sub));
    
    vshift = 0;
    phi_sp = -1j*[0:pi:(num_sp-1)*pi];
    
    
    %% gv calc
    
    % full gradient waveform
    gz = [];
    for step = 1:num_sp-1
        if(mod(step,2))
            gz_ud = [rfsub_wait    gap   (-1)^step*grad_bp     gap  ...
                inv_wait gap (-1)^(step+1)*grad_bp    gap];
            gz = [gz gz_ud];
        else
            gz_ud = [rfsub_wait    gap   (-1)^(step)*grad_bp     gap  ...
                inv_wait gap (-1)^(step+1)*grad_bp    gap];
            gz = [gz gz_ud];
        end
        
    end
    gz = [gz rfsub_wait];
    
    gz_flip = -abs(gz);
    gz_off = zeros(size(gz));
    
    % Gradient samples and time fector
    NGz = length(gz);
    tGz = (0:dtGz:dtGz*(NGz-1));       %In milliseconds
    
    inv_start = length([rfsub_wait gap grad_wait gap zeros(1,length(inv_wait)/2)]);
    inv_dist  = length([inv_wait gap grad_wait gap rfsub_wait gap grad_wait gap]);
    
    kv_locs = ceil(length(rfsub_wait)/2 + [0:(num_sp-1)]*(inv_dist));
    
    %[moment_vec,t] = grad_moment(grad,dt,order,inv_start,inv_dist,b1scale)
    [m0,t] = grad_moment(gz,dtGz,0,inv_start,inv_dist,1);
    [m1,~] = grad_moment(gz,dtGz,1,inv_start,inv_dist,1);
    
    vfov=max(abs((1./diff(m1(kv_locs)))/42.58e6));
    
end



%FVE_bugfix this is to fix the bug with improper fv encoding
[~,ind]=sort((m1(kv_locs)));
hpnew=hpscale;
hpnew=hpscale(ind);

kv_sp=(m1(kv_locs));
figure(250);stem(linspace(min(kv_sp),max(kv_sp),length(hpscale)),(abs(hpnew)))


np_new=num_sp;
rf_weight = dzrf(np_new,tbw,'inv',ftype,d1,d2);
hpscale_new = np_new/sum(rf_weight)*rf_weight;

ind2=abs(2*(floor(np_new/2)-(0:np_new-1)))+1;
ind2(1:floor(length(ind2)/2))=ind2(1:floor(length(ind2)/2))-1;
hpscale_new=hpscale_new(ind);


Fnew=interp1(linspace(min(kv_sp),max(kv_sp),np_new),hpscale_new...
    ,linspace(min(kv_sp),max(kv_sp),np_new*200),'pchip');
hold on; plot(linspace(min(kv_sp),max(kv_sp),np_new*200),abs(Fnew));

Fpick=interp1((linspace(min(kv_sp),max(kv_sp),np_new*200)),...
    ((Fnew)),kv_sp,'linear');
if(strcmp(ftype,'max'))
    %  Fpick(1)=hpscale(1);
else if(strcmp(ftype,'min'))
        %   Fpick(end)=hpscale(end);
    end
end

hold on; stem(kv_sp,abs(Fpick));

hpscale=Fpick;

% rf_weight = dzrf(np_new,4,'inv','max',0.01,0.001);
% hpscale = np_new/sum(rf_weight)*rf_weight;
% ind2=abs(2*(floor(np_new/2)-(0:np_new-1)))+1;
% ind2(1:floor(length(ind2)/2))=ind2(1:floor(length(ind2)/2))-1;
%hpscale(ind2)=hpscale;


%hpscale=F(max(m1(kv_locs))-m1(kv_locs));


%hpscale(ind)=hpscale;

% endFVE_bugfix this is to fix the bug with improper fv encoding
% full b1 waveform
vfov=(1./diff(m1(kv_locs)))/42.58e6
delta_kv=(diff(m1(kv_locs)))

b1 = [];
for itr=1:num_sp-1
    b1_ud = [hpscale(itr)*rf_sub.*exp(phi_sp(itr))/1.0      gap    grad_wait  gap  ...
        inv.*exp(1i*phi(itr)).*exp(phi_sp(itr))  gap   grad_wait gap];
    b1 = [b1   b1_ud];
end
b1 = [b1 hpscale(end)*rf_sub.*exp(phi_sp(end))]*1e-2;

rho = abs(b1);
phi = angle(b1);

%% Plotting

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

