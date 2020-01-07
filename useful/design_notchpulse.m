
clear all; 

flip_angle = pi/2;
TBW = 8;
gamma = 4.257;

sim_type = 0;


%%% slice_dz = TBW/(4.257*sum(gz_notch)*4/1000)   %% unit: cm
dz_target = 5; %1.6;
max_gz = 0.12; %0.32;
gz_a = [ 0:0.0601:max_gz];
gz_d = fliplr(gz_a);

sum_gz_core  = TBW / (gamma *4/1000 * dz_target );  % 469.8144
samp_gz_core = round( 0.5*( sum_gz_core - 4* sum(gz_a) ) /max_gz );

gz_one = max_gz*ones( 1, samp_gz_core);

% samp_tmp = 650;
% max_gz = 0.358;
% gz_a = [0:0.0601:max_gz];
% gz_d = fliplr(gz_a);
% gz_one = max_gz*ones(1,samp_tmp);
gz = [gz_a gz_one gz_d gz_a gz_one gz_d];
sum(gz)

samples = length(gz);


b = wsinc( TBW, samples);
%b = dzrf(samples,TBW, 'ex','pm');


bb = b/sum(b) * sin( 45*pi/180 );

% bb( (length(bb))/2 ) = bb( (length(bb))/2 ) - sin(45*pi/180 );
bb = [ bb(1: length(bb)/2)  -sin(45*pi/180)  bb(length(bb)/2+1: end) ];



% get alpha polynomial from beta polynomial 
aa = B2A(bb);

%aa = aa(1:end-1);  % last element is too high. Check it out later

% get RF pulse
%rf_sinc = ab2rf( aa, bb(1:end-1));
rf_sinc = ab2rf( aa, bb);

% rf_sinc = rf_sinc(2:end-1);
% rf_sinc = [ rf_sinc(2: ( length(rf_sinc)+1)/2 -1 )  rf_sinc( (length(rf_sinc)+1)/2 + 1: end-2) ];
rf_sinc=[ rf_sinc(1:( length(rf_sinc)-1)/2)  rf_sinc( (length(rf_sinc)-1)/2+2:end) ];



rf_sinc = flip_angle * real(rf_sinc)/ sum(real(rf_sinc) );
%rf_sinc = transpose( verse_tae( gz, rf_sinc ) );

rfs_sinc = rfscaleg( rf_sinc, samples*4/1000 );


samp_hard = 134;
rf_hard = ones(1,samp_hard);
rf_hard = rf_hard/samp_hard;
rf_hard = pi/2*rf_hard;

rfs_hard = rfscaleg(rf_hard,samp_hard*4/1000);


rf_notch = [rf_sinc(1,1:samples/2) -rf_hard rf_sinc(1,samples/2+1:samples)];
gz_notch = [gz(1,1:samples/2) zeros(1,samp_hard) gz(1,samples/2+1:samples)];

rf_notch_s = rfscaleg(rf_notch,length(rf_notch)*4/1000);

time_line = [4:4:4*length(rf_notch)]/1000;  %% [ms]
max_b1 = max(abs( rf_notch_s(:) ))
slice_dz = TBW/(4.257*sum(gz_notch)*4/1000)


% simulate the slice profile
zlist = -5:0.1:5;  % cm  (slice thickness is 1~cm)
dt = 4*10^-6;

if (sim_type==1)

    for f = 1:101 %-200:2:200
        f
        for k = 1:length(zlist)
            M = simpulse(rf_notch_s,gz_notch,dt,[0;0;1],zlist(k),4*f-204);
            Mz(k, f) = M(3);  
            Mr(k, f) = M(1) + j*M(2);
        end
    end

    save M_notch_new Mz Mr;


else
    %clear M,Mr;
    for k = 1:length(zlist)
        M(:,k) = simpulse(rf_notch_s,gz_notch,dt,[0;0;1],zlist(k),0);
    end
    Mr = M(1,:) + j*M(2,:);

    figure;
    subplot(311); plot(time_line,rf_notch_s);
    axis([min(time_line),max(time_line),min(rf_notch_s(:)),max(rf_notch_s(:))+0.01]);
    title('RF')
    subplot(312); plot(time_line,gz_notch);
    axis([min(time_line),max(time_line),min(gz_notch(:)),max(gz_notch(:))+0.2]);
    title('Gz')

    % slice profile
    subplot(313)
    %plot(zlist,M(1,:),'r-',zlist,M(2,:),'b-',zlist,M(3,:),'g-');
    plot(zlist,abs(Mr),'r-',zlist,M(3,:),'k-');grid on;
    axis( [ -5, 5, -0.5 1.2 ] );
    legend('|Mr|', 'Mz');
    title('Mz/Mo')


end

% save to HRF file
% saveHRF('notch_peak1100.hrf', rf_notch, 0.0017, length(rf_notch_s)*4, 16.006, 0,0, [],0,0, [],0,0, gz_notch,max(abs(gz_notch)),0);





