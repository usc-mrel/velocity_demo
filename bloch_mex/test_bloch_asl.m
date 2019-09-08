% Simple script to simulate the flow signal in SPASL
clear all;

tp=0.003;
t1=1250;
t2=40;
t2b=200;
t1b=1650;
delt=375;
bt=1150;
flow=2.18;
lambda=0.9;
lab_type=1;
eff=0.5;
mode_imag=0;
mode_asl=0;

opyres = 98;
TR=3.1;
TE=1.29;

A=306;B=79;
load xtime.mat;
im_dur = TR*opyres+10;


clear signal;
fas=10:10:70;
for ii=1:1:length(fas)
    fa=fas(ii);
    for numrr=1:3
t = (im_st_tstamp(1:31,1)*1e-3+ceil(im_dur/2)).';

[sig_spasl_C1, sig_spasl_C] = spaslSimImag([], im_st_tstamp(1:numrr:31,1)*1e-3, im_end_tstamp(1:numrr:31,1)*1e-3, t(7:numrr:31), t(1:numrr:end), flow, delt, bt, ...
    t1b, t1, t2, fa, TR, TE, eff, mode_imag, mode_asl, lab_type, A, B);

t = (im_st_tstamp(32:62,1)*1e-3+ceil(im_dur/2));
[sig_spasl_T1, sig_spasl_T] = spaslSimImag(lab_tstamp(32:end,1)*1e-3, im_st_tstamp(32:numrr:end,1)*1e-3, im_end_tstamp(32:numrr:end,1)*1e-3, t(7:numrr:31), t(1:numrr:end), flow, delt, bt, ...
    t1b, t1, t2, fa, TR, TE, eff, mode_imag, mode_asl, lab_type, A, B);

close all; plot(t(1:numrr:end)*1e-3,(sig_spasl_C.'-sig_spasl_T))
plot(t(1:numrr:end)*1e-3,sig_spasl_C,'r',t(1:numrr:end)*1e-3,sig_spasl_T);axis square;
signal(numrr,ii)=mean(sig_spasl_C1)-mean(sig_spasl_T1);
    end
end
%%
close all
load thickdata_global_seg.mat;
load xtime.mat;
Nrrs=1;
seg=[2:6];
data_set=4;
t_im_st=im_st_tstamp(1:31,1)*1e-3;
t_im_end=im_end_tstamp(1:31,1)*1e-3;
t = (im_st_tstamp(1:31,1)*1e-3+ceil(im_dur/2)).';
xx=mbf_data{data_set}.signalc(:,1:31);
mbf_data{data_set}.signalc(:,2:4:31)=xx(:,3:4:end);
mbf_data{data_set}.signalc(:,3:4:31)=xx(:,2:4:end);
modelfun_spasl_C = @(b,x) spaslSimImag([], t_im_st', t_im_end', ...
    x, x, flow, delt, bt, t1b, b(4), t2, b(3), TR, TE, eff, mode_imag, mode_asl, lab_type, b(2), b(1));
beta0 = [1;1;50;1000];
    params=lsqcurvefit(modelfun_spasl_C,beta0,t(Nrrs:end),nanmean(mbf_data{data_set}.signalc(seg,Nrrs:31),1),[0 0 10 0 ],[2000 2000 50 2000]);
    
A=params(2);B=params(1);fa=params(3);t1=params(4);
Nrrs=1
figure(3);plot(t(Nrrs:31)*1e-3,modelfun_spasl_C(params,t),t(Nrrs:31)*1e-3,nanmean(mbf_data{data_set}.signalc(seg,Nrrs:31),1)); axis([0 31 150 250]); hold on;

Nrrs=10;
xx=mbf_data{data_set}.signalt(:,1:31);
mbf_data{data_set}.signalt(:,2:4:31)=xx(:,3:4:end);
mbf_data{data_set}.signalt(:,3:4:31)=xx(:,2:4:end);
t = (im_st_tstamp(32:62,1)*1e-3+ceil(im_dur/2));  
t_im_st=im_st_tstamp(32:62,1)*1e-3;
t_im_end=im_end_tstamp(32:62,1)*1e-3;
lab_tstamp2=lab_tstamp(32:end,1)*1e-3;

clear b;
modelfun_spasl_T = @(b,x) spaslSimImag(lab_tstamp2, t_im_st, t_im_end, ...
    t, x, b(1), b(2), b(3), t1b, t1, t2, fa, TR, TE, eff, mode_imag, mode_asl, lab_type, A, B);
beta0 = [2;400;1250];
    params2=lsqcurvefit(modelfun_spasl_T,beta0,t(Nrrs:31),nanmean(mbf_data{data_set}.signalt(seg,Nrrs:31),1)',[0.5 100 100],[5 1000 1500]);    

    Nrrs=1;
plot(t(Nrrs:31)*1e-3,modelfun_spasl_T(params2,t),t(Nrrs:31)*1e-3,nanmean(mbf_data{data_set}.signalt(seg,Nrrs:31),1)) 
axis([0 31 150 250])

legend('control fitted','control acquired','tagged fitted','tagged acquired');

%% 
plot(t(Nrrs:31)*1e-3,nanmean(mbf_data{data_set}.signalc(seg,Nrrs:31),1),t(Nrrs:31)*1e-3,nanmean(mbf_data{data_set}.signalt(seg,Nrrs:31),1)) 