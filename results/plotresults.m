% L curve - myocardial signal vs. blood signal

%D = dir('./old/resultsmat*');
D = dir('./old_09_24_final_vFOV_140cms/resultsmat*'); % with penalty term |S(v_b)|>0.02
%D = dir('./new_penalty/resultsmat*'); % with penalty term max|S(v_b)|-0.02

myo_signal = zeros(1,length(D)); blood_signal=myo_signal;
lambda=myo_signal; rfweight=zeros(length(D),9);

% max_m  = [0.33037 0.07907 0.054481     0.073169     0.057775  0.070208  0.043657]
% min_m  = [0       0       0            0            0         0         0]
%myo_signal = [0.062004 0.016599 0.0058135  0.0080541 0.0059509 0.0069086 0.0059901];
% std_m  = [0.062004 0.01817  0.0081087  0.011946  0.0087137 0.010927 0.0073466]
% max_b  = [1.9846   1.9795   1.9868     1.9877    1.9863    1.9871    1.9872 1.9869]
% min_b  = [0        0         0         0         0         0         0]
%blood_signal = [1.169    1.0549   1.051      0.9828    1.0386    1.0682    1.0651 1.011];
% std_b  = [0.52136  0.53744 0.55025     0.57843   0.55695   0.56264 0.57416]
% lambda = [0.1      0.4     0.7         1         1.3       1.6    2 10];
% original: 
% myo: 0.10637 0 0.012758 0.014603
% blood: 1.9806 1.2191e-05 1.0218 0.58779

figure(1); clf;
for itr=1:length(D)
    itr
    %load(['./old/' D(itr).name]);
    load(['./old_09_24_final_vFOV_140cms/' D(itr).name ]); % max/mean
%     load(['./new_penalty/' D(itr).name]);
     myo_signal(itr) = results.max_myo_signal;
     myo_signal2(itr) = results.mean_myo_signal;
     myo_signal_std(itr) = results.std_myo_signal;
     blood_signal(itr) = results.mean_blood_signal;
     blood_signal_std(itr) = results.std_blood_signal;
     lambda(itr) = results.lambda;
     rfweight(itr,:) = results.rfweight;
    
    figure(1); 
    %plot(myo_signal(itr),blood_signal(itr),'*','color',co(itr,:)); hold on;
    %errorbar(myo_signal(itr),blood_signal(itr),-blood_signal_std(itr),blood_signal_std(itr),-myo_signal_std(itr),myo_signal_std(itr),'*','color',co(itr,:)); hold on;
    text(myo_signal(itr) + 0.008,blood_signal(itr),['\lambda=' num2str(lambda(itr))],'FontSize',16); hold on;
    if (itr<3)
            text(myo_signal2(itr) + 0.008,blood_signal(itr),['\lambda=' num2str(lambda(itr))],'FontSize',16); hold on;
    end
    % xlim([-0.003 0.08]);
   % ylim([0.65 1.5]);
end

%%

plot(myo_signal,blood_signal,'*-k',myo_signal2,blood_signal,'*-b'); 
xlabel('ASL signal (myocardial velocities)'); ylabel('ASL signal (coronary velocities)');
legend('max ASL signal (myocardial velocities)','mean ASL signal (myocardial velocities)');
   xlim([-0.003 0.5]);
    ylim([0.65 1.5]);
    set(gca,'FontSize',16);