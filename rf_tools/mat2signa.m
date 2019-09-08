%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab_to_signa - this function takes as input an rf pulse variable, which
% is in units of radians, the time of the pulse pw in seconds, the time-band
% width product and the optimal refocus interval.  Two files are written,
% rf.rho and rf.txt
% function [void] = mat2signa (rf,pw,nom_flip,gy,gz,nomsl)

% march 27 2003: modified to output a file in format for krishna's
% 3T flouro sequence.

% ***** Changes *******
% Feb. 4 1999: created flag for complex valued rf, included ability to 
%              write external file for theta board for complex rf.
%
% *********************

function [void] = mat2signa(rf,pw,nom_flip,gy,gz,nomsl)

% some constants are set
gamma = 26752; %radians per second per gauss

% calculate some useful stuff
dt = pw/length(rf);
b1 = rf/(gamma*dt);

% calculate the parameters
area = 1/(length(rf)*abs(max(rf)))*abs(sum(rf))
abswidth = 1/(length(rf)*max(abs(rf)))*sum(abs(rf))
effwidth = 1/(length(rf)*abs(max(rf))^2)*sum(abs(rf.^2))
dtycyc = 1/length(rf)*length(find(abs(rf)>0.2236*max(abs(rf))))
maxpw = dtycyc
max_b1 = abs(max(b1))
max_int_b1=abs(sum(b1))*dt*1000
max_int_b1_sq = abs(sum(b1.*conj(b1)))*dt*1000
max_b1_rms=sqrt(1/(length(b1))*abs(sum(b1.*b1)))
nom_pw = round(pw*1E6)
a_gys = max(abs(gy));
a_gzs = max(abs(gz));

% write the text data file
fid0=fopen('./testrf.dat','w');
fprintf(fid0,'1 \t\t#Spectral-Spatial\n');
fprintf(fid0,'%d\t\t#res\n',length(rf));
fprintf(fid0,'%d\t\t#pw\n',nom_pw);
fprintf(fid0,'%f\t#nom_flip\n',nom_flip);
fprintf(fid0,'%f\t#abswidth\n',abswidth);
fprintf(fid0,'%f\t#effwidth\n',effwidth);
fprintf(fid0,'%f\t#area\n',area);
fprintf(fid0,'%f\t#dtycyc\n',dtycyc);
fprintf(fid0,'%f\t#maxpw\n',maxpw);
fprintf(fid0,'%f\t#max_b1\n',max_b1);
fprintf(fid0,'%f\t#max_int_b1_sqr\n',max_int_b1_sq);
fprintf(fid0,'%f\t#max_rms_b1\n',max_b1_rms);
fprintf(fid0,'%f\t#a_gys\n',a_gys);
fprintf(fid0,'%f\t#a_gzs\n',a_gzs);
fprintf(fid0,'%f\t#nom_thk(mm)\n',nomsl);
fprintf(fid0,'0\t\t#g_pow\n');
fprintf(fid0,'0\t\t#g_pos_pow\n');
fprintf(fid0,'0\t\t#g_neg_pow\n');
fprintf(fid0,'0\t\t#g_abs\n');
fprintf(fid0,'0\t\t#g_dgdt\n');
fprintf(fid0,'0\t\t#g_pwm\n');
fprintf(fid0,'0\t\t#g_pwm_abs\n');
fclose(fid0);




