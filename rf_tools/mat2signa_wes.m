
% function [void] = mat2signa_wes(rf,dt,res,nom_flip,g,nomsl,filename)

% rf:       (G)
% dt:       (ms)
% tbw:      (#)
% res:      (#)
% nom_flip: (degree)
% g:        (G/cm)
% nomsl:    (mm)
% filename: ('...')

% ***** Changes *******
% Aug. 14 2007: Wesley modified the inputs and unit of them.
% Mar. 27 2003: modified to output a file in format for krishna's 3T flouro sequence.
% Feb. 4 1999: created flag for complex valued rf, included ability to write external file for theta board for complex rf.

function [void] = mat2signa_wes(rf,dt,res,tbw,nom_flip,g,nomsl,filename)

% calculate the parameters
area = 1/(length(rf)*abs(max(rf)))*abs(sum(rf));
abswidth = 1/(length(rf)*max(abs(rf)))*sum(abs(rf));
effwidth = 1/(length(rf)*abs(max(rf))^2)*sum(abs(rf.^2));
dtycyc = 1/length(rf)*length(find(abs(rf)>0.2236*max(abs(rf))));
maxpw = dtycyc;
max_b1 = abs(max(rf));
max_int_b1=abs(sum(rf))*dt;
max_int_b1_sq = abs(sum(rf.*conj(rf)))*dt;
max_b1_rms=sqrt(1/(length(rf))*abs(sum(rf.*rf)));
nom_pw = round(res*dt*1000);
a_gzs = max(abs(g));
specspat = 0;
nom_bw = 0;
isodelay = 0;
area_rf1 = 0;
nsubpulse = 0;
% write the text data file
fid0=fopen(filename,'w');

fprintf(fid0,'%d\t\t#res\n',length(rf));
fprintf(fid0,'%d\t\t#pw\n',nom_pw);
fprintf(fid0,'%d\t\t#tbw\n',tbw);
fprintf(fid0,'%f\t#nom_flip\n',nom_flip);
fprintf(fid0,'%f\t#abswidth\n',abswidth);
fprintf(fid0,'%f\t#effwidth\n',effwidth);
fprintf(fid0,'%f\t#area\n',area);
fprintf(fid0,'%f\t#dtycyc\n',dtycyc);
fprintf(fid0,'%f\t#maxpw\n',maxpw);
fprintf(fid0,'%f\t#max_b1\n',max_b1);
fprintf(fid0,'%f\t#max_int_b1\n',max_int_b1);
fprintf(fid0,'%f\t#max_int_b1_sqr\n',max_int_b1_sq);
fprintf(fid0,'%f\t#max_rms_b1\n',max_b1_rms);
fprintf(fid0,'%d\t#specspat\n',specspat);
fprintf(fid0,'%f\t#a_gzs\n',a_gzs);
fprintf(fid0,'%f\t#nom_thk(mm)\n',nomsl);
fprintf(fid0,'%f\t#nom_bw(mm)\n',nom_bw);
fprintf(fid0,'%f\t#isodelay(mm)\n',isodelay);
fprintf(fid0,'%f\t#area_rf1(mm)\n',area_rf1);
fprintf(fid0,'%d\t#area_rf1(mm)\n',nsubpulse);

fclose(fid0);
















