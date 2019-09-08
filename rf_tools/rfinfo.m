%
%	function [rfi] = rfinfo(rf,T,flip)
%
%	Function calculates the RF pulse structure information
%	for the given RF pulse.  The info is printed out, and
%	also put in the array rfi.  NAN is used in rfi where there
%	should be a pointer put in RF_PULSE.
%
%	INPUT:
%		rf	- rf pulse (arb units).
%		T	- pulse duration in ms.
%		flip	- flip angle in radians.
%
%
%	OUTPUT:
%		Array with the rf pulse structure information.
%
%	B. Hargreaves  Aug 2001
%
%
%	See also:  rfsim, rad2gauss

function [rfi] = rfinfo(rf,T,flip)

rfg = rad2gauss(rf,T,flip);
rfn = rfg/max(rfg);

gamma = 4258;
peakB1Hz = 622;
peakB1G	= peakB1Hz / gamma;

% Crude bandwidth calculation.
rfp = [rfn];
rfp(10*length(rfn)) = 0;
rfbw = fftshift(fft(fftshift(rfp)));
rfbw = abs(rfbw)/max(rfbw);

plotc(rfbw);
title('RF bandwidth');

rfhi = find(abs(rfbw) > .5);
bw = (max(rfhi)-min(rfhi))/(10*T)*1000;

rfi(1) = NaN;	% pointer to actual duration in ms
rfi(2) = NaN;	% pointer to actual fractional amplitude 
rfi(3) = sum(rfn)/length(rfn);
rfi(4) = sum(rfn.*conj(rfn))/length(rfn);
rfi(5) = rfi(3);
rfi(6) = length(find(abs(rfn)>.23))/length(rfn); 	% *100; 
rfi(7) = length(find(abs(rfn)>0))/length(rfn);  	% *100; 
rfi(8) = 1;	% number of copies... change if necessary!
rfi(9) = max(abs(rfg));
rfi(10) = sum(rfg.*rfg)/length(rfg)*T;
rfi(11) = sqrt(sum(rfg.*rfg)/length(rfg));
rfi(12) = flip*180/pi;
rfi(13) = NaN;	% pointer to desired flip.
rfi(14) = 1000*T;
rfi(15) = bw;
rfi(16) = NaN;
rfi(17) = NaN;

tt = sprintf('{');
disp(tt);
tt = sprintf('(int *)&rf_dur,		/* actual pulse duration. */');
disp(tt);
tt = sprintf('(FLOAT *)&a_rf,		/* actual pulse amplitude (ffs) */');
disp(tt);
tt = sprintf('%8.5f,		/* absolute width */',rfi(3));
disp(tt);
tt = sprintf('%8.5f,		/* effective width */',rfi(4));
disp(tt);
tt = sprintf('%8.5f,		/* absolute width */',rfi(5));
disp(tt);
tt = sprintf('%8.2f,		/* duty cycle 23 percent*/',rfi(6));
disp(tt);
tt = sprintf('%8.2f,		/* duty cycle */',rfi(7));
disp(tt);
tt = sprintf('%8d,		/* number of copies of pulse */',rfi(8));
disp(tt);
tt = sprintf('%8.5f,		/* max B1 (G) */',rfi(9));
disp(tt);
tt = sprintf('%8.5f,		/* int B1^2 (GG-ms) */',rfi(10));
disp(tt);
tt = sprintf('%8.5f,		/* rms B1 (G) */',rfi(11));
disp(tt);
tt = sprintf('%8.3f,		/* nom fa (deg) */',rfi(12));
disp(tt);
tt = sprintf('(FLOAT *)&flip_rf, 	/* actual flip angle */');
disp(tt);
tt = sprintf('%8.2f,		/* nom pulse width */',rfi(14));
disp(tt);
tt = sprintf('%8.2f,		/* nom pulse bandwidth */',rfi(15));
disp(tt);
tt = sprintf('PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,		/* activity */');
disp(tt);
tt = sprintf('  0,			/* flag for pulse used in TG setting.*/');
disp(tt);
tt = sprintf('%8.2f,		/* isodelay */',rfi(14)/2);

disp(tt);
tt = sprintf('  1.0,			/* duty-cycle scale factor */');
disp(tt);

tt = sprintf('(int *)&res_rf,		/* actual pulse resolution */');
disp(tt);
tt = sprintf('TRUE			/* rfpulse uses external gradient waveform file */');
disp(tt);
tt = sprintf('}');
disp(tt);




