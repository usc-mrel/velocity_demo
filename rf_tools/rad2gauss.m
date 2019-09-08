%[B1gauss] = rad2gauss(rf,T,flip)
%
% 	Scales an RF waveform to it's value in gauss.  
%
%	T is the pulse duration, in ms.
%	flip is the desired flip angle (radians).
%
%	See also pulsetip.m
%	         scalerf.m
%

% ======================== CVS Log Messages ========================
% $Log: rad2gauss.m,v $
% Revision 1.2  2002/03/28 00:50:11  bah
% Added log to source file
%
%
%
% ================================================================== 

function [B1gauss] = rad2gauss(rf,T,flip)
% HPD change
% s = sum(rf);
s = sum(abs(rf));

B1gauss = rf * (max(size(rf))*flip) / (4258*2*pi*T/1000*s);


