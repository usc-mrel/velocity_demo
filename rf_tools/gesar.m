%function [HeadSAR, BodySAR]=gesar(rf,dur,mass,TR)
%
%	Calculates SAR for head coil and body coil
%	using GE's empirical formulas in terms of
%	equivalent pulses.  See Epic 4.x manual for
%	a description.
%
%	rf - pulse profile, in gauss.
%	dur - pulse duration in ms.
%	mass - mass in kg.
%	TR - repetition time in ms.
%
%	SAR = [HeadSAR BodySAR]
%

%	By Brian Hargreaves, 	Sept, 1997.

function [HeadSAR, BodySAR]=gesar(rf,dur,mass,TR)

% Calculate # "equivalent (or standard) pulses"

s = max(size(rf));
gamma = 4258;				% Hz/G.
intbsq =  sum (rf .* conj(rf))/s*dur;
numeq = real(intbsq) / (500/gamma)^2;


% Calculate SAR of equivalent pulse

sareqh = .0016*(mass)^0.835; 
sareqb = .0026*(mass)^1.48;

% Put it all together...

HeadSAR = sareqh * numeq / (TR/1000) / mass;
BodySAR = sareqb * numeq / (TR/1000) / mass;

if (nargout ==1)
  HeadSAR = [HeadSAR BodySAR];
end;





