% plotc(x,y) 
%
% plots the complex vector y as a function of x.
% plots real (blue --) ,imag (red :)  and magnitude (yellow -).
% returns the row as a vector. 

function plotc(x,y); 

if (nargin <2)
	y = x;
	x = 1:max(size(y));
	% color... plot(x,real(y),'b--',x,imag(y),'r:',x,abs(y),'y-');
	p = plot(x,real(y),'k--',x,imag(y),'k:',x,abs(y),'k-');
	set(p, 'LineWidth', 3);
        grid on; xlim([min(x) max(x)]);
	legend('Real', 'Imaginary', 'Absolute');


else
	% color... plot(x,real(y),'b--',x,imag(y),'r:',x,abs(y),'y-');
	p = plot(x,real(y),'k--',x,imag(y),'k:',x,abs(y),'k-');
        set(p, 'LineWidth', 3);
	grid on; xlim([min(x) max(x)]);
        legend('Real', 'Imaginary', 'Absolute');

end;

