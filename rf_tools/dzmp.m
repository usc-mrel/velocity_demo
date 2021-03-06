% dzmp -- design a minimum phase filter by first designing a linear phase
%         filter and factoring it into two terms

%  written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University

function h = dzmp(n,tb,d1,d2)

% n2 = 2*n-1;
n2 = 2*n;
di = 0.5*dinf(2*d1,0.5*d2.*d2);
w = di/tb;
f = [0 (1-w)*(tb/2) (1+w)*(tb/2) (n/2)]/(n/2);
m = [1 1 0 0];
w = [1 2*d1/(0.5*d2*d2)];
f(find(f>1))=1-0.05;
%f= [0 0.2 0.9 1];
% if(length(find(f>1))>0)
%    f(3)=1-0.05;
%    f(2)=f(2)-0.05;
% end
if exist('remez'),     % finally, the m script
  hl = eval('remez(n2-1,f,m,w)');
else
  disp('No PM routine found!');
end;
h = fmp(hl);

