function y = csplinx(t,x,nt);

y = spline(t,real(x),nt)+i*spline(t,imag(x),nt);

