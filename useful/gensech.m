% int gensech(float beta, float mu, int pwsech, float dtus, float cexp)
% {
%   int i;
%   double t, dt,dur,sech,mid,phase,dum;
%
%   dt = 1.e-6 * dtus; 		   /* convert to seconds */
%   dur = 1.e-6 * pwsech; 		   /* convert to seconds */
%   sechres = dur/dt + 0.5;
%
% /* Calculate waveforms */
%   mid = 0.5*(sechres-1);
%   for (i=0;i<sechres;i++) {
%     t = (i-mid) * dt;
%     sech = 2.0 / (exp(beta*t) + exp(-1.0*beta*t));
%     sechrho[i] = EVENIZE(MAX_PG_WAMP * sech * pow(cos(t*PI/dur),(double)cexp));
%     phase = -1.0 * mu * log(sech);
%     sechtheta[i] = EVENIZE(MAX_PG_WAMP * (2. * modf(phase/(2.*PI)+0.5,&dum) - 1.));
%     /* -1 matches sign of phase on signa's sech pulse */
%   }
%
%   return SUCCESS;
%
% }


function [rf_im,grf,bw] =gensech(beta,  mu,  pwsech,  dtus,  cexp, ssinvthick, b1_scale)
gslewmax    = 15;             %G/cm/ms (14.5 conservative)
gmax        = 4;              %G/cm (3.9 convservative)
b1max       = 0.115*b1_scale;
GAM=4257.59; % Hz/G
EVENIZE=@(x) 2*x/2 ;
dt = 1.e-6 * dtus; 		   % convert to seconds */
dur = 1.e-6 * pwsech;      % convert to seconds */
sechres = dur/dt + 0.5;
mid = 0.5*(sechres-1);
for ii=1:sechres
    t = (ii-mid) * dt;
    sech2 = 2.0 / (exp(beta*t) + exp(-1.0*beta*t));
    sechrho(ii) = (1 * sech2 * (cosd(t*pi/dur).^cexp));
    phase = -1.0 * mu * log(sech2);
    sechtheta(ii) = (1 * (2. * (phase/(2.*pi)+0.5) - 1.));
    df(ii)=-1*mu*beta*tanh(beta*t);
    At(ii)= b1max*sech(beta*t);
end
rf_im=sechrho.*exp(pi*1j*(sechtheta));
bw=beta*mu/pi;
a_gz = -1*10.*beta*mu/(pi*GAM*ssinvthick);
grf = const2trap(a_gz * ones(1, length(rf_im)));
rf_im = b1max*[zeros(1,(length(grf)-length(rf_im))/2) rf_im zeros(1,(length(grf)-length(rf_im))/2)];

    function [g, varargout] =  const2trap(g, varargin)
        
        if nargin > 1
            rf = varargin{1};
        else
            rf = [];
        end
        
        % Find max gradient and sign of g
        [amp, ind] = max(abs(g));
        sign_g = sign(g(ind));
        
        % Calculate time for ramp
        tramp = abs(amp)/gslewmax;
        tramp = tramp + dt*1000 - mod(tramp,dt*1000);
        % Recalculate slew rate
        gslew = amp/tramp;
        gramp = sign_g*gslew*(0:dt*1000:tramp-dt*1000);
        rampres = length(gramp);
        
        % Append Ramps
        g = [gramp g fliplr(gramp)];
        
        if nargout > 1
            rf = [zeros(1, rampres) rf zeros(1, rampres)];
            varargout{1} = rf;
        end
        
    end
end
%

% float sechbeta = 260.0 with {0,,,VIS, "Beta of sech pulse (1/s)",};
% float sechmu = 52 with {0,,,VIS, "Mu of sech pulse",};
% float sechdt = 4.0 with {0,,,VIS, "Time resolution of sech pulse (us)",};
% float sechcexp = 0.65 with {0,,,VIS, "exponent for cosine window of sech",};
% float sechpw = 15.0 with {0,,,VIS, "Width of sech pulse (ms)",};
% int sechres = 0 with {0,,,VIS, "resolution for sech pulse",};