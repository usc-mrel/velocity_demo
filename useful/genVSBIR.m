function [b1, gz, dt] = genVSBIR(zeta, kappa, omega, Vmax, ...
    varargin)
%genVSBIR generates velocity selective adiabatic BIR pulses
%   [b1, gz, dt] = genVSBIR(zeta, kappa, omega) generates a adiabatic BIR-4
%   pulse, b1, given parameters zeta(z), kappa(k) and omega(w) in kHz with
%   amplitude modulation of:
%
%   A(t) = tanh(z*(1-t/Tseg)) odd segments
%   A(t) = tanh(z*(t/Tseg)) even segments
%
%   and frequency modulation of:
%
%   phi(t) = -wmax*Tseg*ln(|cos(k*t/Tseg)|/(k(tan(k))) odd segments
%   phi(t) = -wmax*Tseg*ln(|cos(k*t/Tseg-1)|/(k(tan(k))) even segments
%
%   along with gradient, gz, to saturation all spins above a cutoff
%   velocity Vmax. By default, the RF timestep is 2us and the gradient
%   timestep is 4 us.
%
%   [b1, gz, dt] = genVSBIR(..., 'type', X, 'dt', [dtRF dtGz], 'Gmax', G,
%       'Smax', s, 'gamma', g)
%   has the optional paramters, 'type', which generates a BIR-X pulse where
%   X is either 4 or 8 and determines the number of segments in the BIR
%   pulse, 'dt', which sets the RF and Gz timestep, 'Gmax', which
%   denotes the maximum gradient amplitude of the MRI scanner, Smax, which
%   denotes the maximum gradient slew rate, and 'gamma', the gyromagnetic
%   ratio.
%
%   Author: Terrence Jao
%   Date:   1-29-2014
%   Ref:
%   [1] Meakin, James A., and Peter Jezzard. ?An Optimized Velocity
%   Selective Arterial Spin Labeling Module with Reduced Eddy Current
%   Sensitivity for Improved Perfusion Quantification.? Magnetic Resonance
%   in Medicine 69, no. 3 (March 1, 2013): 832?838. doi:10.1002/mrm.24302.
%
%   Moment Derivation
%   [2] Glover, G. H. ?Handbook of MRI Pulse Sequences, M. A. Bernstein,
%   K. F. King and X. J. Zhou. Elsevier Academic Press, 2004,
%   ISBN: 0-12-092861-2.? NMR in Biomedicine 18, no. 3 (May 2005): 285?287.
%   doi:10.1002/nbm.947.

%--------------------------------------------------------------------------
%Input Parser
%--------------------------------------------------------------------------
p = inputParser;
p.addRequired('zeta', @isscalar);
p.addRequired('kappa', @isscalar);
p.addRequired('omega', @isscalar);
p.addRequired('Vmax', @isscalar);
p.addParamValue('sym', 1, @(x) x == 1 || x==0);
p.addParamValue('type', 4, @(x) x==4 || x==8);
p.addParamValue('dt', [2 4], @(x) length(x) == 2);
p.addParamValue('Gmax', 2, @isscalar);          %Max Grad 4 G/cm
p.addParamValue('Smax', 15, @isscalar);         %Max Slew 15 G/cm/ms
p.addParamValue('B1max', 0.16, @isscalar);      %Max B1 ampliutde 0.16G
p.addParamValue('gamma',.026752, @isscalar);    %Gamma radians/us/G
p.addParamValue('ramp_time', 500, @isscalar);   %Ramp Time us
p.addParamValue('grad_wait', 0, @isscalar);     %Gradient Wait
p.addParamValue('RF_width', 2000, @isscalar);   %RF lobe width us
p.parse(zeta, kappa, omega, Vmax, varargin{:});
sym = p.Results.sym;
type = p.Results.type;
dt = p.Results.dt;
gmax = p.Results.Gmax;
smax = p.Results.Smax;
b1max = p.Results.B1max;
gamma = p.Results.gamma;
r = p.Results.ramp_time;
grad_wait = p.Results.grad_wait;
Tseg = p.Results.RF_width;

dtRF = dt(1);   dtGz = dt(2);
omega = pi/45*1e-3*omega;              %conversion from kHz to M-rad/s
smax = 1e-3*smax;                %conversion from G/cm/ms to G/cm/us
%--------------------------------------------------------------------------
%Generate RF and Gradient
%--------------------------------------------------------------------------

dtRatio = dtGz/dtRF;
grad_wait_rf = grad_wait/dt(1);   %grad_wait is in us, must divide by dt to get length of grad_wait
grad_wait_gz = grad_wait/dt(2);

[AO, phiO, AE, phiE] = genSeg(zeta, kappa, omega, dtRF, Tseg);
Gzlobe = genGz(Vmax, gmax, smax, gamma, r, dtGz, grad_wait, sym, type, Tseg);

% lenGz = modsub(length(Gzlobe)*dtRatio, dtRF);
% lenBIR = modsub(length(AO)/dtRatio, dtGz);

lenGz = length(Gzlobe)*dtRatio;
lenBIR = length(AO)/dtRatio;
% if do_grad_wait == 1

%Asymetric BIR8
if sym == 0
    switch type
        case 4
            A = b1max*[AO zeros(1,lenGz) zeros(1,grad_wait_rf) AE AO zeros(1,lenGz) zeros(1,grad_wait_rf) AE];
            phi = [phiO zeros(1,lenGz) zeros(1,grad_wait_rf) phiE phiO zeros(1,lenGz) zeros(1,grad_wait_rf) phiE];
            gz = [zeros(1,lenBIR) Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) Gzlobe zeros(1,grad_wait_gz) ...
                zeros(1,lenBIR)];
        case 8
            A = b1max*[AO AE AO   zeros(1,lenGz) zeros(1,grad_wait_rf) AE AO    zeros(1,grad_wait_rf) zeros(1,2*lenGz) zeros(1,grad_wait_rf) ...
                AE AO  zeros(1,grad_wait_rf) zeros(1,lenGz)  AE];
            phi = [phiO phiE phiO   zeros(1,lenGz) zeros(1,grad_wait_rf) phiE phiO    zeros(1,grad_wait_rf) zeros(1,2*lenGz) zeros(1,grad_wait_rf) ...
                phiE phiO  zeros(1,grad_wait_rf) zeros(1,lenGz) phiE];
            gz = [zeros(1,3*lenBIR)  Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR)  zeros(1,grad_wait_gz) -Gzlobe Gzlobe zeros(1,grad_wait_gz) ...
                zeros(1,2*lenBIR) zeros(1,grad_wait_gz) -Gzlobe  zeros(1,lenBIR)];
        otherwise
            A = b1max*[AO zeros(1,lenGz) AE AO zeros(1,lenGz) AE];
            phi = [phiO zeros(1,lenGz) phiE phiO zeros(1,lenGz) phiE];
            gz = [zeros(1,lenBIR) Gzlobe zeros(1,2*lenBIR) Gzlobe ...
                zeros(1,lenBIR)];
    end
    %Symetric BIR8
else
    switch type
        case 4
%             A = b1max*[AO zeros(1,lenGz) zeros(1,grad_wait_rf) AE AO zeros(1,lenGz) zeros(1,grad_wait_rf) AE];
%             phi = [phiO zeros(1,lenGz) zeros(1,grad_wait_rf) phiE phiO zeros(1,lenGz) zeros(1,grad_wait_rf) phiE];
%             gz = [zeros(1,lenBIR) Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) Gzlobe zeros(1,grad_wait_gz) ...
%                 zeros(1,lenBIR)];
%             A = b1max*[AO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) AE AO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) AE];
%             phi = [phiO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) phiE phiO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) phiE];
%             gz = [zeros(1,lenBIR) Gzlobe zeros(1,grad_wait_rf) -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) ...
%                 -Gzlobe zeros(1,grad_wait_gz) Gzlobe zeros(1,grad_wait_rf) zeros(1,lenBIR)];
            A = b1max*[AO zeros(1,2*lenGz) zeros(1,grad_wait_rf) AE AO zeros(1,2*lenGz) zeros(1,grad_wait_rf) AE];
            phi = [phiO zeros(1,2*lenGz) zeros(1,grad_wait_rf) phiE phiO zeros(1,2*lenGz) zeros(1,grad_wait_rf) phiE];
            gz = [zeros(1,lenBIR) Gzlobe -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) ...
                -Gzlobe Gzlobe zeros(1,grad_wait_rf) zeros(1,lenBIR)];

        case 8
%             A = b1max*[AO zeros(1,lenGz) zeros(1,grad_wait_rf) AE AO zeros(1,lenGz) zeros(1,grad_wait_rf) ...
%                 AE AO zeros(1,lenGz) zeros(1,grad_wait_rf)  AE AO  zeros(1,lenGz) zeros(1,grad_wait_rf)  AE];
%             phi = [phiO  zeros(1,lenGz) zeros(1,grad_wait_rf)  phiE phiO   zeros(1,lenGz) zeros(1,grad_wait_rf) ...
%                 phiE phiO  zeros(1,lenGz) zeros(1,grad_wait_rf) phiE phiO   zeros(1,lenGz) zeros(1,grad_wait_rf) phiE];
%             gz = [zeros(1,lenBIR)  Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) -Gzlobe  zeros(1,grad_wait_gz) zeros(1,2*lenBIR)...
%                 -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) Gzlobe  zeros(1,grad_wait_gz)  zeros(1,lenBIR) ];
%             A = b1max*[AO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) AE AO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) ...
%                 AE AO zeros(1,2*lenGz) zeros(1,2*grad_wait_rf)  AE AO  zeros(1,2*lenGz) zeros(1,2*grad_wait_rf)  AE];
%             phi = [phiO  zeros(1,2*lenGz) zeros(1,2*grad_wait_rf)  phiE phiO   zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) ...
%                 phiE phiO  zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) phiE phiO   zeros(1,2*lenGz) zeros(1,2*grad_wait_rf) phiE];
%             gz = [zeros(1,lenBIR)  Gzlobe zeros(1,grad_wait_gz) -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) -Gzlobe  zeros(1,grad_wait_gz) ...
%                 Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR)...
%                 Gzlobe zeros(1,grad_wait_gz) -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) -Gzlobe zeros(1,grad_wait_gz) Gzlobe  zeros(1,grad_wait_gz)  zeros(1,lenBIR) ];
%                  
            
             A = b1max*[AO zeros(1,2*lenGz) zeros(1,grad_wait_rf) AE AO zeros(1,2*lenGz) zeros(1,grad_wait_rf) ...
                AE AO zeros(1,2*lenGz) zeros(1,grad_wait_rf)  AE AO  zeros(1,2*lenGz) zeros(1,grad_wait_rf)  AE];
            phi = [phiO  zeros(1,2*lenGz) zeros(1,grad_wait_rf)  phiE phiO   zeros(1,2*lenGz) zeros(1,grad_wait_rf) ...
                phiE phiO  zeros(1,2*lenGz) zeros(1,grad_wait_rf) phiE phiO   zeros(1,2*lenGz) zeros(1,grad_wait_rf) phiE];
            gz = [zeros(1,lenBIR)  Gzlobe -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) -Gzlobe  ...
                Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR)...
                Gzlobe -Gzlobe zeros(1,grad_wait_gz) zeros(1,2*lenBIR) -Gzlobe Gzlobe  zeros(1,grad_wait_gz)  zeros(1,lenBIR) ];
                 
        
        otherwise
            A = b1max*[AO zeros(1,lenGz) AE AO zeros(1,lenGz) AE];
            phi = [phiO zeros(1,lenGz) phiE phiO zeros(1,lenGz) phiE];
            gz = [zeros(1,lenBIR) Gzlobe zeros(1,2*lenBIR) Gzlobe ...
                zeros(1,lenBIR)];
    end
end
% else
%     switch type
%         case 4
%             A = b1max*[AO zeros(1,lenGz) AE AO zeros(1,lenGz) AE];
%             phi = [phiO zeros(1,lenGz) phiE phiO zeros(1,lenGz) phiE];
%             gz = [zeros(1,lenBIR) Gzlobe zeros(1,2*lenBIR) Gzlobe ...
%                 zeros(1,lenBIR)];
%         case 8
%             A = b1max*[AO AE AO zeros(1,lenGz) AE AO zeros(1,2*lenGz) ...
%                 AE AO zeros(1,lenGz) AE];
%             phi = [phiO phiE phiO zeros(1,lenGz) phiE phiO zeros(1,2*lenGz) ...
%                 phiE phiO zeros(1,lenGz) phiE];
%             gz = [zeros(1,3*lenBIR) Gzlobe zeros(1,2*lenBIR) -Gzlobe Gzlobe ...
%                 zeros(1,2*lenBIR) -Gzlobe zeros(1,lenBIR)];
%         otherwise
%             A = b1max*[AO zeros(1,lenGz) AE AO zeros(1,lenGz) AE];
%             phi = [phiO zeros(1,lenGz) phiE phiO zeros(1,lenGz) phiE];
%             gz = [zeros(1,lenBIR) Gzlobe zeros(1,2*lenBIR) Gzlobe ...
%                 zeros(1,lenBIR)];
%     end
%
% end
%----Output---%
b1 = A.*exp(1i.*phi);

end

function [AO, phiO, AE, phiE] = genSeg(zeta, kappa, omega, dt, Tseg)
% genSeg generates the odd and even segments of the BIR-X pulse
%Tseg = 2000;
%t = 0:dt:modsub(Tseg,dt)-dt;
t = 0:dt:modsub(Tseg,dt);

AO = tanh(zeta*(1-t/Tseg));
AE = tanh(zeta*(t/Tseg));

phiO = -omega*Tseg*log(abs(cos(kappa*t/Tseg))/(kappa*tan(kappa)));
phiE = -omega*Tseg*log(abs(cos(kappa*((t/Tseg)-1)))/(kappa*tan(kappa)));

end

function Gz = genGz(Vmax, gmax, smax, gamma, r, dt, grad_wait, sym, type, Tseg)
% genGz generates the gradient lobes to achieve a cutoff at Vmax
%    _________
%   /|       |\
%  / |       | \
% /__|_______|__\
%
% M_a   M_p   M_d     First moment of Dephaser = Rephaser
%
% m = integral{t*f(t)}
% Ascending Ramp:
% M_a = integral{t*smax*t} |0 to gmax/smax
% M_a = smax*t^3/3 | 0 to gmax/smax = gmax^3/(3*smax^2)
%
% Plateau:
% M_p = integral{t*gmax}|0 to T
% M_p = T^2*gmax/2
%
% Descending Ramp:
% M_d = integral{t*(-smaxt + gmax)} | 0 to gmax/smax
% M_d = -smax*t^3/3 + gmax*t^2/2    | 0 to gmax/smax
% M_d = -gmax^3/(3*smax^2) + gmax^3/(2*smax^2) = gmax^3/(6*smax^2)
%
% First moment of translated waveform:
% m0*delT + m1 where m0 is zeroth moment and m1 is first moment
%
% M = M_a + M_f + M_d
% let r = gmax/smax, where r is ramp time
% Mlobe = gmax*(r^2 + 3rT/2 + T^2/2)
%
% Zeroth Moment: 2*G*r

%Tseg = 2000;

Vmax = 1e-6*Vmax;           %Conversion from cm/sec to cm/us
M = pi/(2*gamma*Vmax);        %Moment of gradient lobe

%r = gmax/smax;             %Ramp time

%We are fixing the ramp time to be half a millisecond. We can change the
%slew rate accordingly to reach different gradient amplitudes in the
%alloted ramp time.

%r = 500;                    %Ramp time

%T = roots([gmax/2 3*gmax*r/2 gmax*r^2-M]);
%T = roots([2*gmax 18*gmax*r 28*gmax*r^2+8*gmax*r*Tseg-M]);

%T=roots([8*gmax 24*gmax*r+8*gmax*Tseg 16*gmax*r^2+8*gmax*r*Tseg-M]);
%T=roots([-4*gmax*Tseg -8*gmax*r*Tseg-M]);

%only ture for BIR-8 Pulse right now

%BELOW IS CORRECT WITHOUT GRAD_WAIT
%T=roots([4*gmax 12*gmax*r+4*gmax*Tseg 8*gmax*r^2+4*gmax*r*Tseg-M]);

if sym==0
    %------------------------------------------------------
    %Calcualte Gradient Moment for Assymetric BIR
    %------------------------------------------------------
    if type == 4
        % BIR4 Calculation
    end
    if type == 8
        % BIR8 Calculation
%         T=roots([4*gmax 12*gmax*r+4*gmax*(Tseg+grad_wait) 8*gmax*r^2+4*gmax*r*(Tseg+grad_wait)-M]);
        T=roots([4*gmax 12*gmax*r+4*gmax*(Tseg+grad_wait) 8*gmax*r^2+4*gmax*r*(Tseg+grad_wait)-M]);
        
        if imag(T(1)) ~= 0
            display('Error!, Time cannot be complex!')
        end
        
        T = max(T);
        
        if T >= 0
            g = gmax;
            s = g/r;
        else
            %If plateau time is less than 0, then we need to lower the gradient
            %ampliutde to obtain the correct first moment. g<gmax.
            T = 0;
            %g = M/(16*r^2+8*r*Tseg);
            g = M/(8*r^2+4*r*(Tseg+grad_wait));
            s = g/r;
        end
    end
    
else
    %------------------------------------------------------
    % --- Calcaulte Gradient Momvent for Symmetric BIR ---%
    if type == 4
        % BIR4 Calculation 
         % This equation is for bipolar formulation with gradient delay after
         % every bipolar gradient set
        T=roots(gmax*[2 6*r 4*r^2-M/gmax]);
        
        % This equation is for bipolar formulation with gradient delay after
        % every gradient
        %T=roots(gmax*[2 6*r+2*grad_wait 4*r^2+2*r*grad_wait-M/gmax]);
        
        if imag(T(1)) ~= 0
            display('Error!, Time cannot be complex!')
        end
        
        T = max(T);
        %fprintf('T = %0.3f\n', T);
        if T >= 0
            g = gmax;
            s = g/r;
        else
        
        
        T = 4;      %should be T = 0, but we have to have something...
        % This equation is for bipolar formulation with gradient delay after
        % every bipolar gradient set
        g = M/((r+T)*(4*r + 2*T));
        
        % This equation is for bipolar formulation with gradient after
        % every gradient
        % g = M/((r+T)*(4*r + 2*T + 2*grad_wait));
        s = g/r;
        
        end
        
        
    end
    if type == 8
        % BIR8 Calculation
        
         % This equation is for bipolar formulation with gradient delay after
         % every bipolar gradient set
          T=roots(gmax*[2 6*r 4*r^2-M/gmax]);
%          T=roots(gmax*[4 12*r 8*r^2-M/gmax]);
         
        
        % This equation is for bipolar formulation with gradient delay after
        % every gradient
        % T=roots(gmax*[4 12*r+4*grad_wait 8*r^2+4*r*grad_wait-M/gmax]);
         
        % This equation is for unipolar formulation: 
        % T=roots(gmax*[8 24*r+6*(Tseg+grad_wait) 16*r^2+6*r*(Tseg+grad_wait)-M/gmax]);
        
        if imag(T(1)) ~= 0
            display('Error!, Time cannot be complex!')
        end
        
        T = max(T)
        %fprintf('T = %0.3f\n', T);
        if T >= 0
            g = gmax;
            s = g/r;
        else
            %If plateau time is less than 0, then we need to lower the gradient
            %ampliutde to obtain the correct first moment. g<gmax.
            
            %Without T
%             T = 0;
%             g = M/(16*r^2+6*r*Tseg+4*r*grad_wait);
            
            % This equation is for bipolar formulation with gradient delay after
            % every bipolar gradient set
            T = 4;
            g = M/((r+T)*(4*r + 2*T));
            
%             T = 4;
%             g = M/((r+T)*(8*r + 4*T));
            
            % This equation is for bipolar formulation with gradient delay after
            % every gradient
            %With a little bit of T = 4 us
%             T = 4;
%             g = M/((r+T)*(8*r + 4*T + 4*grad_wait));
            
            % This equation is for unipolar formulation:
            %g = M/(4*(r+T)*(4*r + 2*T + 1.5*Tseg + grad_wait));
            
            s = g/r;
        end
    end
    
    
end
t_d = 0:dt:modsub(g/s,dt)-dt;   %time stamp of dephaser
N = modsub(T/dt,dt);            %number of samples in plateau

Gz = [t_d*s g*ones(1,N) fliplr(t_d*s)];
%keyboard;
end

function z = modsub(x, y)
z = x-mod(x,y);
end



