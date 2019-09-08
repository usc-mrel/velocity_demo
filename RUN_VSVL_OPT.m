%RUN_VS_OPT designs VS BIR4 or BIR8 with optimized at a user defined range
% of B0 and B1.
%
%   Author: Vanessa, adapted from Terrence Jao
%   Date: ...
%   References:

%% Pulse Timing Calculation
clear all; close all;
% -------------------------------------------------------------------------
% Do we want to plot?
ploton = 0;
% Do we want to display messages?
verbose = 1;        %Turn off verbose for slightly faster runtime

%--------------------------------------------------------------------------
%Parameters of BGS Optimization Script - Tweak Here!
%--------------------------------------------------------------------------


% ---- Tissue Paramters ---- %
T1 = 10000; %1.932;         % T1 Blood at 3T (x)
T2 = 10000; %0.275;         % T2 blood 3T (s)

% ---- Pulse Paramters ---- %
Vc = 5; %10;            % Cutoff Velocity (cm/s)
b1nom = .15;        % b1 scale
ramp_time = 500;    % length of ramps (us)
% rf_width = 1500;    % length of RF subpulse (us)
grad_wait = 0;      % Delay between gradient and RF (us)
dtRF = 2*51;           % time step RF (us)
dtG = 2*51;            % time step Gradient (us)
%type = 8;           % 4 - BIR4 or 8 - BIR8


% ---- VS Optimization Paramters ---- %
z0 = 0;
Nv = 51;            % Number of velocities to simulate
Nf = 3;           % Number of frequencies to simulate
Nb1 = 5;           % Number of b1 scales to simulate

delf = 50;         % off-resonance range to simulate (-delf, + delf) (Hz);

v = linspace(-4*Vc, 4*Vc, Nv);
df = linspace(-delf,delf,Nf);
%df = 0;
b1scale = linspace(.8,1.3,Nb1);
%b1scale = 1;

%--------------------------------------------------------------------------
%Running the Optimization
%--------------------------------------------------------------------------
 
phi_sp = zeros(1,9); %ones(1,9)*pi/2;
phi_inv = [0 0 pi pi pi 0 0 pi pi pi 0 0 0 pi pi 0];
delay = 0;

%load('J_2018_08_12') % <- Generate initial guesses!?
%idx = find(J == min(J(:)))
%[I1,I2,I3,I4,I5,I6,I7,I8,I9] = ind2sub(size(J),idx);
for itr = 5 %1:size(I1)
  %  vals = [I1(itr),I2(itr),I3(itr),I4(itr),I5(itr),I6(itr),I7(itr),I8(itr),I9(itr)]
    
  %  phi_sp = -pi + pi/2*(vals-1);
    %Initial Guess for Tm (remember you may get stuck at local minima)
    %tm = [2.6520    0.0422    1.2048   -2.7217    0.0096...
    %      0.3470    0.2897    0.2939   -0.0877]; %[phi_sp]; % phi_inv delay];
    %tm = [0.7321    0.5176    1.8270    1.7932    0.1276    0.4268    0.3358    0.1894    0.0074];
    
    %% Initial guess (based on Nb1 = 5, b1scale=0.8-1.3, & Nf=1, df=0;
    % phi_sp =  [-0.2160   -0.4066   -1.9293 0.0047    1.8630   -2.3865 ...
    %             2.4916    0.7299   -0.8933];
    % phi_inv =   [ -0.4273   -0.5228    2.0189 0.9676    0.5584    0.3689 ...
    %         -0.8347    1.6182    0.8704 -1.4080   -1.9820    0.5488 ...
    %         0.0155    2.0742    1.5645 -2.0424];
    % %% Initial guess #2, based on Nb1 & df = 5
    % phi_sp = [1.6715   -2.8425   -0.4132 2.1737    1.8228   -3.1320 ...
    %             1.1338    1.0290   -0.0454];
    % phi_inv = [0.1136    0.4312    2.7771 1.5072    3.1188    1.9540 ...
    %         0.1451    2.2849    1.3915 -1.2471   -1.8339    0.4087 ...
    %         -0.0712    2.7871    2.7855 -0.4849];
    %%
    tm = [phi_sp phi_inv delay];
    %tm = ones(1,16);
    fvalprev = inf;
    init_vals = tm;
    
    model = @(x) vsvlimFun(x, b1nom, Vc, ramp_time, grad_wait, ...
        dtRF, dtG, df, b1scale, v, T1, T2);
    options = optimoptions(@fmincon,'MaxIterations',10000, ...
                    'MaxFunctionEvaluations',10000, ...
                    'OptimalityTolerance',1e-8, ...
                    'ConstraintTolerance',1e-8)
    
    lb = -pi*1.01*ones(1,length(tm));
    ub = pi*1.01*ones(1,length(tm));
    
    %% eq:
    % alpha(t) = gamma*int_t^T(t-s)*v_0*G(s)ds
    % v_0 = 0, ±FOVv, +2FOVv, ...
    % T = pulse duration, G(s) = bipolar gradient pair?
    
    fprintf('-------------------------------\n');
    fprintf('Outer Iteration %d\n', itr);
    fprintf('-------------------------------\n');
    
%     id1 = pi-1;
%     id2 = pi-1;
%     Aeq = [1 0 0 0 0  0 0 0 0 id1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%            1 0 0 0 0  0 0 0 0 0 id2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%            0 1 0 0 0  0 0 0 0 0 0 id1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%            0 1 0 0 0  0 0 0 0 0 0 0 id2 0 0 0 0 0 0 0 0 0 0 0 0 0;
%            0 0 1 0 0  0 0 0 0 0 0 0 0 id1 0 0 0 0 0 0 0 0 0 0 0 0;
%            0 0 1 0 0  0 0 0 0 0 0 0 0 0 id2 0 0 0 0 0 0 0 0 0 0 0;
%            0 0 0 1 0  0 0 0 0 0 0 0 0 0 0 id1 0 0 0 0 0 0 0 0 0 0;
%            0 0 0 1 0  0 0 0 0 0 0 0 0 0 0 0 id2 0 0 0 0 0 0 0 0 0;
%            0 0 0 0 1  0 0 0 0 0 0 0 0 0 0 0 0 id1 0 0 0 0 0 0 0 0;
%            0 0 0 0 1  0 0 0 0 0 0 0 0 0 0 0 0 0 id2 0 0 0 0 0 0 0;
%            0 0 0 0 0  1 0 0 0 0 0 0 0 0 0 0 0 0 0 id1 0 0 0 0 0 0;
%            0 0 0 0 0  1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id2 0 0 0 0 0;
%            0 0 0 0 0  0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id1 0 0 0 0;
%            0 0 0 0 0  0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id2 0 0 0;
%            0 0 0 0 0  0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id1 0 0;
%            0 0 0 0 0  0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id2 0];
%      Beq = zeros(1,16);
           
    [est, fval]= fmincon(model,init_vals,[],[],[],[],lb,ub,[],options);
    est_opts(itr,:) = est;
    fval_opts(itr,:) = fval;
    %%
    
    %Stopping Criteria
    if abs(fvalprev-fval) < 1e-5
        break;
    else
        fvalprev = fval;
    end
    
    if verbose; fprintf('\nOptimization Complete\n'); end;
    
end

% % %% from NAM
% n = 8;
% x_param4 = zeros(n^8,9, 'double');
% x = zeros(1,9);
% for idx1 = 1:n
%     x(1)=-pi + pi/4*(idx1-1);
%     for idx2 = 1:n
%         x(2)=-pi + pi/4*(idx2-1);
%         for idx3 = 1:n
%             x(3)=-pi + pi/4*(idx3-1);
%             for idx4 = 1
%                 x(4)=0;
%                 for idx5 = 1:n
%                     x(5)=-pi + pi/4*(idx5-1);
%                     for idx6 = 1:n
%                         x(6)=-pi + pi/4*(idx6-1);
%                         for idx7 = 1:n
%                             x(7)=-pi + pi/4*(idx7-1);
%                             for idx8 = 1:n
%                                 x(8)=-pi + pi/4*(idx8-1);
%                                 for idx9 = 1:n
%                                     x(9)=-pi + pi/4*(idx9-1);
%                                     x_param4((idx1-1)*n^8 + (idx2-1)*n^7 +  ...
%                                         (idx3-1)*n^6 + (idx4-1)*n^5 +  ...
%                                         (idx5- 1)*n^4 + (idx6- 1)*n^3 + ...
%                                         (idx7- 1)*n^2 + (idx8- 1)*n^1 + idx9,:) = x;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% %%
% 
% %profile on
% tic;
% parfor idx = 1:size(x_param4,1);
%    J(idx) = vsvlimFun([x_param4(idx,:) phi_inv], b1nom, Vc, ramp_time, grad_wait, ...
%             dtRF*1e-3, dtG, df, b1scale, v, T1, T2);
% end
% toc;



%% ------------------------------------------------------------------------
%Running Simulation with Different Parameters
%--------------------------------------------------------------------------
% x = [phi_sp phi_inv];  
% for ii=1:5
%     x(1)=-pi + pi/2*(ii-1);
%     for jj=1:5
%         x(2)=-pi + pi/2*(jj-1);
%         for kk=1:5
%             x(3)=-pi + pi/2*(kk-1);
%             for ll=1
%                 x(4)=0;
%                 for mm=1:5
%                     x(5)=-pi + pi/2*(mm-1);
%                     for nn=1:9
%                         x(6)=-pi + pi/4*(nn-1);
%                         for oo=1:9
%                             x(7)=-pi + pi/4*(oo-1);
%                             for pp=1:5
%                                 x(8)=-pi + pi/2*(pp-1);
%                                 for qq=1:5
%                                     x(9)=-pi + pi/2*(qq-1);   
%                                     J(ii,jj,kk,ll,mm,nn,oo,pp,qq) = ...
%                                     vsvlimFun(x, b1nom, Vc, ramp_time, grad_wait, ...
%                                     dtRF*1e-3, dtG, df, b1scale, v, T1, T2);           
%                                 end 
%                             end
%                         end                       
%                     end
%                 end
%             end
%         end
%         disp([ ' jj= ' num2str(kk) ' completed! ']);
%     end
%     disp([ ' ii= ' num2str(ii) ' completed! ']);
% 
% end
% 
% %%
% %load('J_2018_08_12') % <- Generate initial guesses!?
% idx = find(J == min(J(:)))
% 
% [I1,I2,I3,I4,I5,I6,I7,I8,I9] = ind2sub(size(J),idx);
% for itr = 1:size(I1)
% vals = [I1(itr),I2(itr),I3(itr),I4(itr),I5(itr),I6(itr),I7(itr),I8(itr),I9(itr)]
% end
