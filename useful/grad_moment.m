function [moment_val,t] = grad_moment(grad,dt,order,inv_start,inv_dist,b1scale)
% Calculate the gradient moments
% [moment_vec,t] = grad_moment(grad,dt,order,inv_start,inv_dist,b1scale);

gamma = 42.58e6; % Hz/T
t = 0:dt:dt*(length(grad)-1); % milliseconds  
% moment_val = grad(1)*(t(1)/1000)^order/10000*dt/1000; % G/cm*(s)^order*1T/10000G*(s)
% 
% for itr=1:length(grad)
%     moment_vec(itr) = (moment_val + grad(itr)*(t(itr)/1000)^order/10000*dt/1000); % G/cm*ms*1s/1000ms*1T/10000G= T/cm*s*(s)^order
%     if (mod(itr,inv_dist)-inv_start)==0
%         moment_vec(itr) = -moment_vec(itr);
%     end
%     moment_val = moment_vec(itr);
% end
% excitation k-space: integrate from current point to the end of the pulse
moment_val = grad(end)*(t(end)/1000)^order/10000*dt/1000; % G/cm*(s)^order*1T/10000G*(s)

for itr=length(grad):-1:1
    moment_vec(itr) = (moment_val + grad(itr)*(t(itr)/1000)^order/10000*dt/1000); % G/cm*ms*1s/1000ms*1T/10000G= T/cm*s*(s)^order
    if (mod(itr,inv_dist)-inv_start)==0
        moment_vec(itr) = -moment_vec(itr);
    end
    moment_val = moment_vec(itr);
end

moment_val = gamma*moment_vec; % Hz/T*T/cm*s*(s)^order = (s)^order/cm

end

