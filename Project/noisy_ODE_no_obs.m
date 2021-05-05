function [xdot] = noisy_ODE(t,x,noise_hist,time_hist,A,B,C,Gamma,Deg,K,L)
%noisy_ODS contains 
%
% Inputs:  time, 
% 
% Outputs: (8N*1) vector containing the derivatives of all the states in x
% 
% Notes:   1. Derivatives are computed assuming a zero-order hold noise
%             input
%

numSt = 4;                  %Number of states 
N = length(x)/numSt;        %Number of agents

xdot = zeros(size(x));

dt = time_hist(2)-time_hist(1);
tvec_indicies = find(time_hist<=t+dt &  time_hist>t-dt);
w = noise_hist(:,tvec_indicies(1));    

%Loop over agents 
for ii = 1:N
    
    indSt = (ii-1)*numSt;
    x_i = x(indSt+1:indSt+numSt);    
    x_i_dot = A*x_i;
    %Loop over neighbors. Add noise to output measurement and derivatives
    for jj = 1:N
        if Gamma(ii,jj) == 1
            indSt_jj = (jj-1)*numSt;
            x_jj = x(indSt_jj+1:indSt_jj+numSt);
            y_jj = C*x_jj + w((jj-1)*2+1:(jj-1)*2+2);
            x_i_dot = x_i_dot + B*K*(y_jj-C*x_i) ;
        end
    end
        
    %Package into overall derivative
    xdot(indSt+1:indSt+numSt) = x_i_dot;
    
end






