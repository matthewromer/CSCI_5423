function [xdot] = noisy_ODE(t,x,noise_hist,time_hist,A,B,C,Gamma,Deg,K,L,K2)
%noisy_ODS contains 
%
% Inputs:  time, 
% 
% Outputs: (8N*1) vector containing the derivatives of all the states in x
% 
% Notes:   1. Derivatives are computed assuming a zero-order hold noise
%             input
%

numSt = 8;                  %Number of states 
N = length(x)/numSt;        %Number of agents

xdot = zeros(size(x));

dt = time_hist(2)-time_hist(1);
tvec_indicies = find(time_hist<=t+dt &  time_hist>t-dt);
w = noise_hist(:,tvec_indicies(1));    

%Loop over agents 
for ii = 1:N
    
    indEst = (ii-1)*numSt/2;
    indSt = (ii-1)*numSt/2 +numSt/2*N;
    
    x_hat_star_i = x(indEst+1:indEst+numSt/2);
    x_i = x(indSt+1:indSt+numSt/2);
    
    x_hat_star_i_dot = (A-L*C)*x_hat_star_i;
    %Loop over neighbors. Add noise to output measurement and compute observer
    %state derivatives
    for jj = 1:N
        if Gamma(ii,jj) == 1
            indSt_jj = (jj-1)*numSt/2 +numSt/2*N;
            x_jj = x(indSt_jj+1:indSt_jj+numSt/2);
            y_jj = C*x_jj + w((jj-1)*2+1:(jj-1)*2+2);
            x_hat_star_i_dot = x_hat_star_i_dot + L* y_jj;
        end
    end
    %Compute state derivatives
    
    x_i_dot = A*x_i + B*K*(x_hat_star_i-Deg(ii,ii)*x_i);

    if ii == 1
        x_i_dot = x_i_dot - B*K2*x_i;
    end
    
    %Package into overall derivative
    xdot(indEst+1:indEst+numSt/2) = x_hat_star_i_dot;
    xdot(indSt+1:indSt+numSt/2) = x_i_dot;
    
end






