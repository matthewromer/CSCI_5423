%Script to implement simulation of spacecraft

%Setup
clear all
close all

%System parameters
N = 4;
n = 4;

%Orbit parameters
mu = 3.986e5;
a = 403 + 6378;
w = sqrt(mu/a^3);

%Clohessy-Wiltshire equation SS model
A = [0     0  1   0;
     0     0  0   1;
     3*w^2 0  0   2*w;
     0     0 -2*w 0];
B = [0 0;
     0 0;
     1 0;
     0 1];
C = [1 0 0 0;
     0 1 0 0];

%Ring network adjacency and degree matrices
Gamma = [0 1 0 1;
         1 0 1 0;
         0 1 0 1;
         1 0 0 1];
Deg = diag([2 2 2 2]) ;

%Simulation time: just under one orbit
T = 0:10:(2*pi/w)*1-100;

%Choose system gains to stabilize each agent's dynamics
K = place(A,B,[-0.0015 -0.00225 -0.00175 -0.0025]);

%"Bad" observer design- stabilize observer but destabilize overall system
Lt = place(A',C',[-0.0002 -0.0001 -0.0021 -0.00011]);   
L = Lt';  

%Set noise level (measurement noise covariance matrix for a single agent's
%output)
Qsingle = [1 1];

%Get results for bad case
[XIndep,posIndep] = satelliteSimRun(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Independently');
 
%"Good" observer design- stabilize observer and overall system
Lt = place(A',C',[-0.2 -0.1 -0.21 -0.11]);             
L = Lt';

%Get results for good case
[XTogether,posTogether] = satelliteSimRun(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Together');

%Get results for case with agent 1 as a leader
[Xleader,posLeader] = satelliteSimRunLeaders(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,'Observer-Based Feedback with Leader');
%Increase noise level and get results
Qsingle = [1000 1000];
[XhighNoise,posHighNoise] = satelliteSimRun(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Together (High Noise)');

%Plot results for comparison of leader and leaderless systems
figure
subplot(2,2,1)
plot(T/60,XTogether(:,17),T/60,XTogether(:,21),T/60,XTogether(:,25),T/60,XTogether(:,29))
xlabel('Time (min)')
ylabel('X Deviation (km)')
title({'Position with Respect to Reference','Orbit, no Leaders'},'fontsize',14)
subplot(2,2,3)
plot(T/60,XTogether(:,18),T/60,XTogether(:,22),T/60,XTogether(:,26),T/60,XTogether(:,30))
xlabel('Time (min)')
ylabel('Y Deviation (km)')
legend('Satellite 1','Satellite 2','Satellite 3','Satellite 4','fontsize',12)
subplot(2,2,2)
plot(T/60,Xleader(:,17),T/60,Xleader(:,21),T/60,Xleader(:,25),T/60,Xleader(:,29))
xlabel('Time (min)')
ylabel('X Deviation (km)')
title({'Position with Respect to Reference','Orbit, with Leaders'},'fontsize',14)
subplot(2,2,4)
plot(T/60,Xleader(:,18),T/60,Xleader(:,22),T/60,Xleader(:,26),T/60,Xleader(:,30))
xlabel('Time (min)')
ylabel('Y Deviation (km)')
saveas(gcf, ['./Results/distToReference'], 'png')

function [X,pos] = satelliteSimRun(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,string)
%satelliteSimRun propagates the satellite dynamics without leaders
%Inputs:  Time vector, (A,B,C) matrices from SS model, controller gain 
%         matrix K, observer gain matrix L, adjacency matrix, degree matrx,
%         orbit rate w, semi-major axis a, descriptor string 
%
%Outputs: State history in rotating frame, state history in inertial frame.
%         Also generates plots of X-Y positions and average deviations
%         from center 

%Seed
rng(6)

%Compute system size
n = size(A,1);
N = size(L,1);
numTsteps = length(T);

%Create initial condition - just disperse positions (400 km covariance)
theta = [400 400 0 0]';
x0 = rand(n*N,1).*[theta;theta;theta;theta]-[theta;theta;theta;theta]/2;


Qtrue = diag([Qsingle Qsingle Qsingle Qsingle]);

%Create measurement noise history for true trajectory
noiseHist = zeros(8,numTsteps);
for ii = 1:numTsteps
   noiseHist(:,ii) =  mvnrnd(zeros(2*n,1),Qtrue)';
end

x0Full = [zeros(n*N,1);x0];

%Numerical state solution for trajectory 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,X] = ode45(@cwEqOde,T,x0Full,options,noiseHist,T,A,B,C,Gamma,Deg,K,L);

%Post-process to get position data and plots
[pos] = postProcessSim(T,X,w,a,N,string);

end

function [X,pos] = satelliteSimRunLeaders(T,A,B,C,K,L,Gamma,Deg,w,a,Qsingle,string)
%satelliteSimRunLeaders propagates the satellite dynamics with leaders
%Inputs:  Time vector, (A,B,C) matrices from SS model, controller gain 
%         matrix K, observer gain matrix L, adjacency matrix, degree matrx,
%         orbit rate w, semi-major axis a, descriptor string 
%
%Outputs: State history in rotating frame, state history in inertial frame.
%         Also generates plots of X-Y positions and average deviations
%         from center

%Seed
rng(6)

%Compute system size
n = size(A,1);
N = size(L,1);
num_tsteps = length(T);

%Select gains for leader
K2 = place(A,B,[-0.007 -0.007 -0.008 -0.009]);

%Create initial condition - just disperse positions (400 km covariance)
theta = [400 400 0 0]';
x0 = rand(n*N,1).*[theta;theta;theta;theta]-[theta;theta;theta;theta]/2;

Qtrue = diag([Qsingle Qsingle Qsingle Qsingle]);

%Create measurement noise history for true trajectory
noise_hist = zeros(8,num_tsteps);
for ii = 1:num_tsteps
   noise_hist(:,ii) =  mvnrnd(zeros(2*n,1),Qtrue)';
end

x0Full = [zeros(n*N,1);x0];

%Numerical state solution for trajectory 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,X] = ode45(@cwOdeLeaders,T,x0Full,options,noise_hist,T,A,B,C,Gamma,Deg,K,L,K2);

%Post-process to get position data and plots
[pos] = postProcessSim(T,X,w,a,N,string);

end

function [xdot] = cwEqOde(t,x,noiseHist,timeHist,A,B,C,Gamma,Deg,K,L)
%cwEqOde contains the equations of motion for a 4-agent
%Clohessy-Wiltshire equations
%
% Inputs:  Time, state vector, noise history vector, time history
%          corresponding to noise history, (A,B,C) matrices for SS model,
%          graph adjacency matrix, graph degree matrix, controller gain 
%          matrix K, observer gain matrix L
% 
% Outputs: (8N*1) vector containing the derivatives of all the states in x
% 
% Notes:   1. Derivatives are computed assuming a zero-order hold noise
%             input
%

numSt = 8;                  %Number of states 
N = length(x)/numSt;        %Number of agents

%Initialize xdot to zero
xdot = zeros(size(x));

%Get noise at this timestep
dt = timeHist(2)-timeHist(1);
tvecIndicies = find(timeHist<=t+dt &  timeHist>t-dt);
w = noiseHist(:,tvecIndicies(1));    

%Loop over agents 
for ii = 1:N
    
    %Extract this agent's state
    indEst = (ii-1)*numSt/2;
    indSt = (ii-1)*numSt/2 +numSt/2*N;
    xHatStar_ii = x(indEst+1:indEst+numSt/2);
    x_i = x(indSt+1:indSt+numSt/2);
    
    %Compute derivative of observer state
    xHatStar_ii_dot = (A-L*C)*xHatStar_ii;
    
    %Loop over neighbors. Add noise to output measurement and compute
    %adjustment to observer state derivatives
    for jj = 1:N
        if Gamma(ii,jj) == 1
            indSt_jj = (jj-1)*numSt/2 +numSt/2*N;
            x_jj = x(indSt_jj+1:indSt_jj+numSt/2);
            y_jj = C*x_jj + w((jj-1)*2+1:(jj-1)*2+2);
            xHatStar_ii_dot = xHatStar_ii_dot + L* y_jj;
        end
    end
    
    %Compute state derivatives
    x_i_dot = A*x_i + B*K*(xHatStar_ii-Deg(ii,ii)*x_i);
    
    %Package into overall derivative
    xdot(indEst+1:indEst+numSt/2) = xHatStar_ii_dot;
    xdot(indSt+1:indSt+numSt/2) = x_i_dot;
    
end

end

function [xdot] = cwOdeLeaders(t,x,noiseHist,timeHist,A,B,C,Gamma,Deg,K,L,K2)
%cwOdeLeaders contains the equations of motion for a 4-agent
%Clohessy-Wiltshire equations with leaders
%
% Inputs:  Time, state vector, noise history vector, time history
%          corresponding to noise history, (A,B,C) matrices for SS model,
%          graph adjacency matrix, graph degree matrix, controller gain 
%          matrix K, observer gain matrix L, leader gain matrix K2
% 
% Outputs: (8N*1) vector containing the derivatives of all the states in x
% 
% Notes:   1. Derivatives are computed assuming a zero-order hold noise
%             input
%          2. Agent 1 is a leader
%

numSt = 8;                  %Number of states 
N = length(x)/numSt;        %Number of agents

%Initialize xdot to zero
xdot = zeros(size(x));

%Get noise at this timestep
dt = timeHist(2)-timeHist(1);
tvecIndicies = find(timeHist<=t+dt &  timeHist>t-dt);
w = noiseHist(:,tvecIndicies(1));    

%Loop over agents 
for ii = 1:N
    
    %Extract this agent's state    
    indEst = (ii-1)*numSt/2;
    indSt = (ii-1)*numSt/2 +numSt/2*N;
    xHatStar_i = x(indEst+1:indEst+numSt/2);
    x_i = x(indSt+1:indSt+numSt/2);
    
    xHatStar_i_dot = (A-L*C)*xHatStar_i;
    %Loop over neighbors. Add noise to output measurement and compute observer
    %state derivatives
    for jj = 1:N
        if Gamma(ii,jj) == 1
            indSt_jj = (jj-1)*numSt/2 +numSt/2*N;
            x_jj = x(indSt_jj+1:indSt_jj+numSt/2);
            y_jj = C*x_jj + w((jj-1)*2+1:(jj-1)*2+2);
            xHatStar_i_dot = xHatStar_i_dot + L* y_jj;
        end
    end
    %Compute state derivatives
    
    x_i_dot = A*x_i + B*K*(xHatStar_i-Deg(ii,ii)*x_i);

    if ii == 1
        x_i_dot = x_i_dot - B*K2*x_i;
    end
    
    %Package into overall derivative
    xdot(indEst+1:indEst+numSt/2) = xHatStar_i_dot;
    xdot(indSt+1:indSt+numSt/2) = x_i_dot;
    
end

end

function [pos] = postProcessSim(T,X,w,a,N,string)
%postProcessSim computes and plots positions in the inertial frame from a
%simulation of the Clohessy-Wiltshire model
%   Inputs:  Time, state, orbit rate, semi-major axis, number of agents,
%            description string
%   Outputs: X-Y position histories for all agents in inertial frame

%Compute positions in X-Y based on simulation data
pos = zeros(2,length(X),N);
for ii = 1:length(X)
    t = T(ii);
    trans = [cos(t*w) -sin(t*w)   %Transformation from rotating to inertial
             sin(t*w)  cos(t*w)]; %coordinates
    
    for jj = 1:N
        pos(:,ii,jj) = trans*[a+X(ii,(jj+3)*4+1);X(ii,(jj+3)*4+2)];
    end
end

%Plot X-Y positions
figure
plot(squeeze(pos(1,:,1)),squeeze(pos(2,:,1)))
title({'Satellite X-Y Positions',string},'fontsize',14)
hold on
plot(squeeze(pos(1,:,2)),squeeze(pos(2,:,2)))
plot(squeeze(pos(1,:,3)),squeeze(pos(2,:,3)))
plot(squeeze(pos(1,:,4)),squeeze(pos(2,:,4)))
xlabel('X Position (km)','fontsize',14)
ylabel('Y Position (km)','fontsize',14)
legend('Satellite 1','Satellite 2','Satellite 3','Satellite 4','fontsize',14)
axis equal
saveas(gcf, ['./Results/XYPos',string], 'png')

%Compute average distance to center of group
averagePos = zeros(2,length(X));
averageDist= zeros(1,length(X));
for ii = 1:length(X)
    averagePos(:,ii) = (pos(:,ii,1) + pos(:,ii,2) + pos(:,ii,3) + pos(:,ii,4))/4;
    
    dist1 = norm(pos(:,ii,1)-averagePos(:,ii));
    dist2 = norm(pos(:,ii,2)-averagePos(:,ii));
    dist3 = norm(pos(:,ii,3)-averagePos(:,ii));
    dist4 = norm(pos(:,ii,4)-averagePos(:,ii));
    
    averageDist(ii) = mean([dist1,dist2,dist3,dist4]);
    
end

%Plot average distance to group center
figure
plot(T/60,averageDist)
title({'Average Disance of Satellites to Group Center',string},'fontsize',14)
xlabel('Time (min)','fontsize',14)
ylabel('Distance (km)','fontsize',14)
saveas(gcf, ['./Results/distToCenter',string], 'png')

end








