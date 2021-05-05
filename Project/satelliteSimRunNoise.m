function [X,pos] = satelliteSimRunNoise(A,B,C,K,L,Gamma,Deg,Lap,w,a,Qsingle,string)
%satelliteSimRun propagates the satellite dynamics 
%Inputs:  (A,B,C) matrices from SS model, controller gain matrix K, observer
%         gain matrix L, adjacency matrix, degree matrx, laplacian matrix,
%         orbit rate w, semi-major axis a, descriptor string 
%
%Outputs: None. Generates plots of X-Y positions and average deviations
%         from center 

rng(6)

n = size(A,1);
N = size(L,1);

%Create full system dynamcis matrices
Atilde = [kron(eye(N),A-L*C), kron(Gamma,L*C);
    kron(eye(N),B*K),kron(eye(N),A)-kron(Deg,B*K)];
Btilde = zeros(N*n*2,1);
Ctilde = [kron(eye(N),zeros(size(C))),kron(eye(N),C)];

%Create SS model object
sys = ss(Atilde,Btilde,Ctilde,0);

%Create initial condition - just disperse positions
theta = [10 10 0 0]'*40;
x0 = rand(n*N,1).*[theta;theta;theta;theta]-[theta;theta;theta;theta]/2;

%Create time vector
T = 0:10:(2*pi/w)*1-100;
num_tsteps = length(T);

Qtrue = diag([Qsingle Qsingle Qsingle Qsingle]);

%Create measurement noise history for true trajectory
noise_hist = zeros(8,num_tsteps);
for ii = 1:num_tsteps
   noise_hist(:,ii) =  mvnrnd(zeros(2*n,1),Qtrue)';
end

x0Full = [zeros(n*N,1);x0];

%Numerical state solution for trajectory 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,X] = ode45(@noisy_ODE,T,x0Full,options,noise_hist,T,A,B,C,Gamma,Deg,K,L);

%Compute positions in X-Y based on simulation data 
pos = zeros(2,length(X),4);
for ii = 1:length(X)
    t = T(ii);
    trans = [cos(t*w) -sin(t*w)   %Transformation from rotating to inertial
             sin(t*w)  cos(t*w)]; %coordinates
    
    for jj = 1:4    
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
saveas(gcf, ['XYPos',string], 'png')

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
saveas(gcf, ['distToCenter',string], 'png')
end
