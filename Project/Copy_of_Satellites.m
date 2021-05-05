%Script to implement simulation of spacecraft

clear all
close all

N = 4;
n = 4;

mu = 3.986e5;
a = 403 + 6378;
w = sqrt(mu/a^3);
A = [0 0 1 0;
    0 0 0 1;
    3*w^2 0 0 2*w;
    0 0 -2*w 0];

B = [0 0;
    0 0;
    1 0;
    0 1];

C = [1 0 0 0;
    0 1 0 0];

Gamma = [0 1 0 1;
         1 0 1 0;
         0 1 0 1;
         1 0 0 1];

Deg = diag([2 2 2 2]) ;
Lap = Deg-Gamma;

K = place(A,B,[-0.0015 -0.00225 -0.00175 -0.0025]);

Qsingle = [1 1];

Lt = place(A',C',[-0.0002 -0.0001 -0.0021 -0.00011]);  %Bad 
L = Lt';    
[XIndep,posIndep] = satelliteSimRunNoise(A,B,C,K,L,Gamma,Deg,Lap,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Independently');
 
Lt = place(A',C',[-0.2 -0.1 -0.21 -0.11]);              %Good
L = Lt';
[XTogether,posTogether] = satelliteSimRunNoise(A,B,C,K,L,Gamma,Deg,Lap,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Together');

Qsingle = [1000 1000];
[XhighNoise,posHighNoise] = satelliteSimRunNoise(A,B,C,K,L,Gamma,Deg,Lap,w,a,Qsingle,'Observer-Based Feedback - Gains Tuned Together (High Noise)');

Qsingle = [1 1];
[Xleader,posLeader] = satelliteSimRunLeaders(A,B,C,K,L,Gamma,Deg,Lap,w,a,Qsingle,'Observer-Based Feedback with Leader');

T = 0:10:(2*pi/w)*1-100;
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
saveas(gcf, ['distToReference'], 'png')


