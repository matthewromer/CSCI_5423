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

K = place(A,B,[-0.001 -0.00125 -0.0015 -0.00175]);

Lt = place(A',C',[-0.2 -0.1 -0.21 -0.11]);              %Good
%Lt = place(A',C',[-0.0002 -0.0001 -0.0021 -0.00011]);  %Bad 
L = Lt';

satelliteSimRun(A,B,C,K,L,Gamma,Deg,Lap,w,a,'Observer-Based Feedback - Gains Tuned Together')

