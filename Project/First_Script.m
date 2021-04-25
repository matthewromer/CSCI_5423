%Script to test some basic details of model

%Define Dynamics
w = 1;
A = [-0.1 1;
     0 0];
B = [0;1];
C = [1 0];

n = length(A);                   %Number of states in agent state-space
p = size(C,1);                   %Number of inputs/outputs in agent s-s

N = 2;                           %Number of nodes in network

% Gamma = [0 1 0;
%          1 0 1;
%          0 1 0]; 
% Deg =   [1  0  0;
%          0  2  0;
%          0  0  1];
% Lap = Deg-Gamma;

Gamma = [0 1;
         1 0];
Deg   = [1 0;
         0 1];
Lap = Deg-Gamma;     
     
ss_single = ss(A,B,C,0);   
QCtl = eye(n)*10;
RCtl = 1;
[K, ~, ~] = lqr(ss_single,QCtl,RCtl);
K = place(A,B,[-1,-2]);

Lt = place(A',C',[-0.01 -0.02]);
L = Lt';

%Create full system dynamcis matrices 
Atilde = [kron(eye(N),A-L*C), kron(Gamma,L*C);
          kron(eye(N),A) + kron(Deg,B*K), kron(eye(N),-B*K)];
Btilde = zeros(N*n*2,1);
Ctilde = [kron(eye(N),zeros(size(C))),kron(eye(N),C)];

ss_1 = ss(Atilde,Btilde,Ctilde,0);

eig_obs = eig(kron(eye(N),A-L*C))
eig_plant = eig(kron(eye(N),A)-kron(Lap,B*K))

[V,Diag] = eig(Atilde);
V
eig(Diag)

x_0 = [0.1 0 0.2 0];

ss_plant_only = ss(kron(eye(N),A)-kron(Lap,B*K),zeros(4,1),[1 0 0 0;0 0 1 0],0);

A_tilde_no_obs = (kron(eye(N),A)-kron(Lap,B*K))
%[V,Diag] = eig(A_tilde_no_obs)