%Script to support discussion of model properties

%Define Dynamics
w = 1;
A = [-0.1 1;
     0 0];
B = [0;1];
C = [1 0];

n = length(A);                   %Number of states in agent state-space
p = size(C,1);                   %Number of inputs/outputs in agent s-s

N = 2;                           %Number of nodes in network

%Define graph. Two-agent toy system 
Gamma = [0 1;
         1 0];
Deg   = [1 0;
         0 1];
Lap = Deg-Gamma;     
     
%Place controller and observer poles
K = place(A,B,[-1,-2])

Lt = place(A',C',[-0.01 -0.02]);
L = Lt'

%Compute eigenvalues of controller and observer individually 
eig_obs = eig(kron(eye(N),A-L*C))
eig_control = eig(kron(eye(N),A)-kron(Lap,B*K))

%Create full system dynamcis matrices 
Atilde = [kron(eye(N),A-L*C), kron(Gamma,L*C);
          kron(eye(N),B*K),kron(eye(n),A)-kron(Deg,B*K)];
Btilde = zeros(N*n*2,1);
Ctilde = [kron(eye(N),zeros(size(C))),kron(eye(N),C)];

%Compute eigenvalues of full system A matrix 
[V,Diag] = eig(Atilde);
diag(Diag)
V

Kval = 1;

%Define set of leaders
G = zeros(N,Kval); 
    
%Create G, First m agents are leaders
for i = 1:Kval
   G(i,i) = 1;
end

%System without observers
Btilde = kron(G,B);
rank(ctrb(kron(eye(N),A)-kron(Lap,B*C),Btilde))

%System with observers
G = [1; 0];
Btilde = [kron(G,zeros(size(B)));kron(G,B)];
rank(ctrb(Atilde,Btilde))

%System with observers, different L value
Lt = place(A',C',[-0.1 -2]);
L = Lt';
Atilde = [kron(eye(N),A-L*C), kron(Gamma,L*C);
          kron(eye(N),B*K),kron(eye(n),A)-kron(Deg,B*K)];
rank(ctrb(Atilde,Btilde))
 
