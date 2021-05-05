mu = 3.986e5;
a = 403 + 6378;
w = sqrt(mu/a^3);
A = [0 0 1 0;
    0 0 0 1;
    3*w^2 0 0 2*w;
    0 0 -2*w 0];

B = [1; 1; 1; 1];
C = eye(4);

T = 0:10:(2*pi/w);
sys = ss(A,B,C,0);
x0 = [50; 0; 0; 0];
[Y,T,X] = initial(sys,x0,T);

leaderPos = zeros(2,length(X));
for ii = 1:length(X)
   t = T(ii); 
   trans = [cos(t*w) -sin(t*w)
            sin(t*w)  cos(t*w)];
    
   leaderPos(:,ii) = trans*[a+X(ii,1);X(ii,2)]; 
       
   
end

figure 
plot(leaderPos(1,:),leaderPos(2,:))
axis equal

figure
plot(T,X(:,1))
