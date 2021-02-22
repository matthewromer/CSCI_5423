function [diTP1] = indivModel(cMat, vMat, i, alpha, omegaMat, g, isLeaderMat, rho)
%indivModel implements the discre-time individual direction-selection model
%from "Effective leadership and decision- making in animal groups on the
%move" by Couzin et. al. (2005).
%
%   Inputs:  Matrix of individual positions (as 2x1 column vectors), matrix
%            of individual headings (as 2x1 column vectors), index of
%            current individual, distance at which individuals are
%            repulsed, weight leaders put on the diresired direction,
%            desired direction for leaders, matrix providing leader status
%            of each individual, distance at which individuals are
%            attracted
%   Outputs: Heading unit vector at t plus 1

numIndivs = size(cMat,2);
maxAngle = 2;               %Radians

%Compute distance to other indivs and distance magnitudes
deltaCmat = cMat - cMat(:,i);
distMat = zeros(1,numIndivs);
for jj = 1:numIndivs
    distMat(jj) = norm(deltaCmat(:,jj));
end
distMat(i) = inf;

diTP1 = [0;0];
if min(distMat) < alpha
    %An individual is close enough that need to prioritize avoidance. 
    %Implement eqn 1
    for jj = 1:numIndivs
        if jj ~= i && distMat(jj) < alpha
            diTP1  = diTP1 - deltaCmat(:,jj)/norm(deltaCmat(:,jj));
        end
    end

elseif min(distMat) < rho
    %No individuals close, prioritize alignment and attraction
    %Implement eqn 2    
    for jj = 1:numIndivs
        if jj ~= i && distMat(jj) < rho
            diTP1  = diTP1 + deltaCmat(:,jj)/norm(deltaCmat(:,jj)) + vMat(:,jj)/norm(vMat(:,jj));
        end
        
    end
else
    %No individuals close whatsoever. Continue on previous path.
    diTP1 = vMat(:,i);
    
end

diTP1 = diTP1/norm(diTP1);

%For leaders, implement eqn 3
if isLeaderMat(i)
    diTP1 = (diTP1 + omegaMat(:,i)*g)/ norm(diTP1 + omegaMat(:,i)*g);
    
end

%Limit the turn angle
diTP1 = limitTurnAngle(vMat(:,i),diTP1,maxAngle);

%Add noise to turn angle 
%NOTE: Using gaussian here, but paper uses circular-wrapped gaussian.
%Difference should be very small for this angle 
theta = normrnd(0,0.01);    

noiseRotMat = [cos(theta) -sin(theta);
               sin(theta)  cos(theta)];

diTP1 = noiseRotMat*diTP1;

end