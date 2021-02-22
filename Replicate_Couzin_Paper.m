%Script to implement model from "Effective leadership and decision-making
%in animal groups on the move"
%
% Author: Matthew Romer
% Date:   2-21-2021

clearvars
close all

%Simulation parameters

%Scalar params
numIndivs = 50;
alpha = 1;
numTsteps = 500;
tStep = 0.2;
numLeaders = 10;
rho = 6;

%Vector params
g = [0;-1];
omegaMat = repmat(0.5,1,numIndivs);
isLeaderMat = zeros(1,numIndivs);
isLeaderMat(1:numLeaders) = ones(1,numLeaders);
sMat = repmat(alpha,1,numIndivs);
tVec = 0:tStep:numTsteps*tStep;

%Initial conditions
cMat = normrnd(0,2,2,numIndivs);
vMat = normrnd(0,0.5,2,numIndivs);
for i = 1:numIndivs
   vMat(:,i) = vMat(:,i)/norm(vMat(:,i));
end


%Preallocate
diTP1Mat = zeros(size(vMat));
cHistMat = zeros(2,numIndivs,numTsteps);
cHistMat(:,:,1) = cMat;
centroidPos = zeros(2,numTsteps);

%Loop over time 
for kk = 1:numTsteps
    
    %Compute directions of motion for this timestep
    for i = 1:numIndivs
        diTP1Mat(:,i) = indivModel(cMat, vMat, i, alpha, omegaMat, g, isLeaderMat,rho);
    end
    
    %Compute updated positions after this timestep
    for i = 1:numIndivs
        cMat(:,i) = cMat(:,i) + sMat(i) * diTP1Mat(:,i) * tStep;
    end
    
    %Compute centroid position at this timestep
    centroidPos(:,kk) = (sum(cMat')/numIndivs)';
    
    %Update history matrix 
    vMat = diTP1Mat;
    cHistMat(:,:,kk+1) = cMat;

end

%Plot results
% figure
% for i = 1:numIndivs
%     scatter(squeeze(cHistMat(1,i,:)),squeeze(cHistMat(2,i,:)))
%     hold on
% end
% legend('Agent 1','Agent 2','Agent 3')
% axis equal

figure
pointsize = 14;      %adjust as needed
pointidx = (1 : numTsteps)/tStep;
scatter(centroidPos(1,:),centroidPos(2,:), pointsize, pointidx);
colormap( jet(numTsteps) )
axis equal
xlabel('X Position')
ylabel('Y Position')
c = colorbar;
c.Label.String = 'Time (sec)';
title('Centroid Postion')


dir = (centroidPos(:,numTsteps)-centroidPos(:,numTsteps-50))/50/0.2
dTheta = acosd(dot(dir,g)/norm(dir)/norm(g))

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

else
    %No individuals close, prioritize alignment and attraction
    %Implement eqn 2    
    for jj = 1:numIndivs
        if jj ~= i && distMat(jj) < rho
            diTP1  = diTP1 + deltaCmat(:,jj)/norm(deltaCmat(:,jj)) + vMat(:,jj)/norm(vMat(:,jj));
        end
        
    end
    
end

diTP1 = diTP1/norm(diTP1);

%For leaders, implement eqn 3
if isLeaderMat(i)
    diTP1 = (diTP1 + omegaMat(:,i)*g)/ norm(diTP1 + omegaMat(:,i)*g);
    
end

%Limit the turn angle
diTP1 = limitTurnAngle(vMat(:,i),diTP1,maxAngle);

%Add noise to turn angle 
theta = normrnd(0,0.01);    %TODO: Pull from circular-wrapped gaussian

noiseRotMat = [cos(theta) -sin(theta);
               sin(theta)  cos(theta)];

diTP1 = noiseRotMat*diTP1;

end

function [limitedDir] = limitTurnAngle(currentDir,targetDir,maxAngle)
%limitTurnAngle updates the target direction based on a maximum turn angle
%   Inputs:  Current direction, initial target direction, maximum turn
%            angle (radians)
%   Outputs: Updated target direction with turn angle limited 
%

%Compute angle between current angle and origin and target angle and origin
currAng = atan2(currentDir(2),currentDir(1));
targAng = atan2(targetDir(2),targetDir(1));

%Compute difference in angles 
deltaAng = acos(dot(currentDir,targetDir)/norm(currentDir)/norm(targetDir));

%If angle is greater than the maximum, change the target direction to one
%that is no more than the maximum angle different from the initial
%direction
if deltaAng > maxAngle
    
    %Determine whether to turn CW or CCW
    q = sin(currAng-targAng);
    if q > 0
        dir = -1;
    else
        dir = 1;
    end
        
    %Form rotation matrix and turn target direction 
    theta = maxAngle*dir;
    mat = [cos(theta) -sin(theta);
           sin(theta)  cos(theta)];
         
    limitedDir = mat*currentDir;
else
    limitedDir = targetDir;
end

end