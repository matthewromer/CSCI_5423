function [angularDev] = flockTrial(numTsteps,omega,alpha,numIndivs,tStep,numLeaders,rho,g)
%flockTrial runs one trial of the flocking logic and reports the angular
%deivation at the end of the trial
%   Inputs:  Number of timesteps, goal vector weight omega, keep-out radius
%            alpha, number of individuals, timestep, number of leaders,
%            local attention radius rho, target vector g
%   Outputs: Angular deviation at the end of the trial (deg)
%

%Create vector parameters
omegaMat = repmat(omega,1,numIndivs);
isLeaderMat = zeros(1,numIndivs);
isLeaderMat(1:numLeaders) = ones(1,numLeaders);
sMat = repmat(alpha,1,numIndivs);

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

dir = (centroidPos(:,numTsteps)-centroidPos(:,numTsteps-50))/50/0.2;
angularDev = acosd(dot(dir,g)/norm(dir)/norm(g));

if isnan(angularDev)
    fprint('here')
end