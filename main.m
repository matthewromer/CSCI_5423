%Script to implement model from "Effective leadership and decision-making
%in animal groups on the move"
%
% Version that uses functions
%
% Author: Matthew Romer
% Date:   2-21-2021

%Setup
clearvars
close all

%Scalar params
numIndivsMat = [10 30 50];
alpha = 1;
numTsteps = 500;
tStep = 0.2;
numLeadersMat = [1 2 4 6 8 10 15 20 30 40 50];
maxLeadersIndex = [6 9 11];
rho = 6;
omega = 0.5;

%Vector params
g = [0;-1];
colorMat = ['k','b','g'];

%Looping params
numTrials = 100;

tic
angularDev = zeros(length(numLeadersMat),numTrials,length(numIndivsMat));
speedMat   = zeros(length(numLeadersMat),numTrials,length(numIndivsMat)); 
for kk = 1:length(numIndivsMat)
    numIndivs = numIndivsMat(kk);
    kk
    for ii = 1:maxLeadersIndex(kk)
        numLeaders = numLeadersMat(ii);
        ii
        for jj = 1:numTrials
            [angularDev(ii,jj,kk),speedMat(ii,jj,kk)] = flockTrial(numTsteps,omega,alpha,numIndivs,tStep,numLeaders,rho,g);
        end
        
    end
end


figure
for kk = 1:length(numIndivsMat)
    angularDevsKk = squeeze(angularDev(:,:,kk));
    angularDevsKk = angularDevsKk(1:maxLeadersIndex(kk),:);
    plot(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),1-mean(angularDevsKk')/90,'color',[colorMat(kk)],'linewidth',2)
    hold on
    
end

ylim([0 1])
title('Replication of Figure 1')
legend('N = 10','N = 30','N = 50')

figure
for kk = 1:length(numIndivsMat)
    speedsKk = squeeze(speedMat(:,:,kk));
    speedsKk = speedsKk(1:maxLeadersIndex(kk),:);
    plot(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),mean(speedsKk'),'color',[colorMat(kk)],'linewidth',2)
    hold on
    
end

ylim([0 1])
title('Average group speed')
legend('N = 10','N = 30','N = 50')

toc