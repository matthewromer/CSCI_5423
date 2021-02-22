%Script to implement model from "Effective leadership and decision-making
%in animal groups on the move"
%
% Author: Matthew Romer
% Date:   2-22-2021

%Setup
clearvars
close all

%Scalar params
numIndivsMat = [10 30 50 100];
alpha = 1;
numTsteps = 500;
tStep = 0.2;
numLeadersMat = [1 2 4 6 8 10 15 20 30 40 50 75 100];
maxLeadersIndex = [6 9 11 13];
rho = 6;
omega = 0.5;
speed = alpha;

%Vector params
g = [0;-1];
colorMat = {'k','b','g',[0.9100    0.4100    0.1700]};

%Looping params
numTrials = 20;

%Preallocation
angularDev = zeros(length(numLeadersMat),numTrials,length(numIndivsMat));
speedMat   = zeros(length(numLeadersMat),numTrials,length(numIndivsMat)); 

tic

%Main Loop: For each number of individuals, iterate over number of leaders
%in group. For each number of leaders, run multiple trials to compute
%angular deviation of group from target direction and average speed of
%group
for kk = 1:length(numIndivsMat)
    numIndivs = numIndivsMat(kk);
    for ii = 1:maxLeadersIndex(kk)
        numLeaders = numLeadersMat(ii);
        
        %Allow user to see current progress in console
        fprintf('%i Individuals, %i leaders \n',numIndivs,numLeaders)

        for jj = 1:numTrials
            [angularDev(ii,jj,kk),speedMat(ii,jj,kk)] = flockTrial(numTsteps,omega,alpha,numIndivs,tStep,numLeaders,rho,g);
        end
        
    end
end


%Plot results. First figure: plot the normalized accuracy 
%Note: Formula for normalized accuracy is a best-guess, paper does not say
%what specific formula is. 
figure
for kk = 1:length(numIndivsMat)
    angularDevsKk = squeeze(angularDev(:,:,kk));
    angularDevsKk = angularDevsKk(1:maxLeadersIndex(kk),:);
    plot(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),1-mean(angularDevsKk')/90,'color',[colorMat{kk}],'linewidth',2)
    hold on
    scatter(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),1-mean(angularDevsKk')/90,20,[colorMat{kk}],'HandleVisibility','off')
    
end

ylim([0 1])
title({'Replication of Figure 1A. Directional Accuracy','as a Function of Leader Proportion'},'FontSize',14)
legend('N = 10','N = 30','N = 50','N = 100','Location','se','FontSize',14)
xlabel('Leader Proportion','FontSize',14)
ylabel('Scaled Directional Accuracy','FontSize',14)

%Plot results. Second figure: plot the normalized group speed
figure
for kk = 1:length(numIndivsMat)
    speedsKk = squeeze(speedMat(:,:,kk));
    speedsKk = speedsKk(1:maxLeadersIndex(kk),:);
    plot(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),mean(speedsKk')/speed,'color',[colorMat{kk}],'linewidth',2)
    hold on
    scatter(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),mean(speedsKk')/speed,20,[colorMat{kk}],'linewidth',2,'HandleVisibility','off')

end

ylim([0 1])
title('Average Group Speed as a Function of Leader Proportion','FontSize',14)
legend('N = 10','N = 30','N = 50','N = 100','Location','se','FontSize',14)
xlabel('Leader Proportion','FontSize',14)
ylabel('Normalized Group Speed','FontSize',14)


toc