%Script to plot stored results
%
% Author: Matthew Romer
% Date:   2-22-2021

%Load stored data
clearvars
close all
load('Results_100_Trials.mat')

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
speed = alpha;

figure
for kk = 1:length(numIndivsMat)
    speedsKk = squeeze(speedMat(:,:,kk));
    speedsKk = speedsKk(1:maxLeadersIndex(kk),:);
    plot(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),mean(speedsKk'),'color',[colorMat{kk}],'linewidth',2)
    hold on
    scatter(numLeadersMat(1:maxLeadersIndex(kk))/numIndivsMat(kk),mean(speedsKk'),20,[colorMat{kk}],'linewidth',2,'HandleVisibility','off')

    
end

ylim([0 1])
title('Average Group Speed as a Function of Leader Proportion','FontSize',14)
legend('N = 10','N = 30','N = 50','N = 100','Location','se','FontSize',14)
xlabel('Leader Proportion','FontSize',14)
ylabel('Normalized Group Speed','FontSize',14)
