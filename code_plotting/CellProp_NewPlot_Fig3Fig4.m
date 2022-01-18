clear
gro =0;
Sparsity_All = [];
Skaggs_All = [];
pfs_All = [];
MeanR_All = [];
MaxR_All = [];
MPos_All = [];
pfsasy_All = [];
theta_All = [];

data_origin = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';

for anm = [2 3 5]
for sess = 2:3

    
for direction = 1:2
    
load([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],...
    'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')    
%         
% load(['/Users/federico/GitHub/Basins_GLM_paper/Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],...
%     'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')   
    
Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 
pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>10 & pf_limits(:,2)<30)); 
theta_All = cat(1,theta_All,theta_scores);

end
end
end 

HH = histcounts(MPos_All,0:40,'Normalization','probability');
HH = smoothdata(HH,'movmean',5);

theta_all = zscore(theta_All);

figure(103)
subplot(2,3,1)
[cc_MaxR,pp_MaxR]=corr(theta_All(:),MaxR_All(:),'rows','complete','type','Spearman')
scatter(theta_All,MaxR_All*20,30,'filled')
% errorbar(gg,mean(MaxR_All)*20,std(MaxR_All)*20,'LineWidth',2)
title(['Max Rate, p = ' num2str(pp_MaxR) ', cc = ' num2str(cc_MaxR)])
lsline 
% xlim([0.5 2.5])
% ylabel('Hz')
% xticks([1 2])
% xticklabels({'Prec','Lock'})
% hold on


subplot(2,3,2)
[cc_MeanR,pp_MeanR]=corr(theta_All(:),MeanR_All(:),'rows','complete','type','Spearman')
scatter(theta_All,MeanR_All*5,30,'filled')
%errorbar(gg,mean(MeanR_All)*5,std(MeanR_All)*5,'LineWidth',2)
title(['Mean Rate, p = ' num2str(pp_MeanR) ', cc = ' num2str(cc_MeanR)])
lsline 
% xlim([0.5 2.5])
% ylabel('Hz')
% xticks([1 2])
% xticklabels({'Prec','Lock'})
% hold on


subplot(2,3,3)
[cc_Sparsity,pp_Sparsity]=corr(theta_All(:),Sparsity_All(:),'rows','complete','type','Spearman')
scatter(theta_All,Sparsity_All,30,'filled')
%errorbar(gg,mean(Sparsity_All),std(Sparsity_All),'LineWidth',2)
title(['Sparsity, p = ' num2str(pp_Sparsity) ', cc = ' num2str(cc_Sparsity)])
lsline 
%xlim([0.5 2.5])

% xticks([1 2])
% xticklabels({'Prec','Lock'})
% hold on


subplot(2,3,4)
[cc_Skaggs,pp_Skaggs]=corr(theta_All(:),Skaggs_All(:),'rows','complete','type','Spearman')
scatter(theta_All,Skaggs_All,30,'filled')
%errorbar(gg,mean(Skaggs_All),std(Skaggs_All),'LineWidth',2)
title(['Skaggs Info, p = ' num2str(pp_Skaggs) ', cc = ' num2str(cc_Skaggs)])
lsline 
% xlim([0.5 2.5])
% ylabel('bits')
% xticks([1 2])
% xticklabels({'Prec','Lock'})
% hold on


subplot(2,3,5)
[cc_pfs,pp_pfs]=corr(theta_All(:),pfs_All(:),'rows','complete','type','Spearman')
scatter(theta_All,pfs_All,30,'filled')
%errorbar(gg,mean(pfs_All),std(pfs_All),'LineWidth',2)
title(['Place Field Size, p = ' num2str(pp_pfs) ', cc = ' num2str(cc_pfs)])
lsline 
% xlim([0.5 2.5])
% ylabel('cm')
% xticks([1 2])
% xticklabels({'Prec','Lock'})
% hold on


subplot(2,3,6)
[cc_MPos,pp_MPos]=corr(MPos_All(:),theta_All(:),'rows','complete','type','Spearman')
scatter(MPos_All,theta_All,30,'filled')%plot((1:40)*2.5,HH,'LineWidth',2)
title(['Place Field Size, p = ' num2str(pp_MPos) ', cc = ' num2str(cc_MPos)])
lsline
% xlim([0 100])
% ylabel('p')
% xlabel('cm')

figure(301)
histogram(pfsasy_All,-10:1:10,'Normalization','probability')
hold on

mean(pfsasy_All)

figure(1000)
histfit(theta_All,35,'normal')

%% Saving

set(gcf,'renderer','Painters')
saveas(gcf,'Theta_score_Fig4','epsc')

%% Field 1 vs Field 2
dir_w = {'UP';'DOWN'};


kk=0;
th_corr_all = [];
for direction = 1:2
for anm = 1:6
for ss = 1:2 
if (anm==6 && ss == 2)
continue
end

kk = kk +1;
load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{direction} '_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
p1 = p_cells;

th1 = theta_scores;

s2 = StatisticsCSD.mvlCirc(:,:,1,1);
s1 = StatisticsCSD.CC_Circ(:,:,1,1);

s1 = mean(s1(:,:)-s2(:,:),2);
ss1 = abs(s1(p_cells));


load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F2.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{direction} '_A' num2str(anm) '_S' num2str(ss) '_F2.mat'])
p2 = p_cells;

th2 = theta_scores;
s2 = StatisticsCSD.mvlCirc(:,:,1,1);
s1 = StatisticsCSD.CC_Circ(:,:,1,1);

s2 = mean(s1(:,:)-s2(:,:),2);
ss2 = abs(s2(p_cells));


[pp,ia,ib] = intersect(p1,p2);

th1(isnan(th1)) = 0;
th2(isnan(th2)) = 0;

th_simil(kk)=corr(th1(ia),th2(ib));

th_corr_all = cat(1,th_corr_all,[th1(ia),th2(ib)]);

end
end
end

% plot(th_simil,'LineWidth',2)
% refline(0,0)
% xlabel('Sessions')
% ylabel('UP and DOWN Theta Score Correlation')
% switching = find((g1(ia)-g2(ib))~=0);
% no_switching = find((g1(ia)-g2(ib))==0);
% mean(ss1(ia))
% mean(ss2(ib))
% 
% mean(ss1(ia(switching)))
% mean(ss2(ib(switching)))
% 
% mean(ss1(ia(no_switching)))
% mean(ss2(ib(no_switching)))



[cc,pp]=corr(th_corr_all(:,1),th_corr_all(:,2),'type','Spearman')


figure(1)
subplot(1,2,1)
scatter(th_corr_all(:,1),th_corr_all(:,2),30,'filled')
lsline
xlabel('Theta Score FIRST FIELD')
ylabel('Theta Score SECOND FIELD')
title(['Same Direction Correlation r = ' num2str(cc) '; p-value = ' num2str(pp)])



%% Direction 1 vs Direction 2
dir_w = {'UP';'DOWN'};


kk=0;
th_corr_all = [];
for direction = 1:1
for anm = 1:6
for ss = 1:2 
    
if (anm==6 && ss == 2)
continue
end    
    
kk = kk +1;
load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{1} '_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
p1 = p_cells;

th1 = theta_scores;

s2 = StatisticsCSD.mvlCirc(:,:,1,1);
s1 = StatisticsCSD.CC_Circ(:,:,1,1);

s1 = mean(s1(:,:)-s2(:,:),2);
ss1 = abs(s1(p_cells));


load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{2} '_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
p2 = p_cells;

th2 = theta_scores;
s2 = StatisticsCSD.mvlCirc(:,:,1,1);
s1 = StatisticsCSD.CC_Circ(:,:,1,1);

s2 = mean(s1(:,:)-s2(:,:),2);
ss2 = abs(s2(p_cells));


[pp,ia,ib] = intersect(p1,p2);

th1(isnan(th1)) = 0;
th2(isnan(th2)) = 0;

th_simil(kk)=corr(th1(ia),th2(ib));

th_corr_all = cat(1,th_corr_all,[th1(ia),th2(ib)]);

end
end
end

% plot(th_simil,'LineWidth',2)
% refline(0,0)
% xlabel('Sessions')
% ylabel('UP and DOWN Theta Score Correlation')
% switching = find((g1(ia)-g2(ib))~=0);
% no_switching = find((g1(ia)-g2(ib))==0);
% mean(ss1(ia))
% mean(ss2(ib))
% 
% mean(ss1(ia(switching)))
% mean(ss2(ib(switching)))
% 
% mean(ss1(ia(no_switching)))
% mean(ss2(ib(no_switching)))

[cc,pp]=corr(th_corr_all(:,1),th_corr_all(:,2),'type','Spearman')


figure(1)
subplot(1,2,2)
scatter(th_corr_all(:,1),th_corr_all(:,2),30,'filled')
lsline
xlabel('Theta Score UP')
ylabel('Theta Score DOWN')
title(['Across Directions Correlation r = ' num2str(cc) '; p-value = ' num2str(pp)])

%% Theta Score difference vs Place field Distance 
clear
gro =0;
kk = 0; 
dir_w = {'UP';'DOWN'};
th_diff_all = [];
pf_distance_all = [];

th_corr_all = [];
for direction = 1:2
for anm = 1:6
for ss = 1:2 
if (anm==6 && ss == 2) || (anm==4 && ss == 2)
continue
end

kk = kk +1;
load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{direction} '_A' num2str(anm) '_S' num2str(ss) '_F1.mat'])
p1 = p_cells;

th1 = theta_scores;


load(['AnimalNG' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S' num2str(ss) '_F2.mat'])
load(['AnimalNG' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{direction} '_A' num2str(anm) '_S' num2str(ss) '_F2.mat'])
p2 = p_cells;

th2 = theta_scores;


[pp,ia,ib] = intersect(p1,p2);
th1 = th1(ia);
th2 = th2(ib);
% th1(isnan(th1)) = 0;
% th2(isnan(th2)) = 0;

 th_simil=abs(th1-th2);
 th_diff_all = cat(1,th_diff_all,th_simil);


 load(['place_field_limits/PF_limit_A' num2str(anm) '_S' num2str(ss)])

    if direction == 1
    pf_limits = pf_limit_upMPF;
    elseif direction == 2
    pf_limits = -pf_limit_downMPF;
    
    end
    
    
pf_distance = abs(pf_limits(pp,2,1) - pf_limits(pp,2,2));
pf_distance_all = cat(1,pf_distance_all,pf_distance);


end
end
end

%%%


[cc,pp]=corr(pf_distance_all(:),th_diff_all(:),'type','Spearman','rows','complete')


figure(10),clf
scatter(pf_distance_all(:),th_diff_all(:),30,'filled')
lsline
ylabel('Theta Score difference')
xlabel('Place Field distance')
title(['Theta Score difference vs Place Field distance Correlation r = ' num2str(cc) '; p-value = ' num2str(pp)])














