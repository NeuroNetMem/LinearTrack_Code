%% COMPUTE THETA SCORES



dir_w = {'UP';'DOWN'};

GLM_dir = '/home/fstella/Data/Matteo_GLMEarly/';
destination_save_1 = '/home/fstella/Data/Matteo_GLMEarly/';
destination_save_2 = '/home/fstella/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/Cluster_Info/';


th_scores_all = [];
for anm = [2 3 5]

 for field = 1:2   
    
     


for ss = 2:3

    for dire = 1:2
    
load([GLM_dir 'Animal' num2str(anm) '/StatsCSDLockingMPF_A' num2str(anm) '_S0' num2str(ss) '_F' num2str(field) '.mat'])
load([GLM_dir 'Animal' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{dire} '_A' num2str(anm) '_S0' num2str(ss) '_F' num2str(field) '.mat'])


figure(101)
clf;
for GG =1:2

t_c = p_cells(groups==GG);
s2 = StatisticsCSD.mvlCirc(t_c,:,dire);
s1 = StatisticsCSD.CC_Circ(t_c,:,dire);


for cont = 1:14
   scatter(ones(size(s1,1),1)*cont,(s1(:,cont)-s2(:,cont)),80,[1-(GG-1) 0 GG-1],'filled')
   hold on 
   %scatter(ones(size(s1,1),1)*cont,s2(:,cont),50,'b','filled')
    
    
end


th_score = mean(s1(:,6:12)-s2(:,6:12),2);
scatter(ones(size(s1,1),1)*15,th_score,80,[1-(GG-1) 0 GG-1],'filled')



end
xlabel('Channels')
ylabel('Precession - Locking')
refline(0,0)



t_c = p_cells;
s2 = StatisticsCSD.mvlCirc(t_c,:,dire,1);
s1 = StatisticsCSD.CC_Circ(t_c,:,dire,1);

theta_scores = nanmean(s1(:,6:12)-s2(:,6:12),2);
th_scores_all = cat(1,th_scores_all,theta_scores);

% figure(102)
% histfit(exp(theta_scores),15,'lognormal')
%histogram(theta_scores,-0.3:0.05:0.3,'Normalization','probability')


%pause()


save([destination_save_1 'Animal' num2str(anm) '/Cluster_GroupingsMPF_' dir_w{dire} '_A' num2str(anm) '_S0' num2str(ss) '_F' num2str(field) '.mat'],'theta_scores','-append')
save([destination_save_2  '/Cluster_GroupingsMPF_' dir_w{dire} '_A' num2str(anm) '_S0' num2str(ss) '_F' num2str(field) '.mat'],'theta_scores','-append')

    end
end
 end
end


%%

figure(102)
histfit(th_scores_all,35,'normal')



%%





