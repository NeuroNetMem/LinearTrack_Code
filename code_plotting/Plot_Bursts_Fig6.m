
close all

%% PANEL - BURST PROBABILITY AGAINST GAMMA

load('Theta_Limits.mat')
G_Lim = th_limits;



typ =1 ;

kk = 0; 


cc_ss_all = [];
cc_sb_all = [];
cc_lb_all = [];
burst_per_all = [];
theta_all = [];
for anm =[1 2 3 4 5 6]
for sess = 1:2
    
 if (anm==6 && sess == 2)
continue
end  
  
    kk = kk+1;
    
   if typ == 1
       
load(['BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga2.mat'])
   elseif typ == 2
    
load(['BurstDati/BurstAnalysisREM_A' num2str(anm) '_S' num2str(sess) '.mat'])

   end
cc_ss_all = cat(3,cc_ss_all,smoothdata(cc_ss,1,'gaussian',3));
cc_sb_all = cat(3,cc_sb_all,smoothdata(cc_sb,1,'gaussian',3));
cc_lb_all = cat(3,cc_lb_all,smoothdata(cc_ab,1,'gaussian',3));
%burst_per_all = burst_per_all + burst_per;
burst_per_all = cat(3,burst_per_all,burst_per);
theta_all = cat(1,theta_all,theta_session);
end


end


th_group = discretize(theta_all,G_Lim);
gr_list = unique(th_group);


figure(20+typ)

sgtitle('Burst Prominence')
subplot(221)
title('Short Bursts')
for gr = 1:numel(G_Lim)-1
errorbar(squeeze(nanmean(burst_per_all(1,:,th_group==gr),3)),squeeze(nanstd(burst_per_all(1,:,th_group==gr),[],3))/sqrt(sum(th_group==gr)),'LineWidth',3)
hold on
end

title('Short Bursts')

legend('Group 1','Group 2','Group 3')
ylabel('Burst Percentage')
xlabel('Gamma Balance: Medium->Slow')

subplot(223)
for gr = 1:numel(G_Lim)-1
errorbar(squeeze(nanmean(burst_per_all(2,:,th_group==gr),3)),squeeze(nanstd(burst_per_all(2,:,th_group==gr),[],3))/sqrt(sum(th_group==gr)),'LineWidth',3)
hold on
end
title('All Bursts')

legend('Group 1','Group 2','Group 3')
ylabel('Burst Percentage')
xlabel('Gamma Balance: Medium->Slow')


if typ == 1
sgtitle('Burst Probability RUN')
elseif typ == 2
sgtitle('Burst Probability REM')    
end


%Stats 1
% 1Vs3 
a1 = squeeze(burst_per_all(2,:,th_group==1));
a3 = squeeze(burst_per_all(2,:,th_group==3));

for g1 = 1:6
[h1(g1),p1(g1)]=ttest2(a1(g1,:),a3(g1,:));
end


% 3Vs3 Short
a1 = squeeze(burst_per_all(1,:,th_group==1));
a3 = squeeze(burst_per_all(1,:,th_group==3));

for g1 = 1:6
[h2(g1),p2(g1)]=ttest2(a3(1,:),a3(g1,:));
end

% 3Vs3 All
a1 = squeeze(burst_per_all(2,:,th_group==1));
a3 = squeeze(burst_per_all(2,:,th_group==3));
for g1 = 1:6
[h3(g1),p3(g1)]=ttest2(a3(1,:),a3(g1,:));
end


figure(20+typ)
subplot(222)
for gr = 1:numel(G_Lim)-1
errorbar(squeeze(nanmean(burst_per_all(1,:,th_group==gr),[2,3])),squeeze(nanstd(burst_per_all(1,:,th_group==gr),[],[2,3]))/sqrt(sum(th_group==gr)*3),'LineWidth',3)
hold on
end
title('Short Bursts - Overall')

subplot(224)
for gr = 1:numel(G_Lim)-1
errorbar(squeeze(nanmean(burst_per_all(2,:,th_group==gr),[2,3])),squeeze(nanstd(burst_per_all(2,:,th_group==gr),[],[2,3]))/sqrt(sum(th_group==gr)*3),'LineWidth',3)

ylim([0.15 0.24])

hold on
end
title('Long Bursts - Overall')

%Stats 1
% 1Vs3 
%Locking vs Precessing 
mean_a1 = mean(squeeze(burst_per_all(2,:,th_group==1)));
mean_a3 = mean(squeeze(burst_per_all(2,:,th_group==3)));

[mean_h1,mean_p1]=ttest2(mean_a1,mean_a3)



% %% Saving
% figure(21)
% set(gcf,'renderer','Painters')
% saveas(gcf,'Fig6_BurstOverall','epsc')


%% PANEL - BURST PROBABILITY AGAINST THETA PHASE

load('Theta_Limits.mat')
G_Lim = th_limits;



typ =1 ;

kk = 0; 


cc_ss_all = [];
cc_sb_all = [];
cc_lb_all = [];
burst_per_all = [];
theta_all = [];
for anm =[1 2 3 4 5 6]

for sess = 1:2
if (anm==6 && sess == 2)
continue
end    
  
    kk = kk+1;
    
   if typ == 1
       
load(['BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga2RR.mat'])
       
   elseif typ == 2
    
load(['BurstDati/BurstAnalysisREM_A' num2str(anm) '_S' num2str(sess) '.mat'])

   end
cc_ss_all = cat(3,cc_ss_all,smoothdata(cc_ss,1,'gaussian',3));
cc_sb_all = cat(3,cc_sb_all,smoothdata(cc_sb,1,'gaussian',3));
cc_lb_all = cat(3,cc_lb_all,smoothdata(cc_ab,1,'gaussian',3));
%burst_per_all = burst_per_all + burst_per;
burst_per_all = cat(3,burst_per_all,burst_per);
theta_all = cat(1,theta_all,theta_session);
end


end


th_group = discretize(theta_all,G_Lim);
gr_list = unique(th_group);


for ga = [1 6]

figure(100*typ+ga)
for gro = 1:numel(G_Lim)-1




cc = nanmean(cc_ss_all(:,:,th_group==gro),3);
cc_s = nanstd(cc_ss_all(:,:,th_group==gro),[],3);

cc = nanmean(cc(:,ga),2);
cc_s = cc_s(:,ga);

cc_s = cc_s./nansum(cc,'all');
cc = cc./nansum(cc,'all');
cc = smoothdata(repmat(cc,3,1),1,'gaussian',5);
cc_s = repmat(cc_s,3,1);
subplot(3,1,1)
for l = 1:size(cc,2)
 %plot(cc(:,l),'Color',[1-l/size(cc,2) 0 l/size(cc,2)],'LineWidth',3) 
  errorbar(cc(:,l),cc_s(:,l)/20,'LineWidth',3) 



hold on
end
%ylim([0 0.04])



cc = nanmean(cc_sb_all(:,:,th_group==gro),3);
cc_s = nanstd(cc_sb_all(:,:,th_group==gro),[],3);

cc = nanmean(cc(:,ga),2);
cc_s = cc_s(:,ga);

cc_s = cc_s./nansum(cc,'all');
cc = cc./nansum(cc,'all');
cc = smoothdata(repmat(cc,3,1),1,'gaussian',5);
cc_s = repmat(cc_s,3,1);
subplot(3,1,2)
for l = 1:size(cc,2)
 %plot(cc(:,l),'Color',[1-l/size(cc,2) 0 l/size(cc,2)],'LineWidth',3) 
  errorbar(cc(:,l),cc_s(:,l)/20,'LineWidth',3)  

hold on
end
%ylim([0 0.04])


cc = nanmean(cc_lb_all(:,:,th_group==gro),3);
cc_s = nanstd(cc_lb_all(:,:,th_group==gro),[],3);

cc = nanmean(cc(:,ga),2);
cc_s = cc_s(:,ga);

cc_s = cc_s./nansum(cc,'all');
cc = cc./nansum(cc,'all');




cc = smoothdata(repmat(cc,3,1),1,'gaussian',5);
cc_s = repmat(cc_s,3,1);
subplot(3,1,3)
for l = 1:size(cc,2)

 %plot(cc(:,l),'Color',[1-l/size(cc,2) 0 l/size(cc,2)],'LineWidth',3) 
  errorbar(cc(:,l),cc_s(:,l)/20,'LineWidth',3) 

hold on
end
%ylim([0 0.04])


% if typ == 1
% sgtitle('Spike Phase Probability RUN')
% elseif typ == 2
% sgtitle('Spike Phase Probability REM')    
% end

end
subplot(3,1,1)
plot([20 20],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
title('Single Spikes')
xlabel('Theta Phase')
ylabel('Spike Prob')


subplot(3,1,2)
plot([20 20],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
title('Short Burst Spikes')
xlabel('Theta Phase')
ylabel('Spike Prob')
subplot(3,1,3)
plot([20 20],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
plot([40 40],[0.005 0.1],'--k')
title('All Burst Spikes (n>3)')
xlabel('Theta Phase')
ylabel('Spike Prob')


if typ == 1 && ga==1
sgtitle({'Spike Distribution RUN','Medium Gamma'})
elseif typ == 1 && ga == 6
sgtitle({'Spike Distribution RUN','Slow Gamma'})

elseif typ == 2 && ga==1
sgtitle({'Spike Distribution REM','Medium Gamma'})
elseif typ == 2 && ga == 6
sgtitle({'Spike Distribution REM','Slow Gamma'})
end

end



 % Single Spikes 1Vs3 MG
a1 = squeeze(cc_ss_all(:,1,th_group==1));
a3 = squeeze(cc_ss_all(:,1,th_group==3));
for c = 1:size(a1,2)
a1(:,c) = a1(:,c)./sum(a1(:,c));
end
for c = 1:size(a3,2)
a3(:,c) = a3(:,c)./sum(a3(:,c));
end
for g1 = 1:20
[h11(g1),p11(g1)]=ttest2(a1(g1,:),a3(g1,:));
end


clear a1
clear a3
% Single Spikes 1Vs3 SG
a1 = squeeze(cc_ss_all(:,6,th_group==1));
a3 = squeeze(cc_ss_all(:,6,th_group==3));
for c = 1:size(a1,2)
a1(:,c) = a1(:,c)./sum(a1(:,c));
end
for c = 1:size(a3,2)
a3(:,c) = a3(:,c)./sum(a3(:,c));
end
for g1 = 1:20
[h12(g1),p12(g1)]=ttest2(a1(g1,:),a3(g1,:));
end


clear a1
clear a3
% Long Bursts 1Vs3 MG
a1 = squeeze(cc_lb_all(:,1,th_group==1));
a3 = squeeze(cc_lb_all(:,1,th_group==3));
for c = 1:size(a1,2)
a1(:,c) = a1(:,c)./sum(a1(:,c));
end
for c = 1:size(a3,2)
a3(:,c) = a3(:,c)./sum(a3(:,c));
end
for g1 = 1:20
[h15(g1),p15(g1)]=ttest2(a1(g1,:),a3(g1,:));
end


clear a1
clear a3
% Long Bursts 1Vs3 SG
a1 = squeeze(cc_lb_all(:,6,th_group==1));
a3 = squeeze(cc_lb_all(:,6,th_group==3));
for c = 1:size(a1,2)
a1(:,c) = a1(:,c)./sum(a1(:,c));
end
for c = 1:size(a3,2)
a3(:,c) = a3(:,c)./sum(a3(:,c));
end
for g1 = 1:20
[h16(g1),p16(g1)]=ttest2(a1(g1,:),a3(g1,:));
end


clear a1
clear a3
% Long Bursts 3Vs3 SG+MG
a1 = squeeze(cc_lb_all(:,1,th_group==3));
a3 = squeeze(cc_lb_all(:,6,th_group==3));
for c = 1:size(a1,2)
a1(:,c) = a1(:,c)./sum(a1(:,c));
end
for c = 1:size(a3,2)
a3(:,c) = a3(:,c)./sum(a3(:,c));
end
for g1 = 1:20
[h19(g1),p19(g1)]=ttest2(a1(7,:),a1(g1,:));
[h20(g1),p20(g1)]=ttest2(a3(8,:),a3(g1,:));
end


% %% Saving
% figure(101)
% set(gcf,'renderer','Painters')
% saveas(gcf,'Fig6_MG','epsc')
% 
% %%Saving
% figure(106)
% set(gcf,'renderer','Painters')
% saveas(gcf,'Fig6_SG','epsc')

%% PANEL PREFERRED PHASE

for typ = 1:2
figure(5000 + typ )
clf; 

for gro=1:numel(G_Lim)-1
Spike_t=[];
for anm =[1 2 3 4 5 6]
for sess = 1:2
    
if (anm==6 && sess == 2)
continue
end  
  
    kk = kk+1;
if typ == 1
load(['BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga2.mat'])
elseif typ == 2
   load(['BurstDati/BurstAnalysisREM_A' num2str(anm) '_S' num2str(sess) '.mat']) 
    
end
th_group = discretize(theta_session,G_Lim);
Spike_t = cat(1,Spike_t,CellPref_Angle(th_group==gro));

end
end

Spike_t(Spike_t<0)=Spike_t(Spike_t<0)+2*pi;


cc=histcounts(Spike_t,linspace(0,2*pi,21));
cc = smoothdata(repmat(cc/sum(cc),1,3),'movmean',3);
plot(cc,'LineWidth',3)
hold on

end
ylim([0 0.12])
plot([20 20],[0 0.1],'--k')
plot([40 40],[0 0.1],'--k')

if typ == 1
title('Preferred Phase Probability RUN')
elseif typ == 2
title('Preferred Phase Probability REM')    
end
xlabel('Theta Phase')
ylabel('Pref Phase Prob')
legend 

end

% %%Saving
% figure(5001)
% set(gcf,'renderer','Painters')
% saveas(gcf,'Fig6_PhasePref','epsc')
% 
% figure(5002)
% set(gcf,'renderer','Painters')
% saveas(gcf,'Fig6_PhasePrefREM','epsc')

%%
%% PANEL - BURST PROBABILITY AGAINST THETA PHASE - COMPARE WITH RANDOM
close all
load('Theta_Limits.mat')
G_Lim = th_limits;



typ =1 ;

kk = 0; 


cc_ss_all = [];
cc_sb_all = [];
cc_lb_all = [];
burst_per_all = [];
theta_all = [];
for anm =[1 2 3 4 5 6]

for sess = 1:2
if (anm==6 && sess == 2)
continue
end    
  
    kk = kk+1;
    

load(['BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga2.mat'])
cc_ss_1 = cc_ss;
cc_sb_1 = cc_sb;
cc_ab_1 = cc_ab;
       
load(['BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga2RR.mat'])
cc_ss_2 = nanmean(cc_ss,4);
cc_sb_2 = nanmean(cc_sb,4);
cc_ab_2 = nanmean(cc_ab,4);

cc_ss_d = cc_ss_1 - cc_ss_2;
cc_sb_d = cc_sb_1 - cc_sb_2;
cc_ab_d = cc_ab_1 - cc_ab_2;

cc_ss_all = cat(3,cc_ss_all,smoothdata(cc_ss_d,1,'gaussian',3));
cc_sb_all = cat(3,cc_sb_all,smoothdata(cc_sb_d,1,'gaussian',3));
cc_lb_all = cat(3,cc_lb_all,smoothdata(cc_ab_d,1,'gaussian',3));
%burst_per_all = burst_per_all + burst_per;

theta_all = cat(1,theta_all,theta_session);
end


end


th_group = discretize(theta_all,G_Lim);
gr_list = unique(th_group);


for ga = [1 6]

figure(100*typ+ga)
for gro = 1:numel(G_Lim)-1




cc = nanmean(cc_ss_all(:,:,th_group==gro),3);
cc_s = nanstd(cc_ss_all(:,:,th_group==gro),[],3);

cc = nanmean(cc(:,ga),2);
cc_s = cc_s(:,ga);

% cc_s = cc_s./nansum(cc,'all');
% cc = cc./nansum(cc,'all');
cc = smoothdata(repmat(cc,3,1),1,'gaussian',5);
cc_s = repmat(cc_s,3,1);
subplot(2,1,1)
for l = 1:size(cc,2)
 %plot(cc(:,l),'Color',[1-l/size(cc,2) 0 l/size(cc,2)],'LineWidth',3) 
  errorbar(cc(:,l),cc_s(:,l)/20,'LineWidth',3) 



hold on
end
%ylim([0 0.04])






cc = nanmean(cc_lb_all(:,:,th_group==gro),3);
cc_s = nanstd(cc_lb_all(:,:,th_group==gro),[],3);

cc = nanmean(cc(:,ga),2);
cc_s = cc_s(:,ga);

% cc_s = cc_s./nansum(cc,'all');
% cc = cc./nansum(cc,'all');




cc = smoothdata(repmat(cc,3,1),1,'gaussian',5);
cc_s = repmat(cc_s,3,1);
subplot(2,1,2)
for l = 1:size(cc,2)

 %plot(cc(:,l),'Color',[1-l/size(cc,2) 0 l/size(cc,2)],'LineWidth',3) 
  errorbar(cc(:,l),cc_s(:,l)/20,'LineWidth',3) 

hold on
end
%ylim([0 0.04])


% if typ == 1
% sgtitle('Spike Phase Probability RUN')
% elseif typ == 2
% sgtitle('Spike Phase Probability REM')    
% end

end
subplot(2,1,1)
plot([20 20],[-1 1]*0.00005,'--k')
plot([40 40],[-1 1]*0.00005,'--k')
plot([40 40],[-1 1]*0.00005,'--k')
title('Single Spikes')
xlabel('Theta Phase')
ylabel('Delta Spike Prob')



subplot(2,1,2)
plot([20 20],[-1 1]*0.00005,'--k')
plot([40 40],[-1 1]*0.00005,'--k')
plot([40 40],[-1 1]*0.00005,'--k')
title('All Burst Spikes (n>3)')
xlabel('Theta Phase')
ylabel('Delta Spike Prob')


if typ == 1 && ga==1
sgtitle({'Spike Distribution RUN','Medium Gamma'})
elseif typ == 1 && ga == 6
sgtitle({'Spike Distribution RUN','Slow Gamma'})

elseif typ == 2 && ga==1
sgtitle({'Spike Distribution REM','Medium Gamma'})
elseif typ == 2 && ga == 6
sgtitle({'Spike Distribution REM','Slow Gamma'})
end

end

% %%
% %%Saving
% figure(101)
%   set(gcf,'renderer','Painters')
% saveas(gcf,'./Figures/Fig6_SuppBurstPrefMG','epsc')
% 
% figure(106)
%   set(gcf,'renderer','Painters')
% saveas(gcf,'./Figures/Fig6_SuppBurstPrefSG','epsc')

%% Significance Testing
% Test for slow gamma 
a1 = squeeze(cc_lb_all(:,6,th_group==1)); % Locking
a3 = squeeze(cc_lb_all(:,6,th_group==3)); % Precessing
[h11,p11]=ttest(a1');
[h12,p12]=ttest(a3');


% Test for medium gamma 
a1 = squeeze(cc_lb_all(:,1,th_group==1)); % Locking
a3 = squeeze(cc_lb_all(:,1,th_group==3)); % Precessing
[h21,p21]=ttest(a1');
[h22,p22]=ttest(a3');

% Test for medium gamma 
a1 = squeeze(cc_ss_all(:,1,th_group==1)); % Locking
a3 = squeeze(cc_ss_all(:,1,th_group==3)); % Precessing




[h31,p31]=ttest(a1');
[h32,p32]=ttest(a3');
