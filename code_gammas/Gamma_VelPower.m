
%% Calculate average power in XD 
clear

%Your Layers
lay_n = [3 6 9;
  2 5 9;
  1 5 9;
  1 5 9;
  2 5 9;
  1 6 9];

for animal = 4
    
pyr = lay_n(animal,1)
rad = lay_n(animal,2)
slm = lay_n(animal,3)

for session = 2
    
    
load(['Animal' num2str(animal) '\'  num2str(session) 'C_Raw_CSD'])
load(['Animal' num2str(animal) '\'  num2str(session) 'Velocity'])
load(['Animal' num2str(animal) '\'  num2str(session) 'XY'])
load(['Animal' num2str(animal) '\'  num2str(session) 'pos_times'])
    
lfp = [];
lfp = Raw_CSD;

vel = Vel.data;
pos = XY.data;

[SGpyr_Power,SGpyr_Phase,SGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[20 45]);
[SGrad_Power,SGrad_Phase,SGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[20 45]);
[SGslm_Power,SGslm_Phase,SGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[20 45]);

[MGpyr_Power,MGpyr_Phase,MGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[60 90]);
[MGrad_Power,MGrad_Phase,MGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[60 90]);
[MGslm_Power,MGslm_Phase,MGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[60 90]);

[FGpyr_Power,FGpyr_Phase,FGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[100 150]);
[FGrad_Power,FGrad_Phase,FGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[100 150]);
[FGslm_Power,FGslm_Phase,FGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[100 150]);

%Reshape lfp timeframe to match velocity
timevector = [];
timevector = 0:1/1000:(length(lfp)/1000);timevector(end)=[];

    gamma = [];
for gamma=1:3
    if gamma == 1 
        Gpyr_Power = SGpyr_Power;
        Grad_Power = SGrad_Power;
        Gslm_Power = SGslm_Power;
    elseif gamma == 2   
        Gpyr_Power = MGpyr_Power;
        Grad_Power = MGrad_Power;
        Gslm_Power = MGslm_Power;
    elseif gamma == 3   
        Gpyr_Power = FGpyr_Power;
        Grad_Power = FGrad_Power;
        Gslm_Power = FGslm_Power;              
    end     
    lay=[];
for lay=1:3
  
switch lay
    
    case 1
        G_Pow=Gpyr_Power; 
    case 2
        G_Pow=Grad_Power; 
    case 3
        G_Pow=Gslm_Power;       

end 


%Average speed in windows of 330ms
clear frq cwt_binned vel_2
count = 1;

window =6;

for idx =1:window:numel(vel)-window;
    idx
vel_2(count) = mean(vel(idx:idx+window));
pos_2(count) = mean(pos(idx:idx+window));
G_power2(count,lay,gamma) = mean(G_Pow(:,timevector >= pos_times(idx) & timevector <= pos_times(idx+window)));

% cwt_binned(count,:) = max(real(clean_cwt(:,timevector >= pos_times(idx) & timevector <= pos_times(idx+window)))');
% cwt_binned_power(count,:) = zscore(abs((cwt_binned(count,:)).^2));

count = count+1;

end

end 
end 
save(['GammaVelPos2_A' num2str(animal) '_S' num2str(session)],'vel_2','pos_2','G_power2')

end

end 



%% Plot Scatters and compute correlations

figure(333),clf
subplot(1,4,1)
scatter(vel_2,g_power_2(:,1))
title('SG Pyr')
xlabel('Binned Velocity')
ylabel('Power')
correlation(1) = corr(vel_2',g_power_2(:,1),'rows','complete')

subplot(1,4,2)
scatter(vel_2,g_power_2(:,2))
xlabel('Binned Velocity')
ylabel('Power')
title('SG Rad')
correlation(2) = corr(vel_2',g_power_2(:,2),'rows','complete')

subplot(1,4,3)
scatter(vel_2,g_power_2(:,3))
xlabel('Binned Velocity')
ylabel('Power')
title('SG SLM')
correlation(3) = corr(vel_2',g_power_2(:,3),'rows','complete')


%% All contact of the probe
clear
clc

animal = 1;
session = 1;

load ([num2str(session) 'C_Raw_CSD'])
% load ([num2str(session) 'C_Raw_LFP'])
load('Velocity.mat')
load('pos_times.mat')

lfp = Raw_CSD;

vel = Vel.data;

%Reshape lfp timeframe to match velocity
timevector = 0:1/1000:(length(lfp)/1000);timevector(end)=[];


for lay_n = 1:14
lay_n

[SGPower{lay_n},SG_Phase,SG_Cycle]=ComputeWavelet_2(Raw_CSD,lay_n,[20 45]);
[MGPower{lay_n},MG_Phase,MG_Cycle]=ComputeWavelet_2(Raw_CSD,lay_n,[60 90]);

end 


%Average speed in windows of 330ms
clear frq cwt_binned vel_2

count = 1;
window = 15;

for lay_n = 1:14
    lay_n
    count = 1;
for idx =1:window:numel(vel)-window;
    
    
vel_2(count) = mean(vel(idx:idx+window));

pow1 =SGPower{lay_n};
pow2 =MGPower{lay_n};

g_power_SG(count,lay_n) = mean(pow1(:,timevector >= pos_times(idx) & timevector <= pos_times(idx+window)));
g_power_MG(count,lay_n) = mean(pow2(:,timevector >= pos_times(idx) & timevector <= pos_times(idx+window)));

correlationSG(lay_n) = corr(vel_2',g_power_SG(:,lay_n),'rows','complete');
correlationMG(lay_n) = corr(vel_2',g_power_MG(:,lay_n),'rows','complete');

count = count+1;

end
end 



% %% Re-order power for velocities
% 
% 
% vel_2 = vel_2';
% 
% % max(zscore(g_power_2)')
% % [max_bin, idx_max] = max(zscore(g_power_2)');
% % [max_bin2, idx_max2] = max((g_power_2)');
% 
% figure(11),clf
% bar(vel)
% 
% %Plot sorted by velocity
% [se es] = sort(vel_2);
% % % [se es] = sort(acc_2);
% 
% 
% figure(30),clf
% contourf((g_power_2(es,:)'),10, 'LineStyle','none')
% colormap(jet)
% 
% xlabel('Time sorted by Speed: 0 cm/s - 40cm/s')
% ylabel('Frq (Hz)')
% set(gca,'ytick',1:5:30,'yticklabels',{'20','30','43','75','105','150','200'})
% axis tight
% 
% pb_smooth = imgaussfilt((g_power_2(es,:)'),2);
% contourf(pb_smooth,50,'linecolor','none')
% 
% 
% 
% % %% Plot with superimposed velocity
% % 
% % figure(12),clf
% % contourf((cwt_binned_power'),100, 'LineStyle','none')
% % colorbar;colormap(jet)
% % hold on
% % plot(vel_2,'k','LineWidth',1.5)
% % ylim([0 35])
% % 
% % xlabel('Time')
% % ylabel('Frq (Hz)')
% % set(gca,'ytick',2:5:30,'yticklabels',{'20','30','43','75','105','150','200'})
% % axis tight
% % 
% % 
% % 
% % % figure(22),clf
% % % contourf((cwt_binned_power'),50, 'LineStyle','none')
% % % colorbar;colormap(jet)
% % % hold on
% % % plot(acc_2,'k','LineWidth',2)
% % % % ylim([0 50])
% 
% 
% % %%
% % for i=1:size(cwt_binned_power,1)
% %     [asd dsa] = max(cwt_binned_power(i,:));
% %     rodvar(i) = dsa; %based on gamma power
% %         [asd dsa] = max(cwt_binned_power(i,:));
% %     rodvar(i) = dsa; %based on gamma power
% % end
% % [ee er] = sort(rodvar);
% % figure(22),clf
% % contourf((cwt_binned_power(er,:)'),100, 'LineStyle','none')
% % %%
% % figure(22),clf
% % % plot(vel_2,'k')
% % % hold on
% % % plot(smooth(vel_2,10),'r')
% % [asd dsa] = max((cwt_binned)');
% % 
% % scatter(smooth(vel_2,10),(ee))
% % xlabel('Speed')
% xlabel('Speed')
% 
% 
