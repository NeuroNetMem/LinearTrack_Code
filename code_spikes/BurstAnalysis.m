lay_n = [3 6 9;
    2 5 9;
    1 6 9;
    1 5 9;
    2 5 9;
    1 6 9];

direction_name = {'UP','DOWN'};

random_control = 1;
n_rand = 5;

lfp_dir = '~/Data/Matteo_Early/';
prop_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/'; 
save_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/'; 

for anm = [2 3 5]
    disp(['Animal ' num2str(anm)])
for sess = 2:3

%% Gamma Power Extraction 
if(~exist([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_CSDAllGamma_forSpCross.mat'],'file'))

load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'C_Raw_CSD.mat'])

lay_a = lay_n(anm,:);

pyr = lay_a(1);
rad = lay_a(2);
slm = lay_a(3);

[SGpyr_Power,SGpyr_Phase,SGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[20 45]);
[SGrad_Power,SGrad_Phase,SGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[20 45]);
[SGslm_Power,SGslm_Phase,SGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[20 45]);
[MG_Power,MG_Phase,MG_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[60 90]);
%[FG_Power,FG_Phase,FG_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[100 180]);



save([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_CSDAllGamma_forSpCross.mat'],'SGpyr_Power','SGrad_Power','SGslm_Power','MG_Power')
else
    load([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_CSDAllGamma_forSpCross.mat'],'SGpyr_Power','SGrad_Power','SGslm_Power','MG_Power')
end


if(~exist([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP.mat'],'file'))
slm = lay_a(3);
load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'C_Raw_LFP.mat'])
[~,Th_Phase] = ComputeWavelet_2(Raw_LFP,slm,[6 10]);

save([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP.mat'],'Th_Phase')

else
  load([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP.mat'],'Th_Phase')  
end

phzr = Th_Phase;
phzr(phzr < 0) = phzr(phzr < 0) + 2 * pi;





%%
for comp = [1 2 3]
switch comp
    
    case 1
        G1_Pow=SGpyr_Power; 
        G2_Pow=MG_Power;   
    case 2
        G1_Pow=SGrad_Power; 
        G2_Pow=MG_Power; 
    case 3
        G1_Pow=SGslm_Power; 
        G2_Pow=MG_Power; 
    case 4
        G1_Pow=SGrad_Power; 
        G2_Pow=MG_Power; 
    case 5
        G1_Pow=SGpyr_Power; 
        G2_Pow=SGrad_Power;   
        
        
end


G1_Pow = zscore(G1_Pow);
G2_Pow = zscore(G2_Pow);

G1_Pow(G1_Pow<-4) = -4;
G1_Pow(G1_Pow>4) = 4;
G2_Pow(G2_Pow<-4) = -4;
G2_Pow(G2_Pow>4) = 4;


G1_Pow = G1_Pow - min(G1_Pow);
G1_Pow = G1_Pow/max(G1_Pow);
G2_Pow = G2_Pow - min(G2_Pow);
G2_Pow = G2_Pow/max(G2_Pow);

GammaPow_U=(G1_Pow-G2_Pow)./(G1_Pow+G2_Pow); 

GammaPow_U = smoothdata(GammaPow_U, 'gaussian',20);


%%
load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'Running'])

load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'XY'])


TT = XY.t;
PP = XY.data;
lost = find(isnan(PP));

PP(lost)=interp1(TT(~isnan(PP)),PP(~isnan(PP)),TT(lost));


PP(isnan(PP))=0;
Speed = smoothdata(abs(diff(PP))*30,'movmedian',20); 
Sp_Th = mean((Speed))*1;


Speed(Speed < Sp_Th) = 0;
Speed(Speed >=Sp_Th) = 1;


CC = bwconncomp(Speed);


GammaPow_R = ones(size(GammaPow_U))*NaN;
for tt = 1:CC.NumObjects
t1 = TT(CC.PixelIdxList{tt}(1))*1000-30;
t2 = TT(CC.PixelIdxList{tt}(end))*1000+30;

GammaPow_R(round(t1):round(t2)) = GammaPow_U(round(t1):round(t2));



end

Trial_Up = ones(size(phzr))*NaN;
for tt = 5:numel(up.start)
   Trial_Up(ceil(up.start(tt)*1000):ceil(up.stop(tt)*1000)) = tt;
end

Trial_Down = ones(size(phzr))*NaN;
for tt = 5:numel(down.start)
   Trial_Down(ceil(down.start(tt)*1000):ceil(down.stop(tt)*1000)) = tt;
end


% GammaPow_R = ones(size(GammaPow_U))*NaN;
% for tt = 1:numel(down.start)
%    GammaPow_R(ceil(down.start(tt)*1000):ceil(down.stop(tt)*1000)) = GammaPow_U(ceil(down.start(tt)*1000):ceil(down.stop(tt)*1000));
% end
% 
% for tt = 1:numel(up.start)
%    GammaPow_R(ceil(up.start(tt)*1000):ceil(up.stop(tt)*1000)) = GammaPow_U(ceil(up.start(tt)*1000):ceil(up.stop(tt)*1000));
% end


%% Calculate Phase of spikes and Burst vs Single Spikes during MG vs SG



 load([prop_dir 'Theta_Limits/Theta_Limits.mat'])
strat_h = th_limits_lin;
n_strat = numel(strat_h)-1;

n_groups = 6;
cc_bin = 20;

Pbound = prctile(GammaPow_R, linspace(0,100,n_groups+1));

%Pbound = linspace(-1,1,5);

if random_control == 0
cc_sb = NaN(cc_bin,numel(Pbound)-1,200);
cc_lb = NaN(cc_bin,numel(Pbound)-1,200);
cc_ss = NaN(cc_bin,numel(Pbound)-1,200);
cc_ab = NaN(cc_bin,numel(Pbound)-1,200);

burst_per = NaN(2,numel(Pbound)-1,200);

elseif random_control == 1
    
cc_sb = NaN(cc_bin,numel(Pbound)-1,200,n_rand);
cc_lb = NaN(cc_bin,numel(Pbound)-1,200,n_rand);
cc_ss = NaN(cc_bin,numel(Pbound)-1,200,n_rand);
cc_ab = NaN(cc_bin,numel(Pbound)-1,200,n_rand);

burst_per = NaN(2,numel(Pbound)-1,200,n_rand);  
end

load([lfp_dir 'Animal' num2str(anm)  '/0' num2str(sess) 'Spikes.mat'])



CellPref_Angle = [];
theta_session = [];
ll = 0; 

ISI_Pop= [];

for dir=1:2
    
    if dir ==1
     %load(['~/Data/Matteo/Animal' num2str(anm) '/AnimalNG' num2str(anm)  '/Cluster_GroupingsUP_A' num2str(anm) '_S' num2str(sess) '.mat'])    
     load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_UP_A' num2str(anm) '_S0' num2str(sess) '_F1.mat'   ])
     
     spikes = spikes_up;
     Trial=Trial_Up;  
       
    elseif dir ==2
     %load(['~/Data/Matteo/Animal' num2str(anm) '/AnimalNG' num2str(anm)  '/Cluster_GroupingsDOWN_A' num2str(anm) '_S' num2str(sess) '.mat'])          
     load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_DOWN_A' num2str(anm) '_S0' num2str(sess) '_F1.mat'   ])
     spikes = spikes_down;
     
     Trial=Trial_Down;
    end 
 
%       if gro == 1
%      dataspikes = spikes(p_cells(fr_phi_corr<-0.2));
%      elseif gro == 2
%      dataspikes = spikes(p_cells(abs(fr_phi_corr)<0.2));    
%      end   
     
%     th_strat = discretize(theta_scores,strat_h);
%     dataspikes = spikes(p_cells(th_strat==gro));
    theta_session = cat(1,theta_session,theta_scores(:));
    dataspikes = spikes(p_cells);

load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_' direction_name{dir} '_A' num2str(anm) '_S0' num2str(sess) '_F1.mat'   ],'pf_limits')
    
%pf_limits = pf_limits(th_strat==gro,:);


% load(['~/Dropbox/Projects_NIJ/Matteo/Nijmegen/CCDati/n_ord_2_A' num2str(anm) 'S' num2str(sess) 'D' num2str(dir) 'G' num2str(gro) '.mat']) 


for ii=1:size(dataspikes,2)
        
    ll = ll+1;
    
        PF_Limits_1 = pf_limits(ii,:);
    
         spikes = dataspikes{1, ii}.t;     

         
         Spike_Pos_1 = interp1(TT,PP,spikes)';
         takespike = find(Spike_Pos_1>=PF_Limits_1(1) & Spike_Pos_1<=PF_Limits_1(3));
         spikes_pf = spikes(takespike);
         
         CellPref_Angle=cat(1,CellPref_Angle,angle(sum(exp(1j.*phzr(round(spikes_pf*1000))))));
           
         
         
         if numel(spikes_pf)<50 || PF_Limits_1(2)>=90 || PF_Limits_1(2)<=10
               continue
           end
         
  
            SpikeVec = spikes_pf;




%Burst vs single calculation

D1 = [diff(SpikeVec(1:2)) ; min(diff(SpikeVec(1:end-1)),diff(SpikeVec(2:end))); diff(SpikeVec(end-1:end))];
D1(D1<0.001) = 1000;

ISI_Pop = [ISI_Pop;D1(:)];

D_Burst = zeros(size(D1));
D_Burst(D1<0.007)=1;

CC = bwconncomp(D_Burst);



for cc = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cc})>3
    D_Burst(CC.PixelIdxList{cc}) = 2;
    
    end
end



kk = 0;
for pp = 1:numel(Pbound)-1 
    
kk = kk+1;
takespike = find(GammaPow_R(round(SpikeVec*1000))>Pbound(pp) & GammaPow_R(round(SpikeVec*1000))<Pbound(pp+1) & D_Burst'==1);
sb = SpikeVec(takespike);

takespike = find(GammaPow_R(round(SpikeVec*1000))>Pbound(pp) & GammaPow_R(round(SpikeVec*1000))<Pbound(pp+1) & D_Burst'==2);
lb = SpikeVec(takespike);

takespike = find(GammaPow_R(round(SpikeVec*1000))>Pbound(pp) & GammaPow_R(round(SpikeVec*1000))<Pbound(pp+1) & D_Burst'==0);
ss = SpikeVec(takespike);

takespike = find(GammaPow_R(round(SpikeVec*1000))>Pbound(pp) & GammaPow_R(round(SpikeVec*1000))<Pbound(pp+1) & D_Burst' >0);
ab = SpikeVec(takespike);




takegamma = find(GammaPow_R>Pbound(pp) & GammaPow_R<Pbound(pp+1));


if random_control == 1

    
for rr = 1:n_rand 
    
    
rp=randperm(numel(takegamma));
sb_r = takegamma(rp(1:numel(sb)));

rp=randperm(numel(takegamma));
lb_r = takegamma(rp(1:numel(lb)));

rp=randperm(numel(takegamma));
ss_r = takegamma(rp(1:numel(ss)));

rp=randperm(numel(takegamma));
ab_r = takegamma(rp(1:numel(ab)));


    sbph = phzr(round(sb_r));
    lbph = phzr(round(lb_r));
    ssph = phzr(round(ss_r)); 
    abph = phzr(round(ab_r)); 

    
    gamma_th = phzr(takegamma); % What is the distribution of theta phases given my choice of gamma power?
     
    
contr = smoothdata(histcounts(gamma_th,linspace(0,2*pi,cc_bin+1)),'gaussian',3);
%contr = ones(cc_bin,1)';    

cc = histcounts(sbph,linspace(0,2*pi,cc_bin+1))./contr; % Spike count normalized to phase probability
cc_sb(:,pp,ll,rr)=cc;

cc = histcounts(lbph,linspace(0,2*pi,cc_bin+1))./contr;
cc_lb(:,pp,ll,rr)=cc;

cc = histcounts(ssph,linspace(0,2*pi,cc_bin+1))./contr;
cc_ss(:,pp,ll,rr)=cc;

cc = histcounts(abph,linspace(0,2*pi,cc_bin+1))./contr;
cc_ab(:,pp,ll,rr)=cc;

burst_per(:,pp,ll,rr) = [numel(sb)/numel(ss) numel(ab)/numel(ss)];    
    
end    
    
    
else

% 
%     phLaser = phzr(round(SpikeVecLaser*1000));
%     phNOLaser =  phzr(round(SpikeVecNOLaser*1000));
    
    sbph = phzr(round(sb*1000));
    lbph = phzr(round(lb*1000));
    ssph = phzr(round(ss*1000)); 
    abph = phzr(round(ab*1000)); 

    
    gamma_th = phzr(takegamma); % What is the distribution of theta phases given my choice of gamma power?
     
    
contr = smoothdata(histcounts(gamma_th,linspace(0,2*pi,cc_bin+1)),'gaussian',3);
%contr = ones(cc_bin,1)';    

cc = histcounts(sbph,linspace(0,2*pi,cc_bin+1))./contr; % Spike count normalized to phase probability
cc_sb(:,pp,ll)=cc;

cc = histcounts(lbph,linspace(0,2*pi,cc_bin+1))./contr;
cc_lb(:,pp,ll)=cc;

cc = histcounts(ssph,linspace(0,2*pi,cc_bin+1))./contr;
cc_ss(:,pp,ll)=cc;

cc = histcounts(abph,linspace(0,2*pi,cc_bin+1))./contr;
cc_ab(:,pp,ll)=cc;

burst_per(:,pp,ll) = [numel(sb)/numel(ss) numel(ab)/numel(ss)];

end




end 


end     
      

end

if random_control == 0
cc_ss(:,:,ll+1:end)=[];
cc_lb(:,:,ll+1:end)=[];
cc_sb(:,:,ll+1:end)=[];
cc_ab(:,:,ll+1:end)=[];

burst_per(:,:,ll+1:end)=[];

elseif random_control == 1
cc_ss(:,:,ll+1:end,:)=[];
cc_lb(:,:,ll+1:end,:)=[];
cc_sb(:,:,ll+1:end,:)=[];
cc_ab(:,:,ll+1:end,:)=[];

burst_per(:,:,ll+1:end,:)=[];

end


if random_control == 0 
save([save_dir 'BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga' num2str(comp)],'cc_ss','cc_sb','cc_lb','cc_ab','burst_per','CellPref_Angle','theta_session','ISI_Pop')
elseif random_control == 1
save([save_dir 'BurstDati/BurstAnalysis_A' num2str(anm) '_S' num2str(sess) '_Ga' num2str(comp) 'RR'],'cc_ss','cc_sb','cc_lb','cc_ab','burst_per','CellPref_Angle','theta_session','ISI_Pop')
end
end
end
end



