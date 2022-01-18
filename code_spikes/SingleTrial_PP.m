lay_n = [3 6 9;
    2 5 9;
    1 6 9;
    1 5 9;
    2 5 9;
    1 6 9];

lfp_dir = '~/Data/Matteo_Early/';
prop_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/'; 
save_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/'; 

n_gammas = 4;
max_gap = 6;
jitter = 0;
j_t = 40;
for spike_type = 1:3
for ph_type = 0:3

for anm = [2 3 5]
    disp(['Animal ' num2str(anm)])
    xm_all = cell(n_gammas,2);
    cm_all = cell(n_gammas,2);
for sess = 2:3


    
    
%% Gamma Power Extraction 

if(~exist([lfp_dir '/Animal' num2str(anm) '/' num2str(sess) 'Filter_CSDAllGamma_forSpCross.mat'],'file'))

load([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'C_Raw_CSD.mat'])

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


if(~exist([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP_forSpCross.mat'],'file'))
slm = lay_a(3);
load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'C_Raw_LFP.mat'])
[~,Th_Phase] = ComputeWavelet_2(Raw_LFP,slm,[6 10]);

save([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP_forSpCross.mat'],'Th_Phase')

else
  load([lfp_dir 'Animal' num2str(anm) '/' num2str(sess) 'Filter_LFP_forSpCross.mat'],'Th_Phase')  
end

phzr = Th_Phase;


phzr(phzr < 0) = phzr(phzr < 0) + 2 * pi;
phzr = phzr+ph_type/2*pi;
phzr(phzr > 2*pi) = phzr(phzr > 2*pi) - 2 * pi;

th_edges = find(diff(phzr)<0);

%%
comp = 1;
switch comp
    
    case 1
        G1_Pow=SGpyr_Power; 
        G2_Pow=MG_Power;   
    case 2
        G1_Pow=MG_Power; 
        G2_Pow=FG_Power;   
    case 3
        G1_Pow=SGpyr_Power; 
        G2_Pow=FG_Power; 
    case 4
        G1_Pow=SGrad_Power; 
        G2_Pow=MG_Power; 
    case 5
        G1_Pow=SGpyr_Power; 
        G2_Pow=SGrad_Power;   
        
        
end


G1_Pow = zscore(G1_Pow);
G2_Pow = zscore(G2_Pow);

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
Sp_Th = mean((Speed))*0.5;


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
% for tt = 1:numel(up.start)
%    GammaPow_R(ceil(up.start(tt)*1000):ceil(up.stop(tt)*1000)) = GammaPow_U(ceil(up.start(tt)*1000):ceil(up.stop(tt)*1000));
% end



%% Calculate Phase of spikes and Burst vs Single Spikes during MG vs SG
close all 

cc_bin = 20;

Pbound = prctile(GammaPow_R, linspace(0,100,n_gammas+1));

%Pbound = linspace(-1,1,5);

load([lfp_dir 'Animal' num2str(anm)  '/0' num2str(sess) 'Spikes.mat'])



for gro = 1:1







for dir =1:2
    
    if dir ==1
     load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_UP_A' num2str(anm) '_S0' num2str(sess) '_F1.mat'   ])
     spikes = spikes_up; 
     pf_limits_1    = pf_limits; 
     pf_limits_2    = pf_limits; 
     dataspikes = spikes(p_cells);
    
     Trial = Trial_Up;
      
       
    elseif dir ==2
       load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_DOWN_A' num2str(anm) '_S0' num2str(sess) '_F1.mat'   ])    
     spikes = spikes_down;  
     pf_limits_1    = pf_limits; 
     pf_limits_2    = pf_limits; 
     dataspikes = spikes(p_cells);
     
     Trial = Trial_Down;
    end 
    

    


    
    cc_all = zeros(201,n_gammas,max_gap,numel(dataspikes));
    cc_tri = {};
    
for ii=1:numel(dataspikes)

    

           SpikeVec = dataspikes{1, ii}.t;    

           

           PF_Limits_1 = pf_limits(ii,:);
           
           
           if numel(SpikeVec)<50 || PF_Limits_1(2)>=90 || PF_Limits_1(2)<=10
               continue
           end
           
           D1 = [diff(SpikeVec(1:2)) ; min(diff(SpikeVec(1:end-1)),diff(SpikeVec(2:end))); diff(SpikeVec(end-1:end))];
D1(D1<0) = 1000;
spikes_bursts_1 = SpikeVec(D1<0.009);
spikes_single_1 = SpikeVec(D1>0.009);


SpikeVec_1 = SpikeVec;



for pp = 1:numel(Pbound)-1 
%for tr = 1:max(Trial)   
   switch spike_type  
       case 1
   d1 = SpikeVec_1;
   d2 = SpikeVec_1;
       case 2
   d1 = spikes_bursts_1;
   d2 = spikes_bursts_1;
       case 3
   d1 = spikes_single_1;
   d2 = spikes_single_1;
   end
    
Spike_Pos_1 = interp1(TT,PP,d1)';
Spike_Pos_2 = interp1(TT,PP,d1)';     
           
% Spike_Pos_1
% PF_Limits_1
% pause()


takespike = find(GammaPow_R(round(d1*1000))>Pbound(pp) & GammaPow_R(round(d1*1000))<Pbound(pp+1) & Spike_Pos_1>=PF_Limits_1(1) & Spike_Pos_1<=PF_Limits_1(3));
%takespike = find(Spike_Pos_1>=PF_Limits_1(1) & Spike_Pos_1<=PF_Limits_1(3));
data_1 = d1(takespike)*1000;

takespike = find(GammaPow_R(round(d2*1000))>Pbound(pp) & GammaPow_R(round(d2*1000))<Pbound(pp+1) & Spike_Pos_2>=PF_Limits_1(1) & Spike_Pos_2<=PF_Limits_1(3));
%takespike = find(Spike_Pos_2>=PF_Limits_1(1) & Spike_Pos_2<=PF_Limits_1(3));
data_2 = d2(takespike)*1000;

Th_Sp1 = phzr(round(data_1))';
Th_Sp2 = phzr(round(data_2))';

Spike_Cycle_1 = discretize(data_1,th_edges);
Spike_Cycle_2 = discretize(data_2,th_edges);

Spike_First_1 = [1; diff(Spike_Cycle_1)>0];
Spike_First_2 = [1; diff(Spike_Cycle_2)>0];

Trial_Sp1 = Trial(round(data_1))';
Trial_Sp2 = Trial(round(data_2))';

% Trial_Sp1(Trial_Sp1~=tr)=NaN;
% Trial_Sp2(Trial_Sp2~=tr)=NaN;

DD1 = pdist2(Spike_Cycle_1,Spike_Cycle_2,@(x,y) x-y);
DD1b = pdist2(Spike_First_1,Spike_First_2,@(x,y) x*y);
DD1c = pdist2(Trial_Sp1,Trial_Sp2,@(x,y) x-y);

if(numel(DD1)==0)
   continue 
end

DD1 = DD1(:);
DD1b = DD1b(:);
DD1c = DD1c(:);

for gap = 1:max_gap
Use_Spikes = find(DD1==gap & DD1b == 1 & DD1c == 0);


DD2 = pdist2(Th_Sp1,Th_Sp2,@(x,y) x-y);
DD2 = DD2(:);
DD2 = DD2(Use_Spikes);

%subplot(numel(Pbound)-1,2,2+(kk-1)*2)
%cc = histcounts(DD,-101:1:101,'Normalization','count');
cc = histcounts(DD2,linspace(-pi,pi,202),'Normalization','count');

cc_all(:,pp,gap,ii) = cc;
%cc_tri{tr,ii}=DD2;

end
end



end 

save([save_dir 'STPPDati/SingleTrialPPDatiA' num2str(anm) 'S' num2str(sess) 'D' num2str(dir) 'Gr' num2str(gro) 'Sp' num2str(spike_type) 'Ph' num2str(ph_type)],'cc_all')



end
end


end


end
end
end