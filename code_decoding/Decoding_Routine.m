clear all

%Your Layers

lay_n = [3 6 9;
  2 5 9;
  1 5 9;
  1 5 9;
  2 5 9;
  1 6 9];

data_dir = '~/Data/Matteo_Early/';
glm_dir = '~/Data/Matteo_GLMEarly/';
prop_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';


for run_groups = [0 1]
for divide_by_session = 0:run_groups
divide_by_session 
 load([prop_dir 'Theta_Limits/Theta_Limits.mat'])
if(run_groups==1)
 thsc_lims = th_limits_prc;
elseif(run_groups==0)
thsc_lims = [th_limits_prc(1) th_limits_prc(end)];
end
n_strat = numel(thsc_lims)-1;

for animal=[2 3 5]

pyr = lay_n(animal,1)
rad = lay_n(animal,2)
slm = lay_n(animal,3)

for session=2:3
    
if(divide_by_session ==1)
 load([prop_dir 'Theta_Limits/Theta_Limits_A' num2str(animal) 'S' num2str(session) '.mat'])
 thsc_lims = th_limits_prc;
end
    
    
    
for run_dir=1:2
    
%% MODIFY FOR UP/DOWN (ALSO CHANGE FILENAME FOR SAVING AT THE BOTTOM)
Fine_Size = 10;
Spike_Size = 3;


Bin_Size=Fine_Size;
TiWi = (Bin_Size*Spike_Size)/10;
Rescale_F=1000/Bin_Size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'C_Raw_CSD'])
load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'C_Raw_LFP'])
%load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'C_Envelope'])
load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'XY'])
load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'Spikes'])
load([data_dir 'Animal' num2str(animal) '/0'  num2str(session) 'Running'])
Pos = XY;

    if(run_dir==1)
     spikes=spikes_up;
    Mom_r=up; 
    
    load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_UP_A' num2str(animal) '_S0' num2str(session) '_F1.mat'   ],'p_cells','theta_scores')  
    
    Ta_Pl=p_cells;

    
    load([glm_dir 'Animal' num2str(animal) '/GammaEnvelopeGLM_SF_A' num2str(animal) '_S' num2str(session) '/GLMPyrGamma_A' ...
    num2str(animal) '_S' num2str(session) '_G' num2str(1) '_T' num2str(3) '_D' num2str(run_dir) '_Sp1.mat'],'PlotVal_Ch')

    MM_Ref = PlotVal_Ch;
    MM_Ref = squeeze(sum(MM_Ref,1));



    elseif(run_dir==2)
     spikes=spikes_down;
    Mom_r=down;    
    load([prop_dir 'Cluster_Info/' '/Cluster_GroupingsMPF_DOWN_A' num2str(animal) '_S0' num2str(session) '_F1.mat'   ],'p_cells','theta_scores')  
    
    Ta_Pl=p_cells;

    
    load([glm_dir 'Animal' num2str(animal) '/GammaEnvelopeGLM_SF_A' num2str(animal) '_S' num2str(session) '/GLMPyrGamma_A' ...
    num2str(animal) '_S' num2str(session) '_G' num2str(1) '_T' num2str(3) '_D' num2str(run_dir) '_Sp1.mat'],'PlotVal_Ch')

    MM_Ref = PlotVal_Ch;
    MM_Ref = squeeze(sum(MM_Ref,1));
    end    




%% GAMMAs

%save_dir = 'D:\processing\m352436\Track9\GammaGLM_G';
 
save_dir = [prop_dir 'Decoding/'];

if ~exist(save_dir, 'dir')
 mkdir(save_dir)
end


if(~exist([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_CSDGam_Deco.mat'],'file'))





    
[SGpyr_Power,SGpyr_Phase,SGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[20 45]);
[SGrad_Power,SGrad_Phase,SGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[20 45]);
[SGslm_Power,SGslm_Phase,SGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[20 45]);
[MG_Power,MG_Phase,MG_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[60 90]);
%[FG_Power,FG_Phase,FG_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[100 180]);
% [FGDG_Power,FGDG_Phase,FGDG_Cycle]=ComputeWavelet_2(Raw_CSD,dg,[100 180]);

save([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_CSDGam_Deco.mat'],'SGpyr_Power','SGrad_Power','SGslm_Power','MG_Power')

else
    
  load([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_CSDGam_Deco.mat'],'SGpyr_Power','SGrad_Power','SGslm_Power','MG_Power')  
    
end
TH_Power_All={};

for th_site = 3:3
    
    
    
    if(~exist([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_LFPThe_Deco.mat'],'file'))
    
switch th_site
    case 1

    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_LFP,pyr,[6 10]);
    case 2
    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_LFP,rad,[6 10]);
    case 3
    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_LFP,slm,[6 10]);

end
save([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_LFPThe_Deco.mat'],'TH_Power_All','TH_Phase_All','TH_Cycle_All')
    else
        load([data_dir 'Animal' num2str(animal) '/' num2str(session) 'Filter_LFPThe_Deco.mat'],'TH_Power_All','TH_Phase_All','TH_Cycle_All')
    end
    
    TH_Power = TH_Power_All{th_site};
    TH_Phase = TH_Phase_All{th_site};
    TH_Cycle = TH_Cycle_All{th_site};
    
%% THETA


TH_Phase_Un=TH_Phase;
TH_Cycle_Un=TH_Cycle;

%% SELECT GAMMA COMPARISON

for comp=[1]
  
PlotVal_Ch=[];    
    
TH_Phase=TH_Phase_Un;
%TH_Cycle=TH_Cycle_Un;
  
TH_Cycle=localMaximum(-abs(TH_Phase-0),round(100*0.75))+1;


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

        
end


%G_Ratio = (G1_Pow-G2_Pow)./(G1_Pow+G2_Pow);

% G_Ratio=G1_Pow./G2_Pow;
% G_Ratio=(G_Ratio-min(G_Ratio))./(max(G_Ratio)-min(G_Ratio));


TH_Phase=TH_Phase(1:Bin_Size:end);
TH_Cycle=floor(TH_Cycle/Bin_Size)+1; 
G1_Pow = G1_Pow(1:Bin_Size:end);
G2_Pow = G2_Pow(1:Bin_Size:end);




%TH_Cycle=localMaximum(TH_Phase,3)+1;

%%


pos_times=XY.t;
XY_data=XY.Data;
PP = XY.Data;
TT = XY.t; 

lost = find(isnan(PP));

PP(lost)=interp1(TT(~isnan(PP)),PP(~isnan(PP)),TT(lost));


PP(isnan(PP))=0;
Speed = smoothdata(abs(diff(PP))*30,'movmedian',20); 
Sp_Th = mean((Speed))*0.5;



pos_tt=(pos_times-pos_times(1))*Rescale_F+1;

pos_int=interp1(pos_tt,XY_data,1:numel(TH_Phase));
speed_int = interp1(pos_tt(1:end-1),Speed,1:numel(TH_Phase));


speed_int(speed_int < Sp_Th) = 0;
speed_int(speed_int >=Sp_Th) = 1;

PosTot=pos_int;
RR=find(isnan(pos_int));
PosTot(RR)=[];

Pos=pos_int;

TH_Phase=TH_Phase(round(pos_times(1)*Rescale_F):end);
TH_Cycle=TH_Cycle-round(pos_times(1)*Rescale_F)+1;

G1_Pow = G1_Pow(round(pos_times(1)*Rescale_F):end);
G2_Pow = G2_Pow(round(pos_times(1)*Rescale_F):end);



%%
SpTime=zeros(numel(TH_Phase),numel(spikes),3);
for spike_type = 1:3
   

for uu=1:numel(spikes)
    
               SpikeVec = spikes{uu}.t;  
               
               if(numel(SpikeVec)>5)
                   
                D1 = [diff(SpikeVec(1:2)) ; min(diff(SpikeVec(1:end-1)),diff(SpikeVec(2:end))); diff(SpikeVec(end-1:end))];
                D1(D1<0) = 1000;
                spikes_bursts = SpikeVec(D1<0.009);
                spikes_single = SpikeVec(D1>0.009);
    
                switch spike_type
                    case 1 
                    SPIKE = SpikeVec;    
                    case 2
                    SPIKE = spikes_single;
                    case 3
                    SPIKE = spikes_bursts;    
                end
          
    TT=round((SPIKE-pos_times(1))*Rescale_F)+1;
    for ti = TT(:)'
    SpTime(ti,uu,spike_type)=SpTime(ti,uu,spike_type)+1;          
    end
               end 
    
end

for ti = 1:size(SpTime,1)-Spike_Size
    SpTime(ti,:,spike_type) = sum(SpTime(ti:ti+Spike_Size-1,:,spike_type),1);

end



end
%%

Take_bin1=[];
for mm=1:numel(Mom_r.start)

Take_bin1=cat(1,Take_bin1,(round((Mom_r.start(mm)-pos_times(1))*Rescale_F):round((Mom_r.stop(mm)-pos_times(1))*Rescale_F))');

end

Take_bin1(Take_bin1==0)=1;

%Add speed filter
Take_bin2 = find(speed_int);

Take_bin1 = intersect(Take_bin1,Take_bin2);
Take_bin1=intersect(Take_bin1,find(~isnan(Pos)));
Take_bin=Take_bin1;

Th_Phase_U=TH_Phase(Take_bin);

Pos_U=Pos(Take_bin);
G1_Pow = zscore(G1_Pow(Take_bin));
G2_Pow = zscore(G2_Pow(Take_bin));

G1_Pow(G1_Pow<-4) = -4;
G1_Pow(G1_Pow>4) = 4;
G2_Pow(G2_Pow<-4) = -4;
G2_Pow(G2_Pow>4) = 4;


G1_Pow = G1_Pow - min(G1_Pow);
G1_Pow = G1_Pow/max(G1_Pow);
G2_Pow = G2_Pow - min(G2_Pow);
G2_Pow = G2_Pow/max(G2_Pow);

GammaPow_U=(G1_Pow-G2_Pow)./(G1_Pow+G2_Pow); 

for GG = 1:n_strat
    take_scores = find(theta_scores>thsc_lims(GG) & theta_scores<thsc_lims(GG+1));
    
    
    
for spike_type = 1:3
SpTime_U=squeeze(SpTime(Take_bin,Ta_Pl(take_scores),spike_type));
MM_Ref_U = MM_Ref(:,Ta_Pl(take_scores));


[BaPos_U,BaPosQ_U,CoPos_U,CoPosQ_U,Ba_De_All,Co_De_All]=Decode_Position_1(SpTime_U,MM_Ref_U,TiWi);





%  




    if(run_groups==1 && divide_by_session == 0)
    save([save_dir 'DecodeRunning_A' num2str(animal) '_S' num2str(session) '_G' num2str(comp) '_T' num2str(th_site) '_D' num2str(run_dir) '_Sp' num2str(spike_type) 'Gr' num2str(GG)],...
        'Pos_U','GammaPow_U','Th_Phase_U','TH_Phase','Take_bin','BaPos_U','BaPosQ_U','CoPos_U','CoPosQ_U','Ba_De_All','Co_De_All')   
    elseif(run_groups==0)

        save([save_dir 'DecodeRunning_A' num2str(animal) '_S' num2str(session) '_G' num2str(comp) '_T' num2str(th_site) '_D' num2str(run_dir) '_Sp' num2str(spike_type)],...
        'Pos_U','GammaPow_U','Th_Phase_U','TH_Phase','Take_bin','BaPos_U','BaPosQ_U','CoPos_U','CoPosQ_U','Ba_De_All','Co_De_All')
    
    elseif(run_groups==1 && divide_by_session == 1)

        save([save_dir 'DecodeRunning_A' num2str(animal) '_S' num2str(session) '_G' num2str(comp) '_T' num2str(th_site) '_D' num2str(run_dir) '_Sp' num2str(spike_type) 'GrInt' num2str(GG)],...
        'Pos_U','GammaPow_U','Th_Phase_U','TH_Phase','Take_bin','BaPos_U','BaPosQ_U','CoPos_U','CoPosQ_U','Ba_De_All','Co_De_All')
    


    end
end
end
end
end

end% run_dir cycle
end % sessions cycle
end % animal cycle

end
end