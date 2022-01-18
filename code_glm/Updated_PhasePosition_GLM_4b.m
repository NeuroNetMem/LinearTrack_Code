
% GLM for Gamma dominated Phase precession 
clear
clc

%Your Layers

lay_n = [3 6 9;
  2 5 9;
  1 5 9;
  1 5 9;
  2 5 9;
  1 6 9;
  2 6 9];

raw_dir = '~/Data/Matteo_Early'; 
process_dir = '~/Data/Matteo_GLMEarly';
for Laser = [0 1]
for animal=5:5

pyr = lay_n(animal,1)
rad = lay_n(animal,2)
slm = lay_n(animal,3)

for session=2:3
for run_dir=1:2
    
%% MODIFY FOR UP/DOWN (ALSO CHANGE FILENAME FOR SAVING AT THE BOTTOM)

Bin_Size=10;
Rescale_F=1000/Bin_Size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Raw_LFP'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Raw_CSD'])
%load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Envelope'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'XY'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'SpikesON'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'SpikesOFF'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'Running'])

Pos = XY;

    if(run_dir==1 && Laser == 0)
     spikes=spikes_upOFF;
    Mom_r=up; 

    elseif(run_dir==2 && Laser == 0)
     spikes=spikes_downOFF;
    Mom_r=down;   
    elseif(run_dir==1 && Laser == 1)
     spikes=spikes_upON;
    Mom_r=up;   
    elseif(run_dir==2 && Laser == 1)
     spikes=spikes_downON;
    Mom_r=down;   
    end    




%% GAMMAs

%save_dir = 'D:\processing\m352436\Track9\GammaGLM_G';
 save_dir = [process_dir '/Animal' num2str(animal) '/GammaGLM_SF_A' num2str(animal) '_S' num2str(session) '/'];


if ~exist(save_dir, 'dir')
 mkdir(save_dir)
end

    
[SGpyr_Power,SGpyr_Phase,SGpyr_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[20 45]);
[SGrad_Power,SGrad_Phase,SGrad_Cycle]=ComputeWavelet_2(Raw_CSD,rad,[20 45]);
[SGslm_Power,SGslm_Phase,SGslm_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[20 45]);
[MG_Power,MG_Phase,MG_Cycle]=ComputeWavelet_2(Raw_CSD,slm,[60 90]);
%[FG_Power,FG_Phase,FG_Cycle]=ComputeWavelet_2(Raw_CSD,pyr,[100 180]);
% [FGDG_Power,FGDG_Phase,FGDG_Cycle]=ComputeWavelet_2(Raw_CSD,dg,[100 180]);

TH_Power_All={};

for th_site = 1:3
switch th_site
    case 1

    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_CSD,pyr,[6 10]);
    case 2
    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_CSD,rad,[6 10]);
    case 3
    [TH_Power_All{th_site},TH_Phase_All{th_site},TH_Cycle_All{th_site}]=ComputeWavelet_2(Raw_CSD,slm,[6 10]);

end
    
    TH_Power = TH_Power_All{th_site};
    TH_Phase = TH_Phase_All{th_site};
    TH_Cycle = TH_Cycle_All{th_site};
    
%% THETA


TH_Phase_Un=TH_Phase;
TH_Cycle_Un=TH_Cycle;




%% SELECT THETA

for comp=1:3
  
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

mean(Speed)

pos_tt=(pos_times-pos_times(1))*Rescale_F+1;

pos_int=interp1(pos_tt,XY_data,1:numel(TH_Phase));
speed_int = interp1(pos_tt(1:end-1),Speed,1:numel(TH_Phase));


speed_int(speed_int < Sp_Th) = 0;
speed_int(speed_int >=Sp_Th) = 1;



PosTot=pos_int;
RR=find(isnan(pos_int));
PosTot(RR)=[];

Pos=pos_int;
Speed = speed_int; 

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
                spikes_bursts = SpikeVec(D1<0.007);
                spikes_single = SpikeVec(D1>0.007);
    
                switch spike_type
                    case 1 
                    SPIKE = SpikeVec;    
                    case 2
                    SPIKE = spikes_single;
                    case 3
                    SPIKE = spikes_bursts;    
                end
          
    TT=round((SPIKE-pos_times(1))*Rescale_F)+1;
    SpTime(TT,uu,spike_type)=1;          
                
               end 
    
end
end
%%

Take_bin1=[];
for mm=1:numel(Mom_r.start)

Take_bin1=cat(1,Take_bin1,(round((Mom_r.start(mm)-pos_times(1))*Rescale_F):round((Mom_r.stop(mm)-pos_times(1))*Rescale_F))');

end
Take_bin1(Take_bin1==0)=1;







%Add speed filter
Take_bin2 = find(Speed);

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


for spike_type = 1:3
SpTime_U=squeeze(SpTime(Take_bin,:,spike_type));


N_Gamma = 6;
%
[ext,alpha,Like,gauPlotVal,gauCentX,gauCentY,gauCentZ] = GLMPhaseThetaPowGamma(SpTime_U,Pos_U,PosTot,Th_Phase_U,GammaPow_U,40,4,2/N_Gamma*2,40,N_Gamma);
cost=repmat(ext,1,40*20*N_Gamma);
PlotVal=permute(reshape(exp(alpha'*gauPlotVal+cost)./(2*cosh(alpha'*gauPlotVal+cost)),size(SpTime_U,2),20,40,N_Gamma),[2 3 1 4]);

PlotVal_Ch=PlotVal;


 
    save([save_dir 'GLMPyrGammaPow_A' num2str(animal) '_S' num2str(session) '_G' num2str(comp) '_T' num2str(th_site) '_D' num2str(run_dir) '_Sp' num2str(spike_type) '_L' num2str(Laser)],'PlotVal_Ch')
    
end

end

end

end

end% run_dir cycle
end % sessions cycle
end % animal cycle






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM for Gamma Envelopes Phase precession!
clear
clc

raw_dir = '~/Data/Matteo_Laser_2'; 
process_dir = '~/Data/Matteo_GLMLaser_2';

%Your Layers

lay_n = [1 3 6;
  3 5 8;
  1 3 6;
  4 6 8;
  6 9 12;
  1 6 9];
for Laser = [0 1]
for animal=5:5

pyr = lay_n(animal,1)
rad = lay_n(animal,2)
slm = lay_n(animal,3)

for session=2:3
for run_dir=1:2
    
    
%% MODIFY FOR UP/DOWN (ALSO CHANGE FILENAME FOR SAVING AT THE BOTTOM)


Bin_Size=10;
Rescale_F=1000/Bin_Size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Raw_LFP'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Raw_CSD'])
%load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'C_Envelope'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'XY'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'SpikesON'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'SpikesOFF'])
load([raw_dir '/Animal' num2str(animal) '/'  num2str(session) 'Running'])


Pos = XY;

   if(run_dir==1 && Laser == 0)
     spikes=spikes_upOFF;
    Mom_r=up; 

    elseif(run_dir==2 && Laser == 0)
     spikes=spikes_downOFF;
    Mom_r=down;   
    elseif(run_dir==1 && Laser == 1)
     spikes=spikes_upON;
    Mom_r=up;   
    elseif(run_dir==2 && Laser == 1)
     spikes=spikes_downON;
    Mom_r=down;   
    end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% GAMMAs

Bin_Size=10;
Rescale_F=1000/Bin_Size;
 save_dir = [process_dir '/Animal' num2str(animal) '/GammaEnvelopeGLM_SF_A' num2str(animal) '_S' num2str(session) '/'];


if ~exist(save_dir, 'dir')
 mkdir(save_dir)
end





%Take the envelope 
for comp = 1:1
    
    switch comp 
        case 1
        Envelope = Raw_CSD;
        
        case 2 
        Envelope = Filter_Envelope_SG;

        case 3
        Envelope = Filter_Envelope_MG;
        
        case 4 
        Envelope = Filter_Envelope_FG;
            
    end 
        
    
for th_site=1:3

    
%% THETA
switch th_site
    case 1
    [TH_Power,TH_Phase,TH_Cycle]=ComputeWavelet_2(Envelope,pyr,[6 10]);
    case 2
    [TH_Power,TH_Phase,TH_Cycle]=ComputeWavelet_2(Envelope,rad,[6 10]);
    case 3
    [TH_Power,TH_Phase,TH_Cycle]=ComputeWavelet_2(Envelope,slm,[6 10]);

end

TH_Phase_Un=TH_Phase;
TH_Cycle_Un=TH_Cycle;

%

    
PlotVal_Ch=[];    
    
TH_Phase=TH_Phase_Un;
TH_Cycle=TH_Cycle_Un;
    

TH_Power=TH_Power(1:Bin_Size:end);
TH_Power = zscore(TH_Power);

TH_Phase=TH_Phase(1:Bin_Size:end);
TH_Cycle=floor(TH_Cycle/Bin_Size)+1;

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
Speed = speed_int; 


TH_Power=TH_Power(round(pos_times(1)*Rescale_F):end);

TH_Phase=TH_Phase(round(pos_times(1)*Rescale_F):end);
TH_Cycle=TH_Cycle-round(pos_times(1)*Rescale_F)+1;



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
    SpTime(TT,uu,spike_type)=1;          
                
               end 
    
end
end
%%

Take_bin1=[];
for mm=1:numel(Mom_r.start)

Take_bin1=cat(1,Take_bin1,(round((Mom_r.start(mm)-pos_times(1))*Rescale_F):round((Mom_r.stop(mm)-pos_times(1))*Rescale_F))');

end
Take_bin1(Take_bin1==0)=1;


%Add speed filter
Take_bin2 = find(Speed);

Take_bin1 = intersect(Take_bin1,Take_bin2);
Take_bin1=intersect(Take_bin1,find(~isnan(Pos)));






Take_bin2=[];

% for mm=1:numel(TH_Cycle)-1
% 
% Take_bin2=cat(1,Take_bin2,(round(TH_Cycle(mm)):round(TH_Cycle(mm+1))-1)');
% 
% end


%Here we take the Envelope periods with amplitude > zscore(mean)
if comp >1 
  Pow_Thr = prctile(TH_Power(Take_bin1),70);
  Take_bin2 = find(TH_Power>Pow_Thr); 
  
else
 Take_bin2=find(TH_Power<Inf);   
end 

Take_bin=intersect(Take_bin1,Take_bin2);

Th_Phase_U=TH_Phase(Take_bin);
Pos_U=Pos(Take_bin);


for spike_type = 1:3

SpTime_U=squeeze(SpTime(Take_bin,:,spike_type));



%
[ext,alpha,Like,gauPlotVal,gauCentX,gauCentY] = GLMPhaseTheta(SpTime_U,Pos_U,PosTot,Th_Phase_U,40,4,40);
cost=repmat(ext,1,40*20);
PlotVal=shiftdim(reshape(exp(alpha'*gauPlotVal+cost)./(2*cosh(alpha'*gauPlotVal+cost)),size(SpTime_U,2),20,40),1);

PlotVal_Ch=cat(4,PlotVal_Ch,PlotVal);


    save([save_dir 'GLMPyrGamma_A' num2str(animal) '_S' num2str(session) '_G' num2str(comp) '_T' num2str(th_site) '_D' num2str(run_dir)  '_Sp' num2str(spike_type) '_L' num2str(Laser)],'PlotVal_Ch')
    

 
end

end
end


end

end % run_dir cycle
end % sessions cycle
end % animal cycle










