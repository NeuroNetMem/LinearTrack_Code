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

for anm = [2 3 5]
    disp(['Animal ' num2str(anm)])
for sess = 2:3

    m_Pow = zeros(64,4);
    s_Pow = zeros(64,4);
    
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
for comp = [1 2 3 4]
switch comp
    
    case 1
        G1_Pow=SGpyr_Power; 
        
    case 2
        G1_Pow=SGrad_Power; 
        
    case 3
        G1_Pow=SGslm_Power; 
        
    case 4
        
        G1_Pow=MG_Power; 
  
        
        
end


% G1_Pow = zscore(G1_Pow);
% G2_Pow = zscore(G2_Pow);
% 
% G1_Pow(G1_Pow<-4) = -4;
% G1_Pow(G1_Pow>4) = 4;
% G2_Pow(G2_Pow<-4) = -4;
% G2_Pow(G2_Pow>4) = 4;
% 
% 
% G1_Pow = G1_Pow - min(G1_Pow);
% G1_Pow = G1_Pow/max(G1_Pow);
% G2_Pow = G2_Pow - min(G2_Pow);
% G2_Pow = G2_Pow/max(G2_Pow);
% 
% GammaPow_U=(G1_Pow-G2_Pow)./(G1_Pow+G2_Pow); 
% 
% GammaPow_U = smoothdata(GammaPow_U, 'gaussian',20);


%%
load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'Running'])

load([lfp_dir 'Animal' num2str(anm) '/0' num2str(sess) 'XY'])

TT = XY.t;
PP = XY.data;
lost = find(isnan(PP));

PP(lost)=interp1(TT(~isnan(PP)),PP(~isnan(PP)),TT(lost));


PP(isnan(PP))=0;
Speed = smoothdata(abs(diff(PP))*30,'movmedian',20); 
Sp_Th = mean((Speed))*5;


Speed(Speed < Sp_Th) = 0;
Speed(Speed >=Sp_Th) = 1;


CC = bwconncomp(Speed);


GammaPow_R = ones(size(G1_Pow))*NaN;
for tt = 1:CC.NumObjects
t1 = TT(CC.PixelIdxList{tt}(1))*1000-30;
t2 = TT(CC.PixelIdxList{tt}(end))*1000+30;

GammaPow_R(round(t1):round(t2)) = G1_Pow(round(t1):round(t2));



end


Gamma_z = ones(size(G1_Pow))*NaN;
Gamma_z(~isnan(GammaPow_R)) = zscore(GammaPow_R(~isnan(GammaPow_R))); 


phzr_b=discretize(phzr,linspace(0,2*pi,65));
for bb = 1:64
    m_Pow(bb,comp)=nanmean(Gamma_z(phzr_b==bb));
    s_Pow(bb,comp)=nanstd(Gamma_z(phzr_b==bb));
end






end

save([prop_dir 'GammaPhase/GammaPhase_A' num2str(anm) '_S' num2str(sess) ],'m_Pow','s_Pow')

end
end



