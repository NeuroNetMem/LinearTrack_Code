%%

close all 

direction_name = {'UP','DOWN'};

aa=0;

GLM_dir = '/home/fstella/Data/Matteo_GLMEarly/';
destination_save = '/home/fstella/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';

gro = 0;

for anm = [2 3 5] 


    
    aa=aa+1;


gg_c = 0;

gg_c = gg_c+1;
   

Place_AA=[];
    



gg_t = 0;

gg_t = gg_t+1;
    M1_All=0;
    M2_All=0;
    M3_All=0;  
    
M11_All=[];
M12_All=[];
M13_All=[];

    jj=0;
    MM_Comp = [];
    MM_Fig_1 = [];
    MM_Fig_2 = [];
    MM_Fig_3 = [];


    
    
    
    M1 = zeros(20,13,10);
    PP_Angle_All = zeros(13,10);
    PP_Angle_Sin = [];
    PP_Place_Sin = [];
    Pref_Phase = [];
    
    
    for direction = [1 2]
    for sess = [2 3]    
    
 
    
    

load([GLM_dir 'Animal' num2str(anm) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(direction) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref = PlotVal_Ch;

 
load([GLM_dir 'Animal' num2str(anm) '/Cluster_GroupingsMPF_'  direction_name{direction} '_A' num2str(anm) '_S0' num2str(sess) '_F' num2str(1) '.mat'])




TC = p_cells;



Plot_Order = squeeze(sum(MM_Ref(:,:,TC),1));


Max_Pos = zeros(numel(TC),1);

Sparsity=zeros(numel(TC),1);

MaxRate=zeros(numel(TC),1);
MeanRate=zeros(numel(TC),1);
Skaggs=zeros(numel(TC),1);


for uu = 1:numel(TC)

M0=Plot_Order(:,uu);
M0 = M0./max(M0);

if direction == 2
M0 = flipud(M0);

end

M0 = smoothdata(M0,'gaussian',5);

%[ma,mi]=max(M0); 

mi = sum((1:numel(M0)).*M0')/sum(M0);

Max_Pos(uu) = mi;

RateL = M0;
OccuL = ones(numel(RateL),1)/numel(RateL);

Sparsity(uu)=(sum((RateL.*OccuL))^2)/sum(OccuL.*(RateL.^2));
MaxRate(uu)=max(RateL);
MeanRate(uu)=mean(RateL);

OccAndSpike=find((RateL./MeanRate(uu))>0.000001);
SkMR=sum(RateL(OccAndSpike).*OccuL(OccAndSpike));

Skaggs(uu)=sum(RateL(OccAndSpike).*OccuL(OccAndSpike).*log(RateL(OccAndSpike)./SkMR));







end




pf_limits=zeros(numel(TC),3);
pf_size=zeros(numel(TC),1);
pf_asym=zeros(numel(TC),1);


for uu=1:numel(TC)




M0=squeeze(Plot_Order(:,uu));



[ma,mi]=max(M0);


jj=jj+1;


ii = localMaximum(-M0(1:mi),3,true);

im = find(M0(ii)<ma*0.3,1,'last');
if(numel(im)==0)
mi_l = 3;

else
mi_l = ii(im);
end



ii = localMaximum(-M0(mi+1:end),3,true);
im = find(M0(ii+mi)<ma*0.3,1,'first');
if(numel(im)==0)
mi_u = numel(M0)-2;

else
mi_u = ii(im)+mi;
end


pl_1=round(linspace(mi_l,mi,7));
pl_1=pl_1(1:6);
pl_2=round(linspace(mi,mi_u,7));
pl_2=pl_2(2:7);

pl_bi=[pl_1 mi pl_2];



pf_limits(uu,:) = [pl_bi(1) mi pl_bi(end)]./40*100;
pf_size(uu)= abs(pl_bi(end)-pl_bi(1));
%pf_asym(uu)= (abs(mi-mi_u)-abs(mi-mi_l))/(abs(mi-mi_u)+abs(mi-mi_l));


pf_asym(uu)= sum(M0(mi_l:mi_u).*(1:pf_size(uu)+1)')/sum(M0(mi_l:mi_u))-(mi+1-mi_l);

% if(pf_asym(uu)>0 && direction==1)
% figure(2001)
% plot(M0)
% pause()
% 
% end



end



save([destination_save 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')
%save(['/Users/federico/GitHub/Basins_GLM_paper/Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')




    end
    
    end
    






end