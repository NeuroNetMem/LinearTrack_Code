close all 


direction_name = {'UP','DOWN'};

glm_dir = '~/Data/Matteo_GLMEarly';
process_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly';

aa=0;
for anm = [2 3 5] 



    anm2 = anm;    

    
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

load([glm_dir '/Animal' num2str(anm2) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(direction) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref = PlotVal_Ch;


load([glm_dir '/Animal' num2str(anm2) '/GammaGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGammaPow_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(direction) '_Sp1.mat'])



     load([glm_dir '/Animal' num2str(anm2) '/Cluster_GroupingsMPF_' direction_name{direction} '_A' num2str(anm2) '_S0' num2str(sess) '_F1.mat'   ])
     

     
TC = p_cells; 
     
Plot_Order = squeeze(sum(MM_Ref(:,:,TC),1));
m_pos = zeros(numel(TC),1);
for uu = 1:numel(TC)

M0=Plot_Order(:,uu);
M0 = M0./max(M0);

if direction == 2
M0 = flipud(M0);

end

M0 = smoothdata(M0,'gaussian',5);

%[ma,mi]=max(M0); 

mi = sum((1:numel(M0)).*M0')/sum(M0);

m_pos(uu) = mi;
Plot_Order(:,uu) = M0;


end



[~,n_ord] = sort(m_pos,'ascend');
Plot_UnOrder = Plot_Order;
Plot_Order = Plot_Order(:,n_ord);


%PLOT ORDERED CELLS


imagesc(Plot_Order')
title(['Animal ' num2str(anm) ' Session ' num2str(sess) ' Direction ' num2str(direction)])
%set(gca,'YDir','normal')
%pause()






pf_limits=zeros(numel(TC),3);




for uu=1:numel(TC)




M0=squeeze(Plot_UnOrder(:,uu));



[ma,mi]=max(M0);


jj=jj+1;


ii = localMaximum(-M0(1:mi),3,true);

im = find(M0(ii)<ma*0.5,1,'last');
if(numel(im)==0)
mi_l = 3;

else
mi_l = ii(im);
end



ii = localMaximum(-M0(mi+1:end),3,true);
im = find(M0(ii+mi)<ma*0.5,1,'first');
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

% if direction == 2
% pl_bi = fliplr(pl_bi);
% 
% end


% TO PLOT PHASE-POSITION FOR EACH CELL
% 
% imagesc(repmat(MM_Ref(:,pl_bi,uu),2,1))
% set(gca,'YDir','normal')
% pause()

pf_limits(uu,:) = [pl_bi(1) mi pl_bi(end)]./40*100;

end

[~,n_ord_2] = sort(mean(pf_limits,2),'ascend');



save([process_dir '/Cluster_Info/' '/Cluster_GroupingsMPF_' direction_name{direction} '_A' num2str(anm2) '_S0' num2str(sess) '_F1.mat'   ],'m_pos','pf_limits','-append')
save([glm_dir '/Animal' num2str(anm2) '/Cluster_GroupingsMPF_' direction_name{direction} '_A' num2str(anm2) '_S0' num2str(sess) '_F1.mat'   ],'m_pos','pf_limits','-append')

    end
    
    end
    






end



