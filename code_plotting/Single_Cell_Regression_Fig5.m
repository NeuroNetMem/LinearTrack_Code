
%%

  set(gcf,'renderer','Painters')
  saveas(gcf,'./Figures/Fig5C','epsc')


%%
clear
close all 

glm_dir = '~/Data/Matteo_GLMEarly/';
prop_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';
save_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';

a_coor = linspace(0,2*pi,21);
a_coor = a_coor(1:20);

strat_h = [-1 1];
n_strat = numel(strat_h)-1;

direction_name = {'UP','DOWN'};

%Choosing spike types: 
%1 = All, 2 = Single Spikes e 3 = Bust Spike
sp_type = 2;

%Taking slow gamma from: 
%1 = PYR, 2 = RAD e 3=SLM

for gamma_pair = [2]



for anm = [2 3 5] 



gg_c = 0;
for gg=1:n_strat
gg_c = gg_c+1;
   

jj=0;
    

for th_site = [3]




    
    
    
    M1 = zeros(20,13,6);
    PP_Angle_All = zeros(13,6);
    PP_Angle_Sin = [];
    PP_Place_Sin = [];
    Pref_Phase = [];
    
    sl_c = [];
    off_ph_c = [];
    th_c = [];
    
    for direction = [1 2]
    for sess = [2 3]    


    
    
    gamma_pair_n = gamma_pair;




load([glm_dir 'Animal' num2str(anm) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(1) '_D' num2str(direction) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref = PlotVal_Ch;

load([glm_dir 'Animal' num2str(anm) '/GammaGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGammaPow_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(gamma_pair_n) '_T' num2str(th_site) '_D' num2str(direction) '_Sp' num2str(sp_type) '.mat'])
n_rep = 1;


load([glm_dir 'Animal' num2str(anm) '/Cluster_GroupingsMPF_' direction_name{direction} '_A' num2str(anm) '_S0' num2str(sess) '_F1.mat' ])


th_strat = discretize(theta_scores,strat_h);

TC = p_cells(th_strat==gg);

%TC=find(groups==gg);
%TC = find(groups);

th_gr = theta_scores(th_strat==gg);


Asy_Cell=zeros(numel(TC),3);
PP_Slope=zeros(3,1);

% figure(8001)
% plot(squeeze(sum(MM_Ref(:,:,TC),1)))
% pause()
ll=0;
for uu=TC'

ll=ll+1;

MM=squeeze(PlotVal_Ch(:,:,uu,:));
MM = repmat(MM,1,1,n_rep);

M0=squeeze(MM_Ref(:,:,uu));
M0=sum(M0,1);
%M0 = M0 - min(M0(4:37));

[ma,mi]=max(M0);

if(mi<4 && mi >37)
continue
end

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

if direction == 2
pl_bi = fliplr(pl_bi);

end

pl_bi_2 = round(linspace(mi_l,mi_u,13));



for rr = [6 5 4 3 2 1]
[sl,off_ph]=fit_linear(MM(:,pl_bi,rr));

sl_c(jj,rr) = -abs(sl);
off_ph_c(jj,rr) = off_ph;
end
th_c(jj)=th_gr(ll);


end
    end
    end
end

figure(100)
for rr = 1:6
    subplot(6,6,(anm-1)*6+rr)
scatter(th_c, sl_c(:,rr),15,'filled')%,'MarkerFaceColor','k','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
xlim([-0.15 0.2])
hold on

th_v = linspace(-0.15,0.2,7);
th_bin = discretize(th_c,th_v);

pl_th = zeros(6,1);

for kk=1:6
pl_th(kk) = mean(sl_c(th_bin==kk,rr));
end

plot(th_v(2:7),pl_th,'LineWidth',3,'Color','r')
xlabel('Theta Scores')
ylabel('Phase-Position Slope')

end
end
end
end


