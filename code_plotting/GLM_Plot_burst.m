close all 



a_coor = linspace(0,2*pi,21);
a_coor = a_coor(1:20);


 
glm_dir = '~/Data/Matteo_GLMEarly/';


load('Theta_Limits/Theta_Limits.mat')
strat_h = th_limits_lin;
n_strat = numel(strat_h)-1;


direction_name = {'UP','DOWN'};


 

for gamma_pair = [2]



for sp_type = 2:3 



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
    
    for anm = [2 3 5]
    for direction = [1 2]
    for sess = [2 3]    


    
    
    gamma_pair_n = gamma_pair;




load([glm_dir 'Animal' num2str(anm) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(1) '_D' num2str(direction) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref = PlotVal_Ch;



%C:\Users\Matteo\Desktop\github_Basin_PP
%/Users/federico/GitHub/Basins_GLM_paper/AnimalNG


load([glm_dir 'Animal' num2str(anm) '/GammaGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGammaPow_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(gamma_pair_n) '_T' num2str(th_site) '_D' num2str(direction) '_Sp' num2str(sp_type) '.mat'])
n_rep = 1;


load([glm_dir 'Animal' num2str(anm) '/Cluster_GroupingsMPF_' direction_name{direction} '_A' num2str(anm) '_S0' num2str(sess) '_F1.mat' ])


th_strat = discretize(theta_scores,strat_h);

TC = p_cells(th_strat==gg);





Asy_Cell=zeros(numel(TC),3);
PP_Slope=zeros(3,1);

% figure(8001)
% plot(squeeze(sum(MM_Ref(:,:,TC),1)))
% pause()

for uu=TC'


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

% mi_l=find(M0(1:mi)<ma*0.2,1,'last')+1;
% 
% if(numel(mi_l)==0)
% mi_l = 3;
% 
% end

ii = localMaximum(-M0(mi+1:end),3,true);
im = find(M0(ii+mi)<ma*0.5,1,'first');
if(numel(im)==0)
mi_u = numel(M0)-2;

else
mi_u = ii(im)+mi;
end



% mi_u=find(M0(mi+1:end)<ma*0.2,1,'first')-1;
% mi_u=mi_u+mi;
% 
% if(numel(mi_u)==0)
% mi_u = numel(M0)-2;
% 
% end



% mi_l=max(1,mi-3);
% mi_u=min(numel(M0),mi+3);

pl_1=round(linspace(mi_l,mi,7));
pl_1=pl_1(1:6);
pl_2=round(linspace(mi,mi_u,7));
pl_2=pl_2(2:7);

pl_bi=[pl_1 mi pl_2];

if direction == 2
pl_bi = fliplr(pl_bi);

end

pl_bi_2 = round(linspace(mi_l,mi_u,13));




% M1=MM(:,mi_l:mi_u,1);
% M2=MM(:,mi_l:mi_u,2);
% M3=MM(:,mi_l:mi_u,3);

PP_Angle_T = zeros(13,1);
ss= 0 ;
for rr =[6 5 4 3 2 1]
ss = ss+1;
    MM(:,:,rr)=MM(:,:,rr)/sum(MM(:,:,rr),'all');
    M1(:,:,rr) = M1(:,:,rr) + MM(:,pl_bi,rr);%/sum(MM(:,pl_bi,rr),'all');

    for tt = 1:13

PP_Angle_T(tt) = angle(sum(MM(:,pl_bi(tt),rr).*exp(1i*a_coor'))/nansum(MM(:,pl_bi(tt),rr)));%a_coor(Ma_Ref);



    end

    if(ss == 1)
    C_Ref = PP_Angle_T(:);
    end
    
PP_Angle_T = PP_Angle_T;% - C_Ref;%PP_Angle_T(7);    
    
PP_Angle_T(PP_Angle_T>pi) = -(2*pi-PP_Angle_T(PP_Angle_T>pi));
PP_Angle_T(PP_Angle_T<-pi) = (2*pi-abs(PP_Angle_T(PP_Angle_T<-pi)));


%PP_Angle_T(1:13) = smoothdata(unwrap(PP_Angle_T(1:13),pi),'gaussian',7);

PP_Angle_T(PP_Angle_T>pi) = -(2*pi-PP_Angle_T(PP_Angle_T>pi));
PP_Angle_T(PP_Angle_T<-pi) = (2*pi-abs(PP_Angle_T(PP_Angle_T<-pi)));

%PP_Angle_All(:,rr) = PP_Angle_All(:,rr) + exp(1i.*PP_Angle_T);
PP_Angle_Sin(:,jj,rr) = exp(1i.*PP_Angle_T);    
PP_Place_Sin(:,jj,rr) = squeeze(sum(MM(:,pl_bi,rr)/sum(MM(:,pl_bi,rr),'all'),1));    
    
    %C_Ref = PP_Angle_T(7);
    
Pref_Phase(:,jj,rr,1) = sum(MM(:,pl_bi(1:5),rr),2)./sum(MM(:,pl_bi(1:5),rr),'all');
Pref_Phase(:,jj,rr,2) = sum(MM(:,pl_bi(6:8),rr),2)./sum(MM(:,pl_bi(6:8),rr),'all');
Pref_Phase(:,jj,rr,3) = sum(MM(:,pl_bi(9:13),rr),2)./sum(MM(:,pl_bi(9:13),rr),'all');

end

%for aa = 1:13
%PP_Angle_Sin(aa,jj,:) = smoothdata(unwrap(angle(PP_Angle_Sin(aa,jj,:))),'gaussian',7);
%PP_Angle_Sin(aa,jj,:) = exp(1i.*PP_Angle_Sin(aa,jj,:));
%end




 
%     figure(2001) 
%     clf; 
%     for zz = 1:10
%     
%     plot(unwrap(angle(PP_Angle_Sin(4:10,jj,zz))),'LineWidth',2,'Color',[1 - zz/10 0 zz/10])
%     hold on
%     %ylim([-pi pi])
%     axis square
%     end
%     pause()


    
end
    end
    end
    end
    
for rr =1:6
    
PP_Angle_All(:,rr) = sum(PP_Angle_Sin(:,:,rr),2);%PP_Angle_All(:,rr) + exp(1i.*PP_Angle_T);
end    
    
    
figure(40)    
    for zz = 6:-1:1
    subplot(4,n_strat,(sp_type-1)*n_strat+gg_c)
    hold on
    p_ang = angle(PP_Angle_All(3:11,zz));
%     p_ang = unwrap(angle(PP_Angle_All(3:11,zz)));
%     if(p_ang(5)>pi)
%     p_ang = p_ang-2*pi;
%     end
%     if(p_ang(5)<-pi)
%     p_ang = p_ang+2*pi;
%     end
    p_ang=unwrap(p_ang,pi);
    if(p_ang(1)<-pi/2)
        p_ang = p_ang+2*pi;
    end
    plot(p_ang,'LineWidth',3,'Color',[1 - zz/6 0 zz/6])
    ylim([-pi*3/2 pi*3/2])
    ylabel('Theta Phase')
    xlabel('Place Field Position')
    axis square
    end
    
    sgtitle('Single Cell Average Regression')
%     
% 
%     figure(30+gg_c)
%     subplot(1,3,1)
%     plot((angle(PP_Angle_All(5,end:-1:3))),'LineWidth',3)
%     ylabel('Delta Phase')
%     xlabel('PYR Slow Gamma to RAd Slow Gamma')
%     title('Early Place Field')
%     hold on
%     subplot(1,3,2)
%     plot((angle(PP_Angle_All(7,end:-1:3))),'LineWidth',3)
%         ylabel('Delta Phase')
%     xlabel('PYR Slow Gamma to RAd Slow Gamma')
%     title('Middle Place Field')
%     hold on
%     subplot(1,3,3)
%     plot((angle(PP_Angle_All(9,end:-1:3))),'LineWidth',3)
%         ylabel('Delta Phase')
%     xlabel('PYR Slow Gamma to RAd Slow Gamma')
%     title('Late Place Field')
%     hold on
    
%     for zz = 1:10
%     subplot(2,5,zz)
%     plot(unwrap(angle(PP_Angle_Sin(4:10,:,zz))),'LineWidth',2)
%     %ylim([-pi pi])
%     axis square
%     end

    
figure(100+gg_c)
aa = 0;
for zz = 6:-1:1
    aa = aa+1;
    subplot(2,6,aa+(sp_type-2)*6)
    imagesc(repmat(M1(:,:,zz),2,1))
    ylabel('Phase')
    xlabel('Place Field Pos')
    set(gca,'YDir','normal')
    axis square
end

sgtitle('Spike Density')

figure(150)
subplot(4,n_strat,(sp_type-1)*n_strat+gg_c)
for zz = 1:6
   M00 = repmat(M1(:,:,zz),2,1);
   M00 = M00./max(M00,[],'all');
    z_thr = 0.75;
   contour(M00,[z_thr z_thr],'LineColor',[1-zz/6 0 zz/6],'LineWidth',2) 
   hold on
    
    
end

sgtitle('Spike Density for Single and Burst Spikes (rows)')

if(gg_c == 3 && sp_type==3)
figure(151)
zz = 1;
M00 = repmat(M1(:,:,zz),2,1);
   M00 = M00./max(M00,[],'all');
    z_thr = 0.85;
   contour(M00,[z_thr z_thr],'LineColor',[1-(zz-1)/5 0 (zz-1)/5],'LineWidth',2,'LineStyle','--') 
   hold on
end
if(gg_c == 3 && sp_type==2)
figure(151)
zz = 1;
M00 = repmat(M1(:,:,zz),2,1);
   M00 = M00./max(M00,[],'all');
    z_thr = 0.85;
   contour(M00,[z_thr z_thr],'LineColor',[1-(zz-1)/5 0 (zz-1)/5],'LineWidth',2) 
   hold on
   
zz = 6;
M00 = repmat(M1(:,:,zz),2,1);
   M00 = M00./max(M00,[],'all');
    z_thr = 0.85;
   contour(M00,[z_thr z_thr],'LineColor',[1-(zz-1)/5 0 (zz-1)/5],'LineWidth',2) 
   hold on
end



xlabel('Position')
ylabel('Phase')

figure(50)

ss = 0 ;
for zz = [6 5 4 3 2 1]
    ss = ss+1;
    
    Ma = M1(:,:,zz);
    
    if(ss == 1)
    Ma_Ref = ones(13,1);
    for tt = 1:13
    Ma_Ref(tt) = angle(sum(Ma(:,tt).*exp(1i*a_coor'))/nansum(Ma(:,tt)));
    %PP_Angle(tt) = a_coor(find(Ma(:,tt)==max(Ma(:,tt))))-Ma_Ref;
    end    

    end
    
    PP_Angle = ones(13,1);
    for tt = 1:13
    PP_Angle(tt) = angle(sum(Ma(:,tt).*exp(1i*a_coor'))/nansum(Ma(:,tt)));
    %PP_Angle(tt) = a_coor(find(Ma(:,tt)==max(Ma(:,tt))))-Ma_Ref;
    end
    %PP_Angle = PP_Angle - Ma_Ref;
    
    PP_Angle(PP_Angle>pi) = -(2*pi-PP_Angle(PP_Angle>pi));
    PP_Angle(PP_Angle<-pi) = (2*pi-abs(PP_Angle(PP_Angle<-pi)));


    %PP_Angle(1:13) = smoothdata(unwrap(PP_Angle(1:13),pi),'gaussian',5);

%     PP_Angle(PP_Angle>pi) = -(2*pi-PP_Angle(PP_Angle>pi));  
%     PP_Angle(PP_Angle<-pi) = (2*pi-abs(PP_Angle(PP_Angle<-pi)));
    
    p_ang=unwrap(PP_Angle(3:11),pi);
    if(p_ang(1)<-pi/2)
        p_ang = p_ang+2*pi;
    end
    
    
    subplot(4,n_strat,(sp_type-1)*n_strat + gg_c)
    hold on
    plot(p_ang,'LineWidth',2,'Color',[1 - zz/6 0 zz/6])
    ylim([-pi*3/2 pi*3/2])
    axis square
    
    



    
end

sgtitle('Spike Density Regression')

figure(3000)
%clf;
z_i = 0;
for zz = [1 6]
z_i = z_i+1;
subplot(2,n_strat,gg_c+(sp_type-2)*n_strat)
%plot(Pref_Phase(:,:,zz))
hold on
%plot(mean(Pref_Phase(:,:,zz),2),'Color',[1-zz/10 0 zz/10],'LineWidth',2)
AA = mean(PP_Place_Sin(:,:,zz),2);
AA_s = std(PP_Place_Sin(:,:,zz),[],2)./sqrt(size(PP_Place_Sin,2));
%plot(AA,'Color',[1-zz/6 0 zz/6],'LineWidth',2)

errorbar(AA,AA_s,'Color',[1-zz/6 0 zz/6],'LineWidth',2)


COM = sum((1:13).*AA');
%scatter(COM,0,50,[1-zz/6 0 zz/6],'filled')
xlabel('In-Field Position')
ylabel('Normalized Firing Rate')
switch gg_c
    case 1
        title('Phase Locking')
            case 2
        title('Intermediate')
            case 3
        title('Phase Precessing')
end

%ylim([0 0.2])


end

[hh,pp]=ttest2(PP_Place_Sin(:,:,1)',PP_Place_Sin(:,:,6)')

end
end
end
end


%%
figure(151)
  set(gcf,'renderer','Painters')
saveas(gcf,'./Figures/Fig6Supp_Contour','epsc')
% figure(50)
%   set(gcf,'renderer','Painters')
% saveas(gcf,'./Figures/Fig5b','epsc')
% figure(150)
%   set(gcf,'renderer','Painters')
% saveas(gcf,'./Figures/Fig5c','epsc')
% for ff = 101:106
%    figure(ff)
%      set(gcf,'renderer','Painters')
% saveas(gcf,['./Figures/Fig5Anm' num2str(ff-100) 'GLM'],'epsc')
% end
