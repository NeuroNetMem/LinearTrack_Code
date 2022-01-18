close all

p_all = [];

cell_range = 10;
b_win = 100;
reg_range_l = [7 9 9];

direct_n = {'UP','DOWN'};
order_by_session = 0;

%G_Lim = linspace(-0.2,0.25,6);
load('Theta_Limits/Theta_Limits.mat')
G_Lim = th_limits_prc;


sl_r = zeros(4,numel(G_Lim)-1);
sl_r_er = zeros(4,numel(G_Lim)-1);
dir_r = zeros(4,numel(G_Lim)-1);
dir_r_er = zeros(4,numel(G_Lim)-1);

for G = [1 2 3]%numel(G_Lim)-1


    reg_range = reg_range_l(G);
    
for ll = 1:4

    
    s_poi = [];  

    
for  A = [1 2 3 4 5 6]
 
Ty = 3;
Ph = 3;



for S = 1:2
    
if (A==6 && S == 2)
continue
end
    
    if(order_by_session==1)
    load(['Theta_Limits/Theta_Limits_A' num2str(A) 'S' num2str(S)  '.mat'])
    G_Lim = th_limits_prc;
    end
    
    
for D = 1:2
 




load(['CCDati/CCDati_4GBothG_A' num2str(A) 'S' num2str(S) 'D' num2str(D) 'Ga' num2str(2) 'Sp' num2str(Ty) 'Ph' num2str(Ph) '.mat'])

load(['Cluster_Info/' '/Cluster_GroupingsMPF_' direct_n{D} '_A' num2str(A) '_S' num2str(S) '_F1.mat'   ])  


n_ord_1 = find(theta_scores>=G_Lim(G) & theta_scores<=G_Lim(G+1));
%n_ord_1 = find(theta_scores>G_Lim(2) & theta_scores<G_Lim(3));
n_ord_2 = find(theta_scores>=G_Lim(G) & theta_scores<=G_Lim(G+1));
%n_ord_2 = find(theta_scores>=G_Lim(1) & theta_scores<=G_Lim(1+1));

% n_ord_1 = find(theta_scores>=-1000 );
% n_ord_2 = find(theta_scores>=-1000 );


m_pos_1 = m_pos;
m_pos_2 = m_pos;


D_Gap = zeros(numel(n_ord_1),numel(n_ord_2));

  
for jj = 1:numel(n_ord_1)

% figure(71)
% clf;



ref_cell = n_ord_1(jj);



% kk=-(jj-cell_range)+1;
% 
% kk = max(0,kk);
kk=0;



for ii = n_ord_2'
    kk=kk+1;
c1 = squeeze(cc_all(:,ll,ref_cell,ii));
 
if(sum(c1)>0)


c1 = c1./sum(c1);



c_bari = c1(101-b_win:101+b_win);
if(sum(c_bari)>0)
p_bari = sum((1:numel(c_bari)).*c_bari')/sum(c_bari)-numel(c_bari)/2-0.5;

%[~,p_bari] = max(smoothdata(c_bari,'movmean',3));

%p_bari = p_bari - numel(c_bari)/2-0.5;

%c1 = c1 - mean(c1,2);

    %subplot(4,2,1+(ll-1)*2)
% figure(71)
% plot(-b_win:+b_win,c1(101-b_win:101+b_win),'LineWidth',3,'Color',[kk/numel(n_ord) 0 1-kk/numel(n_ord)])
% hold on
% scatter(p_bari,0.025,80,[kk/20 0 1-kk/20],'filled')
% xline(0,'LineWidth',2)
% xlabel('Delta t (1bin = 2ms)')
% ylabel('Spike Probability')
% ylim([0 0.03])
% pause()

% subplot(4,2,2+(ll-1)*2)
% scatter(p_bari,kk,80,[kk/20 0 1-kk/20],'filled')
% hold on
% xlim([-10 10])


k_dist = round(m_pos_2(ii)-m_pos_1(ref_cell));

s_poi = cat(1,s_poi,[p_bari,k_dist]);


D_Gap(ref_cell,ii) = p_bari;



end
end
end




end

% sco = sum(D_Gap,2);
% [~,p] = sort(sco,'descend');
% p_all = [p_all,p];
% 
% 
% 
% figure(101)
% imagesc(D_Gap(p,p))
% pause()


end
end
end


x = s_poi(:,1)/100*2*pi;
y = s_poi(:,2)*2.5;
colo = s_poi(:,2)/(cell_range*2+1);


y_take = find(abs(y)<40 & abs(x)<pi/2);
x_fit = x(y_take);
y_fit = y(y_take);

figure(60+G)
sgtitle(['Group ' num2str(G)])
subplot(3,2,ll)
scatter(y,x,60,[s_poi(:,2)/(cell_range*2+1) zeros(size(s_poi,1),1) 1-s_poi(:,2)/(cell_range*2+1)],'filled')
ylabel('Average Phase Difference (rad)')
xlabel('Place Field Distance (cm)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
hold on

xf = -b_win:b_win;
xf = xf/100*2*pi;

xline(0,'LineWidth',2);
yline(0,'LineWidth',2);

dt_m = zeros(reg_range*2+1,1);
c_k = 0;
for k = -reg_range:reg_range
    c_k = c_k+1;
    
pos = find(ismember(s_poi(:,2),[k:k]));
    
dt_m(c_k)=mean(s_poi(pos,1))/100*2*pi;

end

%dt_m = smoothdata(dt_m,'movmean',5);

plot((-reg_range:reg_range)*2.5,dt_m,'k','LineWidth',3)

cf = fit(y_fit,x_fit,'poly1');

cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = cf_coeff(1);
b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;

dir_r(ll,G)=a;
dir_r_er(ll,G)=a_uncert;


xlim([-40 40])
ylim([xf(1)/2 xf(end)/2])
%pause()

reg_l = smoothdata(dt_m,'movmean',7);

figure(91)
subplot(3,1,G)
plot((-reg_range:reg_range)*2.5,reg_l,'LineWidth',3,'Color',[1-ll/6 0 ll/6])
hold on
ylim([-1 1])

RR = polyfit(1:reg_range*2+1,reg_l',1);

cf = fit((1:reg_range*2+1)',reg_l,'poly1');

cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = cf_coeff(1);
b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;

sl_r(ll,G)=a;
sl_r_er(ll,G)=a_uncert;



figure(92)
subplot(3,2,(G-1)*2+1)
plot((-reg_range:reg_range)*2.5,reg_l-reg_l(1),'LineWidth',3,'Color',[1-ll/6 0 ll/6])
hold on
ylim([-1 1])
xlabel('Place Field Delta')
ylabel('Theta Phase Delta')

% figure(90)
% 
% plot((-reg_range:reg_range)*2.5,reg_l,'LineWidth',3,'Color',cmap(floor(250/numel(G_Lim))*G,:))
% hold on
% ylim([-1 1])
end

% for ll = 1:4
%     subplot(4,2,2+(ll-1)*2)
%     
%     
% end





end

figure(92)
subplot(3,2,[2 4 6])
errorbar(sl_r,sl_r_er,'LineWidth',3)
xlabel('Gamma (Medium -> Slow)')
ylabel('Slope')
% 
% figure(92)
% set(gcf,'renderer','Painters')
% saveas(gcf,'./Figures/Fig7_recapSpikeOrder','epsc')
%%
figure(61)
set(gcf,'renderer','Painters')
saveas(gcf,'./Figures/Fig7_spikeOrderG1','epsc')

figure(63)
set(gcf,'renderer','Painters')
saveas(gcf,'./Figures/Fig7_spikeOrderG3','epsc')
