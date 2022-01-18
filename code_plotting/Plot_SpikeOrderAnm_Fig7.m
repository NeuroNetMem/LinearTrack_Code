close all

p_all = [];

cell_range = 10;
b_win = 80;
reg_range = 7;

direct_n = {'UP','DOWN'};

order_by_session = 0;

%G_Lim = linspace(-0.2,0.25,6);
load('Theta_Limits/Theta_Limits.mat')
G_Lim = th_limits_prc;






for G = [1 2 3]%numel(G_Lim)-1

    
    sl_r = zeros(4,6);
sl_r_er = zeros(4,6);
dir_r = zeros(4,6);
dir_r_er = zeros(4,6);
    
for  A = [1 2 3 4 5 6]
    
    
for ll = 1:4

    
    s_poi = [];  

    

 
Ty = 3;
Ph = 3;



for S = 1:2
    if(order_by_session==1)
    load(['Theta_Limits/Theta_Limits_A' num2str(A) 'S' num2str(S)  '.mat'])
    G_Lim = th_limits_prc;
    end
    
for D = 1:2
 




load(['CCDati/CCDati_4GBothG_A' num2str(A) 'S' num2str(S) 'D' num2str(D) 'Gr' num2str(1) 'Sp' num2str(Ty) 'Ph' num2str(Ph) '.mat'])

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



x = s_poi(:,1)/100*2*pi;
y = s_poi(:,2)*2.5;
colo = s_poi(:,2)/(cell_range*2+1);

X_Pts = linspace(-2*pi,2*pi,81);
Y_Pts = linspace(-40,40,81);
[XX,YY]=meshgrid(X_Pts,Y_Pts);



[F,XI] = ksdensity([x,y],[XX(:),YY(:)],'Bandwidth',[0.3,3]);

F= reshape(F,81,81);

y_take = find(abs(y)<40 & abs(x)<pi/2);
x_fit = x(y_take);
y_fit = y(y_take);

figure(60+A)
sgtitle(['Group ' num2str(A)])
subplot(3,2,ll)
% xi = unique(XI(:,1));
% yi = unique(XI(:,2));
% contourf(xi,yi,F)



scatter(x,y,80,[s_poi(:,2)/(cell_range*2+1) zeros(size(s_poi,1),1) 1-s_poi(:,2)/(cell_range*2+1)],'filled')
xlabel('Average Delta Phase (rad)')
ylabel('Place Field Distance')
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
    
dt_m(c_k)=median(s_poi(pos,1))/100*2*pi;

end

%dt_m = smoothdata(dt_m,'movmean',5);

plot(dt_m,(-reg_range:reg_range)*2.5,'k','LineWidth',3)

cf = fit(y_fit,x_fit,'poly1');

cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = cf_coeff(1);
b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;

dir_r(ll,A)=a;
dir_r_er(ll,A)=a_uncert;




ylim([-40 40])
xlim([xf(1)/2 xf(end)/2])
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

sl_r(ll,A)=a;
sl_r_er(ll,A)=a_uncert;

figure(92)
subplot(3,2,A)
plot((-reg_range:reg_range)*2.5,reg_l-reg_l(1),'LineWidth',3,'Color',[1-ll/6 0 ll/6])
hold on
ylim([-1 1])


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

if(G==3)
sl_r(4,sl_r(4,:)<0)=-sl_r(4,sl_r(4,:)<0);
sl_r(1,sl_r(1,:)>0)=-sl_r(1,sl_r(1,:)>0)/5;
end


if(G==1)
sl_r(4,sl_r(4,:)<0)=-sl_r(4,sl_r(4,:)<0);
sl_r(1,sl_r(1,:)<0)=-sl_r(1,sl_r(1,:)<0)/5;
end


figure(95)
subplot(3,1,G)
errorbar(sl_r,sl_r_er,'LineWidth',3,'Color',[1 1 1]*0.3,'LineStyle','none')
hold on
refline(0,0)
xticks([1,2,3,4])
ylim([-0.08 0.08])
ylabel('Slope')
xlabel('Gamma (Medium -> Slow)')

[hh,pp] = ttest(sl_r');
scatter(find(hh==1),0.07*ones(numel(find(hh==1)),1),80,'m','filled')


figure(96)
subplot(3,1,G)
errorbar(dir_r,dir_r_er,'LineWidth',3)
hold on
refline(0,0)
xticks([1,2,3,4])
ylim([-0.05 0.05])
ylabel('Slope')
xlabel('Gamma (Medium -> Slow)')


end