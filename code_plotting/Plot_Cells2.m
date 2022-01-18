
%%

%%Saving
set(gcf,'renderer','Painters')
saveas(gcf,'Figure8','epsc')
% saveas(gcf,'suppl_A1_S1_D1_C41','epsc')



%%

% sigmoid parameters 
a1 = 10; % steepness (how fast is the transition between cold and hot colors) 
a2 = 0.55; % inflection (where is the transition happening, fraction of max value of the matrix)

colormap default 
C1 = colormap(); 
x1 = (1:size(C1,1))-1;
x1 = x1./max(x1);
x2 = 1./(1 + exp(-a1.*(x1-a2)));

figure(1)
subplot(2,1,1)
plot(linspace(0,1,256),x2)
xlabel('Original Value')
ylabel('New Value')

for ch=1:3
C1(:,ch)=interp1(x1,C1(:,ch),x2);
end


close all 


a_coor = linspace(0,2*pi,21);
a_coor = a_coor(1:20);

direction_name = {'UP','DOWN'};

aa=0;

for gg=1:1
field = gg;

x_all = [];
y_all = [];
c_all = [];
th_all = [];
z_all = [];

for anm = [6] 
    anm2 = anm;
 
    aa=aa+1; 

    for sess = [1]  
    for direction = [2]
    

        load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(direction) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref = PlotVal_Ch;
 
%load(['/Users/federico/GitHub/Basins_GLM_paper/AnimalNG' num2str(anm2) '/Cluster_Groupings' direction_name{direction} '_A' num2str(anm2) '_S' num2str(sess) '.mat'  ])


load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/Cluster_GroupingsMPF_'  direction_name{direction} '_A' num2str(anm2) '_S' num2str(sess) '_F' num2str(field) '.mat'])


th_s = theta_scores;


TC = p_cells;







figure(100+gg)
kk= 0;
for uu=TC'
kk = kk+1;
M0=squeeze(MM_Ref(:,:,uu));

if(direction == 2)
M0 = fliplr(M0);

end


FR=sum(M0,1)';
FR = FR./max(FR);
FR_I=interp1(1:numel(FR),FR,1.5:1:numel(FR)-0.5);

%FR_I = diff(FR);

PP_Angle = zeros(size(M0,2),1);
PP_Ref = 0;
    for tt = 1:size(M0,2)
        PP_Ref = PP_Ref + sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt));
        PP_Angle(tt) = angle(sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt)));%a_coor(Ma_Ref);



    end
 
 PP_Angle = PP_Angle-angle(PP_Ref);   
PP_Angle(PP_Angle<-pi) = PP_Angle(PP_Angle<-pi)+2*pi;     
PP_Angle(PP_Angle>pi) = PP_Angle(PP_Angle>pi)-2*pi; 

 PP_Delta = diff(smoothdata(unwrap(PP_Angle),'movmean',5));   
    
 subplot(311)
 imagesc(repmat(M0,2,1))
%   colorbar
% colormap(C1)
 set(gca,'YDir','normal')
 
 subplot(312)
 plot(FR_I,'LineWidth',2)
 hold on
 plot(unwrap(PP_Delta),'LineWidth',2)
 refline(0,0)
 ylim([-1 1])
 xlabel('Position')
 ylabel('FR and DeltaPhi')
 hold off 
 
 subplot(313)

 
 x = FR_I; 
y = PP_Delta';
z = zeros(size(x));
t = 1:numel(x);
lineColor = t;  

if (th_s(kk))>-2   % Here you can select to see cells based on their score

surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 8);
refline(0,0) 
xlabel('Firing Rate')
ylabel('DeltaPhi')
sgtitle(['Animal ' num2str(anm) ' Session ' num2str(sess) ' Direction ' num2str(direction) ' Cell ' num2str(uu) '; TS = ' num2str(th_s(kk))])
pause()
cla;

end

end



    end
    end
    
    
    
    


end




end


%% TWO DIRECTIONS PLOTTING

%%Saving
set(gcf,'renderer','Painters')
% saveas(gcf,'PL_A6_S1_C','epsc')
% saveas(gcf,'PPPL_A5_S2_D1_C61','epsc')
  saveas(gcf,'UPDOWN_A2_S1_C66','epsc')

%% TWO DIRECTIONS PLOTTING


% sigmoid parameters 
a1 = 10; % steepness (how fast is the transition between cold and hot colors) 
a2 = 0.25; % inflection (where is the transition happening, fraction of max value of the matrix)

colormap default 
C1 = colormap(); 
x1 = (1:size(C1,1))-1;
x1 = x1./max(x1);
x2 = 1./(1 + exp(-a1.*(x1-a2)));

figure(1)
subplot(2,1,1)
plot(linspace(0,1,256),x2)
xlabel('Original Value')
ylabel('New Value')

for ch=1:3
C1(:,ch)=interp1(x1,C1(:,ch),x2);
end


close all 

a_coor = linspace(0,2*pi,21);
a_coor = a_coor(1:20);

direction_name = {'UP','DOWN'};

aa=0;

for gg=1:1
field = gg;

x_all = [];
y_all = [];
c_all = [];
th_all = [];
z_all = [];

for anm = [6] 
    anm2 = anm;

    aa=aa+1;

for direction = [1 2]
for sess = [1 2]    
    

load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(1) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref_Up = PlotVal_Ch;


load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/Cluster_GroupingsMPF_'  direction_name{1} '_A' num2str(anm2) '_S' num2str(sess) '_F' num2str(field) '.mat'])

p_cells_Up = p_cells;
th_s_Up = theta_scores;

load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/GammaEnvelopeGLM_SF_A' num2str(anm) '_S' num2str(sess) '/GLMPyrGamma_A' ...
    num2str(anm) '_S' num2str(sess) '_G' num2str(1) '_T' num2str(3) '_D' num2str(2) '_Sp1.mat'],'PlotVal_Ch')

MM_Ref_Do = PlotVal_Ch;


load(['C:/Users/Matteo/Desktop/github_Basin_PP/AnimalNG' num2str(anm2) '/Cluster_GroupingsMPF_'  direction_name{2} '_A' num2str(anm2) '_S' num2str(sess) '_F' num2str(field) '.mat'])

p_cells_Do = p_cells;
th_s_Do = theta_scores;





TC = intersect(p_cells_Up,p_cells_Do);







figure(100+gg)
kk= 0;
for uu=TC'
kk = kk+1;
M0=squeeze(MM_Ref_Up(:,:,uu));

if(direction == 2)
M0 = fliplr(M0);

end


FR=sum(M0,1)';
FR = FR./max(FR);
FR_I=interp1(1:numel(FR),FR,1.5:1:numel(FR)-0.5);

%FR_I = diff(FR);

PP_Angle = zeros(size(M0,2),1);
PP_Ref = 0;
    for tt = 1:size(M0,2)
        PP_Ref = PP_Ref + sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt));
        PP_Angle(tt) = angle(sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt)));%a_coor(Ma_Ref);



    end
 
 PP_Angle = PP_Angle-angle(PP_Ref);   
PP_Angle(PP_Angle<-pi) = PP_Angle(PP_Angle<-pi)+2*pi;     
PP_Angle(PP_Angle>pi) = PP_Angle(PP_Angle>pi)-2*pi; 

 PP_Delta = diff(smoothdata(unwrap(PP_Angle),'movmean',5));   
    
 subplot(321)
 imagesc(repmat(M0,2,1))
 colorbar
colormap(C1)

 set(gca,'YDir','normal')
  title('UP')
 subplot(323)
 plot(FR_I,'LineWidth',2)
 hold on
 plot(unwrap(PP_Delta),'LineWidth',2)
 refline(0,0)
 xlabel('Position')
 ylabel('FR and DeltaPhi')
 
 hold off 
 
 subplot(325)

 
 x = FR_I; 
y = PP_Delta';
z = zeros(size(x));
t = 1:numel(x);
lineColor = t;  

if (th_s_Up(kk))>-2   % Here you can select to see cells based on their score
cla;
surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 8);
refline(0,0) 
xlabel('Firing Rate')
ylabel('DeltaPhi')
sgtitle(['Animal ' num2str(anm) ' Session ' num2str(sess) ' Direction ' num2str(direction) ' Cell ' num2str(uu) ])



end





M0=squeeze(MM_Ref_Do(:,:,uu));
M0 = fliplr(M0);




FR=sum(M0,1)';
FR = FR./max(FR);
FR_I=interp1(1:numel(FR),FR,1.5:1:numel(FR)-0.5);

%FR_I = diff(FR);

PP_Angle = zeros(size(M0,2),1);
PP_Ref = 0;
    for tt = 1:size(M0,2)
        PP_Ref = PP_Ref + sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt));
        PP_Angle(tt) = angle(sum(M0(:,tt).*exp(1i*a_coor'))/nansum(M0(:,tt)));%a_coor(Ma_Ref);



    end
 
 PP_Angle = PP_Angle-angle(PP_Ref);   
PP_Angle(PP_Angle<-pi) = PP_Angle(PP_Angle<-pi)+2*pi;     
PP_Angle(PP_Angle>pi) = PP_Angle(PP_Angle>pi)-2*pi; 

 PP_Delta = diff(smoothdata(unwrap(PP_Angle),'movmean',5));   
    
 subplot(322)
 imagesc(repmat(M0,2,1))
colormap(C1)
 colorbar
 
 set(gca,'YDir','normal')
 title('DOWN')
 
 subplot(324)
 plot(FR_I,'LineWidth',2)
 hold on
 plot(unwrap(PP_Delta),'LineWidth',2)
 refline(0,0)
 xlabel('Position')
 ylabel('FR and DeltaPhi')
 
 hold off 
 
 subplot(326)

 
 x = FR_I; 
y = PP_Delta';
z = zeros(size(x));
t = 1:numel(x);
lineColor = t;  

if (th_s_Do(kk))>-2   % Here you can select to see cells based on their score
cla;
surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 8);
refline(0,0) 
xlabel('Firing Rate')
ylabel('DeltaPhi')
sgtitle(['Animal ' num2str(anm) ' Session ' num2str(sess) ' Cell ' num2str(uu) ])
pause()


end



    
end
    
    
end
    


end




end



end
