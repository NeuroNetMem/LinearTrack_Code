close all
deco_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';
proc_dir = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';
for A = [2 3 5]

    disp(['Animal ' num2str(A)])
    
for gr = [0 1 2 3]
    disp(['Group ' num2str(gr)])
    
for sp = 1:3
for ph_type = 0:3
for ph_half = 0:2    

pc_all = [];
t_gamma_all = [];
Phase_All = [];
PosD_All = [];

Ph_All = [];
De_All = [];


for S = 2:3
    for D = 1:2
        if(gr>0)
load([proc_dir '/Decoding/DecodeRunning_A' num2str(A) '_S' num2str(S) '_G1_T3_D' num2str(D) '_Sp' num2str(sp) 'Gr' num2str(gr) '.mat'])
        else
load([proc_dir '/Decoding/DecodeRunning_A' num2str(A) '_S' num2str(S) '_G1_T3_D' num2str(D) '_Sp' num2str(sp) '.mat'])
        end
        
phzr = TH_Phase;


phzr(phzr < 0) = phzr(phzr < 0) + 2 * pi;
phzr = phzr+ph_type/2*pi;
phzr(phzr > 2*pi) = phzr(phzr > 2*pi) - 2 * pi;

th_edges = find(diff(phzr)<0)+1;        
    
Cycle_N = discretize(Take_bin,th_edges);

Ty_De_All = normalize(Ba_De_All,1,'norm',1);

if(D==2)
Ty_De_All = flipud(Ty_De_All);
Pos_U = 100 - Pos_U;
end

Pos_U = Pos_U+2.5;



GP = zeros(numel(GammaPow_U),3);

   for gamma = 1:3
           if(gr>0)
    load([proc_dir '/Decoding/DecodeRunning_A' num2str(A) '_S' num2str(S) '_G' num2str(gamma) '_T3_D' num2str(D) '_Sp' num2str(sp) 'Gr' num2str(gr) '.mat'],'GammaPow_U')
           else
    load([proc_dir '/Decoding/DecodeRunning_A' num2str(A) '_S' num2str(S) '_G' num2str(gamma)  '_T3_D' num2str(D) '_Sp' num2str(sp) '.mat'],'GammaPow_U')
           end
           
   GP(:,gamma)=GammaPow_U;        
           
          
           
   end


ll=0;


for cy = 1:max(Cycle_N)

        

    
        t_b = find(Cycle_N==cy);
        cy_ph = phzr(t_b);
    t_gamma = zeros(1,3);
    
    for gamma = 1:3
    t_gamma(gamma) = nanmean(GP(t_b,gamma));
    end
    
    
   %t_gamma = discretize(nanmax(GammaPow_U(t_b)),g_lim );

   t_gamma_1 = nanmean(GammaPow_U(t_b));
   t_gamma_2 = nanmax(GammaPow_U(t_b));
    
%     if(ismember(t_gamma,[G G+1]))
%         continue
%     end

   
       n_lim = 5;
   if(ph_half==1)
   t_ph = find(cy_ph<pi);
   t_b = t_b(t_ph);
   n_lim = 3;
   elseif(ph_half==2)
   t_ph = find(cy_ph>pi);
   t_b = t_b(t_ph);
   n_lim = 3;
   end
    
    
    
    
    th_prob = Ty_De_All(:,t_b);
    ss = find(std(th_prob,[],1)>1e-10);     
        
   if(numel(ss)>n_lim)
       ll = ll + 1;
        aa = mod(ll,30);
   
        
        
        if(aa==0)
   aa = 30;
        end  

   
   [pos_score,pos_slope,pos_mean,seq_score,seq_slope,off_score,off_pos,r_line]=decoding_scores(th_prob,Pos_U(t_b)/2.5);
   
% %    
%         subplot(5,6,aa)
%     imagesc(th_prob)
%     set(gca,'YDir','normal')
%     hold on
%     plot(Pos_U(t_b)/2.5,'LineWidth',3,'Color','m')
%     plot(r_line,'LineWidth',3,'Color','g')
%     
%     hold off
%     title([num2str(cy) ' ' num2str(pos_score) ' ' num2str(seq_score) ' ' num2str(seq_slope) ])
%     caxis([0 0.1])
% %     h= refline(0,Pos_U(t_b(1))/2.5);
% %     h.LineWidth = 3;
%    if(aa==30)
%     pause()
%    end
   
   pc_all = [pc_all;pos_score,pos_slope,pos_mean,seq_score,seq_slope,off_score,off_pos];
   t_gamma_all = [t_gamma_all;t_gamma];
   
   
   end
end



    end
end



save([proc_dir 'ThSeqDati/ThSeqDati_A' num2str(A) 'Gr' num2str(gr) 'Sp' num2str(sp) 'Ph' num2str(ph_type) 'H' num2str(ph_half)],'pc_all','t_gamma_all')




% figure(101)
% for tt = 1:4
%     subplot(2,2,tt)
%     if(ismember(tt,[1 3]))
% histogram(pc_all(:,tt),0:0.02:1,'Normalization','probability')
% hold on
%     elseif(ismember(tt,[2 4]))
%     histogram(pc_all(:,tt),0:0.02:5,'Normalization','probability')
%     end
% end




end
end

end

% pause()

end
end
