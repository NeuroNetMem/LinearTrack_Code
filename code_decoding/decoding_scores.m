function [pos_score,pos_slope,pos_mean,seq_score,seq_slope,off_score,off_pos,r_line_seq] = decoding_scores(MM,Pos)

ss = find(std(MM,[],1)>1e-10);

pos_score = 0;
for bin = ss(:)'
    t_b = round(Pos(bin)-1:Pos(bin)+1);
    t_b = t_b(t_b<=size(MM,1) & t_b>0);
    pos_score=pos_score+nansum(MM(t_b,bin),1);
    
    
end
pos_score = pos_score/numel(ss);
pos_slope = (Pos(end)-Pos(1))/numel(Pos);
pos_mean = mean(Pos);





seq_score = 0;
off_score = 0;
off_pos = 0;
x = 1:size(MM,2);
m_p = ceil(size(MM,2)/2);
r_line_seq = zeros(size(x));

[XX,YY]=meshgrid(1:0.2:2.4,3:size(MM,1)-3);

XX = XX(:);
YY = YY(:);


for prm = 1:numel(XX)

    sl_t = XX(prm);
    off_t = YY(prm);

r_line = (x-m_p)*sl_t+off_t;
li_diff = nanmin(abs(r_line - Pos));




use_b=intersect(find(r_line > 0 & r_line <= size(MM,1)),ss);

if(numel(use_b)<3 || li_diff>3)
seq_score_prm = 0;
seq_slope_prm = 0;

else

seq_score_prm = 0;
for bin = use_b(:)'

    
    t_b = round(r_line(bin)-2:r_line(bin)+2);
    t_b = t_b(t_b<=size(MM,1) & t_b>0);
    seq_score_prm=seq_score_prm+nansum(MM(t_b,bin),1);
    

end

seq_score_prm = seq_score_prm/numel(use_b);
seq_slope_prm = sl_t;


end


if(seq_score_prm>off_score && sl_t<0.6)
off_score = seq_score_prm;
off_pos = mean(r_line);
    
end



if(seq_score_prm>seq_score)
seq_score = seq_score_prm;
seq_slope = seq_slope_prm;
r_line_seq = r_line;
end



end
% 
% [~,prm_m] = min(cost);
% 
% slope = XX(prm_m);
% off_ph = YY(prm_m);
% 
% 
% p_line = (x-m_p)*slope+off_ph;
% 
% figure(5001)
% imagesc(MMR)
% hold on
% plot(x,p_line)




