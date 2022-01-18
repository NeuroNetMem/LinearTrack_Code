function [slope,off_ph] = fit_linear(MM)

MM = MM./sum(MM,'all');

MMR = repmat(MM,5,1);


x = 1:size(MMR,2);
m_p = ceil(size(MMR,2)/2);


[XX,YY]=meshgrid(-4:0.2:4,-10+50:10+50);

XX = XX(:);
YY = YY(:);

cost = zeros(numel(XX),1);
for prm = 1:numel(XX)

    sl_t = XX(prm);
    off_t = YY(prm);

r_line = (x-m_p)*sl_t+off_t;

for y = 1:numel(x)
    
    t_bin = round(r_line(y)-10:r_line(y)+10);
    
    cost(prm) = cost(prm) + sum(MMR(t_bin,y).*(abs((-10:10)')));
    
    
    

end

end

[~,prm_m] = min(cost);

slope = XX(prm_m);
off_ph = YY(prm_m);


p_line = (x-m_p)*slope+off_ph;
% 
% figure(5001)
% imagesc(MMR)
% hold on
% plot(x,p_line)




