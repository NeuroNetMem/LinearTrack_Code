th_all = [];
direct_n = {'UP','DOWN'};

theta_origin = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/Cluster_Info/';
destination_save = '~/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/'; 

for  A = [2 3 5]
 




for S = 2:3
th_sess = [];
for D = 1:2


load([theta_origin 'Cluster_GroupingsMPF_' direct_n{D} '_A' num2str(A) '_S0' num2str(S) '_F1.mat'   ])  



th_all = cat(1,th_all,theta_scores(:));
th_sess = cat(1,th_sess,theta_scores(:));
% load(['~/Dropbox/Projects_NIJ/Matteo/Nijmegen/Cluster_Info/' '/Cluster_GroupingsMPF_' direct_n{D} '_A' num2str(A) '_S' num2str(S) '_F2.mat'   ])  
% 
% th_all = cat(1,th_all,theta_scores(:));


end


th_limits_prc=prctile(th_sess,[0 33 66 100]);
th_limits_lin=linspace(min(th_sess),max(th_sess),4);

save(['./Theta_Limits/Theta_Limits_A' num2str(A) 'S' num2str(S) ],'th_limits_prc','th_limits_lin')

end
end


histogram(th_all,20)




th_limits_prc=prctile(th_all,[0 33 66 100]);
th_limits_lin=linspace(min(th_all),max(th_all),4);

save([destination_save 'Theta_Limits/Theta_Limits'],'th_limits_prc','th_limits_lin')