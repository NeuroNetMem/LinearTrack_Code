function [BaPos,BaPos_Q,CoPos,CoPos_Q,Ba_De_All,Co_De_All] = Decode_Position_1(Spikes,Reference,TiWi)
TiWi = TiWi * 0.05;

MM_Ref = Reference;
MM_Spi = Spikes; 

 AA_1 = 1-pdist2(MM_Ref,MM_Spi,'correlation');   
      
 
 ProbDe_1 = zeros(size(MM_Ref,1),size(MM_Spi,1))*NaN;
 ProbDe_Un_1 = zeros(size(MM_Ref,1),size(MM_Spi,1))*NaN;
 for tt = 1:size(MM_Spi,1)
     for bb = 1:size(MM_Ref,1)
        ProbF = 0;
        
        
        %ProbF = ProbF - sum(TiWi*MM_Ref(bb,MM_Spi(tt,:)==0),2);
        
        aa_f =  sum(MM_Spi(tt,:));
        if(aa_f>2)
        
        aa = find(MM_Spi(tt,:)>-1);
        for mm = aa(:)'    
        Rate = TiWi*MM_Ref(bb,mm);
        Obs = MM_Spi(tt,mm);
        ProbFact=(log((Rate)^(Obs)*exp(-Rate)/factorial(Obs)));
        ProbF = ProbF + ProbFact;
        end
        end
     ProbDe_1(bb,tt)=exp(ProbF);
     end
     ProbDe_Un_1(:,tt) = ProbDe_1(:,tt);
     ProbDe_1(:,tt) = ProbDe_1(:,tt)./sum(ProbDe_1(:,tt));
     
 end

 [BaPos_Q,BaPos] = max(ProbDe_1,[],1);
 [CoPos_Q,CoPos] = max(AA_1,[],1);
 
 Ba_De_All = ProbDe_Un_1;
 Co_De_All = AA_1;
 
 
 