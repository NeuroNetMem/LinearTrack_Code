function [ext,alpha,Like,gauPlotVal,gauCentX,gauCentY,gauCentZ] = GLMPhaseThetaPowGamma(SpBinOcc,Pos,PosTot,ThetaPh,GammaPow,NGau,StGau,StGauPow,NBin,N_Gamma)

posx=Pos(:);
posy=ThetaPh(:);
posz=GammaPow(:);

posxT=PosTot(:);
posyT=[-pi pi];
poszT=[-1 1];


learn=[100 10 1 0.1 0.01];

NLearn=numel(learn);
LikeP=zeros(NLearn,1);


NUnits=size(SpBinOcc,2);

TimeLength=size(SpBinOcc,1)-1;


UnitSp=1*(SpBinOcc(1:end,:)>0);
UnitSp(UnitSp==0)=-1;

cost=repmat(rand(NUnits,1),1,TimeLength);


Gau = @(xd,yd,zd,six,siy,siz) exp(-(xd.^2)/(2*six^2)).*exp(-(yd.^2)/(siy^2)).*exp(-(zd.^2)/(siz^2));


gauCent1=linspace(min(posxT),max(posxT),NGau+1);
gauCent1=gauCent1(1:end-1)+(gauCent1(2)-gauCent1(1))/2;
gauCent2=linspace(min(posyT),max(posyT),11);
gauCent2=gauCent2(1:end-1)+(gauCent2(2)-gauCent2(1))/2;
gauCent3=linspace(min(poszT),max(poszT),N_Gamma+1);
gauCent3=gauCent3(1:end-1)+(gauCent3(2)-gauCent3(1))/2;

[gauCentXP,gauCentYP,gauCentZP]=meshgrid(gauCent1,gauCent2,gauCent3);

gauCentX=gauCentXP(:);
gauCentY=gauCentYP(:);
gauCentZ=gauCentZP(:);

gauDistX=repmat(gauCentX,1,numel(posx)-1)-repmat(posx(1:end-1)',numel(gauCentX),1);
gauDistY=abs(repmat(gauCentY,1,numel(posy)-1)-repmat(posy(1:end-1)',numel(gauCentY),1));
gauDistY(gauDistY>pi)=2*pi-gauDistY(gauDistY>pi);

gauDistZ=repmat(gauCentZ,1,numel(posz)-1)-repmat(posz(1:end-1)',numel(gauCentZ),1);

gauVal=Gau(gauDistX,gauDistY,gauDistZ,StGau,2*pi/6,StGauPow);

alpha=zeros(numel(gauCentX),NUnits);



PlotCent1=linspace(min(posxT),max(posxT),NBin+1);
PlotCent1=PlotCent1(1:end-1)+(PlotCent1(2)-PlotCent1(1))/2;
PlotCent2=linspace(min(posyT),max(posyT),21);
PlotCent2=PlotCent2(1:end-1)+(PlotCent2(2)-PlotCent2(1))/2;
PlotCent3=linspace(min(poszT),max(poszT),N_Gamma+1);
PlotCent3=PlotCent3(1:end-1)+(PlotCent3(2)-PlotCent3(1))/2;

[PlotCentXP,PlotCentYP,PlotCentZP]=meshgrid(PlotCent1,PlotCent2,PlotCent3);

PlotCentX=PlotCentXP(:);
PlotCentY=PlotCentYP(:);
PlotCentZ=PlotCentZP(:);

gauPlotDistX=repmat(gauCentX,1,numel(PlotCentX))-repmat(PlotCentX,1,numel(gauCentX))';
gauPlotDistY=abs(repmat(gauCentY,1,numel(PlotCentY))-repmat(PlotCentY,1,numel(gauCentY))');
gauPlotDistY(gauPlotDistY>pi)=2*pi-gauPlotDistY(gauPlotDistY>pi);

gauPlotDistZ=repmat(gauCentZ,1,numel(PlotCentZ))-repmat(PlotCentZ,1,numel(gauCentZ))';

gauPlotVal=Gau(gauPlotDistX,gauPlotDistY,gauPlotDistZ,StGau,2*pi/6,StGauPow);






%Decision=1;
l=0;

M=1;
lear=learn;
while(l<500)
  Mp=M;
    %*Adapt;
    
    dAlphaLike=(gauVal*UnitSp(2:end,:)-(tanh(+alpha'*gauVal+cost)*gauVal')')./TimeLength;
    

    dCostLike=sum(UnitSp(2:end,:)'-tanh(+alpha'*gauVal+cost),2)./TimeLength;

for k=1:NLearn

alphaN=alpha+lear(k)*dAlphaLike;

costN=cost+repmat(lear(k)*dCostLike,1,TimeLength);

LikeP(k)=trace((+alphaN'*gauVal+costN)*UnitSp(2:end,:))-sum(sum(log(2*cosh(+alphaN'*gauVal+costN)),2));



end

[LL,M]=max(LikeP);

if(l>20)
alphaN=alpha+lear(M)*dAlphaLike;
end



costN=cost+repmat(lear(M)*dCostLike,1,TimeLength);




if(l>20)
alpha=alphaN;
end
cost=costN;

l=l+1;




if(mod(l,1)==0)

learp2=linspace(lear(M)/4,lear(M)*1.5,NLearn);
lear=learp2;
disp(lear);
disp([exp(LL/(NUnits*TimeLength)) l]);


end



end


Like=trace((+alpha'*gauVal+cost)*UnitSp(2:end,:))-sum(sum(log(2*cosh(+alpha'*gauVal+cost)),2));
ext=cost(:,1);

% 
% PlotVal=reshape(alpha'*gauPlotVal+cost(1),100,100);
% figure(3)
% imagesc(PlotVal);
end