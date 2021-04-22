function mask=genMask4(radius,Nd,thickness,rhole,step)
if(~exist('step','var'))
    step=1;
end
tmp=zeros(2*radius+1,1);
tmp([flip(0:-step:-radius) 0:step:radius]+radius+1)=1;
stepmask=tmp*tmp';
[X,Y]=meshgrid(-radius:radius,-radius:radius);
mask=zeros(2*radius+1,2*radius+1,Nd);
for i=1:Nd
    mask(:,:,i)=(subfunc(X,Y,(i-1)*pi/Nd)<thickness).*stepmask;
end
hole=sum(mask,3)>1;
holeRange=radius+1+[-rhole:rhole];
hole(holeRange,holeRange)=1;
for i=1:Nd
    mask(:,:,i)=mask(:,:,i).*(1-hole);
end
hole=hole.*stepmask;
mask=cat(3,hole,mask);
end



function v=subfunc(x,y,theta)
theta=-theta;
v=sqrt(x.^2+y.^2-(x+tan(theta)*y).^2/(tan(theta)^2+1));
end