r=512,c=r;                         %物面采样数
Uo=zeros(c,r);                      %预设物
d=30;a=10;                        %光栅常数和缝宽
for n=1:d:c                        %循环生成物(二维光栅)
   Uo(n:n+a,:)=1;
end
for m=1:d:r
   Uo (:,m:m+a)=1;
end
Uo=Uo(1:c,1:r);
figure,imshow(Uo,[])                 %显示物分布
lamda=6328*10^(-10);k=2*pi/lamda;    %赋值波长,单位:米,波矢
f=0.004; Lo=0.001                   %赋值透镜的焦距,物面的尺寸Lo,单位:米
D1=0.00005                        %赋值滤波片直径,单位:米
D2=0.00005                        %赋值滤波片宽度,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算物光传递到透镜的衍射过程(S-FFT)
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c); %赋值物面的坐标
[xo,yo]=meshgrid(xo,yo);             %生成物面的坐标网格
do=0.0041;                         %物面到透镜的距离do,单位:米
L=r*lamda*do/Lo                    %衍射光在透镜前表面上的尺寸L,单位:米
xl=linspace(-L/2,L/2,r);yl=linspace(-L/2,L/2,c); %赋值透镜前表面的坐标
[xl,yl]=meshgrid(xl,yl);               %生成透镜前表面的坐标网格
F0=exp(j*k*do)/(j*lamda*do)*exp(j*k/2/do*(xl.^2+yl.^2));
F=exp(j*k/2/do*(xo.^2+yo.^2));
FU=(Lo*Lo/r/r).*fftshift(fft2(Uo.*F)); 
U1=F0.*FU;                         %透镜前表面上的光场复振幅分布
I1=U1.*conj(U1);                     %透镜前表面上的光强分布
figure,imshow(I1,[]), colormap(pink),title('透镜上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算通过透镜后的光场
U1yp=U1.*exp(-j*k.*(xl.^2+yl.^2)/2/f);   %计算通过透镜后的光场
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算通过透镜后的光场到达后焦面的过程(S-FFT)
dlf=f,
Lyp=r*lamda*dlf/L,                   %给出后焦面的尺寸,单位:米
xf=linspace(-Lyp/2,Lyp/2,r);yf=linspace(-Lyp/2,Lyp/2,c); %给出后焦面的坐标
[xf,yf]=meshgrid(xf,yf);               %生成后焦面的坐标网格
F0=exp(j*k*dlf)/(j*lamda*dlf)*exp(j*k/2/dlf*(xf.^2+yf.^2));
F=exp(j*k/2/dlf*(xl.^2+yl.^2));
Uf=(L*L/r/r).*fft2(U1yp.*F);Uf=Uf.*F0;  %计算后焦面上的光场分布
I2=Uf.*conj(Uf);                     %后焦面上的光强分布
figure,imshow(I2,[0,max(I2(:))/100]), colormap(pink),title('后焦面上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面生成多边形滤波器
[cc,rr]=getline(gcf,'closed'); %在图中选择多边形，并将坐标赋值给cc和rr
H=zeros(c,r);                            %预设滤波器H
H=roipoly(H,cc,rr);                      %生成滤波器H,使选定区域内的值等于1
figure,imshow(H,[])
Ufyp=H.*Uf;
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算通过焦面后的光场到达观察面的过程
dfi=do*f/(do-f)-f;                     %满足物像公式
Li=r*lamda*dfi/Lyp,                      %给像面的尺寸,单位:米
xi=linspace(-Li/2,Li/2,r);yi=linspace(-Li/2,Li/2,c); %给出像面的坐标
[xi,yi]=meshgrid(xi,yi);             %生成像面的坐标网格
F0=exp(j*k*dfi)/(j*lamda*dfi)*exp(j*k/2/dfi*(xi.^2+yi.^2));
F=exp(j*k/2/dfi*(xf.^2+yf.^2));
% 计算再现光场及其光强分布
Ui=fft2(Ufyp.*F);Ui=Ui.*F0;
Ii=Ui.*conj(Ui);                          %滤波后像面上的光强分布
figure,imshow(Ii,[]),title('再现像'),colormap(gray)
figure,plot(Ii(round(c/2),:))