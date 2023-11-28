Uo=imread('JTC1.bmp');                    %读入待识别目标图像和参考图像
Uo=double(Uo (:,:,1));                    %调取第一层，转换为双精度
Uo(Uo<=100)=0;Uo(Uo~=0)=1;                %将物二值化
figure,imshow(Uo,[])
[c,r]=size(Uo);                           %读取物面采样数
lamda=6328*10^(-10);k=2*pi/lamda;         %赋值波长,单位:米,波矢
f=0.004; Lo=0.001                         %赋值透镜的焦距,物面的尺寸Lo,单位:米
D=0.00005                                 %赋值滤波片直径,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = 
%计算物光传递到透镜前表面的衍射过程(S-FFT)
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c); %赋值物面的坐标
[xo,yo]=meshgrid(xo,yo);                  %生成物面的坐标网格
do=0.0041;                                %物面到透镜的距离do,单位:米
L=r*lamda*do/Lo                           %衍射光在透镜前表面上的尺寸L,单位:米
xl=linspace(-L/2,L/2,r);yl=linspace(-L/2,L/2,c); %赋值透镜前表面的坐标
[xl,yl]=meshgrid(xl,yl);                  %生成透镜前表面的坐标网格
F0=exp(j*k*do)/(j*lamda*do)*exp(j*k/2/do*(xl.^2+yl.^2));
F=exp(j*k/2/do*(xo.^2+yo.^2));
FU=fftshift(fft2(Uo.*F));
U1=F0.*FU;                                %透镜前表面上的光场复振幅分布
I1=U1.*conj(U1);                          %透镜前表面上的光强分布
figure,imshow(I1,[]), colormap(pink),title('透镜上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算通过透镜后的光场
U1yp=U1.*exp(-j*k.*(xl.^2+yl.^2)/2/f);   %计算通过透镜后的光场
%= = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算通过透镜后的光场到达焦面的衍射过程(S-FFT)
dlf=f,
Lyp=r*lamda*dlf/L,                        %给出焦面的尺寸,单位:米
xf=linspace(-Lyp/2,Lyp/2,r);yf=linspace(-Lyp/2,Lyp/2,c); %给出焦面的坐标
[xf,yf]=meshgrid(xf,yf);                  %生成焦面的坐标网格
F0=exp(j*k*dlf)/(j*lamda*dlf)*exp(j*k/2/dlf*(xf.^2+yf.^2));
F=exp(j*k/2/dlf*(xl.^2+yl.^2));
Uf=fft2(U1yp.*F);Uf=Uf.*F0;               %焦面上的光场
I2=Uf.*conj(Uf);                          %焦面上的光强分布
figure,imshow(I2,[0,max(I2(:))/100]), colormap(pink),title('焦面上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = 
Uf=I2;                                    %光场经平方律转换器转换成功率谱分布
%= = = = = = = = = = = = = = = = = = = = = = = = 
%计算通过焦面后的光场到达像面的衍射过程(S-FFT)
dfi=do*f/(do-f)-f;
Li=r*lamda*dfi/Lyp,                      %给出像面的尺寸,单位:米
xi=linspace(-Li/2,Li/2,r);yi=linspace(-Li/2,Li/2,c); %给出像面的坐标
[xi,yi]=meshgrid(xi,yi);                 %生成像面的坐标网格
F0=exp(j*k*dfi)/(j*lamda*dfi)*exp(j*k/2/dfi*(xi.^2+yi.^2));
F=exp(j*k/2/dfi*(xf.^2+yf.^2));
Ui=fftshift(fft2(Uf.*F));Ui=Ui.*F0;      %像面上的再现光场
Ii=Ui.*conj(Ui);                         %像面上的光强分布
Ii=flipud(fliplr(Ii/max(Ii(:))));        %翻转图像并归一化
figure,imshow(Ii,[0,max(Ii(:))/100]),title('再现像光强分布'),colormap(pink)
[xi,yi]=ginput(1);                       %用交互式手动输入相关点坐标
figure,surfl(Ii),shading interp,colormap(gray)
figure,plot(Ii(round(yi),:)),title('剖线')