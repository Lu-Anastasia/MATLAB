Uo=imread('xuehua.jpg');              %读入作为物的图像
Uo=double(Uo (:,:,1));                %调取第一层,转换为双精度
figure,imshow(Uo,[])
Uo=exp(j.*Uo/255);                  %生成纯相位物光场
[c,r]=size(Uo);                      %读取物面采样数
lamda=6328*10^(-10);k=2*pi/lamda;    %赋值波长,单位:米,波矢
f=0.004; Lo=0.001                   %赋值显微物镜的焦距,物面的尺寸,单位:米
D=0.00005                         %赋值滤波片直径,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算物光传递到透镜前表明的衍射过程(S-FFT)
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c); %赋值物面的坐标
[xo,yo]=meshgrid(xo,yo);             %生成物面的坐标网格
do=0.0041;                         %物面到透镜的距离do,单位:米
L=r*lamda*do/Lo                   %计算衍射光在透镜前表面上的尺寸,单位:米
xl=linspace(-L/2,L/2,r);yl=linspace(-L/2,L/2,c); %赋值透镜前表面的坐标
[xl,yl]=meshgrid(xl,yl);               %生成透镜前表面的坐标网格
F0=exp(j*k*do)/(j*lamda*do)*exp(j*k/2/do*(xl.^2+yl.^2));
F=exp(j*k/2/do*(xo.^2+yo.^2));
FU=(Lo*Lo/r/r).*fftshift(fft2(Uo.*F));
U1=F0.*FU;                        %透镜前表面上的光场复振幅分布
I1=U1.*conj(U1);                    %透镜前表面上的光强分布
figure,imshow(I1,[]), colormap(pink),title('透镜上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算通过透镜后的光场
U1yp=U1.*exp(-j*k.*(xl.^2+yl.^2)/2/f);  %计算通过透镜后的光场
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算通过透镜后的光场到达焦面的过程(S-FFT)
dlf=f,
Lyp=r*lamda*dlf/L,                  %给出焦面的再现尺寸,单位:米
xf=linspace(-Lyp/2,Lyp/2,r);yf=linspace(-Lyp/2,Lyp/2,c); %给出焦面的坐标
[xf,yf]=meshgrid(xf,yf);               %生成焦面的坐标网格
F0=exp(j*k*dlf)/(j*lamda*dlf)*exp(j*k/2/dlf*(xf.^2+yf.^2));
F=exp(j*k/2/dlf*(xl.^2+yl.^2));
Uf=(L*L/r/r).*fft2(U1yp.*F);Uf=Uf.*F0; %焦面上的光场分布
I2=Uf.*conj(Uf);                    %焦面上的光强分布
figure,imshow(log(I2),[]), colormap(pink),title('焦面上的光强分布')
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%生成滤波器
DD=round(D*r/Lyp);                 %赋值滤波器直径,单位:像素
H=zeros(c,r);                        %生成滤波器
for n=1:c
   for m=1:r
      if (n-c/2-1).^2+(m-r/2-1).^2<=(DD/2).^2;
      H(n,m)=1;
      end
   end
end
figure,imshow(H,[]);title('理想低通滤波器')
H_p=exp(j*pi/2).*H;                  %正相衬滤波器
H_n=exp(j*3*pi/2).*H;                %负相衬滤波器
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%计算经过“+”、“-”相衬后焦面处的光场
Uf_p=H_p.*Uf+(1-H).*Uf;              %正相衬后焦面处的光场
Uf_n=H_n.*Uf+(1-H).*Uf;              %负相衬后焦面处的光场
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算通过焦面后的光场到达观察面的过程(S-FFT)
dfi=do*f/(do-f)-f;
Li=r*lamda*dfi/Lyp,                   %给出像面的尺寸,单位:米
xi=linspace(-Li/2,Li/2,r);yi=linspace(-Li/2,Li/2,c); %给出像面的坐标
[xi,yi]=meshgrid(xi,yi);                 %生成像面的坐标网格
F0=exp(j*k*dfi)/(j*lamda*dfi)*exp(j*k/2/dfi*(xi.^2+yi.^2));
F=exp(j*k/2/dfi*(xf.^2+yf.^2));
% 计算再现像
Ui=(Lyp*Lyp/r/r).*fft2(Uf.*F);Ui=Ui.*F0;  
Ii=Ui.*conj(Ui);                       %无相衬时再现像的光强分布
figure,imshow(Ii,[]),title('无相衬时的再现像'),colormap(gray)
Ui_p=(Lyp*Lyp/r/r).*fft2(Uf_p.*F);Ui_p=Ui_p.*F0;  
Ii_p=Ui_p.*conj(Ui_p);                 %正相衬时再现像的光强分布
figure,imshow(Ii_p,[]),title('正相衬时的再现像'),colormap(gray)
Ui_n=(Lyp*Lyp/r/r).*fft2(Uf_n.*F);Ui_n=Ui_n.*F0;  
Ii_n=Ui_n.*conj(Ui_n);                 %负相衬时再现像的光强分布
figure,imshow(Ii_n,[]),title('负相衬时的再现像'),colormap(gray)
%= = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面绘制剖线比较相衬效果
Ii_p_yp=flipud(fliplr(Ii_p));             %对再现像作翻转
Ii_n_yp=flipud(fliplr(Ii_n));             %对再现像作翻转
figure,subplot(2,1,1),plot(1-2*angle(Uo(round(r/2),:)),'--k'),axis([1 r -1.1 pi]) 
hold on,plot(1+2*angle(Uo(round(r/2),:)),'-+r'),title('理论结果')
subplot(2,1,2),plot(Ii_n_yp(round(r/2)-2,91:477),'--k'),axis([1 386 0 2.5e-3])
hold on,plot(Ii_p_yp(round(r/2)-2,91:477),'-+r'),title('实验结果')