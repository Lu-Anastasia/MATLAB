%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%构建物面
r=512,c=r;                           %给出物面上的抽样数
Uo=zeros(c,r);                       %预设平整物面
Uo(c/2-c/4+1:c/2+c/4, r/2-r/4+1:r/2+r/4)=1; %生成矩形平面（振动前的物体）
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成全息图记录过程的计算
lamda=6328*10^(-10);k=2*pi/lamda;    %赋值波长、波数
zo=0.3086;                           %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                         %赋值衍射面(物)的尺寸,单位:米
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c);
[xo,yo]=meshgrid(xo,yo);
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(Uo);fF1=fft2(F1);
Fuf=fa.*fF1; 
O=F0.*fftshift(ifft2(Fuf));           %在全息记录面上的物光场分布
I=O.*conj(O);                         %全息记录面上的光强分布
O=O./max(max(sqrt(I)));               %调节光束比
figure,imshow(I,[]),colormap(pink),title('振动前物光的衍射图')
%下面加入参考光
alpha=pi/2.00;                        %参考光与x轴间的夹角
beita=pi/2.0175;                      %参考光与y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=O+R;                            %参、物光干涉
I=inter.*conj(inter);                 %干涉得到全息图
figure,imshow(I,[]),title('振动前的全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%计算物体振动后在全息记录面上的物光场分布
%构建振动物体的形变
thita1=pi/2;thita2=pi/3;               %入射角和反射角
amplitude=zeros(c,r);                  %预设振幅
amplitude(c/2-c/4+1:c/2+c/4,r/2-r/4+1:r/2+r/4)=abs(peaks(256)-min(min(peaks(256))));%用peaks函数赋值振幅
amplitude=amplitude./max(max(amplitude)).*5.*lamda; %将振幅的最大值调整到波长的5倍
phai0=pi;omiga=5;n=50;   %赋值振动的初相、角频率，一个周期中的抽样次数
t=linspace(0,2.*pi/omiga,n);            %赋值抽样时间
Icontinue=zeros(c,r);                   %预设连续曝光全息图
for ii=1:n
Uoyp=Uo.*exp(j.*k.*amplitude.*cos(omiga.*t(ii)+phai0).*(cos(thita1)+cos(thita2)));%振动后的物光波
%用T-FFT算法完成衍射计算，得到全息图
fayp=fft2(Uoyp);Fufyp=fayp.*fF1; 
Oyp=F0.*fftshift(ifft2(Fufyp));         %在全息记录面上的物光场分布
Iyp=Oyp.*conj(Oyp);                     %全息记录面上的光强分布
Oyp=Oyp./max(max(sqrt(Iyp)));           %调节光束比
interyp=Oyp+R;                          %参、物光干涉
Iyp=interyp.*conj(interyp);             %干涉得到全息图
Icontinue=Icontinue+Iyp;                %连续(多次)曝光
end
figure,imshow(Icontinue,[]),title('时间平均全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面用S-FFT算法再现全息像
zi=0.3086                              %给出像面的位置,单位:米
Li=r*lamda*zi/Lo                       %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,r);y=linspace(-Li/2,Li/2,c);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));  %用T-FFT算法得到的全息图尺寸与物面一致
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 取再现照明光为参考光R，计算没有振动物体的再现像
holo1=Lo/r*Lo/c*fftshift(fft2(I.*F.*R)); holo1=holo1.*F0;
Ii1=holo1.*conj(holo1);
figure,imshow(Ii1,[0,max(max(Ii1))./1]),colormap(pink),title('没有振动物体的再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 取再现照明光为参考光R，计算时间平均全息的再现像
holo2=Lo/r*Lo/c*fftshift(fft2(Icontinue.*F.*R)); holo2=holo2.*F0;
Ii2=holo2.*conj(holo2);
figure,imshow(Ii2,[0,max(max(Ii2))./1]),colormap(pink),title('时间平均全息的再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 下面计算贝塞尔函数对再现像调制分布图
a=k.*amplitude(c/2-c/4+1:c/2+c/4, r/2-r/4+1:r/2+r/4).*(cos(thita1)+cos(thita2));
J=besselj(0,a);
J2=J.*J;
figure,imshow(J2,[]),colormap(pink),title('贝塞尔函数调制分布')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =