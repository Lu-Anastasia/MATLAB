r=512,c=r;                       %给出物面上的抽样数
Uo=zeros(r,c);                   %预设物光场
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4)=2; %设置物光场振幅
Uo(r/2-round(r/6)+1:r/2+round(r/6),c/2-round(c/6)+1:c/2+round(c/6))=4; 
Uo(r/2-round(r/12)+1:r/2+round(r/12),c/2-round(c/12)+1:c/2+round(c/12))=6;
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4)=...
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4).*exp(j*peaks(256)*2); %设置物光的相位
lamda=6328*10^(-10);k=2*pi/lamda; %赋值波长、波数
zo=0.3086;                        %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                      %赋值衍射面(物)的尺寸,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成全息图记录过程的计算
xo=linspace(-Lo/2,Lo/2,c);yo=linspace(-Lo/2,Lo/2,r);
[xo,yo]=meshgrid(xo,yo);
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(Uo); fF1=fft2(F1);
Fuf=fa.*fF1; 
O=F0.*fftshift(ifft2(Fuf));      %在全息记录面上的光场分布
I=O.*conj(O);                    %全息记录面上物光强分布
figure,imshow(I,[]),colormap(pink),title('衍射图')
%下面加入平面参考光
alpha=pi/2.02;beita=pi/2.00;     %参考光与x轴、y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=O./mean(sqrt(I(:)))+R;     %调节光束比，并使参、物光干涉
IH=inter.*conj(inter);           %干涉得到全息图
figure,imshow(IH,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面可控放大率的方法重构再现像
M=1.8+eps                        %赋值放大倍数,加eps处理M=1的情况
zp=M.*zo/(M-1)                   %球面照明光的半径
%生成使像缩放的球面参考光
RF=exp(j*k*zp).*exp(j*k.*(xo.^2+yo.^2)/2/zp);  %球面照明光波
zi=M.*zo                         %新的成像距离
F0=exp(j*k*zi)/(j*lamda*zo);
F1=exp(j*k/2/zi.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa=fft2(IH.*conj(R).*RF);        %在全息图上乘R*,再以球面光垂直照射
Fuf=fa.*fF1; 
Ui=F0.*fftshift(ifft2(Fuf));     %在观察屏上的光场分布
Ii=Ui.*conj(Ui);                 %观察屏上的光强
figure,imshow(Ii,[]),colormap(gray),title('缩放的-1级重构像')
ph=angle(Ui);                    %计算重建相位
figure,imshow(ph,[]),title('重建相位')
figure,plot(ph(:,257))           %绘制相位剖线