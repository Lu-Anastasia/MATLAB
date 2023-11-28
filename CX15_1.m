%构建物体
r=512,c=r;                       %给出物面上的抽样数
Uo=zeros(c,r);                   %预设平整物面
Uo(c/2-c/4:c/2+c/4, r/2-r/4:r/2+r/4)=1; %生成矩形（变形前的物体）
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%构建变形后的物体
Uoyp=Uo.*exp(j.*peaks(r).*2);    %在物体上叠加相位，模拟变形后的物体
lamda=6328*10^(-10);k=2*pi/lamda;%赋值波长、波数
zo=0.3086;                       %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                     %赋值衍射面(物)的尺寸,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成物体变形前全息图记录过程的计算
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c);
[xo,yo]=meshgrid(xo,yo);
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(Uo); fF1=fft2(F1);
Fuf=fa.*fF1; 
O=F0.*fftshift(ifft2(Fuf));      %在全息记录面上的物光场分布
I=O.*conj(O);                    %全息记录面上的光强分布
O=O./max(max(sqrt(I)));          %调节光束比
figure,imshow(I,[]),colormap(pink),title('变形前物光的衍射图')
%下面加入参考光
alpha=pi/2.00;                   %参考光与x轴间的夹角
beita=pi/2.0175;                 %参考光与y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=O+R;                       %参、物光干涉
IH= inter.*conj(inter);          %干涉得到的全息图
figure,imshow(IH,[]),title('变形前物光的全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%计算物体变形后在全息记录面上的物光场分布
fayp=fft2(Uoyp);
Fufyp=fayp.*fF1; 
Oyp=F0.*fftshift(ifft2(Fufyp));  %在全息记录面上的物光场分布
Iyp=Oyp.*conj(Oyp);              %全息记录面上的衍射光强分布
figure,imshow(Iyp,[]),colormap(pink),title('变形物光的衍射图')
Oyp=Oyp./max(max(sqrt(Iyp)));    %调节光束比
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用透镜辅助成像
f=0.10                           %透镜焦距，单位：米
zi=zo*f/(zo-f)                   %成像相距
lens=exp(-j*k.*(xo.^2+yo.^2)/2/f);%透镜的相位变换系数
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算再现
%用S-FFT算法再现全息像
Li=r*lamda*zi/Lo                 %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,r);y=linspace(-Li/2,Li/2,c);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));  %用T-FFT算法得到的全息图尺寸与物面一致
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 先取再现照明光为参考光R
UR=Lo/r*Lo/c*fftshift(fft2(IH.*F.*R.*lens));UR=UR.*F0; %再现光场
Ii=UR.*conj(UR);
figure,imshow(Ii,[0,max(max(Ii))./1]),colormap(pink),title('参考光照明得到的再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 再取再现照明光为变形后的物光
UO=Lo/r*Lo/c*fftshift(fft2(IH.*F.*Oyp.*lens));UO=UO.*F0; %再现光场
Iiyp=UO.*conj(UO);
figure,imshow(Iiyp,[0,max(max(Iiyp))./1]),colormap(pink),title('变形物光照明得到的再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% 最后计算两再现光场的干涉检测
II=(UR+UO).*conj(UR+UO);
figure,imshow(II,[0,max(max(II))]),colormap(pink),title('同时照明得到的干涉条纹')