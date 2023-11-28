Uo=imread('guang.bmp');          %调入作为物的图像
Uo=double(Uo(:,:,1));            %取第一层，并转为双精度
figure(1),imshow(Uo,[]),title('原始像')
[r,c]=size(Uo);
lamda=6328*10^(-10);k=2*pi/lamda;%赋值波长、波数
zo=0.20;                         %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                     %赋值衍射面(物)的尺寸,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成全息图记录过程的计算
xo=linspace(-Lo/2,Lo/2,c);yo=linspace(-Lo/2,Lo/2,r);
[xo,yo]=meshgrid(xo,yo);
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(Uo); fF1=fft2(F1);
Fuf=fa.*fF1; 
O=F0.*fftshift(ifft2(Fuf));      %在全息记录面上的光场分布
I=O.*conj(O);                    %全息记录面上的光强分布
figure(2),imshow(I,[]),colormap(pink),title('衍射图')
%下面加入参考光
alpha=pi/2.00;                   %参考光与x轴间的夹角
beita=pi/2.02;                   %参考光与y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=O./max(max(sqrt(I)))+R;    %调节光束比，并使参、物光干涉
II=inter.*conj(inter);           %干涉得到全息图
figure(3),imshow(II,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%对全息图作处理
%先在空域中计算平均值
s1=5                             %设定滤波窗口的大小(像素)
B=ones(s1,s1);                   %设置用于均值滤波的mark
meanII=filter2(B,II)/s1/s1;      %用均值滤波计算全息图在(x,y)附近的平均值
II_meanII=II-meanII;             %全息图与平均值相减,只留下干涉项的近似值
figure(4),subplot(4,1,1),plot(II(:,257)),title('全息图') %绘制全息图的剖线
subplot(4,1,2),plot(meanII(:,257)),title('空域平均') %绘制全息图平均值的剖线
subplot(4,1,4), plot(I(:,257)),title('衍射光强') %绘制衍射光强的剖线
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%再在频域中作高通滤波
s2=25                             %设定高通滤波器的窗口大小(像素)
H=ones(r,c);                      %预设高通滤波器
H(round(r/2)+1-s2:round(r/2)+s2,round(c/2)+1-s2:round(c/2)+s2)=0;%理想高通滤波器
IIFF=fftshift(fft2(II));          %计算全息图的频谱
HFF=H.*IIFF;                      %滤波
figure(5),imshow(abs(IIFF),[0,100])%显示全息图的频谱
figure(6),imshow(H,[])            %显示理想高通滤波器
IIyp=ifft2(HFF);                  %计算滤波后的全息图
figure(4),subplot(4,1,3),plot(abs(IIyp(:,257))),title('频域滤波') %绘制全息图的剖线
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面用S-FFT算法再现全息像
zi=0.20;                           %全息图到观察面的距离,单位:米
Li=r*lamda*zi/Lo                   %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,c);y=linspace(-Li/2,Li/2,r);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));     %用T-FFT算法得到的全息图尺寸与物面一致
% 先计算未处理的再现像
% 取再现照明光垂直入射C=1
holo1=Lo/r*Lo/c*fftshift(fft2(II.*F*1));holo1=holo1.*F0;
Ii1=holo1.*conj(holo1);
figure,imshow(Ii1,[0,max(max(Ii1))./1]),colormap(pink),title('未处理的再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 再计算空域作平均值处理后的再现像
holo2=Lo/r*Lo/c*fftshift(fft2(II_meanII.*F*1));holo2=holo2.*F0;
Ii2=holo2.*conj(holo2);
figure,imshow(Ii2,[0,max(max(Ii2))./1]),colormap(pink),title('空域平均再现像')
% 最后计算频域滤波后的再现像
holo3=Lo/r*Lo/c*fft2(IIyp.*F*1);holo3=holo3.*F0;
Ii3=holo3.*conj(holo3);
figure,imshow(Ii3,[0,max(max(Ii3))./1]),colormap(pink),title('频域滤波再现像')