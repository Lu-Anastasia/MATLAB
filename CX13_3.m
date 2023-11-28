Uo=imread('guang.bmp');          %调入作为物的图像
Uo=double(Uo(:,:,1));            %取第一层，并转为双精度
figure,imshow(Uo,[]),title('原始像')
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
O=F0.*fftshift(ifft2(Fuf));      %在全息记录面上的物光场分布
O=O./max(max(abs(O)));           %将物光的振幅调整到1
I=O.*conj(O);                    %全息记录面上的光强分布
figure,imshow(I,[]),colormap(pink),title('衍射光强')
%下面加入参考光
alpha=pi/2.00;                   %参考光与x轴间的夹角
beita=pi/2.02;                   %参考光与y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
r2xy=R.*conj(R);                 %参考光强
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算参、物光在全息记录面上的干涉,得到全息图
IH=(O+R).*conj(O+R);             %干涉得到全息图IH
figure,imshow(IH,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面用S-FFT算法再现全息像
zi=0.20;                         %全息图到观察面的距离,单位:米
Li=r*lamda*zi/Lo                 %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,c);y=linspace(-Li/2,Li/2,r);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));   %用T-FFT算法得到的全息图尺寸与物面一致
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 先计算不改变成像位置的情况
% 取再现照明光为垂直入射平行光C=1
holo1=Lo/r*Lo/c*fftshift(fft2(IH.*F.*1));holo1=holo1.*F0;
Ii1=holo1.*conj(holo1);
figure,imshow(Ii1,[0,max(max(Ii1))./1]),colormap(pink),title('照明光为垂直入射平行光')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 再计算将共轭像移动到观察屏中心
% 取再现照明光为参考光的共轭光C=R*
holo2=Lo/r*Lo/c*fftshift(fft2(IH.*F.*conj(R)));holo2=holo2.*F0;
Ii2=holo2.*conj(holo2);
figure,imshow(Ii2,[0,max(max(Ii2))./1]),colormap(pink),title('照明光为参考光的共轭光')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 最后计算将原始像移动到观察屏中心
% 取再现照明光为参考光C=R
holo3=Lo/r*Lo/c*fftshift(fft2(IH.*F.*R));holo3=holo3.*F0;
Ii3=holo3.*conj(holo3);
figure,imshow(Ii3,[0,max(max(Ii3))./1]),colormap(pink),title('照明光为参考光')