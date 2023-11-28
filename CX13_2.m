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
%再计算任意一次相移法的参考光
derta=pi/3                       %参考光的相移量δ,非2π整数倍
Ryp=R.*exp(j*derta);
%最后计算四次相移法的参考光
R1=R.*exp(j*2*pi*0/4);           %参考光1
R2=R.*exp(j*2*pi*1/4);           %参考光2
R3=R.*exp(j*2*pi*2/4);           %参考光3
R4=R.*exp(j*2*pi*3/4);           %参考光4
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面计算参、物光在全息记录面上的干涉,得到全息图
%先计算任意一次相移法的全息图
IH=(O+R).*conj(O+R);             %干涉得到全息图IH
IHyp=(O+Ryp).*conj(O+Ryp);       %干涉得到全息图IH'
figure,imshow(IH,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%再计算四次相移法的全息图
I1=(O+R1).*conj(O+R1);           %干涉得到全息图I1
I2=(O+R2).*conj(O+R2);           %干涉得到全息图I2
I3=(O+R3).*conj(O+R3);           %干涉得到全息图I3
I4=(O+R4).*conj(O+R4);           %干涉得到全息图I4
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面用S-FFT算法再现全息像
zi=0.20;                         %全息图到观察面的距离,单位:米
Li=r*lamda*zi/Lo                 %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,c);y=linspace(-Li/2,Li/2,r);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));   %用T-FFT算法得到的全息图尺寸与物面一致
% 取再现照明光垂直入射C=1
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 先计算直接消除法的再现像
holo1=Lo/r*Lo/c*fftshift(fft2((IH-I-r2xy).*F*1));holo1=holo1.*F0;
Ii1=holo1.*conj(holo1);
figure,imshow(Ii1,[0,max(max(Ii1))./1]),colormap(pink),title('直接消除法再现像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 再计算一次任意消除法的再现像
holo2=Lo/r*Lo/c*fftshift(fft2((IH-IHyp).*F*1));holo2=holo2.*F0;
Ii2=holo2.*conj(holo2);
figure,imshow(Ii2,[0,max(max(Ii2))./1]),colormap(pink),title('一次任意消除法再现像')
% 最后计算四次相移法的再现像
holo3=Lo/r*Lo/c*fftshift(fft2((I1.*R1+I2.*R2+I3.*R3+I4.*R4).*F*1));holo3=holo3.*F0;
Ii3=holo3.*conj(holo3);
figure,imshow(Ii3,[0,max(max(Ii3))./1]),colormap(pink),title('四次相移法再现像')
holo3yp=Lo/r*Lo/c*fftshift(fft2(conj(I1.*R1+I2.*R2+I3.*R3+I4.*R4).*F*1));holo3=holo3.*F0;
Ii3yp=holo3yp.*conj(holo3yp);
figure,imshow(Ii3yp,[0,max(max(Ii3yp))./1]),colormap(pink),title('四次相移法共轭光的再现像')