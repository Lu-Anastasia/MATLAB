Uo=imread('guang.bmp');          %调入作为物的图像
Uo=double(Uo(:,:,1));            %取第一层，并转为双精度
figure,imshow(Uo,[]),title('原始像')
[r,c]=size(Uo);
lamda=6328*10^(-10);k=2*pi/lamda;%赋值波长、波数
zo=0.3086;                       %物到全息记录面的距离,单位:米
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
I=O.*conj(O);                    %全息记录面上物光强分布
figure,imshow(I,[]),colormap(pink),title('衍射图')
%下面加入平面参考光
alpha=pi/2.02;                   %参考光与x轴间的夹角
beita=pi/2.00;                   %参考光与y轴间的夹角
R=exp(j*k*(xo*cos(alpha)+yo*cos(beita))); %参考光
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=O./max(max(sqrt(I)))+R;    %调节光束比，并使参、物光干涉
IH=inter.*conj(inter);           %干涉得到全息图
figure,imshow(IH,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面用T-FFT算法再现全息像
%先演示用不缩放的方法重构再现像
zi=zo;                           %全息图到观察面的距离
fa2=fft2(IH.*conj(R).*1);        %在全息图上乘R*,再以平行光垂直照射
Fuf2=fa2.*fF1; 
Ui1=F0.*fftshift(ifft2(Fuf2));   %在观察面上的光场分布
Ii1=Ui1.*conj(Ui1);              %观察屏上的光强
figure,imshow(Ii1,[]),colormap(gray),title('未缩放的-1级重构像')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%再演示可控放大率的方法重构再现像
M=1.5+eps                       %赋值放大倍数,加eps处理M=1的情况
zp=M.*zo/(M-1)                  %球面照明光的半径
%生成使像缩放的球面参考光
RF=exp(j*k*zp).*exp(j*k.*(xo.^2+yo.^2)/2/zp); %球面照明光波
zi=M.*zo                       %新的成像距离
F0=exp(j*k*zi)/(j*lamda*zo);
F1=exp(j*k/2/zi.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa2=fft2(IH.*conj(R).*RF);      %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1; 
Ui2=F0.*fftshift(ifft2(Fuf2));  %在观察屏上的光场分布
Ii2=Ui2.*conj(Ui2);             %观察屏上的光强
figure,imshow(Ii2,[]),colormap(gray),title('缩放的-1级重构像')