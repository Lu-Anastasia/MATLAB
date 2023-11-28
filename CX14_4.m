a=imread('o1.bmp');a=double(a(:,:,1));  %读取全息图
b=imread('o0.bmp');b=double(b(:,:,1));  %读取参考光强
a=a(:,18:779); b=b(:,18:779);           %去除图像边缘
figure,imshow(a,[]),figure,imshow(b,[])
IH=a-b;                                 %得到全息图
figure,imshow(IH,[])
[c,r]=size(IH);
zr=0.1096;                              %参考光到CCD的距离
zo=0.0355;                              %细胞到CCD的距离
Lox=r*10*10^(-6);Loy=c*10*10^(-6);      %全息图的尺寸,单位:米
lamda=6328*10^(-10); k=2*pi/lamda;      %波长和波数
xo=linspace(-Lox/2,Lox/2,r);yo=linspace(-Loy/2,Loy/2,c);%全息图的坐标
[xo,yo]=meshgrid(xo,yo);                %全息图的坐标网格
%生成球面离轴参考光
xr=0.0001;                              %参考点光源的x坐标
yr=0.0007;                              %参考点光源的y坐标
R=exp(j*k*zr).*exp(j*k.*((xo-xr).^2+(yo-yr).^2)/2/zr); %球面参考光波
M=5+eps                                 %赋值放大率
zp=M.*zo/(M-1)                          %球面照明光的半径
RF=exp(j*k*zp).*exp(j*k.*(xo.^2+yo.^2)/2/zp); %球面照明光波
zi=M.*zo                                %新的成像距离
F0=exp(j*k*zi)/(j*lamda*zo);
F1=exp(j*k/2/zi.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa2=fft2(IH.*conj(R).*RF); %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1;
Ui=F0.*fftshift(ifft2(Fuf2));           %在观察屏上的光场分布
Ii=Ui.*conj(Ui);                        %观察屏上的光强
figure,imshow(Ii,[]),colormap(pink),title('缩放的-1级重构像')