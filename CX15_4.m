I1=imread('b1.bmp');I1=double(I1(:,:,1)); %调入第一幅数字全息图
I2=imread('b4.bmp');I2=double(I2(:,:,1)); %调入第二幅数字全息图
IH=I1+I2;                                 %合成双曝光数字全息图
meanIH=filter2(ones(7,7),IH,'same')./49;  %计算双曝光全息图的平均值
meanI1=filter2(ones(7,7),I1,'same')./49;  %计算第一幅全息图的平均值
meanI2=filter2(ones(7,7),I2,'same')./49;  %计算第二幅全息图的平均值
IH=IH-meanIH;                             %减全息图平均值，抑制零级
I1=I1-meanI1;                             %减全息图平均值，抑制零级
I2=I2-meanI2;                             %减全息图平均值，抑制零级
figure,imshow(IH,[])
[c,r]=size(IH);
lamda=6328*10^(-10);k=2*pi/lamda;         %赋值波长、波数
Lox=r*3.2*10^(-6);Loy=c*3.2*10^(-6);      %给出全息图的尺寸,单位:米
%下面计算物光传递到观察屏的衍射成像过程
xo=linspace(-Lox/2,Lox/2,r);yo=linspace(-Loy/2,Loy/2,c);%全息图的坐标
[xo,yo]=meshgrid(xo,yo);                  %全息图的坐标网格
%构建离轴平面参考光，将再现像移至图像中心
alpha=pi/2.045;                           %与x轴的夹角
beita=pi/1.9655;                          %与y轴的夹角
R=exp(j.*k*(xo*cos(alpha)+yo*cos(beita))); %离轴平面参考光
zo=0.120;                                 %物到全息面的距离,单位:米
M=1/3+eps                                 %赋值再现像的放大率
%生成使像缩放的球面参考光
zp=M.*zo/(M-1)                            %计算球面照明光的半径
RF=exp(j*k*zp).*exp(j*k.*(xo.^2+yo.^2)/2/zp); %球面照明光波
zi=M.*zo                                  %新的成像距离
F0=exp(j*k*zi)/(j*lamda*zo);
F1=exp(j*k/2/zi.*(xo.^2+yo.^2));
fF1=fft2(F1);
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%先计算双曝光全息图的再现像
fa2=fft2(IH.*conj(R).*RF);               %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1; 
Ui=F0.*fftshift(ifft2(Fuf2));            %在观察屏上的光场分布
Ii=Ui.*conj(Ui);                         %观察屏上的光强分布
Ii=medfilt2(Ii,[5,5]);                   %中值滤波抑制噪声
figure,imshow(Ii,[0,max(max(Ii))/50]),colormap(pink),title('双曝光全息再现像')
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%计算第一幅全息图的再现像
fa2=fft2(I1.*conj(R).*RF);               %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1; 
Ui1=F0.*fftshift(ifft2(Fuf2));           %在观察屏上的光场分布
Ii1=Ui1.*conj(Ui1);                      %观察屏上的光强分布
Ii1=medfilt2(Ii1,[5,5]);                 %中值滤波抑制噪声
figure,imshow(Ii1,[0,max(max(Ii1))/50]),colormap(pink),title('第一幅全息图再现像')
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%计算第二幅全息图的再现像 
fa2=fft2(I2.*conj(R).*RF);              %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1; 
Ui2=F0.*fftshift(ifft2(Fuf2));          %在观察屏上的光场分布
Ii2=Ui2.*conj(Ui2);                     %观察屏上的光强分布
Ii2=medfilt2(Ii2,[5,5]);                %中值滤波抑制噪声
figure,imshow(Ii2,[0,max(max(Ii2))/50]),colormap(pink),title('第二幅全息图再现像')
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%计算两幅全息图再现光场间的干涉
U=Ui2./Ui1;
I=cos(angle(U));                       %两再现光场相互干涉
I=medfilt2(I,[5,5]);                   %中值滤波抑制噪声
figure,imshow(I,[]),colormap(gray),title('两全息再现光场形成的干涉条纹')