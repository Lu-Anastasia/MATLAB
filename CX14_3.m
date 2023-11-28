IH=imread('ho.bmp');IH=double(IH(:,:,2));%调入数字全息图
meanIH=filter2(ones(7,7),IH,'same')./49; %计算全息图平均值
IH=IH-meanIH;                            %减全息图平均值，抑制零级
figure,imshow(IH,[])
[c,r]=size(IH);
lamda=5328*10^(-10);k=2*pi/lamda;        %赋值波长、波数
Lox=r*3.2*10^(-6);Loy=c*3.2*10^(-6);     %给出全息图的尺寸,单位:米
%计算物光传递到观察屏的衍射过程
xo=linspace(-Lox/2,Lox/2,r);yo=linspace(-Loy/2,Loy/2,c);%全息图的坐标
[xo,yo]=meshgrid(xo,yo);                %全息图的坐标网格
%构建离轴平面参考光
alpha=pi/1.962;                         %与x轴的夹角
beita=pi/1.9855;                        %与y轴的夹角
R=exp(j.*k*(xo*cos(alpha)+yo*cos(beita))); %离轴平面参考光
zo=1.7000;                              %物到全息面的距离,单位:米
M=1/14+eps                              %赋值再现像的放大率
%生成使像缩放的球面参考光
zp=M.*zo/(M-1)                          %计算球面照明光的半径
RF=exp(j*k*zp).*exp(j*k.*(xo.^2+yo.^2)/2/zp); %球面照明光波
zi=M.*zo                                %新的成像距离
F0=exp(j*k*zi)/(j*lamda*zo);
F1=exp(j*k/2/zi.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa2=fft2(IH.*conj(R).*RF);              %在全息图上乘R*,再以球面光垂直照射
Fuf2=fa2.*fF1; 
Ui=F0.*fftshift(ifft2(Fuf2));           %在观察屏上的光场分布
Ii=Ui.*conj(Ui);                        %观察屏上的光强
figure,imshow(Ii,[0,max(max(Ii))/1000]),colormap(gray),title('缩放的-1级重构像')