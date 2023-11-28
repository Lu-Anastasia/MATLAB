r=512,c=r;                            %给出衍射面上的抽样数
a=zeros(r,c);                         %预设衍射面
a(r/2-r/4:r/2+r/4,c/2-c/4:c/2+c/4)=1; %生成衍射孔
lamda=6328*10^(-10);k=2*pi/lamda;     %赋值波长、波数
L0=5*10^(-3);                         %赋值衍射面尺寸，单位:米
d=0.1;                                %赋值观察屏到衍射面的距离,单位:米 
x0=linspace(-L0/2,L0/2,c);            %生成衍射面x轴坐标
y0=linspace(-L0/2,L0/2,r);            %生成衍射面y轴坐标
[x0,y0]=meshgrid(x0,y0);              %生成衍射面的二维坐标网格
%下面开始用式（4-13）计算衍射积分
F0=exp(j*k*d)/(j*lamda*d);            %赋值exp(ikd)/(iλd)
F1=exp(j*k/2/d*(x0.^2+y0.^2));        %赋值exp[ik(x02+y02)/2d]
fa=fft2(a);                           %完成光场U0(x0,y0)的傅里叶变换
fF1=fft2(F1);                         %完成exp[ik(x02+y02)/2d]的傅里叶变换
Fuf=fa.*fF1;                          %完成频谱相乘
U=F0.*fftshift(ifft2(Fuf));           %得到观察屏上的光场分布U(x,y)
I=U.*conj(U);                         %计算观察屏上的光强分布
figure, imshow(I,[0,max(max(I))]),colormap(gray)
figure,imshow(max(max(I))-I,[]),colormap(gray)