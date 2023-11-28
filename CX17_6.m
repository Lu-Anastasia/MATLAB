Uo=imread('guang.bmp');           %调入作为物的图像
Uo=double(Uo(:,:,1));             %取第一层，并转为双精度
figure,imshow(Uo,[]),title('原始像')
[r,c]=size(Uo);
lamda=6328*10^(-10);k=2*pi/lamda; %赋值波长、波数
zo=0.3086;                        %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                      %赋值衍射面(物)的尺寸,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成衍射（物）光场的计算
xo=linspace(-Lo/2,Lo/2,c);yo=linspace(-Lo/2,Lo/2,r);
[xo,yo]=meshgrid(xo,yo);
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(Uo); fF1=fft2(F1);
Fuf=fa.*fF1; 
O=F0.*fftshift(ifft2(Fuf));       %在全息记录面上的光场分布
alpha=pi/2.025;                   %参考光的倾斜角度
beita=0;
R=exp(j*k*(xo*cos(alpha)));       %参考光场
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%绘制计算干涉全息图
U=O./max(abs(O(:)));
inter=(U+R).*conj(U+R);           %连续型计算干涉全息图
figure,imshow(inter,[]),colormap(gray),title('连续型全息图')
qq=asin(abs(U))./pi;
out=inter-cos(qq.*pi);
figure,plot(inter(257, 1:258))
hold on,plot(qq(257, 1:258),'k')
hold on,plot(out(257,1:258),'r')
CGH=out;
CGH(CGH>=0)=1;                    %二值化计算干涉全息图
CGH(CGH~=1)=0;
figure,imshow(CGH,[]),colormap(pink),title('二值化全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = =
%用S-FFT算法完成全息图再现全息像
zi=0.3086;                        %全息图到观察面的距离,单位:米
Li=r*lamda*zi/Lo                  %给出像面的尺寸,单位:米
x=linspace(-Li/2,Li/2,c);y=linspace(-Li/2,Li/2,r);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2));    %用T-FFT算法得到的全息图尺寸与物面一致
holo1=Lo/r*Lo/c*fftshift(fft2(inter.*F*1)); holo1=holo1.*F0; % 再现照明光垂直入射
Ii1=holo1.*conj(holo1);           %计算再现像的光强
figure,imshow(Ii1,[0,max(max(Ii1))./2]),title('连续型CGH再现像')
holo2=Lo/r*Lo/c*fftshift(fft2(CGH.*F*1)); holo2=holo2.*F0; % 再现照明光垂直入射
Ii2=holo2.*conj(holo2);          %计算再现像的光强
figure,imshow(Ii2,[0,max(max(Ii2))./10]),title('二值化CGH再现像')