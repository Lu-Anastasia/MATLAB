lamda=6328*(10^(-10));         %波长
k=2*pi/lamda;                  %波数
alpha=pi/2.005;                %光与X轴的夹角
beita=pi/2.005;                %光与Y轴的夹角
L=0.004                        %观察面的尺寸，单位：米
x=linspace(-L/2,L/2,512);y=x;
[x,y]=meshgrid(x,y);
U=exp(j.*k.*(x.*cos(alpha)+y.*cos(beita))); %构建入射平行光场
ph=k.*(x.*cos(alpha)+y.*cos(beita)); %直接计算实际相位
figure,surfl(ph),shading interp,colormap(gray)
phyp=angle(U);                %计算光场的相位（包裹相位）
figure,imshow(phyp,[])
figure,plot(ph(257,:) , '--')
hold on,plot(phyp(257,:),'r')
diff=U+1;                     %观察面上入射平行光与垂直照射平行光的干涉
I=diff.*conj(diff);           %观察面上的光强
figure,imshow(I,[])
UFuv=fftshift(fft2(U));       %计算光场的频谱
figure,imshow(abs(UFuv),[0,max(max(abs(UFuv)))./50]) %光场的频谱
IFuv=fftshift(fft2(I));       %计算干涉条纹的频谱
figure,imshow(abs(IFuv),[0,max(max(abs(IFuv)))./50]) %干涉条纹的频谱