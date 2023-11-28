lamda=6328e-10;                      %波长，单位：米
k=2*pi/lamda;                        %波数
x0=0.001;                             %点光源的x坐标，单位：米
y0=0.001;                            %点光源的y坐标，单位：米
z=0.3;                            %观察面到点光源的垂直距离，单位：米
L=0.005;                              %观察面的尺寸，单位：米
x=linspace(-L/2,L/2,512);y=x;        %构建x坐标和y坐标
[x,y]=meshgrid(x,y);                 %构建二维坐标网格
U1=exp(j*k*z).*exp(j*k.*((x-x0).^2+(y-y0).^2)/2/z);  %发散球面光波
ph1=k.*((x-x0).^2+(y-y0).^2)/2/z;    %发散球面波的实际相位
figure,surfl(ph1),shading interp,colormap(gray) 
phyp1=angle(U1);                     %发散球面波的包裹相位
figure,imshow(phyp1,[])
U2=exp(-j*k*z).*exp(-j*k.*((x-x0).^2+(y-y0).^2)/2/z); %会聚球面光波
ph2=-k.*((x-x0).^2+(y-y0).^2)/2/z;   %会聚球面波的实际相位
figure,surfl(ph2),shading interp,colormap(gray) 
phyp2=angle(U2);                     %会聚球面波的包裹相位
figure,imshow(phyp2,[]) 
figure, plot(ph2(257,:),'--')        %实际相位的剖线
hold on                              %保持当前图像
plot(phyp2(257,:),'r')               %包裹相位的剖线
diff1=U1+1;                          %观察面上发散球面光与垂直照射平行光的干涉
I1=diff1.*conj(diff1);               %观察面上的光强
figure,imshow(I1,[0,max(max(I1))])
diff2=U2+1;                          %观察面上会聚球面光与垂直照射平行光的干涉
I2=diff2.*conj(diff2);               %观察面上的光强
figure,imshow(I2,[0,max(max(I2))])