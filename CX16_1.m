r=512,c=r;                             %给出物面上的抽样数
Uo=zeros(r,c);                         %预设物光场
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4)=2;  %设置物光场振幅
Uo(r/2-round(r/6)+1:r/2+round(r/6),c/2-round(c/6)+1:c/2+round(c/6))=4; 
Uo(r/2-round(r/12)+1:r/2+round(r/12),c/2-round(c/12)+1:c/2+round(c/12))=6;
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4)=...
Uo(r/2-r/4+1:r/2+r/4,c/2-c/4+1:c/2+c/4).*exp(j*peaks(256)*2); %设置物光的相位
lamda=6328*10^(-10);k=2*pi/lamda;      %赋值波长,单位:米,波矢
D=0.06;f=0.4;                          %赋值透镜的孔径、焦距,单位:米
figure,imshow(Uo.*conj(Uo),[])         %显示物光光强分布
%下面计算物光传递到透镜的衍射过程
L0=0.005                               %赋值物面的尺寸L0,单位:米
x0=linspace(-L0/2,L0/2,c);y0=linspace(-L0/2,L0/2,r); %赋值物面的坐标
[x0,y0]=meshgrid(x0,y0);
d1=1.2;                                %物面到透镜的距离d1,单位:米
L=r*lamda*d1/L0                  %衍射光在透镜前表面上的尺寸L,单位:米
p=linspace(-L/2,L/2,c);q=linspace(-L/2,L/2,r); %赋值透镜前表面的坐标
[p,q]=meshgrid(p,q);
F00=exp(j*k*d1)/(j*lamda*d1)*exp(j*k/2/d1*(p.^2+q.^2));
Fpq=exp(j*k/2/d1*(x0.^2+y0.^2));
a= Uo.*Fpq;FUpq=fft2(a); Ffpq=fftshift(FUpq);
Fufpq=F00.*Ffpq;                      %透镜前表面上的光场复振幅分布
%下面计算通过透镜后的光场
DD=round(D*r/L);                      %计算孔径对应的采样数
pxy=zeros(r,c);                       %生成孔径函数
for n=1:r
   for m=1:c
      if (n-r/2).^2+(m-c/2).^2<=(DD/2).^2;
      pxy(n,m)=1;
      end
   end
end
Fufpqyp=Fufpq.*pxy.*exp(-j*k.*(p.^2+q.^2)/2/f); %计算通过透镜后的光场
%下面计算从透镜到像面的衍射过程
d2=d1*f/(d1-f)                        %由物、像公式给出像距d2,单位:米
Lyp=r*lamda*d2/L,                     %给出像面的尺寸,单位:米
x=linspace(-Lyp/2,Lyp/2,c);y=linspace(-Lyp/2,Lyp/2,r); %给出像面的坐标
[x,y]=meshgrid(x,y);
F0=exp(j*k*d2)/(j*lamda*d2)*exp(j*k/2/d2*(x.^2+y.^2));
F=exp(j*k/2/d2*(p.^2+q.^2));
re_image=fft2(Fufpqyp.*F);re_image=re_image.*F0;%重构再现物光场
figure,imshow(re_image.*conj(re_image),[]),title('再现像')
%下面生成再现物光场与参考光的干涉条纹（像面全息图）
meanA=mean(abs(re_image(:))).*4;     %为调整参、物光的光束比作准备
alpha=pi/2.00; beita=pi/2.02;        %参考光与x轴、y轴间的夹角
R=meanA.*exp(j*k*(x*cos(alpha)+y*cos(beita))); %准备生成四个参考光
R1=R;                         %参考光1
R2=R.*exp(j*pi/2);                   %参考光2
R3=R.*exp(j*2*pi/2);                 %参考光3
R4=R.*exp(j*3*pi/2);                 %参考光4
I1=(re_image+R1).*conj(re_image+R1); %计算像面全息图1
I2=(re_image+R2).*conj(re_image+R2); %计算像面全息图2
I3=(re_image+R3).*conj(re_image+R3); %计算像面全息图3
I4=(re_image+R4).*conj(re_image+R4); %计算像面全息图4
figure,imshow(I1,[])                 %显示像面全息图1
ph=atan((I4-I2)./(I1-I3));           %计算重建相位
figure,imshow(ph,[])                 %显示重建相位
figure,plot(ph(r/2-r/4+1:r/2+r/4,c/2)) %绘制相位剖线
FI=fft2(I1);                         %计算全息图1的频谱,准备用傅里叶变换法重建相位
figure,imshow(log(abs(FI)),[])       %显示全息图1的频谱
FIyp=zeros(r,c);                     %预设滤波后的频谱
FIyp(258:393,174:323)=FI(258:393,174:323); %滤波后的频谱
OR=ifft2(FIyp);                      %作逆傅里叶变换重构OR*
phyp=atan(imag(OR)./real(OR));       %计算重建相位
figure,imshow(phyp,[])               %显示重建相位
figure,plot(phyp (r/2-r/4+1:r/2+r/4,c/2)) %绘制相位剖线