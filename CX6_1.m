Uo=imread('分辨率板_2.bmp');               %读入作为物的图像
Uo=double(Uo (:,:,1));                     %调取第一层，转换为双精度  RGB取一层
lamda=6328*10^(-10);k=2*pi/lamda;          %赋值波长,单位:米,波矢
D=0.01;                                    %赋值透镜的孔径,单位:米
f=0.4;                                     %赋值透镜的焦距,单位:米
figure,imshow(Uo,[])
[c,r]=size(Uo);                            %读取物面采样数
%下面计算物光传递到透镜的衍射过程
L0=0.005                                   %赋值物面的尺寸L0,单位:米
x0=linspace(-L0/2,L0/2,r);y0=linspace(-L0/2,L0/2,c); %赋值物面的坐标
[x0,y0]=meshgrid(x0,y0);
d1=1.2;                                    %物面到透镜的距离d1,单位:米
L=r*lamda*d1/L0                            %衍射光在透镜前表面上的尺寸L,单位:米 用的S-FFT上的观察屏尺寸公式
p=linspace(-L/2,L/2,r);q=linspace(-L/2,L/2,c); %赋值透镜前表面的坐标
[p,q]=meshgrid(p,q);
F00=exp(j*k*d1)/(j*lamda*d1)*exp(j*k/2/d1*(p.^2+q.^2));
Fpq=exp(j*k/2/d1*(x0.^2+y0.^2));
a= Uo.*Fpq;
FUpq=fft2(a);                              %做FFT变换
Ffpq=fftshift(FUpq);
Fufpq=F00.*Ffpq;                           %透镜前表面上的光场复振幅分布
I=Fufpq.*conj(Fufpq);                      %透镜前表面上的光强分布
figure,imshow(I,[]), colormap(pink),title('透镜上的光强分布')
%下面计算通过透镜后的光场
DD=round(D*r/L);                           %计算孔径对应的采样数
pxy=zeros(c,r);                            %生成孔径函数，这步只是建了个都是0的矩阵，下面这一大串都是，造了个圆形的透镜
for n=1:c
   for m=1:r
      if (n-c/2).^2+(m-r/2).^2<=(DD/2).^2;
      pxy(n,m)=1;
      end
   end
end
figure,imshow(pxy,[]);title('孔径函数')
Fufpqyp=Fufpq.*pxy.*exp(-j*k.*(p.^2+q.^2)/2/f); %计算通过透镜后的光场
%下面计算从透镜到观察面的衍射过程
d2=d1*f/(d1-f)                             %由物、像公式给出像距d2,单位:米
Lyp=r*lamda*d2/L,                          %给出观察面（像面）的尺寸,单位:米
x=linspace(-Lyp/2,Lyp/2,r);y=linspace(-Lyp/2,Lyp/2,c); %给出观察面的坐标
[x,y]=meshgrid(x,y);
F0=exp(j*k*d2)/(j*lamda*d2)*exp(j*k/2/d2*(x.^2+y.^2));
F=exp(j*k/2/d2*(p.^2+q.^2));
% 计算再现像
re_image=fft2(Fufpqyp.*F);re_image=re_image.*F0;
if Lyp<0                                   %成虚像时倒像
   re_image=flipud(re_image);re_image=fliplr(re_image);%左右、上下翻转
end
figure,imshow(re_image.*conj(re_image),[]),title('再现像')