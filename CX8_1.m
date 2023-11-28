Uo=imread('分辨率板_2.bmp');               %读入作为物的图像
Uo=double(Uo(:,:,1));                      %调取第一层，转换为双精度
[c,r]=size(Uo);                            %读取物面采样数
lamda=6328*10^(-10);k=2*pi/lamda;          %赋值波长,单位:米,波矢
D=0.01;                                    %赋值透镜的孔径,单位:米
f=0.4;                                     %赋值透镜的焦距,单位:米
figure,imshow(Uo,[])
Lo=0.005;                                  %物面尺寸,单位:米
do=1.2;                                    %物到透镜的距离,亦即物瞳距,单位:米,可以改变
di=do*f/(do-f);                            %透镜到观察屏的距离,单位:米
cf=D/2/lamda/di;                           %截止频率(瑞利观点)
Li=Lo*di/do;                               %像面尺寸,单位:米
kethi=linspace(-1./2./Li,1./2./Li,r).*r;nenta=linspace(-1./2./Li,1./2./Li,c).*c;
[kethi,nenta]=meshgrid(kethi,nenta);       %像的二维频谱网格
H=zeros(c,r);                              %预设传相干递函数
for n=1:c
   for m=1:r
      if kethi(n,m).^2+nenta(n,m).^2<=cf.^2;
      H(n,m)=1;
      end
   end
end
figure,surfl(H),shading interp,colormap(gray);title('相干传递函数CTF ')
Gg=fftshift(fft2(Uo));                     %理想像的频谱
Gic=Gg.*H;                                 %相干照明下像的频谱
Uic=ifft2(Gic);                            %相干照明下像的光场分布
Iic=Uic.*conj(Uic);                        %相干照明下像的光强分布
figure,imshow(Iic,[]),title('相干照明下像的光强分布')
%下面开始完成非相干照明下的成像计算
%先用算法1计算OTF
h=fftshift(fft2(H));                       %用(8-30)计算相干照明下脉冲相应
HH=abs(fftshift(fft2(h.*conj(h))));        %计算CTF的自相关运算
OTF1=HH./max(max(HH));                     %对自相关运算结果作归一化
figure,surfl(OTF1),shading interp,colormap(gray);title('算法1得到的光学传递函数')
%再用算法2计算OTF
[phai,rou]=cart2pol(kethi,nenta);          %将二维频域从直角坐标转为极坐标
OTF2=zeros(c,r);                           %预设光学传递函数
for n=1:c                                  %用式(8-31)循环赋值光学传递函数
   for m=1:r
      if rou(n,m)<=2.*cf
      OTF2(n,m)=2*(acos(rou(n,m)/2/cf)-rou(n,m)/2/cf.*sqrt(1-(rou(n,m)/2/cf).^2))/pi;
      end
   end
end
figure,surfl(OTF2),shading interp,colormap(gray);title('算法2得到的光学传递函数')
%下面给出两种算法的成像结果
Gii1=Gg.*OTF1;                             %非相干照明下像的频谱（用OTF1）
Iii1=abs(ifft2(Gii1));                     %非相干照明下的成像（用OTF1）
figure,imshow(Iii1,[]),title('用OTF1得到的非相干照明成像'),colormap(gray)
Gii2=Gg.*OTF2;                             %非相干照明下像的频谱（用OTF2）
Iii2=abs(ifft2(Gii2));                     %非相干照明下的成像（用OTF2）
figure,imshow(Iii2,[]),title('用OTF2得到的非相干照明成像'),colormap(gray)