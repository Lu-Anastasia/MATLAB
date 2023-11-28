Uo=imread('分辨率板_2.bmp');               %读入作为物的图像
Uo=double(Uo(:,:,1));                      %调取第一层，转换为双精度
lamda=6328*10^(-10);k=2*pi/lamda;          %赋值波长,单位:米,波矢
D=0.01;                                    %赋值透镜的孔径,单位:米
f=0.4;                                     %赋值透镜的焦距,单位:米
figure,imshow(Uo,[])
[c,r]=size(Uo);                            %读取物面采样数
Lo=0.005;                                  %物面尺寸,单位:米
do=1.2;                                    %物到透镜的距离,亦即物瞳距,单位:米,可以改变
cutoff_frequency=D/2/lamda/do;             %截止频率(阿贝观点)
kethi=linspace(-1./2./Lo,1./2./Lo,r).*r; nenta=linspace(-1./2./Lo,1./2./Lo,c).*c;
[kethi,nenta]=meshgrid(kethi,nenta);       %物的二维频谱网格
H=zeros(c,r);                              %预设传递函数
for n=1:c
   for m=1:r
      if kethi(n,m).^2+nenta(n,m).^2<= cutoff_frequency.^2;
         H(n,m)=1;
      end
   end
end
figure,imshow(H,[]);title('传递函数')
Gg=fftshift(fft2(Uo));                      %物的频谱
Gi=Gg.*H;                                   %像的频谱
Ui=ifft2(Gi);                               %像的光场分布
Ii=Ui.*conj(Ui);                            %像的光强分布
figure,imshow(Ii,[]),title('像的光强分布')