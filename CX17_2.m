alpha=600000;                          %设置载频
Uo=imread('guang1.jpg');               %读入图像
Uo=imresize(Uo,[256,256]);             %调整图像大小
Uo=double(Uo (:,:,1));                 %读取图像第一层为物
figure;imshow(Uo,[]);                  %显示输入图像
[r,c]=size(Uo);                        %计算矩阵大小
l0x=0.005;l0y=0.005;                   %设置物面尺寸
x=linspace(-l0x/2,l0x/2,c);y=linspace(-l0y,l0y,r);
[x,y]=meshgrid(x,y);                   %设置物面二维网格
ef=0.5;                                %设置随机纯相位噪声的系数
FUo=fftshift(fft2(Uo.*exp(j.*rands(r,c).*pi.*ef))); %傅里叶变换计算频谱并移频
A=abs(FUo);                            %计算频谱的振幅
A=A./max(A(:));                        %归一化振幅
phi=angle(FUo);                        %计算频谱的相位分布
H=0.5*[1+A.*cos(2*pi*alpha.*x-phi)];   %进行博奇编码
H=round(H./max(H(:)).*255);            %将全息图量化为256个灰阶
rU=ifft2(H);                           %逆傅里叶变换得到再现光场
figure;imshow(H,[]);                   %显示博奇编码得到的全息图
figure;imshow(rU.*conj(rU),[0,0.01]);  %显示再现的图像