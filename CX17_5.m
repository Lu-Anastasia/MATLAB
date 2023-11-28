Uo=imread('guang.bmp');       %读入图像
Uo=imresize(Uo,[256,256]);    %将图像调整到256×256像素
Uo=double(Uo(:,:,1));
figure;imshow(Uo,[]);         %显示输入图像
[r,c]=size(Uo);               %返回矩阵大小
ef=0.5;                       %设置纯相位随机噪声系数
FUo=fftshift(fft2(Uo.*exp(j.*rands(r,c).*pi.*ef)));   %计算傅里叶变换
phi=angle(FUo);               %计算相位分布
H=mod(phi,2*pi);              %计算相位取模数2π的余数
H1=round(H/max(H(:))*255);    %将余数量化为255个灰度级，绘制全息图
CGH=exp(j.*H1/40.58);  %漂白制作成纯相位全息图（理想情况下相位还原到0~2π）
rU=ifft2(CGH);                %进行逆傅里叶变换得到再现光场
figure;imshow(H1,[]);         %显示编码得到的全息图
figure;imshow(rU.*conj(rU),[0,0.0001]); %显示再现像