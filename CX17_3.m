Uo=imread('guang1.jpg');         %读入图像
Uo=double(Uo(:,:,1))./255;       %将图像的灰阶二值化为0和1
figure,imshow(Uo,[])             %显示作为物的图像
[r,c]=size(Uo);
ef=0.5;                          %设置随机相位噪声系数
phai=rands(r,c).*pi.*ef;         %生成随机相位
FUo=fftshift(fft2(Uo.*exp(j.*phai))); %作FFT变换
A=abs(FUo);                      %光场的振幅分布
figure,imshow(A,[])
real1=real(FUo)./max(A(:)).*255; %将实部量化到256灰阶
imag1=imag(FUo)./max(A(:)).*255; %将虚部量化到256灰阶
%下面计算李氏连续灰阶全息图
CGH=zeros(r*4,c*4);              %预设CGH矩阵
for n=1:r
   for m=1:c
   cgh=zeros(4,4);
      if real(FUo(n,m))>=0&imag(FUo(n,m))>=0         %若实部大于0且虚部也大于0
        cgh(:,1)=real1(n,m);cgh(:,2)=imag1(n,m);     %赋值第1列和第2列
        elseif real(FUo(n,m))>=0&imag(FUo(n,m))<0    %若实部大于0而虚部小于0
        cgh(:,1)=real1(n,m);cgh(:,4)=abs(imag1(n,m));%赋值第1列和第4列
        elseif real(FUo(n,m))<0&imag(FUo(n,m))>=0    %若实部小于0而虚部大于0
        cgh(:,3)=abs(real1(n,m));cgh(:,2)=imag1(n,m);%赋值第3列和第2列
        elseif real(FUo(n,m))<0&imag(FUo(n,m))<0     %若实部小于0且虚部也小于0
        cgh(:,3)=abs(real1(n,m));cgh(:,4)=abs(imag1(n,m)); %赋值第3列和第4列
      end
    CGH((n-1)*4+1:n*4,(m-1)*4+1:m*4)=cgh; %将生成的抽样单元放到计算全息图中
    end
end
figure,imshow(CGH,[])            %显示计算全息图
rU=fftshift(ifft2(CGH));         %逆傅里叶变换得到再现光场
rI= rU.*conj(rU);                %计算再现像
figure,imshow(rI,[0,max(max(rI))/10000])%显示再现像