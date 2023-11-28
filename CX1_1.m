fxy=cos(peaks(256).*2+pi)+1;            %构建连续带限函数
[rr,cc]=size(fxy);                      %计算连续函数的大小
figure,imshow(fxy,[])                   %显示连续函数
F=fftshift(fft2(fxy));                  %计算连续函数的频谱
figure,plot(abs(F(round(rr/2)+1,:))),   %观察带宽
figure,plot(abs(F(:,round(cc/2)+1))),   %观察带宽
figure,surfl(abs(F)),shading interp,colormap(gray);   %频谱3D图
combxy=zeros(rr,cc);                    %开始生成comb函数
X=4;Y=4;                                %抽样间隔
for n=1:Y:rr
   for m=1:X:cc
   combxy(n,m)=1;
   end
end
figure,imshow(combxy,[]);               %显示comb函数
C=fftshift(fft2(combxy));               %计算comb函数的频谱
figure,surfl(abs(C)),shading interp,colormap(gray); % 频谱3D图
gxy=zeros(rr,cc);                       %开始生成抽样函数
gxy=fxy.*combxy;                        %生成抽样函数
figure,imshow(gxy,[]);                  %显示抽样函数
Gs=fftshift(fft2(gxy));                 %计算抽样函数的频率
figure,surfl(abs(Gs)),shading interp,colormap(gray); %频谱3D图
figure,plot(abs(Gs(:,cc/2+1))),         %观察频谱是否有重叠
By=round(rr/2/Y);Bx=round(cc/2/X);      %二维矩函数滤波器的宽度
H=zeros(rr,cc);                         %开始生成二维矩函数滤波器
H(round(rr/2)+1-By:round(rr/2)+1+By-1,round(cc/2)+1-Bx:round(cc/2)+1+Bx-1)=1;   
figure,imshow(H,[])                     %显示二维矩函数滤波器
Gsyp=H.*Gs;                             %滤波计算原函数频谱
figure,surfl(abs(Gsyp)),shading interp,colormap(gray);
gxyyp=X*Y.*abs(ifft2(Gsyp));            %逆傅里叶变换计算原函数
figure,imshow(gxyyp,[])                 %显示还原的原函数