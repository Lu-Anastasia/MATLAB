I4=imread('cell21.bmp');I4=double(I4(7:518,10:521,1));  %读取第一幅像面全息图
I3=imread('cell22.bmp');I3=double(I3(7:518,10:521,1));  %读取第二幅像面全息图
I2=imread('cell23.bmp');I2=double(I2(7:518,10:521,1));  %读取第三幅像面全息图
I1=imread('cell24.bmp');I1=double(I1(7:518,10:521,1));  %读取第四幅像面全息图
I1=medfilt2(I1,[3,3]);                            %用中值滤波对全息图降噪
I2=medfilt2(I2,[3,3]);                            %用中值滤波对全息图降噪
I3=medfilt2(I3,[3,3]);                            %用中值滤波对全息图降噪
I4=medfilt2(I4,[3,3]);                            %用中值滤波对全息图降噪
C1=I4-I2;                                         %全息图相减
C2=I1-I3;                                         %全息图相减
ph=atan(C1./(C2+eps));                            %计算包裹相位
figure;imshow(ph,[]);title('包裹相位图');         %显示重建相位
a=ph;                                             %将包裹相位ph赋值给a
%下面开始进行最小二乘解包裹运算
[M,N]=size(a);                                    %计算二维包裹相位的大小（行、列数）
dx=zeros(M,N);dy=zeros(M,N);                      %预设包裹相位沿x方向和y方向的梯度
m=1:M-1; 
dx(m,:)=a(m+1,:)-a(m,:);                          %计算包裹相位沿x方向的梯度
dx=dx-pi*round(dx/pi);                            %去除梯度中的跳跃
n=1:N-1;
dy(:,n)=a(:,n+1)-a(:,n);                          %计算包裹相位沿y方向的梯度
dy=dy-pi*round(dy/pi);                            %去除梯度中的跳跃
p=zeros(M,N);p1=zeros(M,N);p2=zeros(M,N); %为计算ρnm作准备
m=2:M;
p1(m,:)=dx(m,:)-dx(m-1,:);                        %计算Δgxnm-Δgx(n-1)m
n=2:N;
p2(:,n)=dy(:,n)-dy(:,n-1);                        %计算Δgynm–Δgyn(m-1)
p=p1+p2;                                          %计算ρnm
p(1,1)=dx(1,1)+dy(1,1);                           %计算ρnm
n=2:N;
p(1,n)=dx(1,n)+dy(1,n)-dy(1,n-1);                 %赋值Neumann边界条件
m=2:M;
p(m,1)=dx(m,1)-dx(m-1,1)+dy(m,1);
pp=dct2(p)+eps;                                   %计算ρnm的DCT
fi=zeros(M,N);
for m=1:M                                         %计算Φnm在DCT域的精确解
   for n=1:N  
      fi(m,n)=pp(m,n)/(2*cos(pi*(m-1)/M)+2*cos(pi*(n-1)/N)-4+eps);
   end
end
fi(1,1)=pp(1,1);                                  %赋值DCT域的Φ11
phs=idct2(fi);                                    %用iDCT计算解包裹相位在空域中的值
figure,surfl(phs),shading interp,colormap(gray)   %显示解包裹相位
%为检验解包裹效果，对解包裹相位作再包裹运算
U=exp(j.*phs);                                    %用解包裹相位构建光场
rewrap_ph=atan(imag(U)./real(U));                 %作再包裹运算
figure,imshow(rewrap_ph,[])                       %显示再包裹结果