ph=peaks(256).*4;              %预设原始相位
mu=0.2;                        %噪声系数
noise=rands(256,256).*mu.*pi;  %噪声
ph=ph+noise;                   %叠加有噪声的预设原始相位
figure,imshow(ph,[]);          %显示原始含噪声连续相位
U=exp(j*ph); 
wrapped_phs =angle(U);         %对应的包裹相位
figure,imshow(wrapped_phs,[]); %显示包裹相位
a=wrapped_phs;                 %将包裹相位wrapped_phs赋值给a
%下面开始进行最小二乘解包裹运算
[M,N]=size(a);                 %计算二维包裹相位的大小（行、列数）
dx=zeros(M,N);dy=zeros(M,N);   %预设包裹相位沿x方向和y方向的梯度
m=1:M-1; 
dx(m,:)=a(m+1,:)-a(m,:);       %计算包裹相位沿x方向的梯度
dx=dx-pi*round(dx/pi);         %去除梯度中的跳跃
n=1:N-1;
dy(:,n)=a(:,n+1)-a(:,n);       %计算包裹相位沿y方向的梯度
dy=dy-pi*round(dy/pi);         %去除梯度中的跳跃
p=zeros(M,N);p1=zeros(M,N);p2=zeros(M,N); %为计算ρnm作准备
m=2:M;
p1(m,:)=dx(m,:)-dx(m-1,:);     %计算Δgxnm-Δgx(n-1)m
n=2:N;
p2(:,n)=dy(:,n)-dy(:,n-1);     %计算Δgynm–Δgyn(m-1)
p=p1+p2;                       %计算ρnm
p(1,1)=dx(1,1)+dy(1,1);        %计算ρnm
n=2:N;
p(1,n)=dx(1,n)+dy(1,n)-dy(1,n-1);%赋值Neumann边界条件
m=2:M;
p(m,1)=dx(m,1)-dx(m-1,1)+dy(m,1);
pp=dct2(p)+eps;                %计算ρnm的DCT
fi=zeros(M,N);
for m=1:M                      %计算Φnm在DCT域的精确解
   for n=1:N  
      fi(m,n)=pp(m,n)/(2*cos(pi*(m-1)/M)+2*cos(pi*(n-1)/N)-4+eps);
   end
end
fi(1,1)=pp(1,1);                %赋值DCT域的Φ11
phs=idct2(fi);                  %用iDCT计算解包裹相位在空域中的值
figure,imshow(phs,[])           %显示解包裹相位