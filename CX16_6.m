ph=peaks(256)*1;                 %预设原始相位
mu=0.5;                          %噪声系数
noise=rands(256,256).*mu.*pi;    %噪声
ph=ph+noise;                     %叠加有噪声的预设原始相位
U=exp(j*ph); ph0=angle(U);       %对应的包裹相位
[r,c]=size(ph0);
figure,imshow(ph0,[])
phs=ph0;                         %预设解包裹后的相位
for n=1:r-1                      %逐行解包裹
   for m=1:c-1                   %逐列解包裹
      if phs(n,m)-phs(n+1,m)>=pi
        phs(n+1,m)=phs(n+1,m)+2*pi;
      end
      if phs(n,m)-phs(n+1,m)<=-pi
        phs(n+1,m)=phs(n+1,m)-2*pi;
      end
   end
end
figure,imshow(phs,[])