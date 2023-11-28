ph=peaks(256).*4+rands(256,256).*0.5.*pi;  %生成原始相位（含噪声）
U=exp(j*ph);ph0=angle(U);                  %生成对应的包裹相位
[r,c]=size(ph0);
figure,imshow(ph0,[])
residue=zeros(r,c);                        %预设残点矩阵
for n=1:r-1                                %开始计算残点
   for m=1:c-1
       C(n,m)=round((ph0(n,m)-ph0(n,m+1))./2./pi)+round((ph0(n+1,m)-ph0(n,m))./2./pi)...
       +round((ph0(n+1,m+1)-ph0(n+1,m))./2./pi)+round((ph0(n,m+1)-ph0(n+1,m+1))./2./pi);
       if C(n,m)==1
          residue(n,m)=1;                  %用“1”标识正残点
       end
       if C(n,m)==-1
          residue(n,m)=-1;                 %用“-1”标识负残点
       end
   end
end
figure,imshow(residue,[])