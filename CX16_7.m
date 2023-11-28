ph=peaks(256)*3;                     %预设原始相位
mu=0.5;                              %噪声系数
noise=rands(256,256).*mu.*pi;        %噪声
ph=ph+noise;                         %叠加有噪声的预设原始相位
U=exp(j*ph); ph0=angle(U);           %对应的包裹相位
[r,c]=size(ph0);
figure,imshow(ph0,[])
unph=ph0;                             %预设解包裹后的相位
unph(1,:)=unwrap_one_d(ph0(1,:),ph0(1,1),pi);%先完成第一行元素的解包裹运算
for n=1:c
    %下面完成第n列元素的解包裹运算，起始相位值为该列第一行上元素的解包裹相位值
    unph(:,n)=unwrap_one_d(ph0(:,n),unph(1,n),pi);
end
figure,imshow(unph,[])