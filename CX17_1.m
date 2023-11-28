Uo=imread('guang.bmp');Uo=double(Uo(:,:,1)); %调入作为物的图像
sizeUor=64; sizeUoc=64;                    %设置图像大小
Uo=imresize(Uo,[sizeUor,sizeUoc]);         %调整图像大小
figure,imshow(Uo,[])                       %显示物
[r,c]=size(Uo);                            %查询图像大小
ef1=0.8                                    %设置噪声系数（为了降低频谱的动态范围）
FUo=fftshift(fft2(Uo.*exp(j.*rands(r,c).*pi.*ef1))); % 对物光场进行傅里叶变换
Am=abs(FUo);                               %计算模的大小
Ph=mod(angle(FUo),2*pi);                   %计算相位的大小并取2π的余数
Ph=Ph./2/pi;                               %余数除2π
w=6                                        %通光孔径的宽度（设为偶数）
wth=round(w/2);                            %通光孔径的一半
s=2*w;                                     %设置编码单元的大小（像素）
CGH=zeros(r*s,c*s);                        %预设计算全息矩阵
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面完成编码的参数运算（罗曼Ⅲ型）
maxAm=max(abs(FUo(:)));                    %计算频谱的模的最大值
ef2=1.5                                    %设置阈值系数
th=maxAm/ef2;                              %设置阈值
Am(Am>th)=th;                              %将大于阈值的模设置为等于阈值
lmn=round(Am/th*w);                        %将模进行归一化,然后量化到整数0,1,2,…,w
pmn=round(Ph*w);                           %将相位量化到0,1,2,…,w
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面完成全息图赋值运算
for m=1:r;                                 %遍历各行
    for n=1:c;                             %遍历各列
        cgh=zeros(s,s);                    %预设计算全息图抽样单元
        if lmn(m,n)==0;                    %对模不为零的点进行编码
            elseif pmn(m,n)<=wth           %相位值小于π/2，不出现模式溢出
            cgh(w+1-lmn(m,n):w+lmn(m,n),w-wth+pmn(m,n)+1:2*w-wth+pmn(m,n))=1;
            elseif pmn(m,n)>wth            %相位值大于π/2，出现模式溢出
            cgh(w+1-lmn(m,n):w+lmn(m,n),w-wth+pmn(m,n)+1:s)=1;
            cgh(w+1-lmn(m,n):w+lmn(m,n),1:pmn(m,n)-wth)=1;
        end
        CGH((m-1)*s+1:m*s,(n-1)*s+1:n*s)=cgh;%将生成的抽样单元放到计算全息图中
    end
end
figure;imshow(CGH,[]);                     %迂回位相编码结果（罗曼Ⅲ型）
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%下面完成再现运算
RU=fftshift(ifft2(CGH));                   %用逆傅里叶变换计算再现光场
RI=RU.*conj(RU);                           %计算再现像的光强分布
figure;imshow(RI,[0,max(RI(:))/1000]),colormap(pink); %显示再现图像