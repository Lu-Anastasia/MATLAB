I4=imread('cell21.bmp');I4=double(I4(7:518,10:521,1));  %读取第一幅像面全息图
I3=imread('cell22.bmp');I3=double(I3(7:518,10:521,1));  %读取第二幅像面全息图
I2=imread('cell23.bmp');I2=double(I2(7:518,10:521,1));  %读取第三幅像面全息图
I1=imread('cell24.bmp');I1=double(I1(7:518,10:521,1));  %读取第四幅像面全息图
I1=medfilt2(I1,[3,3]);                                  %用中值滤波对全息图降噪
I2=medfilt2(I2,[3,3]);                                  %用中值滤波对全息图降噪
I3=medfilt2(I3,[3,3]);                                  %用中值滤波对全息图降噪
I4=medfilt2(I4,[3,3]);                                  %用中值滤波对全息图降噪
C1=I4-I2;                                               %全息图相减
C2=I1-I3;                                               %全息图相减
ph=atan(C1./(C2+eps));                                  %计算包裹相位
figure;imshow(ph,[]);title('包裹相位图');               %显示重建相位
[r,c]=size(ph);
unph=ph;                                                %预设解包裹后的相位
unph(1,:)=unwrap_one_d(ph(1,:),ph(1,1),pi/2);%先完成第一行元素的解包裹运算
for n=1:c
    %下面完成第n列元素的解包裹运算，起始相位值为第一行上元素的解包裹相位值
    unph(:,n)=unwrap_one_d(ph(:,n),unph(1,n),pi/2);
end
figure,surfl(unph),shading interp,colormap(gray)        %显示解包裹相位
figure,imshow(unph,[])