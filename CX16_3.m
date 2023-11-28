I4=imread('cell21.bmp');I4=double(I4(7:518,10:521,1));  %读取第一幅像面全息图
I3=imread('cell22.bmp');I3=double(I3(7:518,10:521,1));  %读取第二幅像面全息图
I2=imread('cell23.bmp');I2=double(I2(7:518,10:521,1));  %读取第三幅像面全息图
I1=imread('cell24.bmp');I1=double(I1(7:518,10:521,1));  %读取第四幅像面全息图
figure,imshow(I1,[])                                    %显示像面全息图
figure,imshow(I2,[])                                    %显示像面全息图
figure,imshow(I3,[])                                    %显示像面全息图
figure,imshow(I4,[])                                    %显示像面全息图
I1=medfilt2(I1,[3,3]);                                  %用中值滤波对全息图降噪
I2=medfilt2(I2,[3,3]);                                  %用中值滤波对全息图降噪
I3=medfilt2(I3,[3,3]);                                  %用中值滤波对全息图降噪
I4=medfilt2(I4,[3,3]);                                  %用中值滤波对全息图降噪
C1=I4-I2;                                               %全息图相减
C2=I1-I3;                                               %全息图相减
ph=atan(C1./(C2+eps));                                  %计算包裹相位
figure;imshow(ph,[]);title('重建相位图');               %显示重建相位
figure,plot(ph(257,:))                                  %绘制重建相位剖线