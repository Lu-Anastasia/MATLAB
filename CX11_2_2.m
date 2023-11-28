I=imread('ho.bmp');                   %读入数字全息图       
II=zeros(2048,2048);                  %预设全息图（补零）
II(257:1792,:)=double(I(:,:,2));      %将读取数字全息图的第二层(绿色)填到补零后的全息图中
figure,imshow(II,[])
FF=abs(fftshift(fft2(II)));           %计算全息图的频谱
figure,imshow(FF,[0,max(max(FF))/200])%显示计算全息图的频谱
lamda=5328*10^(-10);k=2*pi/lamda;     %赋值波长、波数
[r,c]=size(II);                       %计算全息图的大小(像素)
Lox=c*3.2*10^(-6)                     %赋值全息图的尺寸(纵向),单位:米
Loy=r*3.2*10^(-6)                     %赋值全息图的尺寸(纵向),单位:米
xo=linspace(-Lox/2,Lox/2,c);          %生成全息图的x坐标
yo=linspace(-Loy/2,Loy/2,r);          %生成全息图的y坐标
[xo,yo]=meshgrid(xo,yo);              %生成全息图的坐标网格
%下面用S-FFT算法重构再现像
zi=1.70;                              %全息记录面到像面的距离,单位:米
Lix=c*lamda*zi/Lox                    %给出像面的尺寸(x方向),单位:米
Liy=r*lamda*zi/Loy                    %给出像面的尺寸(y方向),单位:米
x=linspace(-Lix/2,Lix/2,c);y=linspace(-Liy/2,Liy/2,r);
[x,y]=meshgrid(x,y);
F0=exp(j*k*zi)/(j*lamda*zi)*exp(j*k/2/zi*(x.^2+y.^2));
F=exp(j*k/2/zi*(xo.^2+yo.^2)); 
% 取再现照明光垂直入射C=1
holo=Lox/c*Loy/r*fftshift(fft2(II.*F*1)); holo=holo.*F0;
Ii=holo.*conj(holo);
figure,imshow(Ii,[0,max(max(Ii))./1000]),colormap(pink),title('S-FFT再现像')