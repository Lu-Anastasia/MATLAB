o=imread('分辨率板_1.bmp');       %调入作为物的图像
o=double(o(:,:,1));               %取第一层，并转为双精度
figure,imshow(o,[]),title('原始像')
[r,c]=size(o);
lamda=6328*10^(-10);k=2*pi/lamda; %赋值波长、波数
zo=0.3086;                        %物到全息记录面的距离,单位:米
Lo=5*10^(-3)                      %赋值衍射面(物)的尺寸,单位:米
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%用T-FFT算法完成全息图记录过程的计算
xo=linspace(-Lo/2,Lo/2,c);yo=linspace(-Lo/2,Lo/2,r);
[xo,yo]=meshgrid(xo,yo);          %生成物面的坐标网格
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fa=fft2(o); fF1=fft2(F1);
Fuf=fa.*fF1; 
OH=F0.*fftshift(ifft2(Fuf));      %在全息记录面上的光场分布
I=OH.*conj(OH);                   %全息记录面上的光强分布
figure,imshow(I,[]),colormap(pink),title('衍射图')
%下面加入参考光
zr=zo                             %发散球面波(参考光）的半径
xr=0.004;yr=-0.004;               %发散球面波(参考光）的中心位置(米)
R=exp(j.*k.*zr).*exp(j.*k.*(((xo-xr).^2+(yo-yr).^2))/2./zr);  %参考球面光波
%下面计算参、物光在全息记录面上的干涉,得到全息图
inter=OH./max(max(sqrt(I)))+R;    %调节光束比，并使参、物光干涉
II=inter.*conj(inter);            %干涉得到全息图
figure,imshow(II,[]),title('全息图')
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%下面计算再现
holo=fftshift(fft2(II)); 
Ii=holo.*conj(holo);
figure,imshow(Ii,[0,max(max(Ii))./20000]),colormap(pink),title('再现像')