Uo=imread('guang.bmp');             %调入作为物的图像
Uo=double(Uo (:,:,1));              %取第一层，并转为双精度
[r,c]=size(Uo);
Uo=ones(r,c)*0.98-Uo/255*0.5;       %将物转换为高透射率射系数体
figure,imshow(Uo,[0,1]),title('物')
lamda=6328*10^(-10);k=2*pi/lamda;   %赋值波长和波数
Lo=5*10^(-3)                        %赋值衍射面(物)的尺寸
xo=linspace(-Lo/2,Lo/2,r);yo=linspace(-Lo/2,Lo/2,c);
[xo,yo]=meshgrid(xo,yo);            %生成衍射面(物)的坐标网格
zo=0.20;                            %全息记录面到衍射面的距离,单位:米
%下面用T-FFT算法完成物面到全息记录面的衍射计算
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa1=fft2(Uo);
Fuf1=fa1.*fF1; 
Uh=F0.*fftshift(ifft2(Fuf1)); 
Ih=Uh.*conj(Uh);
figure,imshow(Ih,[0,max(max(Ih))/1]),title('全息图')
%下面用T-FFT算法完成全息面到观察面的衍射计算(重构再现像)
zi=0.20                            %赋值再现距离(可以调整)
F0i=exp(j*k*zi)/(j*lamda*zi);
F1i=exp(j*k/2/zi.*(xo.^2+yo.^2));  %T-FFT算法，物面、全息图和再现像尺寸相同
fF1i=fft2(F1i);
fIh=fft2(Ih); 
FufIh=fIh.*fF1i; 
Ui=F0i.*fftshift(ifft2(FufIh)); 
Ii=Ui.*conj(Ui);
figure,imshow(Ii,[0,max(max(Ii))/1]),title('再现像')