Uo=imread('1.bmp');                     %读入数字全息图
Uo =double(Uo (:,:,1));                 %读取数字全息图的第一层(红色)
figure,imshow(Uo,[])
lamda=6328*10^(-10);k=2*pi/lamda;       %赋值波长、波数
[r,c]=size(Uo);                         %计算全息图的大小(像素)
Lox=c*4.65*10^(-6)                      %赋值全息图的尺寸(纵向),单位:米
Loy=r*4.65*10^(-6)                      %赋值全息图的尺寸(纵向),单位:米
xo=linspace(-Lox/2,Lox/2,c);            %生成全息图的x坐标
yo=linspace(-Loy/2,Loy/2,r);            %生成全息图的y坐标
[xo,yo]=meshgrid(xo,yo);                %生成全息图的坐标网格
zo=0.052;                               %衍射面(物)到全息记录面的距离,单位:米
%下面用T-FFT算法重构再现像
F0=exp(j*k*zo)/(j*lamda*zo);
F1=exp(j*k/2/zo.*(xo.^2+yo.^2));
fF1=fft2(F1);
fa1=fft2(Uo);
Fuf1=fa1.*fF1; 
Ui=F0.*fftshift(ifft2(Fuf1));          %在观察平面上的重构光场
Ii=Ui.*conj(Ui);                       %再现像的光强
figure,imshow(Ii,[0,max(max(Ii))/1]),colormap(gray),title('-1级像')