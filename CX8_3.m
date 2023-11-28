x=linspace(0,0.0001*pi,512);y=linspace(0,0.0001*pi,512);
[x,y]=meshgrid(x,y);
b=0.00002; lamda=6328*10^(-10);di=0.02;
tx=cos(2*pi*x./b);
[c,r]=size(tx);
figure,imshow(tx,[])
a=1.5*lamda*di/b
Lo=0.0001*pi;
cutoff_frequency=a/lamda/di;
kethi=linspace(-1./2./Lo,1./2./Lo,r).*r;nenta=linspace(-1./2./Lo,1./2./Lo,c).*c;
[kethi,nenta]=meshgrid(kethi,nenta); 
CTF=zeros(c,r);
for n=1:c
   for m=1:r
      if kethi(n,m).^2+nenta(n,m).^2<= cutoff_frequency.^2;
      CTF(n,m)=1;
      end
   end
end
figure,imshow(CTF,[]);title('相干传递函数')
h=fftshift(fft2(CTF));
OTF=abs(fftshift(fft2(h.*conj(h))));
OTF=OTF./max(max(OTF));
figure,surfl(OTF),shading interp,colormap(gray);title('光学传递函数')
Gg=fftshift(fft2(tx));                             %相干照明时物的频谱
figure,imshow(abs(Gg),[0,max(max(abs(Gg)))/10])    %显示频谱
figure,plot(abs(Gg(257,:)))                        %频谱的剖线
Gic=Gg.*CTF;
Uic=ifft2(Gic);
Iic=Uic.*conj(Uic);
figure,imshow(Iic,[0,1]),title('相干照明'),colormap(gray)
Gg=fftshift(fft2(tx.*conj(tx)));  %非相干照明时物的频谱(物为复振幅的平方)
figure,imshow(abs(Gg),[0,max(max(abs(Gg)))/10])    %显示频谱
figure,plot(abs(Gg(257,:)))                        %频谱的剖线
Gii=Gg.*OTF; 
Iii=abs(ifft2(Gii));
figure,imshow(Iii,[0,1]),title('非相干照明'),colormap(gray)
figure, plot(Iii(257,:))
hold on,plot(Iic(257,:),'r')