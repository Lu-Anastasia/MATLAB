II=imread('rw630mm.jpg');              %读入数字全息图
II=double(II(:,:,1));                  %读取数字全息图的第一层(红色)
figure,imshow(II,[])
holo=fftshift(fft2(II));               %计算再现光场
Ii=holo.*conj(holo);                   %计算再现像光强
figure,imshow(Ii,[0,max(max(Ii))./ 200000]),colormap(pink),title('再现像')