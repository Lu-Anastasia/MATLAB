function varargout=unwrap_one_d(w_ph,s_ph,th)
%本函数完成一行（或列）包裹相位的解包裹运算
%格式 un_ph= unwrap_one_d(w_ph,s_ph,th)；
%输入的三个变量:w_ph为包裹相位，s_ph（为起始相位值，而th为域值，
%若相位值包裹在(-π,π]，域值th取π；若相位值包裹在(-π/2,π/2]，th取π/2
%输出变量un_ph为解包裹相位值
numb=length(w_ph);            %计算包裹相位的大小（元素个数）
un_ph=w_ph;                   %将包裹相位值赋值给un_ph——预设解包裹相位
un_ph(1)=s_ph;                %将起点相位赋值给un_ph的第一个元素（不再参加解包裹）
for n=2:numb                  %从第二个元素开始解包裹运算
   delta=w_ph(n)-w_ph(n-1);   %相邻两个包裹相位作比较（相减）
   if (abs(delta)<=th)        %差值不大于域值
       un_ph(n)=un_ph(n-1)+delta;  %解包裹相位等于前一个解包裹相位值加差值
   elseif (delta<0)           %差值大于域值，且差值为负
       un_ph(n)=un_ph(n-1)+delta+2*th;%前一个包解裹相位值加差值再加2倍域值
   else                       %差值大于域值，且差值为正
       un_ph(n)=un_ph(n-1)+delta-2*th; %前一个解包裹相位值加差值再减2倍域值
   end
end
varargout{1}=un_ph;           %将解包裹相位值赋值给输出变量