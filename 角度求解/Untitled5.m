x1 = 1;y1 = 0;
x2 = 0;y2 =0;
x3 = 5;y3 = 5;
a2 = (x1-x2).^2+(y1-y2).^2;
b2 = (x3-x2).^2+(y3-y2).^2;
c2 = (x1-x3).^2+(y1-y3).^2;
a = sqrt(a2);
b = sqrt(b2);
c = sqrt(c2);
pos = (a2+b2-c2)/(2*a*b);    %求出余弦值
angle = acos(pos);         %余弦值装换为弧度值
realangle = angle*180/pi;   %弧度值转换为角度值
disp(realangle);