function [y2] =RK(h,y1,G,Fy,m,v)
%ËÄ½×Áú¸ñ¿âËþ

k1=Dclb(y1,G,Fy,m,v);

%t2=t+0.5*h;
yy2=y1+0.5*h*k1;
k2=Dclb(yy2,G,Fy,m,v);

yy3=y1+0.5*h*k2;
k3=Dclb(yy3,G,Fy,m,v);

%t3=t+h;
yy4=y1+h*k3;
k4=Dclb(yy4,G,Fy,m,v);

y2=y1+h*(k1+2*k2+2*k3+k4)/6;

end

