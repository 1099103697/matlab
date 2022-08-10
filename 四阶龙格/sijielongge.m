clc;clear all;close all;
h=0.01
x=0:h:1
n=length(x)
x(1)=0
y(1)=0;
for i=1:n-1
    k1=h*fun(x(i),y(i));
    k2=h*fun(x(i)+h/2,y(i)+k1*1/2);
    k3=h*fun(x(i)+h/2,y(i)+k2*1/2);
    k4=h*fun(x(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end
plot(x,y)