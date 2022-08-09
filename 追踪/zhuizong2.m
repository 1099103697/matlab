clear;

clc

h = 0.0001;%时间步长

k = 1;

L = 120; Vs = 90; Vm = 450;

t(1) = 0;  x(1) = 0;  y(1) = 0;%初始值

while y <= L

    x(k+1) = x(k) + Vm*h / sqrt(1+((L-y(k))/(Vs*t(k)-x(k)))^2);

    y(k+1) = y(k) + Vm*h / sqrt(1+((Vs*t(k)-x(k))/(L-y(k)))^2);

    t(k+1) = h*k;

    k = k+1;

end

plot(x,y,x(1):0.05:x(end),L)

t = t(end),x = x(end),y = y(end)