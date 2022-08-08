
%by dynamic
%see also http://www.matlabsky.com
%2008.12.3
%
clear all
clc
clear
tt=0.1;
sm=0.6*tt;
st=0.42*tt;
x(1)=0;y(1)=0;z(1)=0;
pmr(:,1)=[x(1);y(1);z(1)];
ptr(:,1)=[25;5;10];
m=3;
q(1)=0;
o(1)=0;
a(1)=0;
for(k=2:600)
ptr(:,k)=[25-0.42*cos(pi/6)*tt*k;5;10+0.42*sin(pi/6)*k*tt];
r(k-1)=sqrt((ptr(1,k-1)-pmr(1,k-1))^2+(ptr(2,k-1)-pmr(2,k-1))^2+(ptr(3,k-1)-pmr(3,k-1))^2);
c=sqrt((ptr(1,k)-pmr(1,k-1))^2+(ptr(2,k)-pmr(2,k-1))^2+(ptr(3,k)-pmr(3,k-1))^2);
b=acos((r(k-1)^2+st^2-c^2)/(2*r(k-1)*st));
dq=acos((r(k-1)^2-st^2+c^2)/(2*r(k-1)*c));
if abs(imag(b))>0
    
b=0.0000001;
end 
if abs(imag(dq))>0
dq=0.0000001;
end 
q(k)=q(k-1)+dq;
o(k)=o(k-1)+m*dq;
a(k)=o(k)-q(k);
c1=r(k-1)*sin(b)/sin(a(k)+b);
c2=r(k-1)*sin(a(k))/sin(a(k)+b);
c3=sqrt((c1-sm)^2+(c2-st)^2+2*(c1-sm)*(c2-st)*cos(a(k)+b));
dq=a(k)-acos(((c1-sm)^2+c3^2-(c2-st)^2)/(2*(c1-sm)*c3));
if abs(imag(dq))>0
dq=0.0000001;
end 
q(k)=q(k-1)+dq;
o(k)=o(k-1)+m*dq;
a(k)=o(k)-q(k);
c1=r(k-1)*sin(b)/sin(a(k)+b);
c2=r(k-1)*sin(a(k))/sin(a(k)+b);
c3=sqrt((c1-sm)^2+(c2-st)^2+2*(c1-sm)*(c2-st)*cos(a(k)+b));
dq=a(k)-acos(((c1-sm)^2+c3^2-(c2-st)^2)/(2*(c1-sm)*c3));
if abs(imag(dq))>0
dq=0.0000001;
end
q(k)=q(k-1)+dq;
o(k)=o(k-1)+m*dq;
a(k)=o(k)-q(k);
c1=r(k-1)*sin(b)/sin(a(k)+b);
c2=r(k-1)*sin(a(k))/sin(a(k)+b);
c3=sqrt((c1-sm)^2+(c2-st)^2+2*(c1-sm)*(c2-st)*cos(a(k)+b));
x1(k)=ptr(1,k-1)+c2/st*(ptr(1,k)-ptr(1,k-1));
y1(k)=ptr(2,k-1)+c2/st*(ptr(2,k)-ptr(2,k-1));
z1(k)=ptr(3,k-1)+c2/st*(ptr(3,k)-ptr(3,k-1));
x(k)=pmr(1,k-1)+sm/c1*(x1(k)-pmr(1,k-1));
y(k)=pmr(2,k-1)+sm/c1*(y1(k)-pmr(2,k-1));
z(k)=pmr(3,k-1)+sm/c1*(z1(k)-pmr(3,k-1));
pmr(:,k)=[x(k);y(k);z(k)];
r(k)=sqrt((ptr(1,k)-pmr(1,k))^2+(ptr(2,k)-pmr(2,k))^2+(ptr(3,k)-pmr(3,k))^2);
if r(k)<0.06;
break;
end;
end
sprintf('遭遇时间：%3.1f',0.1*k);
figure(1);
plot3(pmr(1,1:k),pmr(2,1:k),pmr(3,1:k),'k',ptr(1,:),ptr(2,:),ptr(3,:));
axis([0 25 0 5 0 25]);
text(x(80),y(80),z(80),'\leftarrow 比例导引');
grid on