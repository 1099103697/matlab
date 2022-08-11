clc;
clear;
m0=15800;
mp=12912;
tall=115;%发动机点火时间
P=260680;
mu=3.986004e14;
R=6378.14e3;
t1=3;
t2=20;%经估算，加速40s后，导弹速度约为0.8Ma
alp=-1/57.3;
c=-mp/tall;%秒耗量
tout=[];
yout=[];
r=[];
height=[];
yout(1,:)=[0 0 0 0 m0];
i=1;
h=0.1;
for t=0:h:115
    tout(i,1)=t;
    r(i,1)=sqrt(yout(i,1)^2+(R+yout(i,2))^2);
    height(i,1)=r(i,1)-R;
    gr=-mu/(r(i,1)^2);
   
    if t>0
       the=atan(yout(i,4)/yout(i,3));
    end
    
    if t<3
       phi=90/57.3;
    elseif t<20
        phi=the+alp;
    elseif t<115
        phi=the;
    else
        phi=the; 
        c=0; 
        P=0;
    end
    
    K1m=h*c;
    K2m=h*c;
    K3m=h*c;
    K4m=h*c;
    
    K1x=h*yout(i,3);
    K2x=h*(yout(i,3)+K1x/2);
    K3x=h*(yout(i,3)+K2x/2);
    K4x=h*(yout(i,3)+K3x);
    
    K1y=h*yout(i,4);
    K2y=h*(yout(i,4)+K1y/2);
    K3y=h*(yout(i,4)+K2y/2);
    K4y=h*(yout(i,4)+K3y); 
    
    K1vx=h*(P/yout(i,5)*cos(phi)+gr/r(i,1)*yout(i,1));
    K2vx=h*(P/(yout(i,5)+K1m/2)*cos(phi)+gr/r(i,1)*(yout(i,1)+K1x/2));
    K3vx=h*(P/(yout(i,5)+K2m/2)*cos(phi)+gr/r(i,1)*(yout(i,1)+K2x/2));
    K4vx=h*(P/(yout(i,5)+K3m)*cos(phi)+gr/r(i,1)*(yout(i,1)+K3x));
    
    K1vy=h*(P/yout(i,5)*sin(phi)+gr/r(i,1)*(R+yout(i,2)));
    K2vy=h*(P/(yout(i,5)+K1m/2)*sin(phi)+gr/r(i,1)*(R+yout(i,2)+K1y/2));
    K3vy=h*(P/(yout(i,5)+K2m/2)*sin(phi)+gr/r(i,1)*(R+yout(i,2)+K2y/2));
    K4vy=h*(P/(yout(i,5)+K3m)*sin(phi)+gr/r(i,1)*(R+yout(i,2)+K3y));  
    
    yout(i+1,1)=yout(i,1)+1/6*(K1x+2*K2x+2*K3x+K4x);
    yout(i+1,2)=yout(i,2)+1/6*(K1y+2*K2y+2*K3y+K4y);
    yout(i+1,3)=yout(i,3)+1/6*(K1vx+2*K2vx+2*K3vx+K4vx);
    yout(i+1,4)=yout(i,4)+1/6*(K1vy+2*K2vy+2*K3vy+K4vy);
    yout(i+1,5)=yout(i,5)+1/6*(K1m+2*K2m+2*K3m+K4m);
    i=i+1;
end

figure(1);
plot(yout(:,1)/1000,yout(:,2)/1000);
grid on;
title('导弹在发射坐标系下的轨迹');
xlabel('x/km');ylabel('y/km');
        
  