clear all

%% 初始条件定义
m=260;%质量
s=0.24;%参考面积
M=2;%马赫数
hi0=5000;%初始高度
inc0=30;%倾角，向下
v0=680;%音速
g0=9.8;%初始重力加速度
R=6370856;%地球半径
a=2;%攻角
clb0=inc0-a;%爬升角
clb1=clb0*pi/180;

inc=inc0;%初始化倾角
clb=clb0;%初始化爬升角
hi=hi0;%初始化高度
v=v0;%初始速度


%% 重力
g=g0*R^2/(R+hi)^2;%实时重力加速度
G=m*g;%重力

%% 空气动力
Cx=[M^2 M 1]*[0.0002 0.0038 0.1582;-0.0022 -0.0132 -0.8520;0.0115 -0.0044 1.9712]*[a^2;a;1];%阻力系数
Cy=[M^2 M 1]*[-0.026;0.0651;0.4913]*a; %升力系数
P=1.225*exp(-0.00015*hi); %密度
Fx=0.5*Cx*P*v^2*s;%阻力
Fy=0.5*Cy*P*v^2*s;%升力

%% 运动学方程
%推力为0，质点弹道无倾角变化
%m*dv=-Fx+G*sin(clb) 
%m*v*d_clb=Fy-G*cos(clb) 
%dx=v*cos(clb)
%dy=v^sin(clb) 
%Wz=d_clb 

h=0.1;%取步长为0.1s
x=0;%初始横坐标为0
y=5000;%初始纵坐标为5000
i=0;%初始化计数器

datav(1,1)=340;
dataclb(1,1)=28;
datax(1,1)=0;
datahi(1,1)=5000;
for t=0:0.1:80
        if hi<0
            break
        end
        i=i+1;
        clb=RK(h,clb,G,Fy,m,v);%更新爬升角
        clb1=clb*pi/180;
        dataclb(1,i)=clb;
        v=v+(G*sin(clb1)-Fx)*h/m;%更新速度
        datav(1,i)=v;
        hi=hi-v*sin(clb1)*h;%更新高度
        datahi(1,i)=hi;
        x=x+v*cos(clb1)*h;%更新横坐标
        datax(1,i)=x;
        P=1.225*exp(-0.00015*hi); %更新密度
        Fx=0.5*Cx*P*v^2*s;%更新阻力
        Fy=0.5*Cy*P*v^2*s;%更新升力
        g=g0*R^2/(R+hi)^2;%更新实时重力加速度
        G=m*g;%更新重力
       % waitbar(t);
end

%% 图像生成

% [m,n]=size(dataclb)
% nn=0.1:0.1:n*0.1;
% figure(1)
% subplot(2,1,1);
% plot(nn,datahi);xlabel('时间/s');ylabel('高度/m');title('高度/时间');grid on;
% figure(1)
% subplot(2,1,2);
% plot(nn,datav);axis([0 inf,0,inf]);xlabel('时间/s');ylabel('速度/（m/s）');title('速度/时间');grid on;
% figure(2)
% plot(nn,dataclb);xlabel('时间/s');ylabel('速度倾角/°');title('速度倾角/时间');grid on;
% figure(3)
% plot(datax,datahi);xlabel('水平运动距离');ylabel('高度');title('水平运动距离/高度');grid on;
