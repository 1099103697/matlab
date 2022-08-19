clear all

%% ��ʼ��������
m=260;%����
s=0.24;%�ο����
M=2;%�����
hi0=5000;%��ʼ�߶�
inc0=30;%��ǣ�����
v0=680;%����
g0=9.8;%��ʼ�������ٶ�
R=6370856;%����뾶
a=2;%����
clb0=inc0-a;%������
clb1=clb0*pi/180;

inc=inc0;%��ʼ�����
clb=clb0;%��ʼ��������
hi=hi0;%��ʼ���߶�
v=v0;%��ʼ�ٶ�


%% ����
g=g0*R^2/(R+hi)^2;%ʵʱ�������ٶ�
G=m*g;%����

%% ��������
Cx=[M^2 M 1]*[0.0002 0.0038 0.1582;-0.0022 -0.0132 -0.8520;0.0115 -0.0044 1.9712]*[a^2;a;1];%����ϵ��
Cy=[M^2 M 1]*[-0.026;0.0651;0.4913]*a; %����ϵ��
P=1.225*exp(-0.00015*hi); %�ܶ�
Fx=0.5*Cx*P*v^2*s;%����
Fy=0.5*Cy*P*v^2*s;%����

%% �˶�ѧ����
%����Ϊ0���ʵ㵯������Ǳ仯
%m*dv=-Fx+G*sin(clb) 
%m*v*d_clb=Fy-G*cos(clb) 
%dx=v*cos(clb)
%dy=v^sin(clb) 
%Wz=d_clb 

h=0.1;%ȡ����Ϊ0.1s
x=0;%��ʼ������Ϊ0
y=5000;%��ʼ������Ϊ5000
i=0;%��ʼ��������

datav(1,1)=340;
dataclb(1,1)=28;
datax(1,1)=0;
datahi(1,1)=5000;
for t=0:0.1:80
        if hi<0
            break
        end
        i=i+1;
        clb=RK(h,clb,G,Fy,m,v);%����������
        clb1=clb*pi/180;
        dataclb(1,i)=clb;
        v=v+(G*sin(clb1)-Fx)*h/m;%�����ٶ�
        datav(1,i)=v;
        hi=hi-v*sin(clb1)*h;%���¸߶�
        datahi(1,i)=hi;
        x=x+v*cos(clb1)*h;%���º�����
        datax(1,i)=x;
        P=1.225*exp(-0.00015*hi); %�����ܶ�
        Fx=0.5*Cx*P*v^2*s;%��������
        Fy=0.5*Cy*P*v^2*s;%��������
        g=g0*R^2/(R+hi)^2;%����ʵʱ�������ٶ�
        G=m*g;%��������
       % waitbar(t);
end

%% ͼ������

% [m,n]=size(dataclb)
% nn=0.1:0.1:n*0.1;
% figure(1)
% subplot(2,1,1);
% plot(nn,datahi);xlabel('ʱ��/s');ylabel('�߶�/m');title('�߶�/ʱ��');grid on;
% figure(1)
% subplot(2,1,2);
% plot(nn,datav);axis([0 inf,0,inf]);xlabel('ʱ��/s');ylabel('�ٶ�/��m/s��');title('�ٶ�/ʱ��');grid on;
% figure(2)
% plot(nn,dataclb);xlabel('ʱ��/s');ylabel('�ٶ����/��');title('�ٶ����/ʱ��');grid on;
% figure(3)
% plot(datax,datahi);xlabel('ˮƽ�˶�����');ylabel('�߶�');title('ˮƽ�˶�����/�߶�');grid on;
