close all
clear all
clc
f1=10;%�ز�1Ƶ��
f2=50;%�ز�2Ƶ��
%fc=1000;%�ز�Ϊ1000Hz
fbaud=600;%�����źŲ����ʣ�1s����600����Ԫ %��һ����Ԫ����ʱ��TB=1/600s
delta_f=abs(f1-f2);
k=fbaud/delta_f;
Kos=10;
fs=fbaud*Kos;%������
Ts=1/fs;%����ʱ����
T_total=5;%�ܵĲ���ʱ����5s
i=fbaud*T_total;%�����ź���Ԫ��
num=round(T_total/Ts);%������������
t1=linspace(0,T_total,i);
t=linspace(0,T_total,num);%0-5֮�����num����
phase_init=0;%���ʼ��λΪ0
d=round(rand(1,i));%�������ֻ����ź�

temp=zeros(1,i);
for index=1:i
    if d(index)==1
        temp(index)=0;
    else
        temp(index)=1;
    end
end

figure(1);
subplot(411);
plot(t1,d);
title('���ֻ����ź�d');
axis([0,T_total,-1,2]);

%(��ͬbit��Ӧ�Ĳ�ͬ��λ��(n),�����ʼ��λΪ0����ʱ������Ƶƫ��
theta=zeros(i,Kos);
%theta1=zeros(i,Kos);
%theta2=zeros(i,Kos);
for index=1:i
    for n=1:Kos
        if d(index)==1
            theta(index,n)=(2*(n-1)*pi)/(k*Kos);%d(i)=1��Ӧ����λ
            %theta2(index,n)=0;
        else
            %theta1(index,n)=0;
            theta (index,n)=(-2*(n-1)*pi)/(k*Kos);%d(i)=0��Ӧ����λ
        end
    end
end

%�ز��ź������
%freq1=cos(phase_init+theta1)+1j*sin(phase_init+theta1);
%freq2=cos(phase_init+theta2)+1j*sin(phase_init+theta2);
%e_fsk=freq1+freq2;
freq=cos(phase_init+theta)+1j*sin(phase_init+theta);
e_fsk=freq;
%ͨ���ŵ�
h=1;
noise=randn(index,Kos);
%r=h*e_fsk+noise;
SNR=10;
r=awgn(e_fsk,SNR);
%snr=var(h*e_fsk)/var(noise);

theta_est_temp=angle(r);
theta_est1=zeros(i,Kos);
theta_est1(1,1)=0;
theta_est2=zeros(i,Kos);
theta_est2(1,1)=0;
for index=1:i
    for n=1:Kos
        theta_est1(index,n)=(2*(n-1)*pi)/(k*Kos);
        
        theta_est2(index,n)=-(2*(n-1)*pi)/(k*Kos);
    end
end
d_est=zeros(1,i);
for index=1:i
    temp1=0;
    temp2=0;
    for n=1:Kos
        x_est1(index,n)=r(index,n)*conj((cos(phase_init+theta_est1(index,n))+1j*sin(phase_init+theta_est1(index,n))));
        x_est2(index,n)=r(index,n)*conj((cos(phase_init+theta_est2(index,n))+1j*sin(phase_init+theta_est2(index,n))));
    end
    p=1/Kos*sum(x_est1(index,2:Kos).*conj(x_est1(index,1:Kos-1)));
    fe=angle(p)/pi;
    r_offset(index,:)=h*e_fsk(index,:).*exp(1j*2*(1:Kos)*pi*fe)+noise(index,:);
    for n=1:Kos
      
        s_est1(index,n)=1/Kos*r_offset(index,n)*conj((cos(phase_init+theta_est1(index,n))+1j*sin(phase_init+theta_est1(index,n))));
        s_est2(index,n)=1/Kos*r_offset(index,n)*conj((cos(phase_init+theta_est2(index,n))+1j*sin(phase_init+theta_est2(index,n))));
    end
    if s_est1(index,n)>s_est2(index,n)
        temp1=temp1+1;
    else
        temp2=temp2+1;
    end
    
end

ber=0;
for index=1:i
    if d_est(index)~=d(index)
        ber=ber+1;
    end
end

ber_rate=ber/i;
