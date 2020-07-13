close all
clear all
clc
f1=10;%载波1频率
f2=50;%载波2频率
%fc=1000;%载波为1000Hz
fbaud=600;%基带信号波特率，1s发送600个码元 %则一个码元持续时间TB=1/600s
delta_f=abs(f1-f2);
k=fbaud/delta_f;
Kos=10;
fs=fbaud*Kos;%采样率
Ts=1/fs;%采样时间间隔
T_total=5;%总的采样时间是5s
i=fbaud*T_total;%基带信号码元数
num=round(T_total/Ts);%真正采样点数
t1=linspace(0,T_total,i);
t=linspace(0,T_total,num);%0-5之间产生num个点
phase_init=0;%设初始相位为0
d=round(rand(1,i));%产生数字基带信号

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
title('数字基带信号d');
axis([0,T_total,-1,2]);

%(不同bit对应的不同相位θ(n),假设初始相位为0，暂时不考虑频偏）
theta=zeros(i,Kos);
%theta1=zeros(i,Kos);
%theta2=zeros(i,Kos);
for index=1:i
    for n=1:Kos
        if d(index)==1
            theta(index,n)=(2*(n-1)*pi)/(k*Kos);%d(i)=1对应的相位
            %theta2(index,n)=0;
        else
            %theta1(index,n)=0;
            theta (index,n)=(-2*(n-1)*pi)/(k*Kos);%d(i)=0对应的相位
        end
    end
end

%载波信号与调制
%freq1=cos(phase_init+theta1)+1j*sin(phase_init+theta1);
%freq2=cos(phase_init+theta2)+1j*sin(phase_init+theta2);
%e_fsk=freq1+freq2;
freq=cos(phase_init+theta)+1j*sin(phase_init+theta);
e_fsk=freq;
%通过信道
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
