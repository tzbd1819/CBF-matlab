clc;
clear all;
% ---------------------------参数初始化------------------------------------
M=32;                         % 阵元数目
L=1024;                      % 信号长度/快拍数 
lam=0.15;                    % 信号波长
d=1/2*lam;                   % 阵元间距
angle=[-pi/3 pi/6 35*pi/180];      % 信号入射角度
snr=-10;                             %信噪比
N=3;                         %信源数
% ---------------------------信号模型--------------------------------------
A=zeros(N,M);
for k=1:N
A(k,:)=exp(-j*2*pi/lam*([0:M-1]*d)*sin(angle(k)));   
end
A=A';                                             % 导向向量
SS=2*(randn(1,L)+j*randn(1,L));
Wn=[0.36 0.38];
[b,a] = fir1(512,Wn,'bandpass');      % 设计带通滤波器
SS1=filter(b,a,SS); 
Wn=[0.38 0.4];
[b,a] = fir1(512,Wn,'bandpass');      % 设计带通滤波器
SS2=filter(b,a,SS); 
Wn=[0.4 0.42];
[b,a] = fir1(512,Wn,'bandpass');      % 设计带通滤波器
SS3=filter(b,a,SS); 
SS=zeros(N,L);   SS=[SS1;SS2;SS3];
RS=SS*SS'/L; 
gls=trace(RS)/N;
R=A*RS*A'+(gls*10^(-snr/10))*eye(M);             %白噪声协方差矩阵为对角阵
% ---------------------------角度估计--------------------------------------
theta=-pi/2:pi/3600:pi/2;       % ULA估计角度变化的范围和频率选择 
for k2=1:length(theta)      % 角度估计
    AA=exp(-j*2*pi/lam*([0:M-1]*d)*sin(theta(k2)))';   % 方向矢量
    WW=AA'*R*AA;
    cbf(k2)=abs(WW);     % 角谱
end
cbf=10*log10(cbf/max(cbf));
% ---------------------------信号模型--------------------------------------
A2=zeros(N,M);
for t=1:N
A2(t,:)=exp(-j*2*pi/lam*([0:M-1]*d)*sin(angle(t)));   
end
A2=A2';                                             % 导向向量                           
RN=zeros(M);
for t1=1:M
    for t2=1:M
        if t1<t2
            xishu=-1;
        else 
            xishu=1;
        end
        RN(t1,t2)=(gls*10^(-snr/10))*0.9^(abs(t1-t2))*exp(xishu*j*pi/4*(t1-t2)^2);
    end
end
R2=A2*RS*A2'+RN;             %加相关噪声矩阵
% ---------------------------角度估计--------------------------------------
theta2=-pi/2:pi/3600:pi/2;       % ULA估计角度变化的范围和频率选择 
for t7=1:length(theta2)      % 角度估计
    AA2=exp(-j*2*pi/lam*([0:M-1]*d)*sin(theta2(t7)))';   % 方向矢量
    WW2=AA2'*R2*AA2;
    cbf2(t7)=abs(WW2);     % 角谱
end
cbf2=10*log10(cbf2/max(cbf2));

plot(theta*180/pi,cbf,'linewidth',1.5);
set(gca, 'Xtick',-90:15:90,'FontName','Times New Roman','FontSize',12)      %X坐标轴刻度数据点位置
grid on;            %加网格线
box on;            %加坐标边框
xlabel('角度/°');
ylabel('归一化波束输出/dB');
title('CBF算法');
axis([-90 90 0 0]);