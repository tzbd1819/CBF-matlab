clc;
clear all;
% ---------------------------������ʼ��------------------------------------
M=32;                         % ��Ԫ��Ŀ
L=1024;                      % �źų���/������ 
lam=0.15;                    % �źŲ���
d=1/2*lam;                   % ��Ԫ���
angle=[-pi/3 pi/6 35*pi/180];      % �ź�����Ƕ�
snr=-10;                             %�����
N=3;                         %��Դ��
% ---------------------------�ź�ģ��--------------------------------------
A=zeros(N,M);
for k=1:N
A(k,:)=exp(-j*2*pi/lam*([0:M-1]*d)*sin(angle(k)));   
end
A=A';                                             % ��������
SS=2*(randn(1,L)+j*randn(1,L));
Wn=[0.36 0.38];
[b,a] = fir1(512,Wn,'bandpass');      % ��ƴ�ͨ�˲���
SS1=filter(b,a,SS); 
Wn=[0.38 0.4];
[b,a] = fir1(512,Wn,'bandpass');      % ��ƴ�ͨ�˲���
SS2=filter(b,a,SS); 
Wn=[0.4 0.42];
[b,a] = fir1(512,Wn,'bandpass');      % ��ƴ�ͨ�˲���
SS3=filter(b,a,SS); 
SS=zeros(N,L);   SS=[SS1;SS2;SS3];
RS=SS*SS'/L; 
gls=trace(RS)/N;
R=A*RS*A'+(gls*10^(-snr/10))*eye(M);             %������Э�������Ϊ�Խ���
% ---------------------------�Ƕȹ���--------------------------------------
theta=-pi/2:pi/3600:pi/2;       % ULA���ƽǶȱ仯�ķ�Χ��Ƶ��ѡ�� 
for k2=1:length(theta)      % �Ƕȹ���
    AA=exp(-j*2*pi/lam*([0:M-1]*d)*sin(theta(k2)))';   % ����ʸ��
    WW=AA'*R*AA;
    cbf(k2)=abs(WW);     % ����
end
cbf=10*log10(cbf/max(cbf));
% ---------------------------�ź�ģ��--------------------------------------
A2=zeros(N,M);
for t=1:N
A2(t,:)=exp(-j*2*pi/lam*([0:M-1]*d)*sin(angle(t)));   
end
A2=A2';                                             % ��������                           
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
R2=A2*RS*A2'+RN;             %�������������
% ---------------------------�Ƕȹ���--------------------------------------
theta2=-pi/2:pi/3600:pi/2;       % ULA���ƽǶȱ仯�ķ�Χ��Ƶ��ѡ�� 
for t7=1:length(theta2)      % �Ƕȹ���
    AA2=exp(-j*2*pi/lam*([0:M-1]*d)*sin(theta2(t7)))';   % ����ʸ��
    WW2=AA2'*R2*AA2;
    cbf2(t7)=abs(WW2);     % ����
end
cbf2=10*log10(cbf2/max(cbf2));

plot(theta*180/pi,cbf,'linewidth',1.5);
set(gca, 'Xtick',-90:15:90,'FontName','Times New Roman','FontSize',12)      %X������̶����ݵ�λ��
grid on;            %��������
box on;            %������߿�
xlabel('�Ƕ�/��');
ylabel('��һ���������/dB');
title('CBF�㷨');
axis([-90 90 0 0]);